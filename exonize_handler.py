from sqlite_utils import *
import gffutils                                        # for creating/loading DBs
import subprocess                                      # for calling gffread
import portion as P                                    # for working with intervals
import os                                              # for working with files
import time                                            # for sleeping between BLAST calls and for timeout on DB creation
import random                                          # for random sleep
import re                                              # regular expressions for genome masking
import tempfile                                        # for creating temporary files
from Bio import SeqIO                                  # for reading FASTA files
from tqdm import tqdm                                  # progress bar
from multiprocessing.pool import ThreadPool            # for parallelization
from Bio.Blast import NCBIXML                          # for parsing BLAST results
from datetime import datetime as dt
import sys


class Exonize(object):
    def __init__(self,
                 gff_file_path,
                 genome_path,
                 specie_identifier,
                 results_db_name='',
                 spec_attribute='ID',
                 save_input_files=False,
                 verbose=True,
                 hard_masking=False,
                 evalue_threshold=1e-2,
                 sleep_max_seconds=5,
                 min_exon_length=30,
                 cds_overlapping_threshold=0.9,
                 masking_perc_threshold=0.8,
                 self_hit_threshold=0.5,
                 batch_number=100,
                 threads=7,
                 timeout_db=160):

        self.genome = None                                              # genome sequence
        self.gene_hierarchy_dict = None                                 # gene hierarchy dictionary (gene -> transcript -> exon)
        self.db_features = None                                         # features in the database
        self.old_filename = None                                        # old filename (if GTF file is provided)
        self.db = None                                                  # database object (gffutils)
        self.specie_identifier = specie_identifier                      # specie identifier
        self.verbose = verbose                                          # verbose mode (True/False)
        self.in_file_path = gff_file_path                               # input file path (GFF/GTF)
        self.genome_path = genome_path                                  # genome path (FASTA)
        self.hard_masking = hard_masking                                # hard masking (True/False)
        self.secs = sleep_max_seconds                                   # max seconds to sleep between BLAST calls
        self.min_exon_len = min_exon_length                             # minimum exon length (bp)
        self.timeout_db = timeout_db                                    # timeout for creating the database (seconds)
        self.evalue = evalue_threshold                                  # e-value threshold for BLAST calls
        self.batch_number = batch_number                                # batch number for BLAST calls
        self.threads = threads                                          # number of threads for BLAST calls
        self.id_spec_attribute = spec_attribute                         # reference for writting intron annotations with gffutils
        self.cds_overlapping_threshold = cds_overlapping_threshold      # CDS overlapping threshold (0-1)
        self.self_hit_threshold = self_hit_threshold                    # self-hit threshold (0-1)
        self.results_db = results_db_name                               # results database name
        self.stop_codons = ["TAG", "TGA", "TAA"]
        self.db_path = f'{self.specie_identifier}_genome_annotations.db'
        if self.results_db == '':
            self.results_db = f'{self.specie_identifier}_results.db'
        self.UTR_features = ['five_prime_UTR', 'three_prime_UTR']
        self.gene_hierarchy_path = f"{self.specie_identifier}_gene_hierarchy.pkl"
        self.feat_of_interest = ['CDS', 'exon', 'intron'] + self.UTR_features
        self.masking_perc_threshold = masking_perc_threshold
        self.save_input_files = save_input_files
        self.logs = []

    @staticmethod
    def reverse_sequence_bool(gene_strand_: str):
        """
        reverse_sequence_bool checks if the gene is in the negative strand and returns True if it is.
        :param gene_strand_: strand
        """
        return gene_strand_ == '-'

    def create_parse_or_update_database(self) -> None:
        """
        create_parse_or_update_database is a function that in the absence of a database it:
        (i)   Provided a GTF file, it converts the file to a gff format.
        (ii)  Creates DB with the gff file by means of the gffutils library.
        Once the db exists it is loaded and the following step is performed:
        (i) Verifies that the database contains intron annotations, if not, it attempts to write them.
        """
        def convert_gtf_to_gff() -> None:
            """
            Convert a GTF file to GFF format using gffread. Flags description:
            -'O': This flag is used to enable the output of the file in GFF3 format.
            -'o': This flag is used to specify the output file name.
            """
            gffread_command = ["gffread", self.old_filename, "-O", "-o", self.in_file_path]
            subprocess.call(gffread_command)

        def create_database() -> None:
            """
            create_database is a function that creates a gffutils database from a GFF3 file.
            - dbfn: path to the database file
            - force: if True, the database will be overwritten if it already exists
            - keep_order: if True, the order of the features in the GFF file will be preserved
            - merge_strategy: if 'create_unique', the database will be created with unique IDs
            - sort_attribute_values: if True, the attribute values will be sorted
            - disable_infer_genes: if True, the function will not attempt to automatically infer gene features
            - disable_infer_transcripts: if True, the function will not attempt to automatically infer transcript features
            """
            try:
                if self.verbose:
                    print("- Creating annotations database", end=" ")
                self.db = gffutils.create_db(self.in_file_path,
                                             dbfn=self.db_path,
                                             force=True,
                                             keep_order=True,
                                             merge_strategy='create_unique',
                                             sort_attribute_values=True,
                                             disable_infer_genes=True,
                                             disable_infer_transcripts=True)
                if self.verbose:
                    print("Done!")
            except ValueError as e:
                print(' ')
                print('---------------------ERROR-------------------------------')
                print(f"Incorrect genome annotations file  {e}")
                print('---------------------------------------------------------')
                sys.exit()

        def search_create_intron_annotations() -> None:
            """
            search_create_intron_annotations is a function that verifies that the gffutils database contains intron
            annotations, if not, it attempts to write them. Some of the db.update() parameters description:
            - make_backup: if True, a backup of the database will be created before updating it.
            - id_spec: dictionary with the following structure: {feature_type: [intron_id]} where the intron_id function
              returns the intron identifier based on the feature attribute (self.id_spec_attribute) present in all annotations
              in the gff file. Common choices are: "ID" or "Parent".
            """
            if 'intron' not in self.db_features:
                print('---------------------WARNING------------------------------')
                print("-The genomic annotations do not contain intron annotations")
                print('----------------------------------------------------------')
                if self.verbose:
                    print(f"- Attempting to write intron annotations in database:", end=" ")
                try:
                    def intron_id(f: dict) -> str:
                        """
                        returns the new intron identifier
                        : param f: annotation dictionary.
                        """
                        return ','.join(f[self.id_spec_attribute])

                    introns_list = list(self.db.create_introns())
                    self.db.update(introns_list, id_spec={'intron': [intron_id]}, make_backup=False)
                except ValueError as e:
                    print("failed to write intron annotations in database, "
                          "please try with a different spec attribute or provide"
                          " a GFF3 file with intron annotations")
                    print(e)
                    sys.exit()
                if self.verbose:
                    print("Done!")

        if not os.path.exists(self.db_path):
            if 'gtf' in self.in_file_path:
                self.old_filename = self.in_file_path
                self.in_file_path = f"{self.old_filename.rsplit('.gtf')[0]}.gff"
                convert_gtf_to_gff()
                print('the GTF file has been converted into a GFF3 file')
                print(f'with filename: {self.in_file_path}')
            create_database()
        if not self.db:
            print("- Reading annotations database:", end=" ")
            self.load_db()
            print("Done!")
        self.db_features = list(self.db.featuretypes())
        search_create_intron_annotations()

    def load_db(self) -> None:
        """
        load_db is a function that loads a gffutils database.
        - dbfn: path to the database file
        - keep_order: This is a parameter that is passed when creating the FeatureDB instance. When keep_order is set to True,
         the order of attributes in the GFF/GTF file will be preserved when they are retrieved from the database.
        """
        try:
            self.db = gffutils.FeatureDB(self.db_path, keep_order=True)
        except ValueError as e:
            print(' ')
            print('---------------------ERROR-------------------------------')
            print(f"Incorrect data base path  {e}")
            print('---------------------------------------------------------')
            sys.exit()

    def read_genome(self) -> None:
        """
        read_genome is a function that reads a FASTA file and stores the masked/unmasked genome sequence in a dictionary.
        The dictionary has the following structure: {chromosome: sequence}
        """
        try:
            tic_genome = time.time()
            if self.verbose:
                print("- Reading genome file:", end=" ")
            parse_genome = SeqIO.parse(open(self.genome_path), 'fasta')
            if self.hard_masking:
                self.genome = {fasta.id: re.sub('[a-z]', 'N', str(fasta.seq)) for fasta in parse_genome}
            else:
                self.genome = {fasta.id: str(fasta.seq) for fasta in parse_genome}
            hms_time = dt.strftime(dt.utcfromtimestamp(time.time() - tic_genome), '%H:%M:%S')
            print(f"Done! [{hms_time}]")
        except (ValueError, FileNotFoundError) as e:
            print(' ')
            print('---------------------ERROR-------------------------------')
            print(f"Incorrect genome file path {e}")
            print('---------------------------------------------------------')
            sys.exit()

    def dump_masking_logs(self) -> None:
        """
        Dumps the logs recorded in the attribute `self.logs` to a file named `exonize_logs.txt`
        """
        if not self.logs:
            self.logs = ['Nothing to report']
        with open(f'exonize_logs.txt', 'w') as f:
            f.write('\n'.join(self.logs))

    def create_gene_hierarchy_dict(self) -> None:
        """
        Constructs a nested dictionary to represent the hierarchical structure and attributes of genes and
        their related mRNA transcripts based on genomic feature data. The created hierarchy is stored in the attribute
        `self.gene_hierarchy_dict` and is also saved as a pickle file.
        Note:
        - GFF coordinates are 1-based. Thus, 1 is subtracted from the start position to convert them to 0-based coordinates.
        - If the gene is in the negative strand the direction of transcription and translation is opposite to the direction
         the DNA sequence is represented meaning that translation starts from the last CDS
        Structure of `self.gene_hierarchy_dict`:
        {gene_id_1: {'coord': gene_coord_1, 'chrom': chromosome_1, 'strand': strand_1,
                    'mRNAs': {mRNA_id_1: {'coord': mRNA_coord_1, 'strand': strand_1,
                             'structure': [{'id': feature_id_1, 'coord': feature_coord_1, 'frame': frame_1,
                             'type': feature_type_1, 'attributes': attribute_dict_1},...]
                             },...}},...}
        """
        self.gene_hierarchy_dict = {}
        for gene in self.db.features_of_type('gene'):
            mrna_transcripts = [mRNA_t for mRNA_t in self.db.children(gene.id, featuretype='mRNA', order_by='start')]
            if mrna_transcripts:
                gene_coord = P.open(gene.start - 1, gene.end)
                mrna_dict = dict(coord=gene_coord, chrom=gene.chrom, strand=gene.strand, mRNAs={})
                for mrna_annot in mrna_transcripts:
                    mrna_coord = P.open(mrna_annot.start - 1, mrna_annot.end)
                    mrna_dict['mRNAs'][mrna_annot.id] = dict(coord=mrna_coord, strand=gene.strand, structure=[])
                    temp_mrna_transcript = []
                    for child in self.db.children(mrna_annot.id, featuretype=self.feat_of_interest, order_by='start'):
                        coord = P.open(child.start - 1, child.end)
                        if coord:
                            temp_mrna_transcript.append(dict(id=child.id,  # ID attribute
                                                             coord=coord,  # ID coordinate starting at 0
                                                             # One of '0', '1' or '2'. '0' indicates that the x base of
                                                             # the feature is the x base of a codon
                                                             frame=child.frame,
                                                             type=child.featuretype,   # feature type name
                                                             attributes=dict(child.attributes))  # feature type name
                                                        )
                    reverse = self.reverse_sequence_bool(gene.strand)
                    mrna_dict['mRNAs'][mrna_annot.id]['structure'] = sort_list_intervals_dict(temp_mrna_transcript, reverse)
                self.gene_hierarchy_dict[gene.id] = mrna_dict
        dump_pkl_file(self.gene_hierarchy_path, self.gene_hierarchy_dict)

    def prepare_data(self) -> None:
        """
        prepare_data is a wrapper function that:
        (i)   creates the database with the genomic annotations (if it does not exist)
        (ii)  reads or creates the gene hierarchy dictionary
        (iii) reads the genome sequence
        (iv)  connects or creates the results database
        """
        self.create_parse_or_update_database()
        if os.path.exists(self.gene_hierarchy_path):
            self.gene_hierarchy_dict = read_pkl_file(self.gene_hierarchy_path)
        else:
            self.create_gene_hierarchy_dict()
        self.read_genome()
        connect_create_results_db(self.results_db, self.timeout_db)
        if self.save_input_files:
            print('---------------------WARNING------------------------------')
            print("-All tblastx input and output files will be saved. "
                  "This may take a large amount of disk space. "
                  "You can disable this option by setting the save_input_files parameter to False.")
            print('----------------------------------------------------------')

    def execute_tblastx(self, query_filename: str, target_filename: str, output_file: str):
        """
        execute_tblastx is a function that executes a tblastx search with the following parameters:
        - tblastx: A BLAST tool that compares the six-frame translations of a nucleotide query sequence
        against the six-frame translations of a nucleotide sequence database.
        - query: query file name
        - subject: subject file name
        - evalue: Expectation value (E) threshold for reporting matches against database sequences.
        - qcov_hsp_perc: Percent query coverage per hsp (high-scoring pair).
        - outfmt: alignment view options - 5: XML output format
        - out: output file name
        """
        tblastx_command = ['tblastx',
                           '-query', query_filename,
                           '-subject', target_filename,
                           '-evalue', str(self.evalue),
                           '-qcov_hsp_perc', str(self.cds_overlapping_threshold * 100),
                           '-outfmt', '5',  # XML output format
                           '-out', output_file]
        subprocess.run(tblastx_command)

    def align_CDS(self, gene_id: str, query_seq: str, hit_seq: str, query_coord: P.Interval) -> dict:
        """
        align_CDS is a function that performs a tblastx search between a query sequence (CDS) and a target sequence (gene).
        :param gene_id: gene identifier
        :param query_seq: query sequence (CDS)
        :param hit_seq: target sequence (gene)
        :param query_coord: query coordinates (CDS) interval
        :return: dict with the following structure: {hsp_id: {'score': '', 'bits': '','evalue': '',...}}
        """
        def tblastx_with_saved_io(identifier_: str, gene_id_: str, hit_seq_: str, query_seq_: str,
                                  query_coord_: P.Interval, gene_coord_: P.Interval):
            output_file = f'output/{identifier_}_output.xml'
            if not os.path.exists(output_file):
                query_filename = f'input/{identifier_}_query.fa'
                target_filename = f'input/{gene_id_}_target.fa'
                if not os.path.exists(target_filename):
                    dump_fasta_file(target_filename, {f"{gene_id_}": hit_seq_})
                dump_fasta_file(query_filename, {identifier_: query_seq_})
                self.execute_tblastx(query_filename, target_filename, output_file)
            with open(output_file, "r") as result_handle:
                blast_records = NCBIXML.parse(result_handle)
                try:
                    temp_ = self.parse_tblastx_output(blast_records, query_coord_, gene_coord_)
                except Exception as e:
                    print(e)
                    sys.exit()
            return temp_

        def execute_tblastx_using_tempfiles(hit_seq_: str, query_seq_: str,
                                            query_coord_: P.Interval, gene_coord_: P.Interval):
            with tempfile.TemporaryDirectory() as tmpdirname:
                query_filename = f'{tmpdirname}/query.fa'
                target_filename = f'{tmpdirname}/target.fa'
                dump_fasta_file(query_filename, {'query': query_seq_})
                dump_fasta_file(target_filename, {'target': hit_seq_})
                output_file = f'{tmpdirname}/output.xml'
                self.execute_tblastx(query_filename, target_filename, output_file)
                with open(output_file, 'r') as result_handle:
                    blast_records = NCBIXML.parse(result_handle)
                    try:
                        temp_ = self.parse_tblastx_output(blast_records, query_coord_, gene_coord_)
                    except Exception as e:
                        print(e)
                        sys.exit()
                return temp_

        chrom = self.gene_hierarchy_dict[gene_id]['chrom']
        gene_coord = self.gene_hierarchy_dict[gene_id]['coord']
        identifier = f'{gene_id}_{chrom}_{str(query_coord.lower)}_{query_coord.upper}'.replace(':', '_')
        if self.save_input_files:
            temp = tblastx_with_saved_io(identifier, gene_id, hit_seq, query_seq, query_coord, gene_coord)
        else:
            temp = execute_tblastx_using_tempfiles(hit_seq, query_seq, query_coord, gene_coord)
        return temp

    def parse_tblastx_output(self, blast_records: dict, q_coord: P.Interval, hit_coord: P.Interval) -> dict:
        """
        the parse_tblastx_output function parses the output of a tblastx search, where a single sequence (CDS)
        has been queried against a single target (gene). Meaning that we only expect to find one BLAST record.
        We only  consider hits that:
            (i)   have an e-value lower than the threshold,
            (ii)  have a minimum alignment length percentage of the query sequence and
            (iii) that do not overlap with the query sequence (self-hit), with a maximum overlap (self.self_hit_threshold)
                  of the query sequence.
        :param blast_records: 'Record' object with blast records
        :param q_coord: query coordinates (CDS) interval
        :param hit_coord: hit coordinates
        :return: dict with the following structure: {target_id {hsp_id: {'score': '', 'bits': '','evalue': '',...}}}
        """
        res_tblastx = {}
        # since we are performing a single query against a single subject, there's only one blast_record
        for blast_record in blast_records:
            if len(blast_record.alignments) == 0:
                continue
            alignment = blast_record.alignments[0]  # Assuming only one alignment per blast_record
            if len([aln for aln in blast_record.alignments]) > 1:
                print("More than one alignment per blast_record")
            for hsp_idx, hsp in enumerate(alignment.hsps):
                blast_target_coord = P.open((hsp.sbjct_start - 1) + hit_coord.lower, hsp.sbjct_end + hit_coord.lower)
                target_query_overlap_percentage = get_overlap_percentage(blast_target_coord, q_coord)
                query_target_overlap_percentage = get_overlap_percentage(q_coord, blast_target_coord)
                if (query_target_overlap_percentage < self.self_hit_threshold  # self-hits are not allowed (max 50% overlap)
                        and target_query_overlap_percentage < self.self_hit_threshold):
                    res_tblastx[hsp_idx] = get_hsp_dict(hsp)
        return res_tblastx

    def get_gene_tuple(self, gene_id: str, has_dup_bin: int) -> tuple:
        """
        get_gene_tuple is a function that given a gene_id, returns a tuple with the following structure:
        (gene_id, chromosome, strand, start_coord, end_coord, 1 if it has a duplication event 0 otherwise)
        """
        gene_coord = self.gene_hierarchy_dict[gene_id]['coord']
        return (gene_id,
                self.gene_hierarchy_dict[gene_id]['chrom'],
                (self.gene_hierarchy_dict[gene_id]['strand']),
                gene_coord.lower,
                gene_coord.upper,
                has_dup_bin)

    def get_candidate_CDS_coords(self, gene_id: str) -> list:
        """
        get_candidate_CDS_coords is a function that given a gene_id, collects all the CDS coordinates with a length
        greater than the minimum exon length (self.min_exon_len) across all transcript.
        If there are overlapping CDS coordinates, they are resolved according to the following criteria:
            - a: if they overlap by more than the overlapping threshold (self.cds_overlapping_threshold) of both CDS
             lengths the shortest CDS is selected.
            - b: both CDSs are selected otherwise
        The rationale behind choosing the shortest CDS is that it will reduce exclusion of short duplications due
        to coverage thresholds.
        :param gene_id: gene identifier
        :return: list of representative CDS coordinates for the gene across all transcripts
        """
        def get_intervals_overlapping_list(intv_list: list, overlap_threshold: float) -> list:
            """
            get_intervals_overlapping_list is a function that given a list of intervals, returns a list of tuples
            with the overlapping pairs described in case a (get_candidate_CDS_coords).
            :param intv_list: list of intervals
            :param overlap_threshold: threshold for resolving overlaps
            :return: list of tuples with pairs of overlapping intervals
            """
            return [(feat_interv, intv_list[idx + 1])
                    for idx, feat_interv in enumerate(intv_list[:-1])
                    if all([get_overlap_percentage(feat_interv, intv_list[idx + 1]) >= overlap_threshold,
                            get_overlap_percentage(intv_list[idx + 1], feat_interv) >= overlap_threshold])]

        def resolve_overlaps_coords_list(coords_list, overlap_threshold):
            """
            resolve_overlaps_coords_list is a function that given a list of coordinates, resolves overlaps according
            to the criteria described above. Since the pairs described in case b (get_candidate_CDS_coords) are not
            considered to be overlapping, they are both included in the search.
            :param coords_list: list of coordinates
            :param overlap_threshold: threshold for resolving overlaps
            :return: list of coordinates without overlaps
            """
            new_list = []
            overlaps_list = get_intervals_overlapping_list(coords_list, overlap_threshold)
            if overlaps_list:
                # We only process the first pair of overlapping intervals since the resolved overlap could also overlap
                # with the next interval in the list.
                intv_a, intv_b = overlaps_list[0]
                for idx, coord in enumerate(coords_list):
                    if coord != intv_a:
                        new_list.append(coord)
                    else:
                        shorter, _ = get_shorter_longer_interv(intv_a, intv_b)
                        new_list.append(shorter)
                        new_list.extend(coords_list[idx + 2:])
                        return resolve_overlaps_coords_list(new_list, overlap_threshold)
            else:
                return coords_list

        CDS_coords_list = list(set([i['coord'] for mrna_id, mrna_annot
                                    in self.gene_hierarchy_dict[gene_id]['mRNAs'].items()
                                    for i in mrna_annot['structure']
                                    if(i['type'] == 'CDS'
                                       and (i['coord'].upper - i['coord'].lower) >= self.min_exon_len)]))
        CDS_coords_list = sorted(CDS_coords_list, key=lambda x: (x.lower, x.upper))
        if CDS_coords_list:
            return resolve_overlaps_coords_list(CDS_coords_list, self.cds_overlapping_threshold)
        return []

    def find_coding_exon_duplicates(self, gene_id: str) -> None:
        """
        find_coding_exon_duplicates is a function that given a gene_id, performs a tblastx for each representative CDS
        (see get_candidate_CDS_coords). If the tblastx search returns hits, they are stored in the "results" database,
        otherwise the gene is recorded as having no duplication event.
        :param gene_id: gene identifier
        """
        def get_sequence_and_check_for_masking(chrom_: str, gene_id_: str, coords_: P.Interval, gene_strand_: str, type_='gene'):
            """
            get_sequence_and_check_for_masking is a function that retrieves a gene/CDS sequence from the genome and checks
            if it is hardmasked. If the sequence is a gene and the percentage of hardmasking is greater than the
            threshold (self.masking_perc_threshold) the CDS duplication search is aborted and the gene is recorded
            as having no duplication event. If the sequence is a CDS and the percentage of hardmasking is greater than
            the threshold, the CDS is not queried. Logs for hardmasked genes and CDS are stored in the logs attribute
            and later dumped into a file.
            :param chrom_: chromosome identifier
            :param gene_id_: gene identifier
            :param coords_: coordinates
            :param gene_strand_: strand
            :param type_: type of sequence (gene or CDS)
            """
            try:
                seq_ = self.genome[chrom_][coords_.lower:coords_.upper]
                if self.reverse_sequence_bool(gene_strand_):
                    seq_ = seq_.reverse_complement()
                masking_perc = round(sequence_masking_percentage(seq_), 2)
                if masking_perc > self.masking_perc_threshold:
                    seq_ = ''
                    if type_ == 'gene':
                        self.logs.append((f'Gene {gene_id_} in chromosome {chrom_} '
                                          f'and coordinates {str(coords_.lower)}, {str(coords_.upper)} is hardmasked.'))
                        insert_gene_ids_table(self.results_db, self.timeout_db, self.get_gene_tuple(gene_id_, 0))
                    if type_ == 'CDS':
                        self.logs.append((f'Gene {gene_id_} - {masking_perc * 100} '
                                          f'% of CDS {cds_coord} located in chromosome {chrom_} is hardmasked.'))
            except KeyError as e_:
                print(f'Either there is missing a chromosome in the genome file '
                      f'or the chromosome identifiers in the GFF3 and FASTA files do not match {e_}')
                sys.exit()
            return seq_

        CDS_blast_dict = {}
        time.sleep(random.randrange(0, self.secs))
        chrom, gene_coord, gene_strand = (self.gene_hierarchy_dict[gene_id]['chrom'],
                                          self.gene_hierarchy_dict[gene_id]['coord'],
                                          self.gene_hierarchy_dict[gene_id]['strand'])
        gene_seq = get_sequence_and_check_for_masking(chrom, gene_id, gene_coord, gene_strand, type_='gene')
        if gene_seq:
            reverse = self.reverse_sequence_bool(gene_strand)
            CDS_coords_list = self.get_candidate_CDS_coords(gene_id)
            CDS_coords_list = sorted(CDS_coords_list, key=lambda x: (x.lower, x.upper), reverse=reverse)
            for cds_coord in CDS_coords_list:
                cds_seq = get_sequence_and_check_for_masking(chrom, gene_id, cds_coord, gene_strand, type_='CDS')
                if cds_seq:
                    tblastx_o = self.align_CDS(gene_id, cds_seq, gene_seq, cds_coord)
                    if temp:
                        CDS_blast_dict[cds_coord] = tblastx_o
            if CDS_blast_dict:
                self.insert_fragments_table(gene_id, CDS_blast_dict)
            else:
                insert_gene_ids_table(self.results_db, self.timeout_db, self.get_gene_tuple(gene_id, 0))

    def insert_fragments_table(self, gene_id: str, blast_cds_dict: dict) -> None:
        """
        insert_fragments_table is a function that given a gene_id and a dictionary with the tblastx results for each
        CDS, inserts in the (i) Genes, and (ii) Fragments table of the results (self.results_db) database.
        These two steps are done in paralell to avoid incomplete data in case of a crash.
        A single entry is recorded in the Genes table, while multiple entries can be recorded in the Fragments table.
        The fragments refer to the tblastx HSPs (high-scoring pairs) that have passed the filtering criteria.
        The fragments tuple has the following structure:
        (gene_id, cds_coord_start, cds_coord_end, hit_query_frame, hit_query_strand, hit_target_frame, hit_target_strand,
         score, bits, evalue, alignment_length, query_start, query_end, target_start, target_end, query_aln_prot_seq,
         target_aln_prot_seq, match (alignment sequence), query_num_stop_codons (count of "*" in target alignment),
         target_num_stop_codons (count of "*" in target alignment))
        :param gene_id: gene identifier
        :param blast_cds_dict: dictionary with the tblastx results. The dictionary has the following structure:
        {cds_coord: {hsp_idx: {'score': '', 'bits': '','evalue': '',...}}}
        """
        tuple_list = [get_fragment_tuple(gene_id, cds_coord, blast_hits, hsp_idx)
                      for cds_coord, blast_hits in blast_cds_dict.items()
                      for hsp_idx, hsp_dict in blast_hits.items()]
        insert_fragments_table_param, insert_gene_table_param = insert_fragments_calls()
        with sqlite3.connect(self.results_db, timeout=self.timeout_db) as db:
            cursor = db.cursor()
            cursor.execute(insert_gene_table_param, self.get_gene_tuple(gene_id, 1))
            cursor.executemany(insert_fragments_table_param, tuple_list)
            db.commit()

    def find_overlapping_annot(self, trans_dict: dict, cds_intv: P.Interval) -> list:
        """
        find_overlapping_annot is a function that given a transcript dictionary and a CDS interval, returns a list of
        tuples with the following structure: (CDS_id, CDS_coord, CDS_frame) for all CDSs that overlap with the query CDS
        interval.
        Note:
        - The set of representative CDS only contains information about the coordinates of the CDSs, not their IDs.Hence,
        if we find hits for a CDS, we need to identify the CDS in the gene transcript. Recall that not all CDS are part
        of all transcripts (e.g. alternative splicing).
        - Following the classification criteria described in 'get_candidate_CDS_coords', we identify the query CDS described in
        case a (get_candidate_CDS_coords), i.e, when they are considerd to overlap.
        """
        return [(i['id'], i['coord'], int(i['frame'])) for i in trans_dict['structure']
                if i['type'] == 'CDS'
                and all([get_overlap_percentage(i['coord'], cds_intv) >= self.cds_overlapping_threshold,
                         get_overlap_percentage(cds_intv, i['coord']) >= self.cds_overlapping_threshold])]

    def identify_full_length_duplications(self) -> None:
        def insert_in_results_db(tuples_fld, tuples_oe, tuples_te):
            instert_full_length_event(self.results_db, self.timeout_db, tuples_fld)
            instert_obligatory_event(self.results_db, self.timeout_db, tuples_oe)
            instert_truncation_event(self.results_db, self.timeout_db, tuples_te)

        rows = query_filtered_full_duplication_events(self.results_db, self.timeout_db)
        tuples_full_length_duplications, tuples_obligatory_events, tuples_truncation_events = [], [], []
        for row in rows:
            fragment_id, gene_id, gene_s, gene_e, cds_s, cds_e, query_s, query_e, target_s, target_e, evalue = row
            cds_intv = P.open(cds_s, cds_e)
            target_intv = P.open(target_s + gene_s, target_e + gene_s)
            for mrna, trans_dict in self.gene_hierarchy_dict[gene_id]['mRNAs'].items():
                found = False
                trans_coord = trans_dict['coord']
                neither, query, target, target_full, target_insertion, target_trunctation, both = 0, 0, 0, 0, 0, 0, 0
                query_CDS, target_CDS = "-", "-"
                annot_target_start, annot_target_end, target_t, query_CDS_frame, target_CDS_frame = None, None, None, None, None
                # ####### QUERY ONLY - FULL LENGTH #######
                # account for allowed shifts in the resolve_overlaps_coords_list function
                query_only = self.find_overlapping_annot(trans_dict, cds_intv)
                if len(query_only) > 1:
                    print(f'overlapping query CDSs: {query_only}')
                    continue
                elif query_only:
                    query_CDS, _, query_CDS_frame = query_only[0]
                    query = 1
                    target_t = 'QUERY_ONLY'
                # ###### CHECK: TARGET REGION NOT IN mRNA #######
                if target_intv.upper < trans_coord.lower or trans_coord.upper < target_intv.lower:
                    if (query + target) == 0:
                        neither = 1
                    tuples_full_length_duplications.append((fragment_id, gene_id, mrna,
                                                            cds_s, cds_e, query_CDS, query_s, query_e,
                                                            "OUT_OF_MRNA",
                                                            target_CDS, annot_target_start, annot_target_end,
                                                            target_s, target_e,
                                                            neither, query, target, both, evalue))
                    continue
                # ####### TARGET ONLY - FULL LENGTH #######
                target_only = self.find_overlapping_annot(trans_dict, target_intv)
                if len(target_only) > 1:
                    print(f'overlapping query CDSs: {target_only}')
                    continue
                # ####### TARGET ONLY #######
                elif target_only:
                    target_CDS, t_CDS_coord, target_CDS_frame = target_only[0]
                    target_full = 1
                    found = True
                    target_t = "FULL"
                    annot_target_start, annot_target_end = t_CDS_coord.lower, t_CDS_coord.upper
                # ####### INSERTION #######
                if target_full == 0:
                    insertion_CDS = filter_structure(trans_dict['structure'], target_intv, 'CDS')
                    if insertion_CDS:
                        target_CDS, t_CDS_coord = insertion_CDS[0]
                        target_insertion = 1
                        target_t = "INS_CDS"
                        found = True
                        annot_target_start, annot_target_end = t_CDS_coord.lower, t_CDS_coord.upper
                    else:
                        insertion_UTR = filter_structure(trans_dict['structure'], target_intv, 'UTR')
                        if insertion_UTR:
                            target_CDS, t_CDS_coord = insertion_UTR[0]
                            target_t = "INS_UTR"
                            found = True
                            annot_target_start, annot_target_end = t_CDS_coord.lower, t_CDS_coord.upper
                        else:
                            insertion_intron = filter_structure(trans_dict['structure'], target_intv, 'intron')
                            if insertion_intron:
                                target_CDS, t_CDS_coord = insertion_intron[0]
                                target_t = "DEACTIVATED"
                                found = True
                                annot_target_start, annot_target_end = t_CDS_coord.lower, t_CDS_coord.upper
                    # ####### TRUNCATION #######
                    if not found:
                        intv_dict = get_interval_dictionary(trans_dict['structure'], target_intv, trans_coord)
                        if intv_dict:
                            target_t = "TRUNC"
                            target_CDS, annot_target_start, annot_target_end = None, None, None
                            for seg_b, value in intv_dict.items():
                                coord_b = value['coord']
                                tuples_truncation_events.append(
                                    (fragment_id, gene_id, mrna, trans_coord.lower, trans_coord.upper,
                                     query_CDS, cds_s, cds_e, query_s, query_e,
                                     target_s, target_e, value['id'], value['type'],
                                     coord_b.lower, coord_b.upper, seg_b.lower, seg_b.upper))
                target = target_full + target_insertion
                # ####### FULL LENGTH DUPL: BOTH #######
                if query + target == 2:
                    both = 1
                    query, target = 0, 0
                    tuples_obligatory_events.append((fragment_id, gene_id, mrna,
                                                    trans_coord.lower, trans_coord.upper,
                                                    cds_s, cds_e,
                                                    query_CDS, query_CDS_frame, query_s, query_e,
                                                    target_CDS, target_CDS_frame, annot_target_start, annot_target_end,
                                                    target_s, target_e, target_t))
                elif query + target == 0:
                    neither = 1
                    target_t = 'NEITHER'
                tuples_full_length_duplications.append((fragment_id, gene_id, mrna,
                                                        cds_s, cds_e, query_CDS, query_s, query_e,
                                                        target_t, target_CDS, annot_target_start, annot_target_end,
                                                        target_s, target_e, neither, query, target, both, evalue))

                insert_in_results_db(tuples_fld, tuples_oe, tuples_te)

    def assign_pair_ids(self, full_matches) -> list:
        fragments = []
        skip_frag = []
        counter = 1
        with tqdm(total=len(full_matches), position=0, leave=True) as progress_bar:
            for frag_a in full_matches:
                frag_id_a, gene_id_a, q_s_a, q_e_a, t_s_a, t_e_a, event_type_a = frag_a
                t_intv_a = P.open(t_s_a, t_e_a)
                q_intv_a = P.open(q_s_a, q_e_a)
                if frag_id_a not in skip_frag:
                    # candidates = query_candidates(self.results_db, self.timeout_db, (gene_id_a, frag_id_a))
                    candidates = [i for i in full_matches if i[1] == gene_id_a and i[0] not in [*skip_frag, frag_id_a]]
                    if candidates:
                        temp_cand = []
                        for frag_b in candidates:
                            frag_id_b, gene_id_b, q_s_b, q_e_b, t_s_b, t_e_b, event_type_b = frag_b
                            t_intv_b = P.open(t_s_b, t_e_b)
                            q_intv_b = P.open(q_s_b, q_e_b)
                            overlapping_pairs = find_overlapping_pairs(q_intv_a, q_intv_b, t_intv_a, t_intv_b)
                            reciprocal_pairs = find_reciprocal_pairs(q_intv_a, q_intv_b, t_intv_a, t_intv_b)
                            if overlapping_pairs or reciprocal_pairs:
                                # Insertion events
                                if "INS_CDS" in event_type_a and "TRUNC" in event_type_b:
                                    # q_1, t_2 and q_2, t_1 have to overlap
                                    if all(perc > 0 for perc in reciprocal_pairs):
                                        temp_cand.append(frag_id_b)
                                    # q_1, q_2 and t_2, t_2 have to overlap
                                    elif all(perc >= self.cds_overlapping_threshold for perc in overlapping_pairs):
                                        temp_cand.append(frag_id_b)
                                # Full events, meaning that the target CDS and query CDS overlap for their greater part
                                elif any(all(perc >= self.cds_overlapping_threshold for perc in pair)  # full dups
                                         for pair in [reciprocal_pairs, overlapping_pairs]):
                                    temp_cand.append(frag_id_b)
                        if temp_cand:
                            skip_frag.extend([frag_id_a, *temp_cand])
                            fragments.extend([(counter, frag) for frag in [frag_id_a, *temp_cand]])
                            counter += 1
                progress_bar.update(1)
        return fragments

    def get_identity_and_sequence_tuples(self) -> list:
        tuples_list = []
        all_fragments = query_fragments(self.results_db, self.timeout_db)
        for fragment in all_fragments:
            (fragment_id, gene_id, gene_start, gene_end, gene_chrom,
             CDS_start, CDS_end, query_start, query_end, target_start, target_end,
             query_strand, target_strand, query_aln_prot_seq, target_aln_prot_seq) = fragment
            query_seq = self.genome[gene_chrom][CDS_start:CDS_end]
            target_seq = self.genome[gene_chrom][gene_start:gene_end]
            query_dna_seq = query_seq[query_start:query_end]
            target_dna_seq = target_seq[target_start:target_end]
            if self.reverse_sequence_bool(query_strand):
                query_dna_seq = reverse_complement(query_dna_seq)
            if self.reverse_sequence_bool(target_strand):
                target_dna_seq = reverse_complement(target_dna_seq)
            if len(query_dna_seq) != len(target_dna_seq):
                raise ValueError(f'{gene_id}: CDS {(CDS_start, CDS_end)} search - sequences must have the same length.')
            dna_identity = round(hamming_distance(query_dna_seq, target_dna_seq), 3)
            prot_identity = round(hamming_distance(query_aln_prot_seq, target_aln_prot_seq), 3)
            tuples_list.append((dna_identity, prot_identity, query_dna_seq, target_dna_seq, fragment_id))
        return tuples_list

    def run_analysis(self) -> None:
        exonize_asci_art()
        self.prepare_data()
        args_list = list(self.gene_hierarchy_dict.keys())
        processed_gene_ids = query_gene_ids_in_res_db(self.results_db, self.timeout_db)
        if processed_gene_ids:
            args_list = [i for i in args_list if i not in processed_gene_ids]
        if args_list:
            batches_list = [i for i in batch(args_list, self.batch_number)]
            tic = time.time()
            gene_n = len(list(self.gene_hierarchy_dict.keys()))
            print(f'- Starting exon duplication search for {len(args_list)}/{gene_n} genes.')
            with tqdm(total=len(args_list), position=0, leave=True) as progress_bar:
                for arg_batch in batches_list:
                    t = ThreadPool(processes=self.threads)
                    t.map(self.find_coding_exon_duplicates, arg_batch)
                    t.close()
                    t.join()
                    progress_bar.update(len(arg_batch))
            hms_time = dt.strftime(dt.utcfromtimestamp(time.time() - tic), '%H:%M:%S')
            print(f' Done! [{hms_time}]')
        else:
            print('All genes already processed if you want to re-run the analysis, '
                  'delete/rename the results DB.')
        tic = time.time()
        insert_percent_query_column_to_fragments(self.results_db, self.timeout_db)
        create_filtered_full_length_events_view(self.results_db, self.timeout_db)
        create_mrna_counts_view(self.results_db, self.timeout_db)
        print('- Classifying events', end=" ")
        tic_ce = time.time()
        self.identify_full_length_duplications()
        hms_time = dt.strftime(dt.utcfromtimestamp(time.time() - tic_ce), '%H:%M:%S')
        print(f' Done! [{hms_time}]')
        create_cumulative_counts_table(self.results_db, self.timeout_db)
        query_concat_categ_pair_list = query_concat_categ_pairs(self.results_db, self.timeout_db)
        reduced_event_types_tuples = generate_unique_events_list(query_concat_categ_pair_list, -1)
        insert_event_categ_full_length_events_cumulative_counts(self.results_db, self.timeout_db, reduced_event_types_tuples)
        identity_and_sequence_tuples = self.get_identity_and_sequence_tuples()
        insert_identity_and_dna_algns_columns(self.results_db, self.timeout_db, identity_and_sequence_tuples)
        full_matches = query_full_events(self.results_db, self.timeout_db)
        print('- Reconciling events')
        fragments = self.assign_pair_ids(full_matches)
        instert_pair_id_column_to_full_length_events_cumulative_counts(self.results_db, self.timeout_db, fragments)
        create_exclusive_pairs_view(self.results_db, self.timeout_db)
        hms_time = dt.strftime(dt.utcfromtimestamp(time.time() - tic), '%H:%M:%S')
        self.dump_masking_logs()
        print(f' Done! [{hms_time}]')
