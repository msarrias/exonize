from .sqlite_utils import *
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
import logging
import datetime


class Exonize(object):
    logger: logging.Logger

    def __init__(self,
                 gff_file_path,
                 genome_path,
                 specie_identifier,
                 results_db_name='',
                 enable_debug=False,
                 hard_masking=False,
                 soft_force=False,
                 hard_force=False,
                 evalue_threshold=1e-2,
                 sleep_max_seconds=5,
                 min_exon_length=30,
                 cds_overlapping_threshold=0.9,
                 masking_perc_threshold=0.8,
                 self_hit_threshold=0.5,
                 batch_number=100,
                 threads=7,
                 timeout_db=160):

        self.configure_logger()
        self._DEBUG_MODE = enable_debug                                 # debug mode (True/False)
        self._SOFT_FORCE = soft_force                                   # (True/False) - will remove results database if it exists
        self._HARD_FORCE = hard_force                                   # (True/False) - will remove results database, genome database and gene hierarchy
        self.genome = None                                              # genome sequence
        self.gene_hierarchy_dict = None                                 # gene hierarchy dictionary (gene -> transcript -> exon)
        self.db_features = None                                         # features in the database
        self.old_filename = None                                        # old filename (if GTF file is provided)
        self.db = None                                                  # database object (gffutils)
        self.specie_identifier = specie_identifier                      # specie identifier
        self.in_file_path = gff_file_path                               # input file path (GFF/GTF)
        self.genome_path = genome_path                                  # genome path (FASTA)
        self.hard_masking = hard_masking                                # hard masking (True/False)
        self.secs = sleep_max_seconds                                   # max seconds to sleep between BLAST calls
        self.min_exon_len = min_exon_length                             # minimum exon length (bp)
        self.timeout_db = timeout_db                                    # timeout for creating the database (seconds)
        self.evalue = evalue_threshold                                  # e-value threshold for BLAST calls
        self.batch_number = batch_number                                # batch number for BLAST calls
        self.threads = threads                                          # number of threads for BLAST calls
        self.cds_overlapping_threshold = cds_overlapping_threshold      # CDS overlapping threshold (0-1)
        self.self_hit_threshold = self_hit_threshold                    # self-hit threshold (0-1)
        self.results_db = results_db_name                               # results database name
        self.masking_perc_threshold = masking_perc_threshold
        self.stop_codons = ["TAG", "TGA", "TAA"]
        self.db_path = f'{self.specie_identifier}_genome_annotations.db'
        self.results_db = self.results_db or f'{self.specie_identifier}_results.db'
        self.UTR_features = ['five_prime_UTR', 'three_prime_UTR']
        self.gene_hierarchy_path = f"{self.specie_identifier}_gene_hierarchy.pkl"
        self.feat_of_interest = ['CDS', 'exon', 'intron'] + self.UTR_features
        self.__neither, self.__query, self.__target, self.__target_full, self.__target_insertion = 0, 0, 0, 0, 0
        self.__both = 0
        self.__annot_target_start, self.__annot_target_end, self.__target_t = None, None, None
        self.__query_CDS, self.__target_CDS = "-", "-"
        self.__found = False
        self.__tuples_full_length_duplications, self.__tuples_insertion_duplications, self.__tuples_truncation_events = [], [], []
        if self._HARD_FORCE:
            self.remove_file_if_exists(self.db_path)
            self.remove_file_if_exists(self.gene_hierarchy_path)
            self.remove_file_if_exists(self.results_db)
        elif self._SOFT_FORCE:
            self.remove_file_if_exists(self.results_db)
        print(self._DEBUG_MODE)

    def configure_logger(self):
        """
        configure_logger is a function that configures the logger.
        INFO level is used for the log file and WARNING and ERROR level for the console.
        """
        log_file_name = f"exonize_log_{datetime.datetime.now():%Y%m%d_%H%M%S}.log"
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        file_handler = logging.FileHandler(log_file_name, mode='w')
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(logging.Formatter('%(message)s'))
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.WARNING)
        console_handler.setFormatter(logging.Formatter('[%(levelname)s]: %(message)s'))
        self.logger.addHandler(file_handler)
        self.logger.addHandler(console_handler)

    @staticmethod
    def dump_pkl_file(out_filepath: str, obj: dict) -> None:
        """
        dump_pkl_file is a function that dumps an object into a pickle file.
        """
        with open(out_filepath, 'wb') as handle:
            pickle.dump(obj, handle)

    @staticmethod
    def dump_fasta_file(out_filepath: str, seq_dict: dict) -> None:
        """
        dump_fasta_file is a function that dumps a dictionary with sequences into a FASTA file.
        """
        with open(out_filepath, "w") as handle:
            for annot_id, annot_seq in seq_dict.items():
                record = SeqRecord(Seq(annot_seq), id=str(annot_id), description='')
                SeqIO.write(record, handle, "fasta")

    @staticmethod
    def read_pkl_file(filepath: str) -> dict:
        """
        read_pkl_file is a function that reads a pickle file and returns the object stored in it.
        """
        with open(filepath, 'rb') as handle:
            read_file = pickle.load(handle)
        return read_file

    @staticmethod
    def reverse_sequence_bool(gene_strand_: str):
        """
        reverse_sequence_bool checks if the gene is in the negative strand and returns True if it is.
        :param gene_strand_: strand
        """
        return gene_strand_ == '-'

    @staticmethod
    def get_shorter_longer_interv(a: P.Interval, b: P.Interval) -> tuple:
        """
        Given two intervals, the function
        get_shorter_longer_interv returns the smaller
        and the larger interval in length.
        """
        return (a, b) if a.upper - a.lower < b.upper - b.lower else (b, a)

    @staticmethod
    def get_overlap_percentage(a: P.Interval, b: P.Interval) -> float:
        """
        Given two intervals, the function
        get_overlap_percentage returns the percentage
        of the overlapping region relative to an interval b.
        """
        intersection = a & b
        return (intersection.upper - intersection.lower) / (b.upper - b.lower) if intersection else 0

    @staticmethod
    def sort_list_intervals_dict(list_dicts: list, reverse=False) -> list:
        """
        sort_list_intervals_dict is a function that sorts a list of dictionaries based on the coordinates of the
        intervals present in the dictionaries. The list is sorted in ascending order by default.
        """
        return sorted(list_dicts, key=lambda x: (x['coord'].lower, x['coord']), reverse=reverse)

    @staticmethod
    def remove_file_if_exists(filepath):
        """
        Removes a file if it exists.
        """
        if os.path.exists(filepath):
            os.remove(filepath)

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
                print("- Creating annotations database", end=" ")
                self.db = gffutils.create_db(self.in_file_path,
                                             dbfn=self.db_path,
                                             force=True,
                                             keep_order=True,
                                             merge_strategy='create_unique',
                                             sort_attribute_values=True,
                                             disable_infer_genes=True,
                                             disable_infer_transcripts=True)
                print("Done!")
            except ValueError as e:
                self.logger.exception(f"Incorrect genome annotations file {e}")
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
                self.logger.warning("The GFF file does not contain intron annotations")
                print(f"- Attempting to write intron annotations in database:", end=" ")
                try:
                    self.db.update(list(self.db.create_introns()), make_backup=False)
                except ValueError as e:
                    self.logger.exception(f"failed to write intron annotations in database. "
                                          f"Please provide a GFF3 file with intron annotations {e}")
                    sys.exit()
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
            self.logger.exception(f"Incorrect data base path {e}")
            sys.exit()

    def read_genome(self) -> None:
        """
        read_genome is a function that reads a FASTA file and stores the masked/unmasked genome sequence in a dictionary.
        The dictionary has the following structure: {chromosome: sequence}
        """
        try:
            tic_genome = time.time()
            print("- Reading genome file:", end=" ")
            parse_genome = SeqIO.parse(open(self.genome_path), 'fasta')
            if self.hard_masking:
                self.genome = {fasta.id: re.sub('[a-z]', 'N', str(fasta.seq)) for fasta in parse_genome}
            else:
                self.genome = {fasta.id: str(fasta.seq) for fasta in parse_genome}
            hms_time = dt.strftime(dt.utcfromtimestamp(time.time() - tic_genome), '%H:%M:%S')
            print(f"Done! [{hms_time}]")
        except (ValueError, FileNotFoundError) as e:
            self.logger.exception(f"Incorrect genome file path {e}")
            sys.exit()

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
                                                             frame=child.frame,  # One of '0', '1' or '2'.
                                                             type=child.featuretype,   # feature type name
                                                             attributes=dict(child.attributes)))  # feature type name
                    reverse = self.reverse_sequence_bool(gene.strand)
                    mrna_dict['mRNAs'][mrna_annot.id]['structure'] = self.sort_list_intervals_dict(temp_mrna_transcript, reverse)
                self.gene_hierarchy_dict[gene.id] = mrna_dict
        self.dump_pkl_file(self.gene_hierarchy_path, self.gene_hierarchy_dict)

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
            self.gene_hierarchy_dict = self.read_pkl_file(self.gene_hierarchy_path)
        else:
            self.create_gene_hierarchy_dict()
        self.read_genome()
        connect_create_results_db(self.results_db, self.timeout_db)
        if self._DEBUG_MODE:
            self.logger.warning("-All tblastx io files will be saved. This may take a large amount of disk space.")

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
        def tblastx_with_saved_io(ident: str, gene_id_: str, hit_seq_: str, query_seq_: str, query_coord_: P.Interval, gene_coord_: P.Interval):
            """
            tblastx_with_saved_io is a function that executes a tblastx search saving input and output files. This
            function is used for debugging purposes. The input and output files are saved in the following paths:
            - input: input/{ident}_query.fa and input/{gene_id_}_target.fa where ident is the identifier of the query
            sequence (CDS) and gene_id_ is the identifier of the target sequence (gene).
            - output: output/{ident}_output.xml where ident is the identifier of the query sequence (CDS).
            """

            output_file = f'output/{ident}_output.xml'
            if not os.path.exists(output_file):
                query_filename = f'input/{ident}_query.fa'
                target_filename = f'input/{gene_id_}_target.fa'
                if not os.path.exists(target_filename):
                    self.dump_fasta_file(target_filename, {f"{gene_id_}": hit_seq_})
                self.dump_fasta_file(query_filename, {ident: query_seq_})
                self.execute_tblastx(query_filename, target_filename, output_file)
            with open(output_file, "r") as result_handle:
                blast_records = NCBIXML.parse(result_handle)
                try:
                    temp_ = self.parse_tblastx_output(blast_records, query_coord_, gene_coord_)
                except Exception as e:
                    self.logger.exception(e)
                    sys.exit()
            return temp_

        def execute_tblastx_using_tempfiles(hit_seq_: str, query_seq_: str, query_coord_: P.Interval, gene_coord_: P.Interval):
            """
            execute_tblastx_using_tempfiles is a function that executes a tblastx search using temporary files.
            :param hit_seq_: target sequence (gene)
            :param query_seq_: query sequence (CDS)
            :param query_coord_: query coordinates (CDS) interval
            :param gene_coord_: target coordinates (gene) interval
            """
            with tempfile.TemporaryDirectory() as tmpdirname:
                query_filename = f'{tmpdirname}/query.fa'
                target_filename = f'{tmpdirname}/target.fa'
                self.dump_fasta_file(query_filename, {'query': query_seq_})
                self.dump_fasta_file(target_filename, {'target': hit_seq_})
                output_file = f'{tmpdirname}/output.xml'
                self.execute_tblastx(query_filename, target_filename, output_file)
                with open(output_file, 'r') as result_handle:
                    blast_records = NCBIXML.parse(result_handle)
                    try:
                        temp_ = self.parse_tblastx_output(blast_records, query_coord_, gene_coord_)
                    except Exception as e:
                        self.logger.exception(e)
                        sys.exit()
                return temp_

        chrom = self.gene_hierarchy_dict[gene_id]['chrom']
        gene_coord = self.gene_hierarchy_dict[gene_id]['coord']
        identifier = f'{gene_id}_{chrom}_{str(query_coord.lower)}_{query_coord.upper}'.replace(':', '_')
        if self._DEBUG_MODE:
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
        def get_hsp_dict(hsp) -> dict:
            return dict(
                score=hsp.score,
                bits=hsp.bits,
                evalue=hsp.expect,
                alignment_len=hsp.align_length * 3,
                hit_frame=hsp.frame,
                query_start=hsp.query_start - 1,
                query_end=hsp.query_end,
                target_start=hsp.sbjct_start - 1,
                target_end=hsp.sbjct_end,
                query_aln_prot_seq=hsp.query,
                target_aln_prot_seq=hsp.sbjct,
                query_num_stop_codons=hsp.query.count('*'),
                target_num_stop_codons=hsp.sbjct.count('*'),
                match=hsp.match)

        res_tblastx = {}
        # since we are performing a single query against a single subject, there's only one blast_record
        for blast_record in blast_records:
            if len(blast_record.alignments) == 0:
                continue
            alignment = blast_record.alignments[0]  # Assuming only one alignment per blast_record
            if len([aln for aln in blast_record.alignments]) > 1:
                print("More than one alignment per blast_record")
            for hsp_idx, hsp_rec in enumerate(alignment.hsps):
                blast_target_coord = P.open((hsp_rec.sbjct_start - 1) + hit_coord.lower, hsp_rec.sbjct_end + hit_coord.lower)
                if (self.get_overlap_percentage(q_coord, blast_target_coord) < self.self_hit_threshold
                        and self.get_overlap_percentage(blast_target_coord, q_coord) < self.self_hit_threshold):
                    res_tblastx[hsp_idx] = get_hsp_dict(hsp_rec)
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
        def get_intervals_overlapping_list(intv_list: list) -> list:
            """
            get_intervals_overlapping_list is a function that given a list of intervals, returns a list of tuples
            with the overlapping pairs described in case a (get_candidate_CDS_coords).
            :param intv_list: list of intervals
            :return: list of tuples with pairs of overlapping intervals
            """
            return [(feat_interv, intv_list[idx + 1])
                    for idx, feat_interv in enumerate(intv_list[:-1])
                    if all([self.get_overlap_percentage(feat_interv, intv_list[idx + 1]) >= self.cds_overlapping_threshold,
                            self.get_overlap_percentage(intv_list[idx + 1], feat_interv) >= self.cds_overlapping_threshold])]

        def resolve_overlaps_coords_list(coords_list):
            """
            resolve_overlaps_coords_list is a function that given a list of coordinates, resolves overlaps according
            to the criteria described above. Since the pairs described in case b (get_candidate_CDS_coords) are not
            considered to be overlapping, they are both included in the search.
            :param coords_list: list of coordinates
            :return: list of coordinates without overlaps
            """
            new_list = []
            overlaps_list = get_intervals_overlapping_list(coords_list)
            if overlaps_list:
                # We only process the first pair of overlapping intervals since the resolved overlap could also overlap
                # with the next interval in the list.
                intv_a, intv_b = overlaps_list[0]
                for idx, coord in enumerate(coords_list):
                    if coord != intv_a:
                        new_list.append(coord)
                    else:
                        shorter, _ = self.get_shorter_longer_interv(intv_a, intv_b)
                        new_list.append(shorter)
                        new_list.extend(coords_list[idx + 2:])
                        return resolve_overlaps_coords_list(new_list)
            else:
                return coords_list

        CDS_coords_list = list(set([i['coord'] for mrna_id, mrna_annot
                                    in self.gene_hierarchy_dict[gene_id]['mRNAs'].items()
                                    for i in mrna_annot['structure']
                                    if(i['type'] == 'CDS'
                                       and (i['coord'].upper - i['coord'].lower) >= self.min_exon_len)]))
        CDS_coords_list = sorted(CDS_coords_list, key=lambda x: (x.lower, x.upper))
        if CDS_coords_list:
            return resolve_overlaps_coords_list(CDS_coords_list)
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

            def sequence_masking_percentage(seq: str) -> float:
                """
                sequence_masking_percentage is a function that given a sequence, returns the percentage of hardmasking
                (N) in the sequence.
                """
                return seq.count('N') / len(seq)
            try:
                seq_ = Seq(self.genome[chrom_][coords_.lower:coords_.upper])
                if self.reverse_sequence_bool(gene_strand_):
                    seq_ = str(seq_.reverse_complement())
                masking_perc = round(sequence_masking_percentage(seq_), 2)
                if masking_perc > self.masking_perc_threshold:
                    seq_ = ''
                    if type_ == 'gene':
                        self.logger.info(f'Gene {gene_id_} in chromosome {chrom_} and coordinates {str(coords_.lower)}, {str(coords_.upper)} is hardmasked.')
                        insert_gene_ids_table(self.results_db, self.timeout_db, self.get_gene_tuple(gene_id_, 0))
                    if type_ == 'CDS':
                        self.logger.info(f'Gene {gene_id_} - {masking_perc * 100} of CDS {cds_coord} located in chromosome {chrom_} is hardmasked.')
                return seq_
            except KeyError as e_:
                self.logger.exception(f'Either there is missing a chromosome in the genome file '
                                      f'or the chromosome identifiers in the GFF3 and FASTA files do not match {e_}')
                sys.exit()

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
                    if tblastx_o:
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
        def get_fragment_tuple(gene_id_: str, cds_coord: P.Interval, blast_hits: dict, hsp_idx: int) -> tuple:
            def reformat_frame_strand(frame: int) -> tuple:
                """
                reformat_frame_strand is a function that converts the frame to a 0-based index and defines a strand variable
                based on the frame sign.
                """
                n_frame = abs(frame) - 1
                n_strand = '-' if frame < 0 else '+'
                return n_frame, n_strand

            hsp_dict = blast_hits[hsp_idx]
            hit_q_frame, hit_t_frame = hsp_dict['hit_frame']
            hit_q_f, hit_q_s = reformat_frame_strand(hit_q_frame)
            hit_t_f, hit_t_s = reformat_frame_strand(hit_t_frame)
            return (gene_id_, cds_coord.lower, cds_coord.upper,
                    hit_q_f, hit_q_s, hit_t_f, hit_t_s,
                    hsp_dict['score'], hsp_dict['bits'], hsp_dict['evalue'],
                    hsp_dict['alignment_len'], hsp_dict['query_start'], hsp_dict['query_end'],
                    hsp_dict['target_start'], hsp_dict['target_end'], hsp_dict['query_aln_prot_seq'],
                    hsp_dict['target_aln_prot_seq'], hsp_dict['match'],
                    hsp_dict['query_num_stop_codons'], hsp_dict['target_num_stop_codons'])

        tuple_list = [get_fragment_tuple(gene_id, cds_coord, blast_hits, hsp_idx)
                      for cds_coord, blast_hits in blast_cds_dict.items() for hsp_idx, hsp_dict in blast_hits.items()]
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
                and all([self.get_overlap_percentage(i['coord'], cds_intv) >= self.cds_overlapping_threshold,
                         self.get_overlap_percentage(cds_intv, i['coord']) >= self.cds_overlapping_threshold])]

    def identify_full_length_duplications(self) -> None:
        """
        identify_full_length_duplications is a function that identifies full-length duplications following our classification
        model. The function iterates over all representative tblastx hits and for each transcript associated with the gene harboring
        the event it identifies the following events:
        - I. Full exon duplication: the match fully (following our coverage criteria) overlaps with a CDS annotation.
        - II. Insertion: the match is found within a larger CDS.
        - III. Deactivation or unnanotated: the match is found in an intron or UTR.
        - IV. Trunctation: the match spans more than one annotation (e.g., CDS, intron, UTR).
        The identify_full_length_duplications function also looks for:
         - I. Obligate events: defined as events where the query CDS and the target CDS are included within the same
          transcript.
        - II. Neither events: defined as events where the query CDS is not found in the transcript and the target
        figures a trunctation or deactivation.
        The identified events are stored in the results database (self.results_db) in the tables:
         - ObligatoryEvents: identifies both, query and target CDSs. A single record is stored per event.
         - TruncationEvents: identifies all covered annotations. The number of records per event will correspond to the
         number of annotations covered by the event.
         - FullLengthDuplications tables: identifies full-length duplications. The number of records per event will correspond
         to the number of transcripts associated with the gene harboring the event.
        """
        def initialize_variables() -> None:
            """
            initializes variables used in the identify_full_length_duplications function
            """
            self.__neither, self.__query, self.__target, self.__target_full, self.__target_insertion = 0, 0, 0, 0, 0
            self.__both = 0
            self.__annot_target_start, self.__annot_target_end, self.__target_t = None, None, None
            self.__query_CDS, self.__target_CDS, self.__query_CDS_frame = "-", "-", "-"
            self.__found = False

        def initialize_list_of_tuples() -> None:
            """
            initializes the list of tuples used to store the identified events in the identify_full_length_duplications function
            """
            self.__tuples_full_length_duplications, self.__tuples_obligatory_events, self.__tuples_truncation_events = [], [], []

        def identify_query(trans_dict_: dict, cds_intv_: P.Interval) -> None:
            """
            identify_query is a function that identifies the tblastx query CDS in the gene transcript. A transcript
            cannot have overlapping CDSs. If the query CDS overlaps with more than one CDS, the program exits.
            """
            query_only_ = self.find_overlapping_annot(trans_dict_, cds_intv_)
            if len(query_only_) > 1:
                self.logger.error(f'Overlapping query CDSs: {query_only_}, please review your GFF3 file')
                sys.exit()
            elif query_only_:
                self.__query_CDS, _, self.__query_CDS_frame = query_only_[0]
                self.__query = 1
                self.__target_t = 'QUERY_ONLY'

        def target_out_of_mRNA(trans_coord_: P.Interval, mrna_: str, row_: list) -> bool:
            """
            target_out_of_mRNA is a function that identifies tblastx hits that are outside the mRNA transcript.
            """
            fragment_id_, gene_id_, gene_s_, _, cds_s_, cds_e_, query_s_, query_e_, target_s_, target_e_, evalue_ = row_
            target_intv_ = P.open(target_s_ + gene_s_, target_e_ + gene_s_)
            if target_intv_.upper < trans_coord_.lower or trans_coord_.upper < target_intv_.lower:
                if (self.__query + self.__target) == 0:
                    self.__neither = 1
                    self.__tuples_full_length_duplications.append((fragment_id_, gene_id_, mrna_,
                                                                   cds_s_, cds_e_, self.__query_CDS, query_s_, query_e_, "OUT_OF_MRNA",
                                                                   self.__target_CDS, self.__annot_target_start, self.__annot_target_end,
                                                                   target_s_, target_e_,
                                                                   self.__neither, self.__query, self.__target, self.__both, evalue_))
                    return True
            return False

        def indetify_full_target(trans_dict_: dict, target_intv_: P.Interval) -> None:
            """
            indetify_full_target is a function that identifies tblastx hits that are full-length duplications as described
            in self.find_overlapping_annot
            """
            target_only_ = self.find_overlapping_annot(trans_dict_, target_intv_)
            if len(target_only_) > 1:
                self.logger.error(f'overlapping target CDSs: {target_only_}')
                sys.exit()
            elif target_only_:
                self.__target_CDS, t_CDS_coord_, self.__target_CDS_frame = target_only_[0]
                self.__target_full = 1
                self.__found = True
                self.__target_t = "FULL"
                self.__annot_target_start, self.__annot_target_end = t_CDS_coord_.lower, t_CDS_coord_.upper

        def indentify_insertion_target(trans_dict_: dict, target_intv_: P.Interval) -> None:
            """
            Identifies tblastx hits that are insertion duplications these can be:
            - INS_CDS: if the insertion is in a CDS
            - INS_UTR: if the insertion is in an untralated region
            - DEACTIVATED: if the insertion is in an intron
            """
            def filter_structure_by_interval_and_type(structure: dict, t_intv_, annot_type: str) -> list:
                return [(i['id'], i['coord']) for i in structure if (i['coord'].contains(t_intv_) and annot_type in i['type'])]

            insertion_CDS_ = filter_structure_by_interval_and_type(trans_dict_['structure'], target_intv_, 'CDS')
            if insertion_CDS_:
                self.__target_CDS, t_CDS_coord_ = insertion_CDS_[0]
                self.__target_insertion = 1
                self.__target_t = "INS_CDS"
                self.__found = True
                self.__annot_target_start, self.__annot_target_end = t_CDS_coord_.lower, t_CDS_coord_.upper
            else:
                insertion_UTR_ = filter_structure_by_interval_and_type(trans_dict_['structure'], target_intv_, 'UTR')
                if insertion_UTR_:
                    self.__target_CDS, t_CDS_coord_ = insertion_UTR_[0]
                    self.__target_t = "INS_UTR"
                    self.__found = True
                    self.__annot_target_start, self.__annot_target_end = t_CDS_coord_.lower, t_CDS_coord_.upper
                else:
                    insertion_intron_ = filter_structure_by_interval_and_type(trans_dict_['structure'], target_intv_, 'intron')
                    if insertion_intron_:
                        self.__target_CDS, t_CDS_coord_ = insertion_intron_[0]
                        self.__target_t = "DEACTIVATED"
                        self.__found = True
                        self.__annot_target_start, self.__annot_target_end = t_CDS_coord_.lower, t_CDS_coord_.upper

        def indentify_truncation_target(trans_dict_: dict, mrna_: str, row_: list) -> None:
            """
            Identifies tblastx hits that are truncation duplications. These are hits that span across more than one annotation
            in the transcript architecture. We record a line per annotation that is truncated.
            """
            fragment_id_, gene_id_, gene_s_, _, cds_s_, cds_e_, query_s_, query_e_, target_s_, target_e_, evalue_ = row_
            target_intv_ = P.open(target_s_ + gene_s_, target_e_ + gene_s_)
            intv_dict_ = get_interval_dictionary(trans_dict_['structure'], target_intv_, trans_dict_['coord'])
            if intv_dict_:
                self.__target_t = "TRUNC"
                self.__target_CDS, self.__annot_target_start, self.__annot_target_end = None, None, None
                for seg_b, value in intv_dict_.items():
                    coord_b = value['coord']
                    self.__tuples_truncation_events.append(
                        (fragment_id_, gene_id_, mrna_, trans_dict_['coord'].lower, trans_dict_['coord'].upper,
                         self.__query_CDS, cds_s_, cds_e_, query_s_, query_e_,
                         target_s_, target_e_, value['id'], value['type'],
                         coord_b.lower, coord_b.upper, seg_b.lower, seg_b.upper))

        def identify_obligate_pair(trans_coord_: P.Interval, mrna_: str, row_: list) -> None:
            """
            Identifies tblastx hits that are obligate pairs. These are hits where the query and target show as CDSs in the
            transcript in question.
            """
            fragment_id_, gene_id_, _, _, cds_s_, cds_e_, query_s_, query_e_, target_s_, target_e_, evalue_ = row_
            self.__both = 1
            self.__query, self.__target = 0, 0
            self.__tuples_obligatory_events.append((fragment_id_, gene_id_, mrna_, trans_coord_.lower, trans_coord_.upper,
                                                    cds_s_, cds_e_, self.__query_CDS, self.__query_CDS_frame,
                                                    query_s_, query_e_, self.__target_CDS, self.__target_CDS_frame,
                                                    self.__annot_target_start, self.__annot_target_end,
                                                    target_s_, target_e_, self.__target_t))

        def identify_neither_pair() -> None:
            """
            Identifies tblastx hits that are neither pairs. These are hits where the query and target do not show as CDSs
            in the transcript in question.
            """
            self.__neither = 1
            self.__target_t = 'NEITHER'

        def insert_full_length_duplication_tuple(mrna_: str, row_: list):
            """
            insert_full_length_duplication_tuple is a function that appends the "row" event to the list of tuples.
            """
            fragment_id_, gene_id_, _, _, cds_s_, cds_e_, query_s_, query_e_, target_s_, target_e_, evalue_ = row_
            self.__tuples_full_length_duplications.append((fragment_id_, gene_id_, mrna_,
                                                           cds_s_, cds_e_, self.__query_CDS, query_s_, query_e_,
                                                           self.__target_t, self.__target_CDS, self.__annot_target_start,
                                                           self.__annot_target_end,
                                                           target_s_, target_e_, self.__neither, self.__query, self.__target,
                                                           self.__both, evalue_))

        def insert_tuples_in_results_db() -> None:
            """
            insert_tuples_in_results_db is a function that inserts the list of tuples collected by the identify_full_length_duplications
            function into the results database.
            """
            instert_full_length_event(self.results_db, self.timeout_db, self.__tuples_full_length_duplications)
            instert_obligatory_event(self.results_db, self.timeout_db, self.__tuples_obligatory_events)
            instert_truncation_event(self.results_db, self.timeout_db, self.__tuples_truncation_events)

        rows = query_filtered_full_duplication_events(self.results_db, self.timeout_db)
        initialize_list_of_tuples()
        for row in rows:
            _, gene_id, gene_s, _, cds_s, cds_e, _, _, target_s, target_e, _ = row
            cds_intv, target_intv = P.open(cds_s, cds_e), P.open(target_s + gene_s, target_e + gene_s)
            for mrna, trans_dict in self.gene_hierarchy_dict[gene_id]['mRNAs'].items():
                initialize_variables()  # we want to classify each transcript independently
                trans_coord = trans_dict['coord']
                # ####### QUERY ONLY - FULL LENGTH #######
                identify_query(trans_dict, cds_intv)
                # ###### CHECK: TARGET REGION NOT IN mRNA #######
                if target_out_of_mRNA(trans_coord, mrna, row):
                    continue
                # ####### TARGET ONLY - FULL LENGTH #######
                indetify_full_target(trans_dict, target_intv)
                if self.__target_full == 0:
                    # ####### INSERTION #######
                    indentify_insertion_target(trans_dict, target_intv)
                    if not self.__found:
                        # ####### TRUNCATION #######
                        indentify_truncation_target(trans_dict, mrna, row)
                self.__target = self.__target_full + self.__target_insertion
                # ####### OBLIGATE PAIR #######
                if self.__query + self.__target == 2:
                    identify_obligate_pair(trans_dict['coord'], mrna, row)
                # ####### NEITHER PAIR #######
                elif self.__query + self.__target == 0:
                    identify_neither_pair()
                insert_full_length_duplication_tuple(mrna, row)
        insert_tuples_in_results_db()

    def assign_pair_ids(self, full_matches: list) -> list:
        def get_shorter_intv_overlapping_percentage(a: P.Interval, b: P.Interval) -> float:
            """
            get_shorter_intv_overlapping_percentage is a function that given two intervals, returns the percentage of
            overlap of the shorter interval with the longer interval.
            """
            shorter, longer = self.get_shorter_longer_interv(a, b)
            return self.get_overlap_percentage(longer, shorter)

        def find_pairs_overlapping_perc(pairs: list) -> list:
            """
            find_pairs_overlapping_perc is a function that given a list of pairs of intervals, returns a list of
            percentages of overlap of the shorter interval with the longer interval.
            """
            return [get_shorter_intv_overlapping_percentage(pair[0], pair[1]) for pair in pairs]

        fragments, skip_frag, pair_id_counter = [], [], 1
        with tqdm(total=len(full_matches), position=0, leave=True, ncols=50) as progress_bar:
            for frag_a in full_matches:
                frag_id_a, gene_id_a, q_s_a, q_e_a, t_s_a, t_e_a, event_type_a = frag_a
                t_intv_a = P.open(t_s_a, t_e_a)
                q_intv_a = P.open(q_s_a, q_e_a)
                if frag_id_a not in skip_frag:
                    candidates = [i for i in full_matches if i[1] == gene_id_a and i[0] not in [*skip_frag, frag_id_a]]
                    if candidates:
                        temp_cand = []
                        for frag_b in candidates:
                            frag_id_b, gene_id_b, q_s_b, q_e_b, t_s_b, t_e_b, event_type_b = frag_b
                            t_intv_b = P.open(t_s_b, t_e_b)
                            q_intv_b = P.open(q_s_b, q_e_b)
                            overlapping_pairs = find_pairs_overlapping_perc([(q_intv_a, q_intv_b), (t_intv_a, t_intv_b)])
                            reciprocal_pairs = find_pairs_overlapping_perc([(t_intv_a, q_intv_b), (t_intv_b, q_intv_a)])
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
                            fragments.extend([(pair_id_counter, frag) for frag in [frag_id_a, *temp_cand]])
                            pair_id_counter += 1
                progress_bar.update(1)
        return fragments

    def get_identity_and_dna_seq_tuples(self) -> list:
        """
        Retrieves DNA sequences from tblastx query and target, computes their DNA and amino acid identity,
        and returns a list of tuples with the structure:
        (DNA_identity, AA_identity, query_dna_seq, target_dna_seq, fragment_id)
        """

        def fetch_dna_sequence(chrom: str, annot_start: int, annot_end: int, pos_annot_s: int, pos_annot_e: int, strand: str):
            """
            Retrieve a subsequence from a genomic region, reverse complementing it if on the negative strand.
            """
            seq = Seq(self.genome[chrom][annot_start:annot_end][pos_annot_s:pos_annot_e])
            return str(seq.reverse_complement()) if self.reverse_sequence_bool(strand) else str(seq)

        def compute_identity(seq_1: str, seq_2: str) -> float:
            """
            Compute the identity between two sequences (seq_1 and seq_2) using the Hamming distance method.
            """
            # Calculate the Hamming distance and return it
            return round(sum(i == j for i, j in zip(seq_1, seq_2)) / len(seq_2), 3)

        def process_fragment(fragment: list) -> tuple:
            """
            process fragment recovers the query/target DNA sequences (since the tblastx alignment is done in amino acids)
            and computes the DNA and amino acid identity. This information is returned in a tuple and later used to update the
            Fragments table.
            """
            (fragment_id, gene_id, gene_start, gene_end, gene_chrom,
             CDS_start, CDS_end, query_start, query_end, target_start, target_end,
             query_strand, target_strand, query_aln_prot_seq, target_aln_prot_seq) = fragment
            query_dna_seq = fetch_dna_sequence(gene_chrom, CDS_start, CDS_end, query_start, query_end, query_strand)
            target_dna_seq = fetch_dna_sequence(gene_chrom, gene_start, gene_end, target_start, target_end,
                                                target_strand)

            if len(query_dna_seq) != len(target_dna_seq):
                raise ValueError(f'{gene_id}: CDS {(CDS_start, CDS_end)} search - sequences must have the same length.')

            return (compute_identity(query_dna_seq, target_dna_seq),
                    compute_identity(query_aln_prot_seq, target_aln_prot_seq),
                    query_dna_seq, target_dna_seq, fragment_id)

        return [process_fragment(fragment) for fragment in query_fragments(self.results_db, self.timeout_db)]

    def run_exonize_pipeline(self) -> None:
        """
        run_exonize_pipeline iterates over all genes in the gene_hierarchy_dict attribute and performs a tblastx search
        for each representative CDS (see get_candidate_CDS_coords).
        The steps are the following:
        - 1. The function checks if the gene has already been processed. If so, the gene is skipped.
        - 2. For each gene, the function iterates over all representative CDSs and performs a tblastx search.
        - 3. The percent_query (hit query coverage) column is inserted in the Fragments table.
        - 4. The Filtered_full_length_events View is created. This view contains all tblastx hits that have passed the
        filtering step.
        - 5. The Mrna_counts View is created. This view contains the number of transcripts associated with each gene.
        - 6. The function creates the Cumulative_counts table. This table contains the cumulative counts of the different
        event types across transcripts.
        - 7. The function collects all the raw concatenated event types (see query_concat_categ_pairs).
        - 8. The function generates the unique event types (see generate_unique_events_list) so that no event type is repeated.
        - 9. The function inserts the unique event types in the Event_categ_full_length_events_cumulative_counts table.
        - 10. The function collects the identity and DNA sequence tuples (see get_identity_and_dna_seq_tuples) and inserts
        them in the Fragments table.
        - 11. The function collects all events in the Full_length_events_cumulative_counts table.
        - 12. The function reconciles the events by assigning a "pair ID" to each event (see assign_pair_ids).
        - 13. The function creates the Exclusive_pairs view. This view contains all the events that follow the mutually exclusive
        category.
        """
        def exonize_asci_art() -> None:
            exonize_ansi_regular = """

            ███████╗██╗  ██╗ ██████╗ ███╗   ██╗██╗███████╗███████╗
            ██╔════╝╚██╗██╔╝██╔═══██╗████╗  ██║██║╚══███╔╝██╔════╝
            █████╗   ╚███╔╝ ██║   ██║██╔██╗ ██║██║  ███╔╝ █████╗
            ██╔══╝   ██╔██╗ ██║   ██║██║╚██╗██║██║ ███╔╝  ██╔══╝
            ███████╗██╔╝ ██╗╚██████╔╝██║ ╚████║██║███████╗███████╗
            ╚══════╝╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚═╝╚══════╝╚══════╝
            """
            print(exonize_ansi_regular)

        def batch(iterable: list, n=1) -> list:
            """
            batch is a function that given a list and a batch size, returns an iterable of lists of size n.
            """
            it_length = len(iterable)
            for indx in range(0, it_length, n):
                yield iterable[indx:min(indx + n, it_length)]

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
            with tqdm(total=len(args_list), position=0, leave=True, ncols=50) as progress_bar:
                for arg_batch in batches_list:
                    t = ThreadPool(processes=self.threads)
                    t.map(self.find_coding_exon_duplicates, arg_batch)
                    t.close()
                    t.join()
                    progress_bar.update(len(arg_batch))
            hms_time = dt.strftime(dt.utcfromtimestamp(time.time() - tic), '%H:%M:%S')
            print(f' Done! [{hms_time}]')
        else:
            print('All genes have been processed. If you want to re-run the analysis, delete/rename the results DB.')
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
        identity_and_sequence_tuples = self.get_identity_and_dna_seq_tuples()
        insert_identity_and_dna_algns_columns(self.results_db, self.timeout_db, identity_and_sequence_tuples)
        full_matches = query_full_events(self.results_db, self.timeout_db)
        print('- Reconciling events')
        fragments = self.assign_pair_ids(full_matches)
        instert_pair_id_column_to_full_length_events_cumulative_counts(self.results_db, self.timeout_db, fragments)
        create_exclusive_pairs_view(self.results_db, self.timeout_db)
        hms_time = dt.strftime(dt.utcfromtimestamp(time.time() - tic), '%H:%M:%S')
        print(f' Done! [{hms_time}]')
