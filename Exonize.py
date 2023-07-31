from sqlite_utils import *
import gffutils                                        # for creating/loading DBs
import subprocess                                      # for calling gffread
import portion as P                                    # for working with intervals
import os                                              # for working with files
import time                                            # for sleeping between BLAST calls and for timeout on DB creation
import random                                          # for random sleep
import re                                              # regular expressions for genome masking
# import tempfile                                        # for creating temporary files
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
                 verbose=True,
                 hard_masking=False,
                 evalue_threshold=1e-1,
                 coverage_threshold=0.95,
                 min_align_len_perc=0.3,
                 sleep_max_seconds=5,
                 min_exon_length=20,
                 batch_number=100,
                 threads=7,
                 timeout_db=160.0):

        self.genome = None                             # genome sequence
        self.gene_hierarchy_dict = None                # gene hierarchy dictionary (gene -> transcript -> exon)
        self.db_features = None                        # features in the database
        self.old_filename = None                       # old filename (if GTF file is provided)
        self.db = None                                 # database object (gffutils)
        self.specie_identifier = specie_identifier     # specie identifier
        self.verbose = verbose                         # verbose mode (True/False)
        self.in_file_path = gff_file_path              # input file path (GFF/GTF)
        self.genome_path = genome_path                 # genome path (FASTA)
        self.hard_masking = hard_masking               # hard masking (True/False)
        self.secs = sleep_max_seconds                  # max seconds to sleep between BLAST calls
        self.min_exon_len = min_exon_length            # minimum exon length (bp)
        self.timeout_db = timeout_db                   # timeout for creating the database (seconds)
        self.evalue = evalue_threshold                 # e-value threshold for BLAST
        self.min_align_len_perc = min_align_len_perc   # minimum alignment length percentage of query length
        self.batch_number = batch_number               # batch number for BLAST calls
        self.threads = threads                         # number of threads for BLAST calls
        self.id_spec_attribute = spec_attribute
        self.coverage_threshold = coverage_threshold
        self.stop_codons = ["TAG", "TGA", "TAA"]
        self.db_path = f'{self.specie_identifier}_genome_annotations.db'
        self.results_db = results_db_name
        if self.results_db == '':
            self.results_db = f'{self.specie_identifier}_results.db'
        self.UTR_features = ['five_prime_UTR', 'three_prime_UTR']
        self.gene_hierarchy_path = f"{self.specie_identifier}_gene_hierarchy.pkl"
        self.feat_of_interest = ['CDS', 'exon', 'intron'] + self.UTR_features

    def create_parse_or_update_database(self) -> None:
        if not os.path.exists(self.db_path):
            if 'gtf' in self.in_file_path:
                self.old_filename = self.in_file_path
                self.in_file_path = f"{self.old_filename.rsplit('.gtf')[0]}.gff"
                self.convert_gtf_to_gff()
                print('the GTF file has been converted into a GFF3 file')
                print(f'with filename: {self.in_file_path}')
            self.create_database()
        if not self.db:
            print("- Reading annotations database:", end=" ")
            self.load_db()
            print("Done!")
        self.db_features = list(self.db.featuretypes())
        self.create_intron_annotations()

    def convert_gtf_to_gff(self) -> None:
        gffread_command = ["gffread", self.old_filename, "-O", "-o", self.in_file_path]
        subprocess.call(gffread_command)

    def create_database(self) -> None:
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

    def load_db(self) -> None:
        try:
            self.db = gffutils.FeatureDB(self.db_path, keep_order=True)
        except ValueError as e:
            print(' ')
            print('---------------------ERROR-------------------------------')
            print(f"Incorrect data base path  {e}")
            print('---------------------------------------------------------')

    def create_intron_annotations(self) -> None:
        """
        run only once otherwise duplicated
        annotations will be created
        """
        if 'intron' not in self.db_features:
            print('---------------------Warning------------------------------')
            print("-The genomic annotations do not contain intron annotations")
            print('----------------------------------------------------------')
            if self.verbose:
                print(f"- Attempting to write intron annotations in database:", end=" ")
            try:
                def intron_id(f) -> str:
                    return ','.join(f[self.id_spec_attribute])

                introns_list = list(self.db.create_introns())
                self.db.update(introns_list, id_spec={'intron': [intron_id]}, make_backup=False)
            except ValueError as e:
                print("failed to write intron annotations in database, "
                      "please try with a different spec attribute or provide a GFF3 file with intron annotations")
                print(e)
            if self.verbose:
                print("Done!")

    def read_genome(self) -> None:
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

    def prepare_data(self) -> None:
        self.create_parse_or_update_database()
        if os.path.exists(self.gene_hierarchy_path):
            self.gene_hierarchy_dict = read_pkl_file(self.gene_hierarchy_path)
        else:
            self.create_gene_hierarchy_dict()
        self.read_genome()
        connect_create_results_db(self.results_db, self.timeout_db)

    def create_gene_hierarchy_dict(self) -> None:
        """
        GFF coordinates are 1-based, so we need to subtract 1 from the start position to convert to 0-based.
        :return:
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
                            temp_mrna_transcript.append({'id': child.id,  # ID attribute
                                                         'coord': coord,  # ID coordinate starting at 0
                                                         # One of '0', '1' or '2'. '0' indicates that the x base of
                                                         # the feature is the x base of a codon
                                                         'frame': child.frame,
                                                         'type': child.featuretype}  # feature type name
                                                        )
                    mrna_dict['mRNAs'][mrna_annot.id]['structure'] = sort_list_intervals_dict(temp_mrna_transcript)
                self.gene_hierarchy_dict[gene.id] = mrna_dict
        dump_pkl_file(self.gene_hierarchy_path, self.gene_hierarchy_dict)

    def run_tblastx(self, gene_id: str, query_seq: str, hit_seq: str, query_coord) -> dict:
        chrom = self.gene_hierarchy_dict[gene_id]['chrom']
        gene_coord = self.gene_hierarchy_dict[gene_id]['coord']
        identifier = f'{gene_id}_{chrom}_{str(query_coord.lower)}_{query_coord.upper}'.replace(':', '_')
        output_file = f'output/{identifier}_output.xml'
        if not os.path.exists(output_file):
            query_filename = f'input/{identifier}_query.fa'
            target_filename = f'input/{gene_id}_target.fa'
            if not os.path.exists(target_filename):
                dump_fasta_file(target_filename, {f"{gene_id}": hit_seq})
            dump_fasta_file(query_filename, {identifier: query_seq})
            tblastx_command = ['tblastx',
                               '-query', query_filename,
                               '-subject', target_filename,
                               '-outfmt', '5',  # XML output format
                               '-out', output_file]
            subprocess.run(tblastx_command)
        with open(output_file, "r") as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            try:
                temp = self.parse_tblastx_output(blast_records, query_seq, hit_seq, query_coord, gene_coord)
            except Exception as e:
                print(e)
                sys.exit()
        return temp

    def parse_tblastx_output(self, blast_records, query_seq: str, hit_seq: str, q_coord, hit_coord) -> dict:
        """
        We only want to consider hits that:
            (i)   have an e-value lower than the threshold,
            (ii)  have a minimum alignment length percentage of the query sequence and
            (iii) that do not overlap with the query sequence (self-hit),
                  with a maximum overlap of 50% of the query sequence.
        :return: dict with the following structure: {target_id {hsp_id:{'score': '', 'bits': '','evalue': '',...}}}
        """
        res_tblastx = {}
        # since we are performing a single query against a single subject, there's only one blast_record
        for blast_record in blast_records:
            if len(blast_record.alignments) == 0:
                # query_len = len(query_seq)
                # print(f"No alignments found {query_id} --- query length: {query_len}")
                continue
            alignment = blast_record.alignments[0]  # Assuming only one alignment per blast_record
            if len([aln for aln in blast_record.alignments]) > 1:
                print("More than one alignment per blast_record")
            for hsp_idx, hsp in enumerate(alignment.hsps):
                blast_target_coord = P.open((hsp.sbjct_start - 1) + hit_coord.lower, hsp.sbjct_end + hit_coord.lower)
                query_aligned_frac = round((hsp.align_length * 3) / blast_record.query_length, 3)
                target_query_overlap_percentage = get_overlap_percentage(blast_target_coord, q_coord)
                query_target_overlap_percentage = get_overlap_percentage(q_coord, blast_target_coord)
                if (hsp.expect < self.evalue
                        and query_aligned_frac > self.min_align_len_perc
                        and query_target_overlap_percentage < 0.5
                        and target_query_overlap_percentage < 0.5):
                    res_tblastx[hsp_idx] = get_hsp_dict(hsp, query_seq, hit_seq)
        return res_tblastx

    def get_gene_tuple(self, gene_id: str, bin_has_dup: int) -> tuple:
        gene_coord = self.gene_hierarchy_dict[gene_id]['coord']
        return (gene_id,
                self.gene_hierarchy_dict[gene_id]['chrom'],
                (self.gene_hierarchy_dict[gene_id]['strand']),
                gene_coord.lower,
                gene_coord.upper,
                bin_has_dup)

    def find_coding_exon_duplicates(self, gene_id: str) -> None:
        """
        The reading frame of the CDS queries is respected in the tblastx search.
        """
        time.sleep(random.randrange(0, self.secs))
        CDS_blast_dict = {}
        chrom = self.gene_hierarchy_dict[gene_id]['chrom']
        gene_coord = self.gene_hierarchy_dict[gene_id]['coord']
        CDS_coords_list = list(set([i['coord']
                                    for mrna_id, mrna_annot in self.gene_hierarchy_dict[gene_id]['mRNAs'].items()
                                    for i in mrna_annot['structure'] if i['type'] == 'CDS']))
        if CDS_coords_list:
            for cds_coord in CDS_coords_list:
                if (cds_coord.upper - cds_coord.lower) >= self.min_exon_len:
                    try:
                        cds_seq = self.genome[chrom][cds_coord.lower:cds_coord.upper]
                        gene_seq = self.genome[chrom][gene_coord.lower:gene_coord.upper]
                    except KeyError as e:
                        print(f'Either there is missing a chromosome in the genome file '
                              f'or the chromosome identifiers in the GFF3 and FASTA files do not match {e}')
                        sys.exit()
                    temp = self.run_tblastx(gene_id, cds_seq, gene_seq, cds_coord)
                    if temp:
                        CDS_blast_dict[cds_coord] = temp
        if CDS_blast_dict:
            self.insert_fragments_table(gene_id, CDS_blast_dict)
        else:
            insert_gene_ids_table(self.results_db, self.timeout_db, self.get_gene_tuple(gene_id, 0))

    def get_fragments_matches_tuples(self, gene_id: str, match_mrna_id: str, event: list, target_coord) -> list:
        cds_mrna_id, cds_id, cds_start, cds_end, query_start, query_end, target_start, target_end, _ = event
        trans_structure = self.gene_hierarchy_dict[gene_id]['mRNAs'][match_mrna_id]
        cds_mrna_coords = self.gene_hierarchy_dict[gene_id]['mRNAs'][cds_mrna_id]['coord']
        return [(gene_id, cds_mrna_id, cds_mrna_coords.lower, cds_mrna_coords.upper,
                 cds_id, cds_start, cds_end, query_start, query_end,
                 target_start, target_end, match_mrna_id,
                 trans_structure['coord'].lower, trans_structure['coord'].upper,
                 annot['id'], annot['type'], annot['coord'].lower, annot['coord'].upper,
                 round(get_overlap_percentage(annot['coord'], target_coord), 3),
                 round(get_overlap_percentage(target_coord, annot['coord']), 3))
                for annot in trans_structure['structure'] if target_coord.overlaps(annot['coord'])]

    def insert_fragments_table(self, gene_id: str, blast_cds_dict: dict) -> None:
        tuple_list = [get_fragment_tuple(gene_id, cds_coord, blast_hits, hsp_idx)
                      for cds_coord, blast_hits in blast_cds_dict.items()
                      for hsp_idx, hsp_dict in blast_hits.items()]

        insert_fragments_table_param, insert_gene_table_param = insert_fragments_calls()
        db = sqlite3.connect(self.results_db, timeout=self.timeout_db)
        cursor = db.cursor()
        cursor.execute(insert_gene_table_param, self.get_gene_tuple(gene_id, 1))
        cursor.executemany(insert_fragments_table_param, tuple_list)
        db.commit()
        db.close()

    def find_full_length_duplications(self):
        rows = query_full_duplication_events(self.results_db, self.timeout_db, self.coverage_threshold)
        tuples_full_length_duplications, tuples_obligatory_events, tuples_truncation_events = [], [], []
        for row in rows:
            fragment_id, gene_id, gene_s, gene_e, cds_s, cds_e, query_s, query_e, target_s, target_e, evalue = row
            cds_intv = P.open(cds_s, cds_e)
            target_intv = P.open(target_s + gene_s, target_e + gene_s)
            for mrna, trans_dict in self.gene_hierarchy_dict[gene_id]['mRNAs'].items():
                found = False
                trans_coord = trans_dict['coord']
                neither, query, target, target_full, target_insertion, both = 0, 0, 0, 0, 0, 0
                # ####### QUERY ONLY - FULL LENGTH #######
                query_only = [i['id'] for i in trans_dict['structure'] if i['type'] == 'CDS' and i['coord'] == cds_intv]
                if len(query_only) > 1:
                    print(f'overlapping query CDSs: {query_only}')
                    continue
                elif query_only:
                    query_CDS = query_only[0]
                    query = 1
                # ###### CHECK: TARGET REGION NOT IN mRNA #######
                if target_intv.lower < trans_coord.lower or trans_coord.upper < target_intv.upper:
                    if (query + target) == 0:
                        neither = 1
                    tuples_full_length_duplications.append((fragment_id, gene_id, mrna,
                                                            cds_s, cds_e, query_s, query_e,
                                                            target_s, target_e,
                                                            neither, query, target, both, evalue))

                    continue
                # ####### TARGET ONLY - FULL LENGTH #######
                target_only = [(i['id'], i['coord']) for i in trans_dict['structure']
                               if (get_overlap_percentage(target_intv, i['coord']) >= self.coverage_threshold
                                   and get_overlap_percentage(i['coord'], target_intv) >= self.coverage_threshold
                                   and i['type'] == "CDS")]
                if len(target_only) > 1:
                    print(f'overlapping query CDSs: {target_only}')
                    continue
                elif target_only:
                    target_full = 1
                    found = True
                    target_t = "FULL"
                    target_CDS, t_CDS_coord = target_only[0]
                # ####### TARGET ONLY - INSERTION #######
                if target_full == 0:
                    insertion_CDS = [(i['id'], i['coord']) for i in trans_dict['structure']
                                     if (i['coord'].contains(target_intv)
                                         and get_overlap_percentage(target_intv, i['coord']) < self.coverage_threshold
                                         and i['type'] == "CDS")]
                    insertion_UTR = [(i['id'], i['coord']) for i in trans_dict['structure']
                                     if (i['coord'].contains(target_intv)
                                         and get_overlap_percentage(target_intv, i['coord']) < self.coverage_threshold
                                         and "UTR" in i['type'])]
                    if insertion_CDS:
                        target_insertion = 1
                        target_t = "INS_CDS"
                        found = True
                        target_CDS, t_CDS_coord = insertion_CDS[0]
                    elif insertion_UTR:
                        target_insertion = 1
                        target_t = "INS_UTR"
                        found = True
                        target_CDS, t_CDS_coord = insertion_UTR[0]
                #                 else:
                #                     print('check here- inserting')
                #                     continue
                target = target_full + target_insertion
                # ####### FULL LENGTH DUPL: BOTH #######
                if query + target == 2:
                    both = 1
                    query, target = 0, 0
                    tuples_obligatory_events.append((fragment_id, gene_id, mrna,
                                                    trans_coord.lower, trans_coord.upper,
                                                    cds_s, cds_e,
                                                    query_CDS, query_s, query_e,
                                                    target_CDS, t_CDS_coord.lower, t_CDS_coord.upper,
                                                    target_s, target_e, target_t))
                elif query + target == 0:
                    neither = 1
                tuples_full_length_duplications.append((fragment_id, gene_id, mrna,
                                                        cds_s, cds_e, query_s, query_e,
                                                        target_s, target_e,
                                                        neither, query, target, both, evalue))
                # ####### FULL LENGTH DUPL: TRUNCATION #######
                if not found:
                    truncation = [(i['id'], i['type'], i['coord']) for i in trans_dict['structure']
                                  if (i['coord'].overlaps(target_intv)
                                      and not i['coord'].contains(target_intv))]
                    reconc_exon_cds = [j for i in [annot for annot in truncation if annot[1] == 'exon']
                                       for j in [annot for annot in truncation if annot[1] == 'CDS']
                                       if (j[-1] & target_intv) == (i[-1] & target_intv)]
                    if reconc_exon_cds:
                        truncation = sorted([*[i for i in truncation if i[1] not in ['CDS', 'exon']],
                                             *reconc_exon_cds],
                                            key=lambda x: (x[-1].lower, x[-1].upper))
                    if truncation:
                        if len(truncation) == 2:
                            a, b = truncation
                            id_a, _, coord_a = a
                            id_b, _, coord_b = b
                            seg_a, seg_b = coord_a & target_intv, coord_b & target_intv
                            tuples_truncation_events.append(
                                (fragment_id, gene_id,
                                 mrna, trans_coord.lower, trans_coord.upper,
                                 query_CDS, cds_s, cds_e, query_s, query_e,
                                 target_s, target_e,
                                 id_a, coord_a.lower, coord_a.upper, seg_a.lower, seg_a.upper,
                                 id_b, coord_b.lower, coord_b.upper, seg_b.lower, seg_b.upper))
                        else:
                            print('look here')
                            print('truncation', gene_id, mrna, truncation, target_intv)
        instert_full_length_event(self.results_db, self.timeout_db, tuples_full_length_duplications)
        instert_obligatory_event(self.results_db, self.timeout_db, tuples_obligatory_events)
        instert_truncation_event(self.results_db, self.timeout_db, tuples_truncation_events)

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
            print('All genes already processed, '
                  'if you want to re-run the analysis, '
                  'delete/rename the results DB.')
        tic = time.time()
        print('- Identifying full length duplications', end=' ')
        self.find_full_length_duplications()
        create_mrna_counts_view(self.results_db, self.timeout_db)
        create_cumulative_counts_view(self.results_db, self.timeout_db)
        create_exclusive_pairs_view(self.results_db, self.timeout_db)
        hms_time = dt.strftime(dt.utcfromtimestamp(time.time() - tic), '%H:%M:%S')
        print(f' Done! [{hms_time}]')
