from utils import *                                                                 # helper functions
import gffutils                                                                     # for creating/loading DBs
import subprocess                                                                   # for calling gffread
import portion as P                                                                 # for working with intervals
import os                                                                           # for working with files
import time                                                                         # for sleeping between BLAST calls and for timeout on DB creation
import random                                                                       # for random sleep
import sqlite3                                                                      # for working with SQLite
import re                                                                           # regular expressions for genome masking
import tempfile                                                                     # for creating temporary files
from Bio import SeqIO                                                               # for reading FASTA files
from tqdm import tqdm                                                               # progress bar
from multiprocessing.pool import ThreadPool                                         # for parallelization
from Bio.Blast import NCBIXML                                                       # for parsing BLAST results
from datetime import datetime as dt                                                 # for timeout


class Exonize(object):
    def __init__(self,
                 gff_file_path,
                 genome_path,
                 specie_identifier,
                 verbose=True,
                 hard_masking=False,
                 evalue_threshold=1e-5,
                 min_align_len_perc=0.9,
                 sleep_max_seconds=5,
                 min_exon_length=10,
                 batch_number=10,
                 threads=6,
                 timeout_db=160.0):

        self.genome = None                                                          # genome sequence
        self.gene_hierarchy_dict = None                                             # gene hierarchy dictionary (gene -> transcript -> exon)
        self.db_features = None                                                     # features in the database
        self.old_filename = None                                                    # old filename (if GTF file is provided)
        self.db = None                                                              # database object (gffutils)
        self.specie_identifier = specie_identifier                                  # specie identifier
        self.verbose = verbose                                                      # verbose mode (True/False)
        self.in_file_path = gff_file_path                                           # input file path (GFF/GTF)
        self.genome_path = genome_path                                              # genome path (FASTA)
        self.hard_masking = hard_masking                                            # hard masking (True/False)
        self.secs = sleep_max_seconds                                               # max seconds to sleep between BLAST calls
        self.min_exon_len = min_exon_length                                         # minimum exon length (bp)
        self.timeout_db = timeout_db                                                # timeout for creating the database (seconds)
        self.evalue = evalue_threshold                                              # e-value threshold for BLAST
        self.min_align_len_perc = min_align_len_perc                                # minimum alignment length percentage of query length
        self.batch_number = batch_number                                            # batch number for BLAST calls
        self.threads = threads                                                      # number of threads for BLAST calls
        self.stop_codons = ["TAG", "TGA", "TAA"]
        self.db_path = f'{self.specie_identifier}_genome_annotations.db'
        self.results_df = f'{self.specie_identifier}_results.db'
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
            if self.verbose:
                print("- Reading annotations database:", end=" ")
            self.load_db()
            if self.verbose:
                print("Done!")
        self.db_features = list(self.db.featuretypes())
        self.create_intron_annotations()

    def convert_gtf_to_gff(self) -> None:
        gffread_command = ["gffread", self.old_filename, "-O", "-o", self.in_file_path]
        subprocess.call(gffread_command)

    def create_database(self) -> None:
        try:
            if self.verbose:
                print("Creating annotations database", end=" ")
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
        except ValueError:
            print("Wrong infile path")

    def load_db(self) -> None:
        try:
            self.db = gffutils.FeatureDB(self.db_path, keep_order=True)
        except ValueError:
            print("Wrong db file path")

    def create_intron_annotations(self) -> None:
        """
        run only once otherwise duplicated
        annotations will be created
        """
        if 'intron' not in self.db_features:
            if self.verbose:
                print(f"Writing intron annotations in database:", end=" ")
            self.db.update(list(self.db.create_introns()),
                           id_spec=dict(intron=[lambda f: ''.join(f"{f[self.specie_identifier]}_{f['Parent']}")]),
                           make_backup=False)
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
            tac_genome = time.time()
            if self.verbose:
                elapsed = tac_genome - tic_genome
                hms_time = dt.strftime(dt.utcfromtimestamp(elapsed), '%H:%M:%S')
                print(f"[{hms_time}]")
        except ValueError:
            print("Wrong genome file path")

    def create_gene_hierarchy_dict(self):
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
                    mrna_dict['mRNAs'][mrna_annot.id] = dict(coord=mrna_coord, strand=gene.strand, structure={})
                    temp_mrna_transcript = {}
                    for child in self.db.children(mrna_annot.id, featuretype=self.feat_of_interest, order_by='start'):
                        coord = P.open(child.start - 1, child.end)
                        if coord:
                            temp_mrna_transcript[coord] = {'id': child.id,  # ID attribute
                                                           'frame': child.frame,  # One of '0', '1' or '2'. '0'
                                                           'type': child.featuretype}  # feature type name
                    mrna_dict['mRNAs'][mrna_annot.id]['structure'] = sort_intervals_dict(temp_mrna_transcript)
                self.gene_hierarchy_dict[gene.id] = mrna_dict
                dump_pkl_file(self.gene_hierarchy_path, self.gene_hierarchy_dict)

    def run_tblastx(self, gene_id: str, mrna_id: str, annot_id: str, query_seq: str, hit_seq: str, query_coord) -> dict:
        with tempfile.TemporaryDirectory() as temp_dir_name:
            query_filename = f'{temp_dir_name}/query.fa'
            target_filename = f'{temp_dir_name}/target.fa'
            output_file = f'{temp_dir_name}/output.xmp'
            dump_fasta_file(query_filename, {annot_id: query_seq})
            dump_fasta_file(target_filename, {f"{gene_id}_{mrna_id}": hit_seq})
            tblastx_command = ['tblastx',
                               '-query', query_filename,
                               '-subject', target_filename,
                               '-outfmt', '5',  # XML output format
                               '-out', output_file]
            subprocess.run(tblastx_command)
            with open(output_file, "r") as result_handle:
                blast_records = NCBIXML.parse(result_handle)
                gene_coord = self.gene_hierarchy_dict[gene_id]['coord']
                temp = self.parse_tblastx_output(blast_records, query_seq, hit_seq, query_coord, gene_coord)
        return temp

    def parse_tblastx_output(self, blast_records, query_seq: str, hit_seq: str, q_coord, gene_coord) -> dict:
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
            alignment = blast_record.alignments[0]  # Assuming only one alignment per blast_record
            for hsp_idx, hsp in enumerate(alignment.hsps):
                blast_target_coord = P.open((hsp.sbjct_start - 1) + gene_coord.lower, hsp.sbjct_end + gene_coord.lower)
                query_aligned_frac = round((hsp.align_length * 3) / blast_record.query_length, 3)
                query_target_overlapping_percentage = get_overlapping_percentage(blast_target_coord, q_coord)
                if (hsp.expect < self.evalue
                        and query_aligned_frac > self.min_align_len_perc
                        and query_target_overlapping_percentage < 0.5):
                    res_tblastx[hsp_idx] = get_hsp_dict(hsp, query_seq, hit_seq)
        return res_tblastx

    def get_gene_tuple(self, gene_id, bin_has_dup):
        gene_coord = self.gene_hierarchy_dict[gene_id]['coord']
        return (gene_id,
                self.gene_hierarchy_dict[gene_id]['chrom'],
                (self.gene_hierarchy_dict[gene_id]['strand']),
                gene_coord.lower,
                gene_coord.upper,
                bin_has_dup)

    def find_coding_exon_duplicates(self, gene_id) -> None:
        """
        The reading frame of the CDS queries is respected in the tblastx search.
        """
        time.sleep(random.randrange(0, self.secs))
        mrna_blast_dict = {}
        chrom = self.gene_hierarchy_dict[gene_id]['chrom']
        for mrna_id, mrna_annot in self.gene_hierarchy_dict[gene_id]['mRNAs'].items():
            mrna_coord = mrna_annot['coord']
            cds_list = {coord: annot for coord, annot in mrna_annot['structure'].items() if annot['type'] == 'CDS'}
            if cds_list:
                blast_mrna_cds_dict = {}
                for coord, annot in cds_list.items():
                    cds_frame = int(annot['frame'])
                    if (coord.upper - coord.lower) >= 30:  # only consider CDS queries with length >= 30
                        cds_seq = self.genome[chrom][(coord.lower+cds_frame):coord.upper]
                        mrna_seq = self.genome[chrom][mrna_coord.lower:mrna_coord.upper]
                        temp = self.run_tblastx(gene_id, mrna_id, annot['id'], cds_seq, mrna_seq, coord)
                        if temp:
                            blast_mrna_cds_dict[annot['id']] = dict(coord=coord,
                                                                    frame=annot['frame'],
                                                                    tblastx_hits=temp)
                if blast_mrna_cds_dict:
                    mrna_blast_dict[mrna_id] = blast_mrna_cds_dict
        if mrna_blast_dict:
            self.insert_fragments_table(gene_id, mrna_blast_dict)
        else:
            self.insert_gene_ids_table(self.get_gene_tuple(gene_id, 0))

    def connect_create_results_db(self):
        db = sqlite3.connect(self.results_df, timeout=self.timeout_db)
        cursor = db.cursor()
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS GeneIDs (
        gene_id VARCHAR(100) PRIMARY KEY,
        gene_chrom VARCHAR(100) NOT NULL,
        gene_strand VARCHAR(1) NOT NULL,
        gene_start INTEGER NOT NULL,
        gene_end INTEGER NOT NULL,
        has_duplicated_exon BINARY(1) NOT NULL,
        UNIQUE(gene_id))
                """)

        cursor.execute("""
        CREATE TABLE IF NOT EXISTS Fragments (
        fragment_id INTEGER PRIMARY KEY AUTOINCREMENT,
        gene_id  VARCHAR(100) NOT NULL REFERENCES GeneIDs(gene_id),
        mrna_id VARCHAR(100) NOT NULL,
        query_id VARCHAR(100) NOT NULL,
        query_frame INTEGER NOT NULL,
        hit_query_frame VARCHAR(100) NOT NULL,
        hit_target_frame VARCHAR(100) NOT NULL,
        score INTEGER NOT NULL,
        bits INTEGER NOT NULL,
        evalue REAL NOT NULL,
        alignment_len INTEGER NOT NULL,
        cds_start INTEGER NOT NULL,
        cds_end INTEGER NOT NULL,               
        query_start INTEGER NOT NULL,
        query_end INTEGER NOT NULL,
        target_start INTEGER NOT NULL,
        target_end INTEGER NOT NULL,
        query_dna_seq VARCHAR NOT NULL,
        target_dna_seq VARCHAR NOT NULL,
        query_aln_prot_seq VARCHAR NOT NULL,
        target_aln_prot_seq VARCHAR NOT NULL,
        match VARCHAR NOT NULL,
        query_num_stop_codons INTEGER NOT NULL,
        target_num_stop_codons INTEGER NOT NULL,
        dna_perc_identity REAL NOT NULL,
        prot_perc_identity REAL NOT NULL
        )
        """)
        cursor.execute(
            """CREATE INDEX IF NOT EXISTS Fragments_idx ON Fragments (gene_id, mrna_id, query_id);""")
        db.commit()
        db.close()

    def query_gene_ids_in_res_db(self):
        db = sqlite3.connect(self.results_df, timeout=self.timeout_db)
        cursor = db.cursor()
        cursor.execute("SELECT gene_id FROM GeneIDs")
        rows = cursor.fetchall()
        db.close()
        return [i[0] for i in rows]

    def insert_gene_ids_table(self, gene_args_tuple: tuple) -> None:
        db = sqlite3.connect(self.results_df, timeout=self.timeout_db)
        cursor = db.cursor()
        insert_gene_table_param = """  
        INSERT INTO GeneIDs (
        gene_id,
        gene_chrom,
        gene_strand,
        gene_start, 
        gene_end,  
        has_duplicated_exon) 
        VALUES (?, ?, ?, ?, ?, ?)
        """
        cursor.execute(insert_gene_table_param, gene_args_tuple)
        db.commit()
        db.close()

    def insert_fragments_table(self, gene_id: str, blast_mrna_cds_dict: dict) -> None:
        tuple_list = [get_fragment_tuple(gene_id, mrna_id, cds_id, cds_dict, hsp_idx)
                      for mrna_id, blast_cds_dict in blast_mrna_cds_dict.items()
                      for cds_id, cds_dict in blast_cds_dict.items()
                      for hsp_idx, hsp_dict in cds_dict['tblastx_hits'].items()]

        insert_fragments_table_param = """
        INSERT INTO Fragments (
        gene_id,
        mrna_id,
        query_id,
        query_frame,
        hit_query_frame,
        hit_target_frame,
        score,
        bits,
        evalue,
        alignment_len,
        cds_start,
        cds_end,
        query_start,
        query_end,
        target_start,
        target_end,
        query_dna_seq,
        target_dna_seq,
        query_aln_prot_seq,
        target_aln_prot_seq,
        match,
        query_num_stop_codons,
        target_num_stop_codons,
        dna_perc_identity,
        prot_perc_identity
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """
        insert_gene_table_param = """
        INSERT INTO GeneIDs
        (gene_id,
        gene_chrom,
        gene_strand,
        gene_start,
        gene_end,
        has_duplicated_exon)
        VALUES (?, ?, ?, ?, ?, ?)
        """
        db = sqlite3.connect(self.results_df, timeout=self.timeout_db)
        cursor = db.cursor()
        cursor.execute(insert_gene_table_param, self.get_gene_tuple(gene_id, 1))
        cursor.executemany(insert_fragments_table_param, tuple_list)
        db.commit()
        db.close()

    def run_analysis(self) -> None:
        self.create_parse_or_update_database()
        if os.path.exists(self.gene_hierarchy_path):
            self.gene_hierarchy_dict = read_pkl_file(self.gene_hierarchy_path)
        else:
            self.create_gene_hierarchy_dict()
        self.read_genome()
        self.connect_create_results_db()
        args_list = list(self.gene_hierarchy_dict.keys())[0:10]
        processed_gene_ids = self.query_gene_ids_in_res_db()
        if processed_gene_ids:
            args_list = [i for i in args_list if i not in processed_gene_ids]
        if args_list:
            batches_list = [i for i in batch(args_list, self.batch_number)]
            tic = time.time()
            if self.verbose:
                print(f'- Starting exon duplication search for {len(args_list)} genes.')
            with tqdm(total=len(args_list)) as progress_bar:
                for arg_batch in batches_list:
                    t = ThreadPool(processes=self.threads)
                    t.map(self.find_coding_exon_duplicates, arg_batch)
                    t.close()
                    t.join()
                    progress_bar.update(len(arg_batch))
            tac = time.time()
            elapsed = tac - tic
            hms_time = dt.strftime(dt.utcfromtimestamp(elapsed), '%H:%M:%S')
            print(f'[{hms_time}]')
