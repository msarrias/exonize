from utils import *                                    # helper functions
import gffutils                                        # for creating/loading DBs
import subprocess                                      # for calling gffread
import portion as P                                    # for working with intervals
import os                                              # for working with files
import time                                            # for sleeping between BLAST calls and for timeout on DB creation
import random                                          # for random sleep
import sqlite3                                         # for working with SQLite
import re                                              # regular expressions for genome masking
# import tempfile                                        # for creating temporary files
from Bio import SeqIO                                  # for reading FASTA files
from tqdm import tqdm                                  # progress bar
from multiprocessing.pool import ThreadPool            # for parallelization
from Bio.Blast import NCBIXML                          # for parsing BLAST results
from datetime import datetime as dt                    # for timeout


class Exonize(object):
    def __init__(self,
                 gff_file_path,
                 genome_path,
                 specie_identifier,
                 spec_attribute='ID',
                 verbose=True,
                 hard_masking=False,
                 evalue_threshold=1e-1,
                 min_align_len_perc=0.3,
                 sleep_max_seconds=5,
                 min_exon_length=20,
                 batch_number=10,
                 threads=6,
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
        self.stop_codons = ["TAG", "TGA", "TAA"]
        self.db_path = f'{self.specie_identifier}_genome_annotations.db'
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
            print('---------------------Warning------------------------------')
            print("-The genomic annotations do not contain intron annotations")
            print('----------------------------------------------------------')
            if self.verbose:
                print(f"- Attempting to write intron annotations in database:", end=" ")
            try:



                def intron_id(f):
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
                    mrna_dict['mRNAs'][mrna_annot.id] = dict(coord=mrna_coord, strand=gene.strand, structure=[])
                    temp_mrna_transcript = []
                    for child in self.db.children(mrna_annot.id, featuretype=self.feat_of_interest, order_by='start'):
                        coord = P.open(child.start - 1, child.end)
                        if coord:
                            temp_mrna_transcript.append({'id': child.id,
                                                         'coord': coord,  # ID attribute
                                                         'frame': child.frame,  # One of '0', '1' or '2'. '0'
                                                         'type': child.featuretype})  # feature type name
                    mrna_dict['mRNAs'][mrna_annot.id]['structure'] = sort_list_intervals_dict(temp_mrna_transcript)
                self.gene_hierarchy_dict[gene.id] = mrna_dict
                dump_pkl_file(self.gene_hierarchy_path, self.gene_hierarchy_dict)

    def run_tblastx(self, gene_id: str, mrna_id: str, annot_id: str,
                    query_seq: str, hit_seq: str, query_coord) -> dict:
        annot_id = annot_id.replace(':', '_')
        identifier = f'{gene_id}_{mrna_id}_{annot_id}'
        output_file = f'output/{identifier}_output.xml'
        if not os.path.exists(output_file):
            query_filename = f'input/{identifier}_query.fa'
            target_filename = f'input/{identifier}_target.fa'
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
            mrna_coord = self.gene_hierarchy_dict[gene_id]['mRNAs'][mrna_id]['coord']
            try:
                temp = self.parse_tblastx_output(blast_records, query_seq, hit_seq,
                                                 query_coord, mrna_coord, annot_id)
            except KeyError:
                print(f"KeyError: {gene_id} {mrna_id} {annot_id}")
        return temp

    def parse_tblastx_output(self, blast_records, query_seq: str, hit_seq: str, q_coord, hit_coord, query_id) -> dict:
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
                query_len = len(query_seq)
                print(f"No alignments found {query_id} --- query length: {query_len}")
                continue
            alignment = blast_record.alignments[0]  # Assuming only one alignment per blast_record
            if len([aln for aln in blast_record.alignments]) > 1:
                print("More than one alignment per blast_record")
            for hsp_idx, hsp in enumerate(alignment.hsps):
                blast_target_coord = P.open((hsp.sbjct_start - 1) + hit_coord.lower, hsp.sbjct_end + hit_coord.lower)
                query_aligned_frac = round((hsp.align_length * 3) / blast_record.query_length, 3)
                query_target_overlap_percentage = get_overlap_percentage(blast_target_coord, q_coord)
                if (hsp.expect < self.evalue
                        and query_aligned_frac > self.min_align_len_perc
                        and query_target_overlap_percentage < 0.5):
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
            cds_dict = sort_key_intervals_dict({i['coord']: {key: m for key, m in i.items() if key != 'coord'}
                                                for i in mrna_annot['structure'] if i['type'] == 'CDS'})
            if cds_dict:
                blast_mrna_cds_dict = {}
                for cds_coord, annot in cds_dict.items():
                    if (cds_coord.upper - cds_coord.lower) >= self.min_exon_len:  # only consider CDS queries with length >= 30
                        cds_seq = self.genome[chrom][cds_coord.lower:cds_coord.upper]
                        mrna_seq = self.genome[chrom][mrna_coord.lower:mrna_coord.upper]
                        temp = self.run_tblastx(gene_id, mrna_id, annot['id'], cds_seq, mrna_seq, cds_coord)
                        if temp:
                            blast_mrna_cds_dict[annot['id']] = dict(coord=cds_coord,
                                                                    frame=annot['frame'],
                                                                    tblastx_hits=temp)
                if blast_mrna_cds_dict:
                    mrna_blast_dict[mrna_id] = dict(coord=mrna_coord, blast_mrna_dict=blast_mrna_cds_dict)
        if mrna_blast_dict:
            self.insert_fragments_table(gene_id, mrna_blast_dict)
        else:
            self.insert_gene_ids_table(self.get_gene_tuple(gene_id, 0))

    def connect_create_results_db(self):
        db = sqlite3.connect(self.results_db, timeout=self.timeout_db)
        cursor = db.cursor()
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS Genes (
        gene_id VARCHAR(100) PRIMARY KEY,
        gene_chrom VARCHAR(100) NOT NULL,
        gene_strand VARCHAR(1) NOT NULL,
        gene_start INTEGER NOT NULL,
        gene_end INTEGER NOT NULL,
        has_duplicated_CDS BINARY(1) NOT NULL,
        UNIQUE(gene_id))
                """)
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS Fragments (
        fragment_id INTEGER PRIMARY KEY AUTOINCREMENT,
        gene_id  VARCHAR(100) NOT NULL REFERENCES Genes(gene_id),
        mrna_id VARCHAR(100) NOT NULL,
        mrna_start INTEGER NOT NULL,
        mrna_end INTEGER NOT NULL,
        CDS_id VARCHAR(100) NOT NULL,
        CDS_frame INTEGER NOT NULL,
        CDS_start INTEGER NOT NULL,
        CDS_end INTEGER NOT NULL,  
        query_frame INTEGER NOT NULL,
        query_strand VARCHAR(1) NOT NULL,
        target_frame INTEGER NOT NULL,
        target_strand VARCHAR(1) NOT NULL,
        score INTEGER NOT NULL,
        bits INTEGER NOT NULL,
        evalue REAL NOT NULL,
        alignment_len INTEGER NOT NULL,
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
        prot_perc_identity REAL NOT NULL)
        """)
        cursor.execute("""CREATE INDEX IF NOT EXISTS Fragments_idx ON Fragments (gene_id, mrna_id, CDS_id);""")
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS Fragments_matches (
        fragment_match_id INTEGER PRIMARY KEY AUTOINCREMENT,
        gene_id VARCHAR(100) NOT NULL REFERENCES Genes(gene_id),
        CDS_mrna_id VARCHAR(100) NOT NULL REFERENCES Fragments(mrna_id),
        CDS_id VARCHAR(100) NOT NULL REFERENCES Fragments(CDS_id),
        CDS_start INTEGER NOT NULL,
        CDS_end INTEGER NOT NULL,
        match_mrna_id VARCHAR(100) NOT NULL,
        match_id VARCHAR(100) NOT NULL,
        feature VARCHAR(100) NOT NULL,
        match_annot_start INTEGER NOT NULL,
        match_annot_end INTEGER NOT NULL,
        overlap_percentage REAL NOT NULL,
        UNIQUE(gene_id, CDS_mrna_id, CDS_id, match_mrna_id, match_id))
                """)
        cursor.execute("""
        CREATE INDEX IF NOT EXISTS Fragments_id_idx ON Fragments_matches (gene_id, CDS_mrna_id, CDS_id);
        """)
        db.commit()
        db.close()

    def get_fragments_matches_tuples(self, gene_id, mrna_id, event, target_coord) -> list:
        cds_mrna_id, cds_id, cds_start, cds_end, target_start, target_end, _ = event
        trans_structure = self.gene_hierarchy_dict[gene_id]['mRNAs'][mrna_id]['structure']
        return [(gene_id, cds_mrna_id, cds_id, cds_start, cds_end, mrna_id,
                 annot['id'], annot['type'], annot['coord'].lower, annot['coord'].upper,
                 round(get_overlap_percentage(target_coord, annot['coord']), 3))
                for annot in trans_structure if target_coord.overlaps(annot['coord'])]

    def insert_fragments_matches_table(self, tuples_list: list) -> None:
        db = sqlite3.connect(self.results_db, timeout=self.timeout_db)
        cursor = db.cursor()
        insert_frag_match_table_param = """
        INSERT OR IGNORE INTO Fragments_matches (
        gene_id,
        CDS_mrna_id,
        CDS_id,
        CDS_start,
        CDS_end,
        match_mrna_id,
        match_id,
        feature,
        match_annot_start,
        match_annot_end,
        overlap_percentage)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ? ,?)
        """
        cursor.executemany(insert_frag_match_table_param, tuples_list)
        db.commit()
        db.close()

    def query_gene_ids_in_res_db(self):
        db = sqlite3.connect(self.results_db, timeout=self.timeout_db)
        cursor = db.cursor()
        cursor.execute("SELECT gene_id FROM Genes")
        rows = cursor.fetchall()
        db.close()
        return [i[0] for i in rows]

    def insert_gene_ids_table(self, gene_args_tuple: tuple) -> None:
        db = sqlite3.connect(self.results_db, timeout=self.timeout_db)
        cursor = db.cursor()
        insert_gene_table_param = """  
        INSERT INTO Genes (
        gene_id,
        gene_chrom,
        gene_strand,
        gene_start, 
        gene_end,  
        has_duplicated_CDS) 
        VALUES (?, ?, ?, ?, ?, ?)
        """
        cursor.execute(insert_gene_table_param, gene_args_tuple)
        db.commit()
        db.close()

    def insert_fragments_table(self, gene_id: str, blast_mrna_cds_dict: dict) -> None:
        tuple_list = [get_fragment_tuple(gene_id, mrna_id, blast_cds_dict['coord'], cds_id, cds_dict, hsp_idx)
                      for mrna_id, blast_cds_dict in blast_mrna_cds_dict.items()
                      for cds_id, cds_dict in blast_cds_dict['blast_mrna_dict'].items()
                      for hsp_idx, hsp_dict in cds_dict['tblastx_hits'].items()]

        insert_fragments_table_param = """
        INSERT INTO Fragments (
        gene_id,
        mrna_id,
        mrna_start,
        mrna_end,
        CDS_id,
        CDS_frame,
        CDS_start,
        CDS_end,
        query_frame,
        query_strand,
        target_frame,
        target_strand,
        score,
        bits,
        evalue,
        alignment_len,
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
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """
        insert_gene_table_param = """
        INSERT INTO Genes
        (gene_id,
        gene_chrom,
        gene_strand,
        gene_start,
        gene_end,
        has_duplicated_CDS)
        VALUES (?, ?, ?, ?, ?, ?)
        """
        db = sqlite3.connect(self.results_db, timeout=self.timeout_db)
        cursor = db.cursor()
        cursor.execute(insert_gene_table_param, self.get_gene_tuple(gene_id, 1))
        cursor.executemany(insert_fragments_table_param, tuple_list)
        db.commit()
        db.close()

    def query_within_gene_events(self, gene_id):
        db = sqlite3.connect(self.results_db, timeout=self.timeout_db)
        cursor = db.cursor()
        fragments_query = """
        SELECT 
        mrna_id,
        CDS_id,
        CDS_start,
        CDS_end,
        target_start,
        target_end,
        MIN(evalue)
        FROM Fragments
        WHERE gene_id==?
        GROUP BY gene_id, mrna_id, CDS_id
        HAVING (ABS(query_start - target_start) > 5 OR ABS(query_end - target_end) > 5)
        """
        cursor.execute(fragments_query, (gene_id,))
        records = cursor.fetchall()
        db.close()
        return records

    def query_genes_with_duplicated_cds(self):
        db = sqlite3.connect(self.results_db, timeout=self.timeout_db)
        cursor = db.cursor()
        cursor.execute("""
        SELECT gene_id 
        FROM Genes 
        WHERE has_duplicated_CDS==1
        """)
        rows = cursor.fetchall()
        db.close()
        return [i[0] for i in rows]

    def identify_events(self):
        genes_with_duplicated_cds = self.query_genes_with_duplicated_cds()
        for gene_id in genes_with_duplicated_cds:
            gene_events = self.query_within_gene_events(gene_id)
            for event in gene_events:
                gene_mrnas_list = list(self.gene_hierarchy_dict[gene_id]['mRNAs'].keys())
                for mrna_id in gene_mrnas_list:
                    target_start, target_end = event[4], event[5]
                    trans_coord = self.gene_hierarchy_dict[gene_id]['mRNAs'][mrna_id]['coord']
                    target_coord = P.open(target_start + trans_coord.lower, target_end + trans_coord.lower)
                    if trans_coord.overlaps(target_coord):
                        fragments_matches_tuples_list = self.get_fragments_matches_tuples(gene_id, mrna_id, event,
                                                                                          target_coord)
                        if fragments_matches_tuples_list:
                            self.insert_fragments_matches_table(fragments_matches_tuples_list)
                        else:
                            print(f"check this case{mrna_id, gene_id, event}")

    def prepare_data(self):
        self.create_parse_or_update_database()
        if os.path.exists(self.gene_hierarchy_path):
            self.gene_hierarchy_dict = read_pkl_file(self.gene_hierarchy_path)
        else:
            self.create_gene_hierarchy_dict()
        self.read_genome()
        self.connect_create_results_db()

    def run_analysis(self) -> None:
        exonize()
        self.prepare_data()
        args_list = list(self.gene_hierarchy_dict.keys())
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
            hms_time = dt.strftime(dt.utcfromtimestamp(time.time() - tic), '%H:%M:%S')
            print(f' Done! [{hms_time}]')
            tic = time.time()
            print('- Identifying region of events', end=' ')
            self.identify_events()
            hms_time = dt.strftime(dt.utcfromtimestamp(time.time() - tic), '%H:%M:%S')
            print(f' Done! [{hms_time}]')
