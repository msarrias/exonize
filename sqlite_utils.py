import sqlite3                                         # for working with SQLite
from utils import *


# #### CREATE TABLES #####
def connect_create_results_db(db_path, timeout_db) -> None:
    db = sqlite3.connect(db_path, timeout=timeout_db)
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
    # cursor.execute("""
    # CREATE TABLE IF NOT EXISTS Fragments_matches (
    # fragment_match_id INTEGER PRIMARY KEY AUTOINCREMENT,
    # fragment_id INTEGER NOT NULL REFERENCES Fragments(fragment_id),
    # gene_id VARCHAR(100) NOT NULL REFERENCES Genes(gene_id),
    # mrna_id VARCHAR(100) NOT NULL,
    # CDS_start INTEGER NOT NULL,
    # CDS_end INTEGER NOT NULL,
    # Neither BINARY(1) NOT NULL,
    # query BINARY(1) NOT NULL,
    # target BINARY(1) NOT NULL,
    # Both BINARY(1) NOT NULL,
    # UNIQUE(fragment_id, gene_id, mrna_id, CDS_start, CDS_end))""")
    # cursor.execute("""
    # CREATE INDEX IF NOT EXISTS Fragments_id_idx ON Fragments_matches (gene_id, mrna_id);
    # """)
    db.commit()
    db.close()


def create_mutually_exclusive_events_table(db_path, timeout_db) -> None:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS MXEs_events (
    event_id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_id VARCHAR(100) NOT NULL REFERENCES Genes(gene_id),
    mrna_A VARCHAR(100) NOT NULL,
    mrna_B VARCHAR(100) NOT NULL,
    CDS_A VARCHAR(100) NOT NULL,
    CDS_A_start INTEGER NOT NULL,
    CDS_A_end INTEGER NOT NULL,
    CDS_B VARCHAR(100) NOT NULL,
    CDS_B_start INTEGER NOT NULL,
    CDS_B_end INTEGER NOT NULL,
    UNIQUE(gene_id, mrna_A, mrna_B))
    """)


def create_obligatory_events_table(db_path, timeout_db) -> None:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("""
            CREATE TABLE IF NOT EXISTS Obligatory_events (
            event_id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_id VARCHAR(100) NOT NULL REFERENCES Genes(gene_id),
            mrna_id VARCHAR(100) NOT NULL,
            mrna_start INTEGER NOT NULL,
            mrna_end INTEGER NOT NULL,
            CDS_A VARCHAR(100) NOT NULL,
            CDS_A_start INTEGER NOT NULL,
            CDS_A_end INTEGER NOT NULL,
            CDS_B VARCHAR(100) NOT NULL,
            CDS_B_start INTEGER NOT NULL,
            CDS_B_end INTEGER NOT NULL,
            UNIQUE(gene_id, mrna_id, CDS_A, CDS_B))
            """)


def query_gene_ids_in_res_db(db_path, timeout_db) -> list:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("SELECT gene_id FROM Genes")
    rows = cursor.fetchall()
    db.close()
    return [i[0] for i in rows]


# #### INSERT INTO TABLES #####
def insert_fragments_matches_table(db_path, timeout_db, tuples_list: list) -> None:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    insert_frag_match_table_param = """
     INSERT OR IGNORE INTO Fragments_matches (
     gene_id,
     CDS_mrna_id,
     CDS_mrna_start,
     CDS_mrna_end,
     CDS_id,
     CDS_start,
     CDS_end,
     query_start,
     query_end,
     target_start,
     target_end,
     match_mrna_id,
     match_mrna_start,
     match_mrna_end,
     match_id,
     feature,
     match_annot_start,
     match_annot_end,
     overlap_percentage_CDS,
     overlap_percentage_match)
     VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ? ,?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
     """
    cursor.executemany(insert_frag_match_table_param, tuples_list)
    db.commit()
    db.close()


def insert_gene_ids_table(db_path, timeout_db, gene_args_tuple: tuple) -> None:
    db = sqlite3.connect(db_path, timeout=timeout_db)
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


def insert_fragments_calls():
    insert_fragments_table_param = """
            INSERT INTO Fragments (
            gene_id,
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
            VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
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
    return insert_fragments_table_param, insert_gene_table_param


# ###### QUERY TABLES ######
def query_within_gene_events(db_path, timeout_db, gene_id: str) -> list:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    fragments_query = """
    SELECT 
    mrna_id,
    CDS_id,
    CDS_start,
    CDS_end,
    query_start,
    query_end,
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


def query_genes_with_duplicated_cds(db_path, timeout_db) -> list:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("""
    SELECT gene_id 
    FROM Genes 
    WHERE has_duplicated_CDS==1
    """)
    rows = cursor.fetchall()
    db.close()
    return [i[0] for i in rows]
