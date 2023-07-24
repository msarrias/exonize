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
    prot_perc_identity REAL NOT NULL,
    UNIQUE(fragment_id, gene_id, CDS_start, CDS_end, query_start, query_end, target_start, target_end))
    """)
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS Full_length_duplications (
    fragment_match_id INTEGER PRIMARY KEY AUTOINCREMENT,
    fragment_id INTEGER NOT NULL REFERENCES Fragments(fragment_id),
    gene_id VARCHAR(100) NOT NULL REFERENCES Genes(gene_id),
    mrna_id VARCHAR(100) NOT NULL,
    CDS_start INTEGER NOT NULL,
    CDS_end INTEGER NOT NULL,
    neither BINARY(1) NOT NULL,
    query BINARY(1) NOT NULL,
    target BINARY(1) NOT NULL,
    both BINARY(1) NOT NULL,
    UNIQUE(fragment_id, gene_id, mrna_id, CDS_start, CDS_end))""")
    cursor.execute("""
    CREATE INDEX IF NOT EXISTS Full_length_duplications_idx ON Full_length_duplications (fragment_id, gene_id, mrna_id);
    """)
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
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS Obligatory_events (
        event_id INTEGER PRIMARY KEY AUTOINCREMENT,
        fragment_id INTEGER NOT NULL REFERENCES Fragments(fragment_id),
        gene_id VARCHAR(100) NOT NULL REFERENCES Genes(gene_id),
        mrna_id VARCHAR(100) NOT NULL,
        mrna_start INTEGER NOT NULL,
        mrna_end INTEGER NOT NULL,
        query_CDS_start INTEGER NOT NULL,
        query_CDS_end INTEGER NOT NULL,
        query_CDS_id VARCHAR(100) NOT NULL,
        query_start INTEGER NOT NULL,
        query_end INTEGER NOT NULL,
        target_CDS_id VARCHAR(100) NOT NULL,
        target_CDS_start INTEGER NOT NULL,
        target_CDS_end INTEGER NOT NULL,
        target_start INTEGER NOT NULL,
        target_end INTEGER NOT NULL,
        UNIQUE(fragment_id, gene_id, mrna_id, query_CDS_id, target_CDS_id))
        """)
    db.commit()
    db.close()


def query_gene_ids_in_res_db(db_path, timeout_db) -> list:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("SELECT gene_id FROM Genes")
    rows = cursor.fetchall()
    db.close()
    return [i[0] for i in rows]


# #### INSERT INTO TABLES #####
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


def instert_full_length_event(db_path, timeout_db, tuples_list):
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    insert_full_length_event_table_param = """
    INSERT INTO Full_length_duplications (
    fragment_id,
    gene_id,
    mrna_id,
    CDS_start,
    CDS_end,
    neither,
    query,
    target,
    both )
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?) 
    """
    cursor.executemany(insert_full_length_event_table_param, tuples_list)
    db.commit()
    db.close()


def instert_obligatory_event(db_path, timeout_db, tuples_list):
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    insert_obl_event_table_param = """
     INSERT OR IGNORE INTO Obligatory_events (
     fragment_id,
     gene_id,
     mrna_id,
     mrna_start,
     mrna_end,
     query_CDS_start,
     query_CDS_end,
     query_CDS_id,
     query_start,
     query_end,
     target_CDS_id,
     target_CDS_start,
     target_CDS_end,
     target_start,
     target_end)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
     """
    cursor.executemany(insert_obl_event_table_param, tuples_list)
    db.commit()
    db.close()


# ###### QUERY TABLES ######
def query_within_gene_events(db_path: str, timeout_db: float) -> list:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    fragments_query = """
    SELECT 
        f.fragment_id,
        f.gene_id,
        g.gene_start,
        g.gene_end,
        f.CDS_start,
        f.CDS_end,
        f.query_start,
        f.query_end ,
        f.target_start ,
        f.target_end 
    FROM Fragments as f
    JOIN Genes AS g ON g.gene_id = f.gene_id
    """
    cursor.execute(fragments_query)
    records = cursor.fetchall()
    db.close()
    return [i for i in records]


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
