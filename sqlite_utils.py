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
    gene_id  VARCHAR(100) NOT NULL,
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
    query_aln_prot_seq VARCHAR NOT NULL,
    target_aln_prot_seq VARCHAR NOT NULL,
    match VARCHAR NOT NULL,
    query_num_stop_codons INTEGER NOT NULL,
    target_num_stop_codons INTEGER NOT NULL,
    FOREIGN KEY (gene_id) REFERENCES Genes(gene_id),
    UNIQUE(fragment_id, gene_id, CDS_start, CDS_end, query_start, query_end, target_start, target_end))
    """)
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS Full_length_duplications (
    fragment_match_id INTEGER PRIMARY KEY AUTOINCREMENT,
    fragment_id INTEGER NOT NULL,
    gene_id VARCHAR(100) NOT NULL,
    mrna_id VARCHAR(100) NOT NULL,
    CDS_start INTEGER NOT NULL,
    CDS_end INTEGER NOT NULL,
    query_id VARCHAR(100) NOT NULL,
    query_start INTEGER NOT NULL,
    query_end INTEGER NOT NULL,
    event_type VARCHAR(100) NOT NULL,
    target_id VARCHAR(100), /* TRUNCTATION envents will take NULL value */
    annot_target_start INTEGER,
    annot_target_end INTEGER,
    target_start INTEGER NOT NULL,
    target_end INTEGER NOT NULL,
    neither BINARY(1) NOT NULL,
    query BINARY(1) NOT NULL,
    target BINARY(1) NOT NULL,
    both BINARY(1) NOT NULL,
    evalue REAL NOT NULL,
    FOREIGN KEY (fragment_id) REFERENCES Fragments(fragment_id),
    FOREIGN KEY (gene_id) REFERENCES Genes(gene_id),
    UNIQUE(fragment_id, gene_id, mrna_id, CDS_start, CDS_end, query_start, query_end, target_start, target_end))""")
    cursor.execute("""
    CREATE INDEX IF NOT EXISTS Full_length_duplications_idx 
    ON Full_length_duplications (fragment_id, gene_id, mrna_id, CDS_start, CDS_end);
    """)
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS Obligatory_events (
    event_id INTEGER PRIMARY KEY AUTOINCREMENT,
    fragment_id INTEGER NOT NULL,
    gene_id VARCHAR(100) NOT NULL,
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
    type VARCHAR(100) NOT NULL,
    FOREIGN KEY (fragment_id) REFERENCES Fragments(fragment_id),
    FOREIGN KEY (gene_id) REFERENCES Genes(gene_id),
    UNIQUE(fragment_id, gene_id, mrna_id, query_CDS_id, target_CDS_id))
    """)
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS Truncation_events (
    event_id INTEGER PRIMARY KEY AUTOINCREMENT,
    fragment_id INTEGER NOT NULL,
    gene_id VARCHAR(100) NOT NULL,
    mrna_id VARCHAR(100) NOT NULL,
    mrna_start INTEGER NOT NULL,
    mrna_end INTEGER NOT NULL,
    query_CDS_id VARCHAR(100) NOT NULL,
    query_CDS_start INTEGER NOT NULL,
    query_CDS_end INTEGER NOT NULL,
    query_start INTEGER NOT NULL,
    query_end INTEGER NOT NULL,
    target_start INTEGER NOT NULL,
    target_end INTEGER NOT NULL,
    id_B VARCHAR(100) NOT NULL,
    type_B VARCHAR(100) NOT NULL,
    B_annot_start INTEGER NOT NULL,
    B_annot_end INTEGER NOT NULL,
    target_B_start INTEGER NOT NULL,
    target_B_end INTEGER NOT NULL,
    FOREIGN KEY (fragment_id) REFERENCES Fragments(fragment_id),
    FOREIGN KEY (gene_id) REFERENCES Genes(gene_id),
    UNIQUE(fragment_id, gene_id, mrna_id, query_CDS_id, id_B))
    """)
    db.commit()
    db.close()


def create_cumulative_counts_table(db_path, timeout_db) -> None:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS Full_length_events_cumulative_counts AS
    SELECT * FROM (
    SELECT
    fn.fragment_id,
    fn.gene_id,
    fn.gene_start,
    fn.gene_end,
    fn.mrna_count,
    fn.CDS_start,
    fn.CDS_end,
    fn.query_start,
    fn.query_end,
    fn.target_start,
    fn.target_end,
    group_concat(fld.event_type) as concat_event_type,
    SUM(fld.both) AS cum_both,
    SUM(fld.query) AS cum_query,
    SUM(fld.target) AS cum_target,
    SUM(fld.neither) AS cum_neither,
    fn.evalue
    FROM (
    SELECT
    fm.fragment_id,
    fm.gene_id,
    g.gene_start,
    g.gene_end,
    gc.mrna_count,
    fm.CDS_start,
    fm.CDS_end,
    fm.query_start,
    fm.query_end,
    fm.target_start,
    fm.target_end,
    fm.evalue
    FROM Fragments as fm
    JOIN Genes AS g
    ON fm.gene_id = g.gene_id
    JOIN Genes_mRNA_counts AS gc
    ON fm.gene_id = gc.gene_id
    ORDER BY fm.fragment_id
    ) AS fn
    JOIN Full_length_duplications AS fld
    ON fn.gene_id = fld.gene_id
    AND fn.fragment_id = fld.fragment_id
    AND fn.CDS_start = fld.CDS_start
    AND fn.CDS_end = fld.CDS_end
    GROUP BY fn.fragment_id
    ORDER BY fn.fragment_id) AS fn2
    """)
    db.commit()
    db.close()


# #### INSERT INTO TABLES #####
def insert_identity_and_dna_algns_columns(db_path, timeout_db, fragments) -> None:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute(""" ALTER TABLE Fragments ADD COLUMN dna_perc_identity REAL;""")
    cursor.execute(""" ALTER TABLE Fragments ADD COLUMN prot_perc_identity REAL;""")
    cursor.execute(""" ALTER TABLE Fragments ADD COLUMN query_dna_seq VARCHAR;""")
    cursor.execute(""" ALTER TABLE Fragments ADD COLUMN target_dna_seq VARCHAR;""")
    cursor.executemany(""" 
    UPDATE Fragments 
    SET dna_perc_identity=?, prot_perc_identity=?, query_dna_seq=?, target_dna_seq=?  WHERE fragment_id=? 
    """, fragments)
    db.commit()
    db.close()


def instert_pair_id_column_to_full_length_events_cumulative_counts(db_path, timeout_db, fragments) -> None:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute(""" ALTER TABLE Full_length_events_cumulative_counts ADD COLUMN pair_id INTEGER;""")
    cursor.executemany(""" UPDATE Full_length_events_cumulative_counts SET pair_id=? WHERE fragment_id=? """, fragments)
    db.commit()
    db.close()


def insert_percent_query_column_to_fragments(db_path, timeout_db) -> None:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("""
    ALTER TABLE Fragments ADD COLUMN percent_query DECIMAL(10, 3);
    """)
    cursor.execute("""
    UPDATE Fragments 
    SET percent_query = ROUND(CAST(int.intersect_end - int.intersect_start AS REAL) / CAST(int.CDS_end - int.CDS_start AS REAL), 3)
    FROM (SELECT 
    Fragment_id,
    MAX(f.CDS_start, (f.query_start + f.CDS_start)) AS intersect_start,
    MIN(f.CDS_end, (f.query_end + f.CDS_start)) AS intersect_end,
    f.CDS_end,
    f.CDS_start
    FROM Fragments AS f
    WHERE f.CDS_end >= (f.query_start + f.CDS_start) 
    AND f.CDS_start <= (f.query_end + f.CDS_start)
    ) AS int
    WHERE Fragments.Fragment_id = int.Fragment_id;
    """)
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


def insert_fragments_calls() -> tuple:
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
    query_aln_prot_seq,
    target_aln_prot_seq,
    match,
    query_num_stop_codons,
    target_num_stop_codons
    )
    VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
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


def instert_full_length_event(db_path, timeout_db, tuples_list) -> None:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    insert_full_length_event_table_param = """
    INSERT INTO Full_length_duplications (
    fragment_id,
    gene_id,
    mrna_id,
    CDS_start,
    CDS_end,
    query_id,
    query_start,
    query_end,
    event_type,
    target_id,
    annot_target_start,
    annot_target_end,
    target_start,
    target_end,
    neither,
    query,
    target,
    both,
    evalue)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) 
    """
    cursor.executemany(insert_full_length_event_table_param, tuples_list)
    db.commit()
    db.close()


def instert_obligatory_event(db_path, timeout_db, tuples_list) -> None:
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
    target_end,
    type)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
     """
    cursor.executemany(insert_obl_event_table_param, tuples_list)
    db.commit()
    db.close()


def instert_truncation_event(db_path, timeout_db, tuples_list) -> None:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    insert_trunc_event_table_param = """
    INSERT OR IGNORE INTO Truncation_events (
    fragment_id,
    gene_id,
    mrna_id,
    mrna_start,
    mrna_end,
    query_CDS_id,
    query_CDS_start,
    query_CDS_end,
    query_start,
    query_end,
    target_start,
    target_end,
    id_B,
    type_B,
    B_annot_start,
    B_annot_end,
    target_B_start,
    target_B_end)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """
    cursor.executemany(insert_trunc_event_table_param, tuples_list)
    db.commit()
    db.close()


# ###### QUERY TABLES ######
def query_gene_ids_in_res_db(db_path, timeout_db) -> list:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("SELECT gene_id FROM Genes")
    rows = cursor.fetchall()
    db.close()
    return [i[0] for i in rows]


def query_filtered_full_duplication_events(db_path: str, timeout_db: float) -> list:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    fragments_query = """
    SELECT 
    fragment_id,
    gene_id,
    gene_start,
    gene_end,
    CDS_start,
    CDS_end,
    query_start,
    query_end,
    target_start,
    target_end,
    evalue
    FROM Filtered_full_length_events
    """
    cursor.execute(fragments_query)
    records = cursor.fetchall()
    db.close()
    return [i for i in records]


def query_tuples_full_length_duplications(db_path, timeout_db) -> list:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("""
    SELECT *
    FROM Full_length_duplications;
    """)
    rows = cursor.fetchall()
    db.close()
    return rows


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


def query_obligatory_events(db_path, timeout_db) -> list:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("""
    SELECT *
    FROM Obligatory_events;
    """)
    rows = cursor.fetchall()
    db.close()
    return rows


def query_truncation_events(db_path, timeout_db) -> list:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("""
    SELECT *
    FROM Truncation_events;
    """)
    rows = cursor.fetchall()
    db.close()
    return rows


def query_exclusive_events(db_path, timeout_db) -> list:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("""
    SELECT *
    FROM Exclusive_pairs;
    """)
    rows = cursor.fetchall()
    db.close()
    return rows


def query_queries_only(db_path, timeout_db) -> list:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("""
    SELECT *
    FROM Full_length_events_cumulative_counts AS fle
    WHERE fle.cum_query==fle.mrna_count;
    """)
    rows = cursor.fetchall()
    db.close()
    return rows


def query_targets_only(db_path, timeout_db) -> list:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("""
    SELECT *
    FROM Full_length_events_cumulative_counts AS fle
    WHERE fle.cum_target==fle.mrna_count;
    """)
    rows = cursor.fetchall()
    db.close()
    return rows


def query_flexible_events(db_path, timeout_db) -> list:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("""
    SELECT *
    FROM Full_length_events_cumulative_counts AS fle
    WHERE (fle.cum_both > 0
    AND (fle.cum_query > 0 OR fle.cum_target > 0));
    """)
    rows = cursor.fetchall()
    db.close()
    return rows


def query_unused_events(db_path, timeout_db) -> list:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("""
    SELECT *
    FROM Full_length_events_cumulative_counts AS fle
    WHERE (fle.cum_neither=fle.mrna_count);
    """)
    rows = cursor.fetchall()
    db.close()
    return rows


def query_optional_events(db_path, timeout_db) -> list:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("""
    SELECT *
    FROM Full_length_events_cumulative_counts AS fle
    WHERE fle.cum_neither>0 AND (fle.cum_query>0 OR fle.cum_target>0 OR fle.cum_both>0);
    """)
    rows = cursor.fetchall()
    db.close()
    return rows


def query_obligate_pairs(db_path, timeout_db) -> list:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("""
    SELECT oe.* 
    FROM Obligatory_events AS oe
    INNER JOIN Full_length_events_cumulative_counts AS fle ON fle.fragment_id=oe.fragment_id
    WHERE fle.cum_both=fle.mrna_count
    """)
    rows = cursor.fetchall()
    db.close()
    return rows


def query_fragments(db_path, timeout_db) -> list:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("""
    SELECT 
    f.fragment_id,
    f.gene_id,
    g.gene_start,
    g.gene_end,
    g.gene_chrom,
    f.CDS_start,
    f.CDS_end,
    f.query_start,
    f.query_end,
    f.target_start,
    f.target_end,
    f.query_strand,
    f.target_strand,
    f.query_aln_prot_seq,
    f.target_aln_prot_seq
    FROM Fragments as f
    INNER JOIN Genes as g ON g.gene_id=f.gene_id
    """)
    rows = cursor.fetchall()
    db.close()
    return rows


def query_candidates(db_path, timeout_db, pars) -> list:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute(
        """
        SELECT 
        f.fragment_id,
        f.gene_id,
        f.CDS_start - f.gene_start as CDS_start,
        f.CDS_end - f.gene_start as CDS_end,
        f.target_start,
        f.target_end
        FROM Full_length_events_cumulative_counts AS f
        WHERE 
        f.gene_id==?
        AND f.fragment_id<>?
        """, pars)
    rows = cursor.fetchall()
    db.close()
    return rows


def query_full_events(db_path, timeout_db) -> list:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    matches_q = """
    SELECT 
    f.fragment_id,
    f.gene_id,
    f.CDS_start - f.gene_start as CDS_start,
    f.CDS_end - f.gene_start as CDS_end,
    f.target_start,
    f.target_end
    FROM Full_length_events_cumulative_counts AS f
    """
    cursor.execute(matches_q)
    full_matches = cursor.fetchall()
    db.close()
    return full_matches


# ### CREATE VIEW ####
def create_filtered_full_length_events_view(db_path, timeout_db) -> None:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("""
    CREATE VIEW IF NOT EXISTS Filtered_full_length_events AS
    WITH filtered_overlapping_queries AS (
    WITH overlaping_fragments AS (
    WITH gene_id_counts AS (
    SELECT f.gene_id, f.CDS_start, f.CDS_end
    FROM Fragments AS f
    WHERE f.percent_query >= 0.9
    GROUP BY f.gene_id, f.CDS_start, f.CDS_end
    HAVING COUNT(*) > 1)
    SELECT
    f.fragment_id, f.gene_id, g.gene_start,
    g.gene_end, f.CDS_start, f.CDS_end,
    f.query_start, f.query_end,
    f.target_start, f.target_end,
    MIN(f.evalue) as evalue
    FROM Fragments as f
    JOIN Genes AS g ON g.gene_id = f.gene_id
    JOIN gene_id_counts AS gc ON gc.gene_id = f.gene_id AND gc.CDS_start = f.CDS_start AND gc.CDS_end = f.CDS_end
    WHERE f.percent_query >= 0.9
    GROUP BY f.gene_id, f.CDS_start, f.CDS_end, f.query_start, f.query_end, f.target_start, f.target_end
    ) 
    SELECT f1.*
    FROM overlaping_fragments AS f1
    LEFT JOIN overlaping_fragments AS f2
    ON f1.gene_id = f2.gene_id
    AND f1.CDS_start = f2.CDS_start
    AND f1.CDS_end = f2.CDS_end
    AND f1.fragment_id <> f2.fragment_id
    AND f1.target_start <= f2.target_end
    AND f1.target_end >= f2.target_start
    WHERE f2.fragment_id IS NULL OR f1.evalue = (
    SELECT MIN(evalue)
    FROM overlaping_fragments
    WHERE gene_id = f1.gene_id
    AND CDS_start = f2.CDS_start
    AND CDS_end = f2.CDS_end
    AND target_start <= f2.target_end
    AND target_end >= f2.target_start
    )
    GROUP BY f1.fragment_id
    ORDER BY f1.fragment_id),
    single_gene_fragments AS (
    WITH single_gene_id_counts AS(
    SELECT f.fragment_id, f.gene_id
    FROM Fragments AS f
    WHERE f.percent_query >= 0.9
    GROUP BY f.gene_id, f.CDS_start, f.CDS_end
    HAVING COUNT(*) = 1)
    SELECT
    f.fragment_id, f.gene_id, g.gene_start, g.gene_end,
    f.CDS_start, f.CDS_end, f.query_start, f.query_end,
    f.target_start, f.target_end, MIN(f.evalue) as evalue
    FROM Fragments as f
    JOIN Genes AS g ON g.gene_id = f.gene_id
    JOIN single_gene_id_counts AS gc ON gc.gene_id = f.gene_id AND gc.fragment_id = f.fragment_id
    WHERE f.percent_query >= 0.9
    GROUP BY f.gene_id, f.CDS_start, f.CDS_end, f.query_start, f.query_end, f.target_start, f.target_end
    ORDER BY f.fragment_id)
    SELECT * FROM (
    SELECT * FROM single_gene_fragments AS sgf
    UNION ALL
    SELECT * FROM filtered_overlapping_queries AS foq)
    ORDER BY fragment_id;
    """)
    db.commit()
    db.close()


def create_mrna_counts_view(db_path, timeout_db) -> None:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("""
    CREATE VIEW IF NOT EXISTS Genes_mRNA_counts AS
    SELECT
    gene_id,
    COUNT(DISTINCT mrna_id) as mrna_count
    FROM Full_length_duplications
    GROUP BY gene_id;
    """)
    db.commit()
    db.close()


def create_exclusive_pairs_view(db_path, timeout_db) -> None:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("""
    CREATE VIEW IF NOT EXISTS Exclusive_pairs AS
    SELECT 
    fm3.fragment_id,
    fm3.gene_id,
    fld.mrna_id,
    fm3.gene_start,
    fm3.gene_end,
    fm3.mrna_count,
    fms.CDS_start,
    fms.CDS_end,
    fms.query_start,
    fms.query_end,
    fms.target_start,
    fms.target_end,
    fms.query_frame,
    fms.query_strand,
    fms.target_frame,
    fms.target_strand,
    fms.score,
    fms.query_dna_seq,
    fms.alignment_len,
    fms.target_dna_seq,
    fms.query_aln_prot_seq,
    fms.target_aln_prot_seq,
    fms.query_num_stop_codons,
    fms.target_num_stop_codons,
    fms.dna_perc_identity,
    fms.prot_perc_identity,
    fld.both,
    fld.query,
    fld.target,
    fld.neither,
    fms.evalue
    FROM (
    SELECT 
    * 
    FROM Full_length_events_cumulative_counts AS fle
    WHERE (fle.cum_query==(fle.mrna_count-fle.cum_target)
    AND fle.cum_target>0
    AND fle.cum_target< fle.mrna_count)) as fm3
    LEFT JOIN Fragments as fms ON
    fm3.gene_id=fms.gene_id
    AND fm3.CDS_start=fms.CDS_start
    AND fm3.CDS_end=fms.CDS_end
    AND fm3.query_start=fms.query_start
    AND fm3.query_end=fms.query_end
    AND fm3.target_start=fms.target_start
    AND fm3.target_end=fms.target_end
    LEFT JOIN Full_length_duplications AS fld
    ON fm3.gene_id = fld.gene_id
    AND fm3.fragment_id = fld.fragment_id
    AND fm3.CDS_start = fld.CDS_start
    AND fm3.CDS_end = fld.CDS_end
    AND fm3.query_start = fld.query_start
    AND fm3.query_end = fld.query_end
    AND fm3.target_start = fld.target_start
    AND fm3.target_end = fld.target_end
    ORDER BY fm3.fragment_id, target, query;
    """)
    db.commit()
    db.close()


def sanity_check(db_path, timeout_db) -> dict:
    db = sqlite3.connect(db_path, timeout=timeout_db)
    cursor = db.cursor()
    cursor.execute("""
    SELECT SUM(obligate_pairs_count) FROM
    (SELECT
        COUNT(*) AS obligate_pairs_count
    FROM (
    select
        f.fragment_id as f_id,
        f.gene_id as g_id,
        f.mrna_id as m_id,
        g.mrna_count as m_count,
        f.query_CDS_id as q_cds_id,
        f.query_CDS_start as q_cds_start,
        f.query_CDS_end as q_cds_end,
        f.query_start as q_start,
        f.query_end as q_end,
        f.target_CDS_id as t_cds_id,
        f.target_CDS_start as t_cds_start,
        f.target_CDS_end as t_cds_end,
        f.target_start as t_start,
        f.target_end as t_end,
        f.type as t,
        fle.cum_both as cum_both
    from Obligatory_events as f
    JOIN Genes_mRNA_counts  as g ON g.gene_id = f.gene_id
    INNER JOIN Full_length_events_cumulative_counts AS fle ON fle.fragment_id==f.fragment_id
    WHERE fle.cum_both==fle.mrna_count
    GROUP BY f.fragment_id
    )
    GROUP BY t);
    """)
    obligate_pairs_count = cursor.fetchone()[0]
    cursor.execute("""
    SELECT SUM(ratio) AS MXEs_full_pairs_count FROM (SELECT
        COUNT(*)/mrna_count as ratio
    FROM Full_length_events_cumulative_counts
    WHERE mrna_count = (cum_query + cum_target)
      AND cum_query = cum_target
    GROUP BY gene_id
    HAVING (COUNT(*) % mrna_count) = 0);
    """)
    MXES_full = cursor.fetchone()[0]

    cursor.execute("""
    SELECT SUM(ratio) AS MXEs_full_pairs_count FROM (SELECT
        COUNT(*)/mrna_count as ratio
    FROM Full_length_events_cumulative_counts
    WHERE mrna_count = (cum_query + cum_target)
      AND cum_query = cum_target
    GROUP BY gene_id
    HAVING (COUNT(*) % mrna_count) = 0);
    """)
    MXES_insertion_missing_match = cursor.fetchone()[0]

    cursor.execute("""
    SELECT COUNT(DISTINCT fld.fragment_id) AS Query_only_counts
    FROM Full_length_events_cumulative_counts AS fld WHERE fld.cum_query==fld.mrna_count;
    """)
    query_only = cursor.fetchone()[0]

    cursor.execute("""SELECT COUNT(DISTINCT fld.fragment_id) AS Target_only_counts
    FROM Full_length_events_cumulative_counts AS fld WHERE fld.cum_target==fld.mrna_count;""")
    target_only = cursor.fetchone()[0]

    cursor.execute("""SELECT COUNT(DISTINCT fld.fragment_id) AS Target_only_counts
    FROM Full_length_events_cumulative_counts AS fld WHERE fld.cum_neither==fld.mrna_count;""")
    neither = cursor.fetchone()[0]

    cursor.execute("""
    SELECT SUM(group_count) FROM (
    SELECT
    cum_query || cum_target || cum_both || cum_neither AS optional_group_category, COUNT(*) AS group_count
    FROM (
    SELECT
        fragment_id,
        CASE WHEN cum_query <> 0 THEN 1 ELSE 0 END AS cum_query,
        CASE WHEN cum_target <> 0 THEN 1 ELSE 0 END AS cum_target,
        CASE WHEN cum_both <> 0 THEN 1 ELSE 0 END AS cum_both,
        CASE WHEN cum_neither <> 0 THEN 1 ELSE 0 END AS cum_neither
    FROM (
    SELECT
        fle.fragment_id,
        fle.cum_query,
        fle.cum_target,
        fle.cum_both,
        fle.cum_neither
    FROM Full_length_events_cumulative_counts AS fle
    WHERE (fle.cum_neither > 0 AND (fle.cum_query > 0 OR fle.cum_target > 0 OR fle.cum_both > 0))
    ) AS subquery
    ) AS grouped_subquery
    GROUP BY optional_group_category);
    """)
    optional = cursor.fetchone()[0]

    cursor.execute("""
    SELECT SUM(group_count) FROM (
    SELECT
        cum_query || cum_target || cum_both || cum_neither AS flexible_group_category, COUNT(*) AS group_count
    FROM (
    SELECT
        fragment_id,
        CASE WHEN cum_query <> 0 THEN 1 ELSE 0 END AS cum_query,
        CASE WHEN cum_target <> 0 THEN 1 ELSE 0 END AS cum_target,
        CASE WHEN cum_both <> 0 THEN 1 ELSE 0 END AS cum_both,
        CASE WHEN cum_neither <> 0 THEN 1 ELSE 0 END AS cum_neither
    FROM (
    SELECT
        fle.fragment_id,
        fle.cum_query,
        fle.cum_target,
        fle.cum_both,
        fle.cum_neither
    FROM Full_length_events_cumulative_counts AS fle
    WHERE (fle.cum_both > 0 AND (fle.cum_query > 0 OR fle.cum_target > 0))
    ) AS subquery
    ) AS grouped_subquery
    GROUP BY flexible_group_category);
    """)
    flexible = cursor.fetchone()[0]
    db.close()

    return {'obligate_pairs_count': obligate_pairs_count,
            'MXES_full': MXES_full,
            'MXES_insertion_missing_match': MXES_insertion_missing_match,
            'query_only': query_only,
            'target_only': target_only,
            'neither': neither,
            'optional': optional,
            'flexible': flexible}
