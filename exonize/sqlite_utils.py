import sqlite3                                         # for working with SQLite


def check_if_table_exists(db_path: str, table_name: str, timeout_db: int) -> bool:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("""SELECT name FROM sqlite_master WHERE type='table';""")
        return any(table_name == col[0] for col in cursor.fetchall())


def check_if_column_in_table_exists(db_path: str, table_name: str, column_name: str, timeout_db: int) -> bool:
    if check_if_table_exists(db_path, table_name, timeout_db):
        with sqlite3.connect(db_path, timeout=timeout_db) as db:
            cursor = db.cursor()
            cursor.execute(f"PRAGMA table_info({table_name})")
            return any(column_name == col[1] for col in cursor.fetchall())
    return False


# #### CREATE TABLES #####
def connect_create_results_db(db_path: str, timeout_db: int) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
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
        CDS_frame VARCHAR NOT NULL,
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
        CREATE TABLE IF NOT EXISTS Obligate_events (
        event_id INTEGER PRIMARY KEY AUTOINCREMENT,
        fragment_id INTEGER NOT NULL,
        gene_id VARCHAR(100) NOT NULL,
        mrna_id VARCHAR(100) NOT NULL,
        mrna_start INTEGER NOT NULL,
        mrna_end INTEGER NOT NULL,
        query_CDS_id VARCHAR(100) NOT NULL,
        query_CDS_frame INTEGER NOT NULL,
        query_CDS_start INTEGER NOT NULL,
        query_CDS_end INTEGER NOT NULL,
        query_start INTEGER NOT NULL,
        query_end INTEGER NOT NULL,
        target_CDS_id VARCHAR(100) NOT NULL,
        target_CDS_frame INTEGER NOT NULL,
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


def create_protein_table(db_path: str, timeout_db: int) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS Proteins (
        gene_id VARCHAR(100) NOT NULL,
        gene_chrom VARCHAR(100) NOT NULL,
        gene_strand VARCHAR(1) NOT NULL,
        gene_start INTEGER NOT NULL,
        gene_end INTEGER NOT NULL,
        transcript_id VARCHAR(100) NOT NULL,
        transcript_start INTEGER NOT NULL,
        transcript_end INTEGER NOT NULL,
        prot_seq VARCHAR NOT NULL,
        primary key (gene_id, transcript_id)
        )
        """)
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS CDSs (
        gene_id VARCHAR(100) NOT NULL,
        transcript_id VARCHAR(100) NOT NULL,
        rank INTEGER NOT NULL,
        cds_id VARCHAR(100) NOT NULL,
        cds_frame INTEGER NOT NULL,
        cds_start INTEGER NOT NULL,
        cds_end INTEGER NOT NULL,
        cds_dna_seq VARCHAR NOT NULL,
        cds_prot_seq VARCHAR NOT NULL,
        primary key (gene_id, transcript_id, rank, cds_id)
        )
        """)
        db.commit()


def create_cumulative_counts_table(db_path: str, timeout_db: int) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
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


# #### INSERT INTO TABLES #####
def insert_identity_and_dna_algns_columns(db_path: str, timeout_db: int, fragments: list) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        if check_if_column_in_table_exists(db_path, 'Fragments', 'dna_perc_identity', timeout_db):
            cursor.execute(""" ALTER TABLE Fragments DROP COLUMN dna_perc_identity;""")
        if check_if_column_in_table_exists(db_path, 'Fragments', 'prot_perc_identity', timeout_db):
            cursor.execute(""" ALTER TABLE Fragments DROP COLUMN prot_perc_identity;""")
        if check_if_column_in_table_exists(db_path, 'Fragments', 'query_dna_seq', timeout_db):
            cursor.execute(""" ALTER TABLE Fragments DROP COLUMN query_dna_seq;""")
        if check_if_column_in_table_exists(db_path, 'Fragments', 'target_dna_seq', timeout_db):
            cursor.execute(""" ALTER TABLE Fragments DROP COLUMN target_dna_seq;""")
        cursor.execute(""" ALTER TABLE Fragments ADD COLUMN dna_perc_identity REAL;""")
        cursor.execute(""" ALTER TABLE Fragments ADD COLUMN prot_perc_identity REAL;""")
        cursor.execute(""" ALTER TABLE Fragments ADD COLUMN query_dna_seq VARCHAR;""")
        cursor.execute(""" ALTER TABLE Fragments ADD COLUMN target_dna_seq VARCHAR;""")
        cursor.executemany(""" 
        UPDATE Fragments 
        SET dna_perc_identity=?, prot_perc_identity=?, query_dna_seq=?, target_dna_seq=?  WHERE fragment_id=? 
        """, fragments)
        db.commit()


def insert_event_categ_full_length_events_cumulative_counts(db_path: str, timeout_db: int, fragments: list) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        if check_if_column_in_table_exists(db_path, 'Full_length_events_cumulative_counts', 'event_type', timeout_db):
            cursor.execute(""" ALTER TABLE Full_length_events_cumulative_counts DROP COLUMN event_type;""")
        cursor.execute(""" ALTER TABLE Full_length_events_cumulative_counts ADD COLUMN event_type VARCHAR(100);""")
        cursor.executemany(""" UPDATE Full_length_events_cumulative_counts SET event_type=? WHERE fragment_id=? """, fragments)
        db.commit()


def instert_pair_id_column_to_full_length_events_cumulative_counts(db_path: str, timeout_db: int, fragments: list) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        if check_if_column_in_table_exists(db_path, 'Full_length_events_cumulative_counts', 'pair_id', timeout_db):
            cursor.execute(""" ALTER TABLE Full_length_events_cumulative_counts DROP COLUMN pair_id;""")
        cursor.execute(""" ALTER TABLE Full_length_events_cumulative_counts ADD COLUMN pair_id INTEGER;""")
        cursor.executemany(""" UPDATE Full_length_events_cumulative_counts SET pair_id=? WHERE fragment_id=? """, fragments)
        db.commit()


def insert_percent_query_column_to_fragments(db_path: str, timeout_db: int) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        if check_if_column_in_table_exists(db_path, 'Fragments', 'percent_query', timeout_db):
            cursor.execute(""" ALTER TABLE Fragments DROP COLUMN percent_query;""")
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


def insert_gene_ids_table(db_path: str, timeout_db: int, gene_args_tuple: tuple) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
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


def insert_fragments_calls() -> tuple:
    insert_fragments_table_param = """
    INSERT INTO Fragments (
    gene_id,
    CDS_start,
    CDS_end,
    CDS_frame,
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
    VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
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


def instert_full_length_event(db_path: str, timeout_db: int, tuples_list: list) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
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


def instert_obligatory_event(db_path: str, timeout_db: int, tuples_list: list) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        insert_obl_event_table_param = """
        INSERT OR IGNORE INTO Obligate_events (
        fragment_id,
        gene_id,
        mrna_id,
        mrna_start,
        mrna_end,
        query_CDS_id,
        query_CDS_frame,
        query_CDS_start,
        query_CDS_end,
        query_start,
        query_end,
        target_CDS_id,
        target_CDS_frame,
        target_CDS_start,
        target_CDS_end,
        target_start,
        target_end,
        type)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
         """
        cursor.executemany(insert_obl_event_table_param, tuples_list)
        db.commit()


def instert_truncation_event(db_path: str, timeout_db: int, tuples_list: list) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
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


def insert_into_proteins(db_path: str, timeout_db: int, gene_args_list_tuple: list) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        insert_gene_table_param = """  
        INSERT INTO Proteins (
        gene_id,
        gene_chrom,
        gene_strand,
        gene_start, 
        gene_end, 
        transcript_id,
        transcript_start,
        transcript_end,
        prot_seq) 
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """
        cursor.executemany(insert_gene_table_param, gene_args_list_tuple)
        db.commit()


def insert_into_CDSs(db_path: str, timeout_db: int, gene_args_tuple_list: list) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        insert_gene_table_param = """  
        INSERT INTO CDSs (
        gene_id,
        transcript_id,
        rank,
        cds_id,
        cds_frame,
        cds_start,
        cds_end, 
        cds_dna_seq, 
        cds_prot_seq) 
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """
        cursor.executemany(insert_gene_table_param, gene_args_tuple_list)
        db.commit()


# ###### QUERY TABLES ######
def query_gene_ids_in_res_db(db_path: str, timeout_db: int) -> list:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("SELECT gene_id FROM Genes")
        return [i[0] for i in cursor.fetchall()]


def query_filtered_full_duplication_events(db_path: str, timeout_db: float) -> list:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
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
        return [i for i in cursor.fetchall()]


def query_tuples_full_length_duplications(db_path: str, timeout_db: int) -> list:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("""
        SELECT *
        FROM Full_length_duplications;
        """)
        return cursor.fetchall()


def query_genes_with_duplicated_cds(db_path: str, timeout_db: int) -> list:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("""
        SELECT gene_id 
        FROM Genes 
        WHERE has_duplicated_CDS==1
        """)
        return [i[0] for i in cursor.fetchall()]


def query_obligate_events(db_path: str, timeout_db: int) -> list:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("""
        SELECT *
        FROM Obligate_events;
        """)
        return cursor.fetchall()


def query_truncation_events(db_path: str, timeout_db: int) -> list:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("""
        SELECT *
        FROM Truncation_events;
        """)
        return cursor.fetchall()


def query_exclusive_events(db_path: str, timeout_db: int) -> list:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("""
        SELECT *
        FROM Exclusive_pairs;
        """)
        return cursor.fetchall()


def query_queries_only(db_path: str, timeout_db: int) -> list:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("""
        SELECT *
        FROM Full_length_events_cumulative_counts AS fle
        WHERE fle.cum_query==fle.mrna_count;
        """)
        return cursor.fetchall()


def query_targets_only(db_path: str, timeout_db: int) -> list:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("""
        SELECT *
        FROM Full_length_events_cumulative_counts AS fle
        WHERE fle.cum_target==fle.mrna_count;
        """)
        return cursor.fetchall()


def query_flexible_events(db_path: str, timeout_db: int) -> list:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("""
        SELECT *
        FROM Full_length_events_cumulative_counts AS fle
        WHERE (fle.cum_both > 0
        AND (fle.cum_query > 0 OR fle.cum_target > 0));
        """)
        return cursor.fetchall()


def query_unused_events(db_path: str, timeout_db: int) -> list:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("""
        SELECT *
        FROM Full_length_events_cumulative_counts AS fle
        WHERE (fle.cum_neither=fle.mrna_count);
        """)
        return cursor.fetchall()


def query_optional_events(db_path: str, timeout_db: int) -> list:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("""
        SELECT *
        FROM Full_length_events_cumulative_counts AS fle
        WHERE fle.cum_neither>0 AND (fle.cum_query>0 OR fle.cum_target>0 OR fle.cum_both>0);
        """)
        return cursor.fetchall()


def query_obligate_pairs(db_path: str, timeout_db: int) -> list:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("""
        SELECT oe.* 
        FROM Obligate_events AS oe
        INNER JOIN Full_length_events_cumulative_counts AS fle ON fle.fragment_id=oe.fragment_id
        WHERE fle.cum_both=fle.mrna_count
        """)
        return cursor.fetchall()


def query_fragments(db_path: str, timeout_db: int) -> list:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
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
        return cursor.fetchall()


def query_concat_categ_pairs(db_path: str, timeout_db: int) -> list:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("""
        SELECT  
        fragment_id,
        concat_event_type
        FROM Full_length_events_cumulative_counts 
        """)
        return cursor.fetchall()


def query_candidates(db_path: str, timeout_db: int, pars: tuple) -> list:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute(
            """
            SELECT 
            f.fragment_id,
            f.gene_id,
            f.CDS_start - f.gene_start as CDS_start,
            f.CDS_end - f.gene_start as CDS_end,
            f.target_start,
            f.target_end,
            f.event_type
            FROM Full_length_events_cumulative_counts AS f
            WHERE 
            f.gene_id==?
            AND f.fragment_id<>?
            """, pars)
        return cursor.fetchall()


def query_full_events(db_path: str, timeout_db: int) -> list:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        matches_q = """
        SELECT 
        f.fragment_id,
        f.gene_id,
        f.CDS_start - f.gene_start as CDS_start,
        f.CDS_end - f.gene_start as CDS_end,
        f.target_start,
        f.target_end,
        f.event_type
        FROM Full_length_events_cumulative_counts AS f
        ORDER BY
            f.gene_id,
            CASE WHEN f.event_type LIKE '%INS_CDS%' THEN 1 ELSE 2 END,
            f.event_type;
        """
        cursor.execute(matches_q)
        return cursor.fetchall()


# ### CREATE VIEW ####
def create_filtered_full_length_events_view(db_path: str, timeout_db: int) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
        cursor = db.cursor()
        cursor.execute("""
        CREATE VIEW IF NOT EXISTS Filtered_full_length_events AS
        WITH
        -- Identify candidate fragments that satisfy our coverage and in-frame criteria.
	    in_frame_candidate_fragments AS (
        	SELECT f.fragment_id, f.gene_id, g.gene_start, g.gene_end, f.CDS_frame, f.query_frame, 
		       g.gene_strand, f.query_strand,
		       f.CDS_start, f.CDS_end, f.query_start, f.query_end, f.target_start, f.target_end, f.evalue
 		FROM Fragments AS f
		JOIN Genes g ON g.gene_id = f.gene_id
 		WHERE f.percent_query >= 0.9 AND f.CDS_frame = f.query_frame
		),
        -- Identify gene_ids+CDSs with more than one dupl fragment
            multi_fragment_genes AS (
	        SELECT cf.gene_id, cf.CDS_start, cf.CDS_end 
        	FROM in_frame_candidate_fragments AS cf
        	GROUP BY cf.gene_id, cf.CDS_start, cf.CDS_end
        	HAVING COUNT(*) > 1
	    ),
        -- Handling multiple fragment genes
	--   Fragments from genes with more than one fragment are kept
            overlapping_fragments AS (
                SELECT cf.*
	    	FROM in_frame_candidate_fragments AS cf
		        -- Joining with Genes and filtered gene_ids
        	JOIN multi_fragment_genes AS mfg ON mfg.gene_id = cf.gene_id
                                                 AND mfg.CDS_start = cf.CDS_start
 					         AND mfg.CDS_end = cf.CDS_end
            ),
            filtered_overlapping_fragments AS (
                SELECT DISTINCT f1.*
            	FROM overlapping_fragments AS f1
            	LEFT JOIN overlapping_fragments AS f2 ON f1.gene_id = f2.gene_id
            	                               AND f1.CDS_start = f2.CDS_start
 					       AND f1.CDS_end = f2.CDS_end
                                               AND f1.fragment_id <> f2.fragment_id
                                               AND f1.target_start <= f2.target_end
                                               AND f1.target_end >= f2.target_start
                -- If step 2 works, introduce step 3, then 4 and so on.
	    	WHERE f2.fragment_id IS NULL -- Keep f1 because it lacks an overlapping fragment
            	      OR f1.fragment_id = (
            	          SELECT fragment_id
            	      	  FROM overlapping_fragments AS ofr
            	      	  WHERE ofr.gene_id = f1.gene_id
		              AND ofr.CDS_start = f2.CDS_start
		              AND ofr.CDS_end = f2.CDS_end
            		      AND ofr.target_start <= f2.target_end
            		      AND ofr.target_end >= f2.target_start
            		  ORDER BY
            		      CASE WHEN gene_strand = query_strand THEN 1 ELSE 2 END,
            		      evalue
            		  LIMIT 1
            	     )
                 ORDER BY f1.fragment_id
            ),
        -- Identify gene_ids+CDSs with exactly one dupl fragment
            single_fragment_genes AS (
                SELECT *
                FROM in_frame_candidate_fragments AS cf
                GROUP BY cf.gene_id, cf.CDS_start, cf.CDS_end
                HAVING COUNT(*) = 1
            ),
        -- Handling single fragment genes
            single_gene_fragments AS (
                SELECT cf.*
		FROM in_frame_candidate_fragments AS cf
        	JOIN single_fragment_genes sfg ON sfg.gene_id = cf.gene_id
 		     AND sfg.fragment_id = cf.fragment_id
            )
        -- Combining the results of single_gene_fragments and filtered_overlapping_fragments
        SELECT * FROM single_gene_fragments
        UNION ALL
        SELECT * FROM filtered_overlapping_fragments
        ORDER BY fragment_id, gene_id, CDS_start, CDS_end, query_start, query_end, target_start, target_end;
        """)
        db.commit()


def create_mrna_counts_view(db_path: str, timeout_db: int) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
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


def create_exclusive_pairs_view(db_path: str, timeout_db: int) -> None:
    with sqlite3.connect(db_path, timeout=timeout_db) as db:
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
