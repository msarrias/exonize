# ------------------------------------------------------------------------
# This module contains the SqliteHandler class, which is used to handle the
# results database.
# ------------------------------------------------------------------------
import sqlite3
import contextlib
from pathlib import Path


class SqliteHandler(object):
    def __init__(
        self,
        results_database_path: Path,
        timeout_database: int,
    ):
        self.results_database_path = results_database_path
        self.timeout_database = timeout_database

    def check_if_table_exists(
        self,
        table_name: str,
    ) -> bool:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute("""SELECT name FROM sqlite_master WHERE type='table';""")
            return any(
                table_name == other_table_name[0]
                for other_table_name in cursor.fetchall()
            )

    def check_if_column_in_table_exists(
        self,
        table_name: str,
        column_name: str,
    ) -> bool:
        if self.check_if_table_exists(table_name=table_name):
            with sqlite3.connect(
                self.results_database_path, timeout=self.timeout_database
            ) as db:
                cursor = db.cursor()
                cursor.execute(f"PRAGMA table_info({table_name})")
                return any(
                    column_name == other_column[1] for other_column in cursor.fetchall()
                )
        return False

    @staticmethod
    def check_if_empty_table(
        db_path: str,
        table_name: str,
    ) -> bool:
        with sqlite3.connect(db_path) as db:
            cursor = db.cursor()
            cursor.execute(
                f"""
                SELECT COUNT(*) FROM {table_name}
                """
            )
            return cursor.fetchone()[0] == 0

    def connect_create_results_database(
        self,
    ) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS Genes (
                    gene_id VARCHAR(100) PRIMARY KEY,
                    gene_chrom VARCHAR(100) NOT NULL,
                    gene_strand VARCHAR(1) NOT NULL,
                    transcript_count INTEGER NOT NULL,
                    gene_start INTEGER NOT NULL,
                    gene_end INTEGER NOT NULL,
                    has_duplicated_cds BINARY(1) NOT NULL
                );
                """
            )
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS Matches (
                    fragment_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    gene_id  VARCHAR(100) NOT NULL,
                    cds_start INTEGER NOT NULL,
                    cds_end INTEGER NOT NULL,
                    cds_frame VARCHAR NOT NULL,
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
                    UNIQUE (
                        fragment_id,
                        gene_id,
                        cds_start,
                        cds_end,
                        query_start,
                        query_end,
                        target_start,
                        target_end
                    )
                );
            """
            )
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS Matches_transcript_classification (
                    fragment_match_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    fragment_id INTEGER NOT NULL,
                    gene_id VARCHAR(100) NOT NULL,
                    transcript_id VARCHAR(100) NOT NULL,
                    cds_start INTEGER NOT NULL,
                    cds_end INTEGER NOT NULL,
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
                    FOREIGN KEY (fragment_id) REFERENCES Matches(fragment_id),
                    FOREIGN KEY (gene_id) REFERENCES Genes(gene_id),
                    UNIQUE (
                        fragment_id,
                        gene_id,
                        transcript_id,
                        cds_start,
                        cds_end,
                        query_start,
                        query_end,
                        target_start,
                        target_end
                    )
                );
                """
            )
            cursor.execute(
                """
                CREATE INDEX IF NOT EXISTS Matches_transcript_classification_idx
                ON Matches_transcript_classification (fragment_id, gene_id, transcript_id, cds_start, cds_end);
                """
            )
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS Obligate_events (
                    event_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    fragment_id INTEGER NOT NULL,
                    gene_id VARCHAR(100) NOT NULL,
                    transcript_id VARCHAR(100) NOT NULL,
                    transcript_start INTEGER NOT NULL,
                    transcript_end INTEGER NOT NULL,
                    query_cds_id VARCHAR(100) NOT NULL,
                    query_cds_frame INTEGER NOT NULL,
                    query_cds_start INTEGER NOT NULL,
                    query_cds_end INTEGER NOT NULL,
                    query_start INTEGER NOT NULL,
                    query_end INTEGER NOT NULL,
                    target_cds_id VARCHAR(100) NOT NULL,
                    target_cds_frame INTEGER NOT NULL,
                    target_cds_start INTEGER NOT NULL,
                    target_cds_end INTEGER NOT NULL,
                    target_start INTEGER NOT NULL,
                    target_end INTEGER NOT NULL,
                    type VARCHAR(100) NOT NULL,
                    FOREIGN KEY (fragment_id) REFERENCES Matches(fragment_id),
                    FOREIGN KEY (gene_id) REFERENCES Genes(gene_id),
                    UNIQUE (
                        fragment_id,
                        gene_id,
                        transcript_id,
                        query_cds_id,
                        target_cds_id
                    )
                );
                """
            )
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS Truncation_events (
                    event_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    fragment_id INTEGER NOT NULL,
                    gene_id VARCHAR(100) NOT NULL,
                    transcript_id VARCHAR(100) NOT NULL,
                    transcript_start INTEGER NOT NULL,
                    transcript_end INTEGER NOT NULL,
                    query_cds_id VARCHAR(100) NOT NULL,
                    query_cds_start INTEGER NOT NULL,
                    query_cds_end INTEGER NOT NULL,
                    query_start INTEGER NOT NULL,
                    query_end INTEGER NOT NULL,
                    target_start INTEGER NOT NULL,
                    target_end INTEGER NOT NULL,
                    id_overlap_annot VARCHAR(100) NOT NULL,
                    type_overlap_annot VARCHAR(100) NOT NULL,
                    overlap_annot_start INTEGER NOT NULL,
                    overlap_annot_end INTEGER NOT NULL,
                    target_overlap_annot_start INTEGER NOT NULL,
                    target_overlap_annot_end INTEGER NOT NULL,
                    FOREIGN KEY (fragment_id) REFERENCES Matches(fragment_id),
                    FOREIGN KEY (gene_id) REFERENCES Genes(gene_id),
                    UNIQUE (
                        fragment_id,
                        gene_id,
                        transcript_id,
                        query_cds_id,
                        id_overlap_annot
                    )
                );
                """
            )
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS Expansions (
                    event_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    gene_id VARCHAR(100),
                    description VARCHAR(100),
                    ref_start INTEGER NOT NULL,
                    ref_end INTEGER NOT NULL,
                    degree INTEGER NOT NULL,
                    cluster_id INTEGER,
                    expansion_id INTEGER NOT NULL,
                    FOREIGN KEY (gene_id) REFERENCES Genes(gene_id),
                    UNIQUE (
                        gene_id,
                        ref_start,
                        ref_end,
                        expansion_id
                    )
                );
                """
            )

    def create_protein_table(self, database_path) -> None:
        with sqlite3.connect(database_path, timeout=self.timeout_database) as db:
            cursor = db.cursor()
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS Gene_Transcript_Proteins (
                    transcript_id VARCHAR(100) PRIMARY KEY,
                    gene_id VARCHAR(100) NOT NULL,
                    gene_chrom VARCHAR(100) NOT NULL,
                    gene_strand VARCHAR(1) NOT NULL,
                    gene_start INTEGER NOT NULL,
                    gene_end INTEGER NOT NULL,
                    transcript_start INTEGER NOT NULL,
                    transcript_end INTEGER NOT NULL,
                    prot_seq VARCHAR NOT NULL
                );
                """
            )
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS Gene_CDS_Proteins (
                    cds_id VARCHAR(100) PRIMARY KEY,
                    gene_id VARCHAR(100) NOT NULL,
                    transcript_id VARCHAR(100) NOT NULL,
                    rank INTEGER NOT NULL,
                    cds_frame INTEGER NOT NULL,
                    cds_start INTEGER NOT NULL,
                    cds_end INTEGER NOT NULL,
                    cds_dna_seq VARCHAR NOT NULL,
                    cds_prot_seq VARCHAR NOT NULL
                    );
                """
            )

    def create_cumulative_counts_table(self) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS Transcript_classification_counts AS
                    SELECT * FROM (
                        SELECT
                            fn.fragment_id,
                            fn.gene_id,
                            fn.gene_start,
                            fn.gene_end,
                            fn.transcript_count,
                            fn.cds_start,
                            fn.cds_end,
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
                                g.transcript_count,
                                fm.cds_start,
                                fm.cds_end,
                                fm.query_start,
                                fm.query_end,
                                fm.target_start,
                                fm.target_end,
                                fm.evalue
                            FROM Matches as fm
                            JOIN Genes AS g ON fm.gene_id = g.gene_id
                            ORDER BY fm.fragment_id
                        ) AS fn
                        JOIN Matches_transcript_classification AS fld ON fn.gene_id = fld.gene_id
                        AND fn.fragment_id = fld.fragment_id
                        AND fn.cds_start = fld.cds_start
                        AND fn.cds_end = fld.cds_end
                        GROUP BY fn.fragment_id
                        ORDER BY fn.fragment_id
                    )
                    AS fn2;
                """
            )

    def create_exclusive_pairs_view(self) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute(
                """
                CREATE VIEW IF NOT EXISTS Exclusive_pairs AS
                    SELECT
                        fm3.fragment_id,
                        fm3.gene_id,
                        fld.transcript_id,
                        fm3.gene_start,
                        fm3.gene_end,
                        fm3.transcript_count,
                        fms.cds_start,
                        fms.cds_end,
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
                        FROM Transcript_classification_counts AS fle
                        WHERE (fle.cum_query==(fle.transcript_count-fle.cum_target)
                        AND fle.cum_target>0
                        AND fle.cum_target< fle.transcript_count)
                    ) as fm3
                    LEFT JOIN Matches as fms ON fm3.gene_id=fms.gene_id
                    AND fm3.cds_start=fms.cds_start
                    AND fm3.cds_end=fms.cds_end
                    AND fm3.query_start=fms.query_start
                    AND fm3.query_end=fms.query_end
                    AND fm3.target_start=fms.target_start
                    AND fm3.target_end=fms.target_end
                    LEFT JOIN Matches_transcript_classification AS fld ON fm3.gene_id = fld.gene_id
                    AND fm3.fragment_id = fld.fragment_id
                    AND fm3.cds_start = fld.cds_start
                    AND fm3.cds_end = fld.cds_end
                    AND fm3.query_start = fld.query_start
                    AND fm3.query_end = fld.query_end
                    AND fm3.target_start = fld.target_start
                    AND fm3.target_end = fld.target_end
                    ORDER BY fm3.fragment_id, target, query;
                """
            )

    def create_filtered_full_length_events_view(
        self,
        query_overlap_threshold: float,
    ) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute(
                f"""
                CREATE VIEW IF NOT EXISTS Full_length_matches AS
                WITH
                -- Identify candidate fragments that satisfy our coverage and in-frame criteria.
                in_frame_candidate_fragments AS (
                    SELECT
                        f.fragment_id,
                        f.gene_id,
                        g.gene_start,
                        g.gene_end,
                        f.cds_frame,
                        f.query_frame,
                        f.target_frame,
                        g.gene_strand,
                        f.query_strand,
                        f.target_strand,
                        f.cds_start,
                        f.cds_end,
                        f.query_start,
                        f.query_end,
                        f.target_start,
                        f.target_end,
                        f.target_aln_prot_seq,
                        f.evalue
                    FROM Matches AS f
                    JOIN Genes g ON g.gene_id = f.gene_id
                    WHERE f.percent_query >= {query_overlap_threshold}
                    AND g.gene_strand = f.target_strand
                    AND f.cds_frame = f.query_frame
                    ),
                    -- Identify gene_ids+cdss with more than one dupl fragment
                    multi_fragment_genes AS (
                    SELECT
                        cf.gene_id,
                        cf.cds_start,
                        cf.cds_end
                    FROM in_frame_candidate_fragments AS cf
                    GROUP BY cf.gene_id, cf.cds_start, cf.cds_end
                    HAVING COUNT(*) > 1
                ),
                    -- Handling multiple fragment genes
                    -- Fragments from genes with more than one fragment are kept
                    overlapping_fragments AS (
                        SELECT
                            cf.*
                        FROM in_frame_candidate_fragments AS cf
                        -- Joining with Genes and filtered gene_ids
                        JOIN multi_fragment_genes AS mfg ON mfg.gene_id = cf.gene_id
                        AND mfg.cds_start = cf.cds_start
                        AND mfg.cds_end = cf.cds_end
                    ),
                    filtered_overlapping_fragments AS (
                        SELECT
                            DISTINCT f1.*
                        FROM overlapping_fragments AS f1
                        LEFT JOIN overlapping_fragments AS f2 ON f1.gene_id = f2.gene_id
                        AND f1.cds_start = f2.cds_start
                        AND f1.cds_end = f2.cds_end
                        AND f1.fragment_id <> f2.fragment_id
                        AND f1.target_start <= f2.target_end
                        AND f1.target_end >= f2.target_start
                        -- If step 2 works, introduce step 3, then 4 and so on.
                        WHERE f2.fragment_id IS NULL -- Keep f1 because it lacks an overlapping fragment
                        OR f1.fragment_id = (
                            SELECT
                                fragment_id
                            FROM overlapping_fragments AS ofr
                            WHERE ofr.gene_id = f1.gene_id
                            AND ofr.cds_start = f2.cds_start
                            AND ofr.cds_end = f2.cds_end
                            AND ofr.target_start <= f2.target_end
                            AND ofr.target_end >= f2.target_start
                            ORDER BY
                                CASE WHEN target_aln_prot_seq NOT LIKE '%*%' THEN 1 ELSE 2 END,
                                evalue
                                LIMIT 1
                            )
                        ORDER BY f1.fragment_id
                    ),
                    -- Identify gene_ids+cdss with exactly one dupl fragment
                    single_fragment_genes AS (
                        SELECT
                            *
                        FROM in_frame_candidate_fragments AS cf
                        GROUP BY cf.gene_id, cf.cds_start, cf.cds_end
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
                SELECT
                    *
                FROM single_gene_fragments
                UNION ALL
                SELECT
                    *
                FROM filtered_overlapping_fragments
                ORDER BY
                    fragment_id,
                    gene_id,
                    cds_start,
                    cds_end,
                    query_start,
                    query_end,
                    target_start,
                    target_end
                ;
            """
            )

    def insert_identity_and_dna_algns_columns(self, list_tuples: list) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            if self.check_if_column_in_table_exists(
                table_name="Matches",
                column_name="dna_perc_identity",
            ):
                cursor.execute(
                    """ DROP VIEW IF EXISTS Exclusive_pairs;"""
                )
                cursor.execute(
                    """ ALTER TABLE Matches DROP COLUMN dna_perc_identity;"""
                )
            if self.check_if_column_in_table_exists(
                table_name="Matches",
                column_name="prot_perc_identity",
            ):
                cursor.execute(
                    """ ALTER TABLE Matches DROP COLUMN prot_perc_identity;"""
                )
            if self.check_if_column_in_table_exists(
                table_name="Matches",
                column_name="query_dna_seq",
            ):
                cursor.execute(
                    """ ALTER TABLE Matches DROP COLUMN query_dna_seq;"""
                )
            if self.check_if_column_in_table_exists(
                table_name="Matches",
                column_name="target_dna_seq",
            ):
                cursor.execute(
                    """ ALTER TABLE Matches DROP COLUMN target_dna_seq;"""
                )
            cursor.execute(
                """ ALTER TABLE Matches ADD COLUMN dna_perc_identity REAL;"""
            )
            cursor.execute(
                """ ALTER TABLE Matches ADD COLUMN prot_perc_identity REAL;"""
            )
            cursor.execute(
                """ ALTER TABLE Matches ADD COLUMN query_dna_seq VARCHAR;"""
            )
            cursor.execute(
                """ ALTER TABLE Matches ADD COLUMN target_dna_seq VARCHAR;"""
            )
            cursor.executemany(
                """
            UPDATE Matches
            SET
                dna_perc_identity=?,
                prot_perc_identity=?,
                query_dna_seq=?,
                target_dna_seq=?
            WHERE fragment_id=?
            """,
                list_tuples,
            )

        self.create_exclusive_pairs_view()

    def insert_event_categ_full_length_events_cumulative_counts(
        self,
        list_tuples: list,
    ) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            if self.check_if_column_in_table_exists(
                table_name="Transcript_classification_counts", column_name="event_type"
            ):
                cursor.execute(
                    """ ALTER TABLE Transcript_classification_counts DROP COLUMN event_type;"""
                )
            cursor.execute(
                """ ALTER TABLE Transcript_classification_counts ADD COLUMN event_type VARCHAR(100);"""
            )
            cursor.executemany(
                """ UPDATE Transcript_classification_counts SET event_type=? WHERE fragment_id=? """,
                list_tuples,
            )

    def insert_event_id_column_to_full_length_events_cumulative_counts(
        self, list_tuples: list
    ) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            if self.check_if_column_in_table_exists(
                table_name="Transcript_classification_counts", column_name="event_id"
            ):
                cursor.execute(
                    """ALTER TABLE Transcript_classification_counts DROP COLUMN event_id;"""
                )
            cursor.execute(
                """ALTER TABLE Transcript_classification_counts ADD COLUMN event_id INTEGER;"""
            )
            cursor.executemany(
                """
                UPDATE Transcript_classification_counts
                SET event_id=?
                WHERE fragment_id=?
                """,
                list_tuples,
            )

    def insert_percent_query_column_to_fragments(self) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            if self.check_if_column_in_table_exists(
                table_name="Matches",
                column_name="percent_query",
            ):
                cursor.execute(
                    """DROP VIEW IF EXISTS Full_length_matches;"""
                )
                cursor.execute(
                    """ALTER TABLE Matches DROP COLUMN percent_query;"""
                )
            cursor.execute(
                """
                ALTER TABLE Matches ADD COLUMN percent_query DECIMAL(10, 3);
                """
            )
            cursor.execute(
                """
                UPDATE Matches
                SET percent_query =
                 ROUND(
                    CAST(int.intersect_end - int.intersect_start AS REAL) /
                    CAST(int.cds_end - int.cds_start AS REAL), 3
                )
                FROM (
                    SELECT
                        Fragment_id,
                        MAX(f.cds_start, (f.query_start + f.cds_start)) AS intersect_start,
                        MIN(f.cds_end, (f.query_end + f.cds_start)) AS intersect_end,
                        f.cds_end,
                        f.cds_start
                    FROM Matches AS f
                    WHERE f.cds_end >= (f.query_start + f.cds_start)
                    AND f.cds_start <= (f.query_end + f.cds_start)
                ) AS int
                WHERE Matches.Fragment_id = int.Fragment_id;
            """
            )

    def insert_gene_ids_table(self, gene_args_tuple: tuple) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            insert_gene_table_param = """
            INSERT INTO Genes (
                gene_id,
                gene_chrom,
                gene_strand,
                transcript_count,
                gene_start,
                gene_end,
                has_duplicated_cds
            )
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """
            cursor.execute(insert_gene_table_param, gene_args_tuple)

    def insert_fragments(
        self,
        gene_args_tuple: tuple,
        fragments_tuples_list: list,
    ) -> None:
        insert_fragments_table_param = """
        INSERT INTO Matches (
            gene_id,
            cds_start,
            cds_end,
            cds_frame,
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
        INSERT INTO Genes (
            gene_id,
            gene_chrom,
            gene_strand,
            transcript_count,
            gene_start,
            gene_end,
            has_duplicated_cds
        )
        VALUES (?, ?, ?, ?, ?, ?, ?)
        """
        with contextlib.closing(
            sqlite3.connect(self.results_database_path, timeout=self.timeout_database)
        ) as db:
            with db:
                with contextlib.closing(db.cursor()) as cursor:
                    cursor.execute(insert_gene_table_param, gene_args_tuple)
                    cursor.executemany(
                        insert_fragments_table_param, fragments_tuples_list
                    )

    def insert_expansion_table(self, list_tuples: list) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute(
                """ SELECT
                gene_id,
                description,
                ref_start,
                ref_end,
                degree,
                cluster_id,
                expansion_id
                FROM Expansions;"""
            )
            records = cursor.fetchall()
            if records:
                list_tuples = [
                    record for record in list_tuples if record not in records
                ]
            insert_gene_table_param = """
            INSERT INTO Expansions (
                gene_id,
                description,
                ref_start,
                ref_end,
                degree,
                cluster_id,
                expansion_id
            )
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """
            cursor.executemany(insert_gene_table_param, list_tuples)

    def insert_full_length_event(self, tuples_list: list) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute("""SELECT * FROM Matches_transcript_classification;""")
            records = cursor.fetchall()
            if records:
                tuples_list = [
                    record
                    for record in tuples_list
                    if record not in [record[1:] for record in records]
                ]
            if tuples_list:
                insert_full_length_event_table_param = """
                INSERT INTO Matches_transcript_classification (
                    fragment_id,
                    gene_id,
                    transcript_id,
                    cds_start,
                    cds_end,
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
                    evalue
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """
                cursor.executemany(insert_full_length_event_table_param, tuples_list)

    def insert_obligate_event(self, tuples_list: list) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            insert_obl_event_table_param = """
            INSERT OR IGNORE INTO Obligate_events (
            fragment_id,
            gene_id,
            transcript_id,
            transcript_start,
            transcript_end,
            query_cds_id,
            query_cds_frame,
            query_cds_start,
            query_cds_end,
            query_start,
            query_end,
            target_cds_id,
            target_cds_frame,
            target_cds_start,
            target_cds_end,
            target_start,
            target_end,
            type)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
             """
            cursor.executemany(insert_obl_event_table_param, tuples_list)

    def insert_truncation_event(self, tuples_list: list) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            insert_trunc_event_table_param = """
            INSERT OR IGNORE INTO Truncation_events (
            fragment_id,
            gene_id,
            transcript_id,
            transcript_start,
            transcript_end,
            query_cds_id,
            query_cds_start,
            query_cds_end,
            query_start,
            query_end,
            target_start,
            target_end,
            id_overlap_annot,
            type_overlap_annot,
            overlap_annot_start,
            overlap_annot_end,
            target_overlap_annot_start,
            target_overlap_annot_end)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """
            cursor.executemany(insert_trunc_event_table_param, tuples_list)

    def insert_into_proteins_table(
        self, database_path: str, gene_args_list_tuple: list
    ) -> None:
        with sqlite3.connect(database_path, timeout=self.timeout_database) as db:
            cursor = db.cursor()
            insert_gene_table_param = """
            INSERT INTO Gene_Transcript_Proteins (
                gene_id,
                gene_chrom,
                gene_strand,
                gene_start,
                gene_end,
                transcript_id,
                transcript_start,
                transcript_end,
                prot_seq
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """
            cursor.executemany(insert_gene_table_param, gene_args_list_tuple)

    def insert_into_cdss_table(
        self, database_path: str, gene_args_tuple_list: list
    ) -> None:
        with sqlite3.connect(database_path, timeout=self.timeout_database) as db:
            cursor = db.cursor()
            insert_gene_table_param = """
            INSERT INTO Gene_CDS_Proteins (
                gene_id,
                transcript_id,
                rank,
                cds_id,
                cds_frame,
                cds_start,
                cds_end,
                cds_dna_seq,
                cds_prot_seq
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """
            cursor.executemany(insert_gene_table_param, gene_args_tuple_list)

    def query_concat_categ_pairs(self) -> list:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute(
                """
            SELECT
                fragment_id,
                concat_event_type
            FROM Transcript_classification_counts;
            """
            )
            return cursor.fetchall()

    def query_full_events(
            self,
            gene_id: str = None
    ) -> list:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            if not gene_id:
                matches_q = """
                SELECT
                    f.fragment_id,
                    f.gene_id,
                    f.cds_start - f.gene_start as cds_start,
                    f.cds_end - f.gene_start as cds_end,
                    f.target_start,
                    f.target_end,
                    f.evalue,
                    f.event_type
                FROM Transcript_classification_counts AS f
                ORDER BY
                    f.gene_id;
                """
            else:
                matches_q = f"""
                SELECT
                f.fragment_id,
                    f.gene_id,
                    f.cds_start - f.gene_start as cds_start,
                    f.cds_end - f.gene_start as cds_end,
                    f.target_start,
                    f.target_end,
                    f.evalue,
                    f.event_type
                FROM Transcript_classification_counts AS f
                WHERE f.gene_id='{gene_id}'
                """
            cursor.execute(matches_q)
            return cursor.fetchall()

    def query_filtered_full_duplication_events(
        self,
    ) -> list:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            fragments_query = """
            SELECT
                fragment_id,
                gene_id,
                gene_start,
                gene_end,
                cds_start,
                cds_end,
                query_start,
                query_end,
                target_start,
                target_end,
                evalue
            FROM Full_length_matches
            """
            cursor.execute(fragments_query)
            return [record for record in cursor.fetchall()]

    def query_fragments(
        self,
    ) -> list:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute(
                """
            SELECT
                f.fragment_id,
                f.gene_id,
                g.gene_start,
                g.gene_end,
                g.gene_chrom,
                f.cds_start,
                f.cds_end,
                f.query_start,
                f.query_end,
                f.target_start,
                f.target_end,
                f.query_strand,
                f.target_strand,
                f.query_aln_prot_seq,
                f.target_aln_prot_seq
            FROM Matches as f
            INNER JOIN Genes as g ON g.gene_id=f.gene_id
            """
            )
            return cursor.fetchall()

    def query_gene_ids_in_results_database(
        self,
    ) -> list:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute("SELECT gene_id FROM Genes")
            return [record[0] for record in cursor.fetchall()]
