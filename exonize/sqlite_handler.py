# ------------------------------------------------------------------------
# This module contains the SqliteHandler class, which is used to handle the
# results database.
# ------------------------------------------------------------------------
import sqlite3
import contextlib
from pathlib import Path
from collections import defaultdict
import portion as P


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

    def clear_results_database(
            self
    ) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute("""
            SELECT
                name,
                type
            FROM sqlite_master
            WHERE type IN ('table', 'view') AND name NOT LIKE 'sqlite_%';
            """)
            items = cursor.fetchall()
            # Drop each table and view except 'Genes'
            for name, type_ in items:
                if name not in ['Genes', 'Matches']:
                    cursor.execute(f"DROP {type_} IF EXISTS {name};")

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
                    GeneID VARCHAR(100) PRIMARY KEY,
                    GeneChrom VARCHAR(100) NOT NULL,
                    GeneStrand VARCHAR(1) NOT NULL,
                    TranscriptCount INTEGER NOT NULL,
                    GeneStart INTEGER NOT NULL,
                    GeneEnd INTEGER NOT NULL,
                    Duplication BINARY(1) DEFAULT 0
                );
                """
            )
            cursor.execute("""CREATE INDEX IF NOT EXISTS Genes_idx ON Genes (GeneID);""")
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS Matches (
                    FragmentID INTEGER PRIMARY KEY AUTOINCREMENT,
                    GeneID  VARCHAR(100) NOT NULL,
                    QueryExonStart INTEGER NOT NULL,
                    QueryExonEnd INTEGER NOT NULL,
                    QueryExonFrame VARCHAR NOT NULL,
                    QueryFrame INTEGER NOT NULL,
                    QueryStrand VARCHAR(1) NOT NULL,
                    TargetFrame INTEGER NOT NULL,
                    TargetStrand VARCHAR(1) NOT NULL,
                    Score INTEGER NOT NULL,
                    Bits INTEGER NOT NULL,
                    Evalue REAL NOT NULL,
                    AlignmentLength INTEGER NOT NULL,
                    QueryStart INTEGER NOT NULL,
                    QueryEnd INTEGER NOT NULL,
                    TargetStart INTEGER NOT NULL,
                    TargetEnd INTEGER NOT NULL,
                    QueryAlnProtSeq VARCHAR NOT NULL,
                    TargetAlnProtSeq VARCHAR NOT NULL,
                    Match VARCHAR NOT NULL,
                    QueryCountStopCodons INTEGER NOT NULL,
                    TargetCountStopCodons INTEGER NOT NULL,
                    FOREIGN KEY (GeneID) REFERENCES Genes(GeneID)
                );
            """
            )
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS Matches_interdependence_classification (
                    MatchID INTEGER NOT NULL,
                    GeneID VARCHAR(100) NOT NULL,
                    TranscriptID VARCHAR(100) NOT NULL,
                    QueryExonStart INTEGER NOT NULL,
                    QueryExonEnd INTEGER NOT NULL,
                    QueryExonID VARCHAR(100) NOT NULL,
                    Mode VARCHAR(100) NOT NULL,
                    TargetID VARCHAR(100), /* TRUNCTATION envents will take NULL value */
                    TargetIDStart INTEGER,
                    TargetIDEnd INTEGER,
                    TargetStart INTEGER NOT NULL,
                    TargetEnd INTEGER NOT NULL,
                    Neither BINARY(1) NOT NULL,
                    Query BINARY(1) NOT NULL,
                    Target BINARY(1) NOT NULL,
                    Both BINARY(1) NOT NULL,
                    FOREIGN KEY (MatchID) REFERENCES Expansions(MatchID),
                    FOREIGN KEY (GeneID) REFERENCES Genes(GeneID),
                    PRIMARY KEY (
                        MatchID,
                        TranscriptID,
                        QueryExonStart,
                        QueryExonEnd
                        )
                );
                """
            )
            cursor.execute(
                """
                CREATE INDEX IF NOT EXISTS Matches_interdependence_classification_idx
                ON Matches_interdependence_classification (match_id, gene_id, transcript_id);
                """
            )
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS Expansions (
                    match_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    gene_id VARCHAR(100),
                    mode VARCHAR(100),
                    start INTEGER NOT NULL,
                    end INTEGER NOT NULL,
                    degree INTEGER NOT NULL,
                    cluster_id INTEGER,
                    expansion_id INTEGER NOT NULL,
                    FOREIGN KEY (gene_id) REFERENCES Genes(gene_id),
                    UNIQUE (
                        gene_id,
                        start,
                        end,
                        expansion_id
                        )
                        );
                """
            )

    def create_matches_interdependence_expansions_counts_table(self):
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS Matches_interdependence_expansions_counts (
                    gene_id VARCHAR(100),
                    expansion_id INTEGER NOT NULL,
                    transcript_count INTEGER NOT NULL,
                    transcript_id VARCHAR(100) NOT NULL,
                    coding_event_count INTEGER NOT NULL,
                    total_coding_event_count INTEGER NOT NULL,
                    missing_coding_event VARCHAR(100) NOT NULL,
                    FOREIGN KEY (gene_id) REFERENCES Genes(gene_id),
                    UNIQUE (
                        gene_id,
                        expansion_id,
                        transcript_id
                    )
                );
                """
            )

    def create_matches_interdependence_counts_table(self) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS Matches_interdependence_counts AS
                    SELECT * FROM (
                        SELECT
                            ROW_NUMBER() OVER (ORDER BY fld.gene_id, e.expansion_id, fld.match_id) AS match_trans_id,
                            fld.match_id,
                            fld.gene_id,
                            e.expansion_id,
                            g.gene_start,
                            g.transcript_count,
                            fld.cds_start,
                            fld.cds_end,
                            fld.target_start,
                            fld.target_end,
                            group_concat(fld.mode) as concat_mode,
                            SUM(fld.both) AS cum_both,
                            SUM(fld.query) AS cum_query,
                            SUM(fld.target) AS cum_target,
                            SUM(fld.neither) AS cum_neither
                        FROM Matches_interdependence_classification AS fld
                        INNER JOIN Genes AS g ON g.gene_id = fld.gene_id
                        INNER JOIN Expansions AS e ON e.gene_id = fld.gene_id AND e.match_id = fld.match_id
                        GROUP BY fld.match_id, fld.cds_start, fld.cds_end
                        ORDER BY fld.gene_id, e.expansion_id
                    );
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
                CREATE TABLE IF NOT EXISTS Matches_full_length AS
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
                        f.evalue,
                        f.dna_identity,
                        f.prot_identity,
                        f.target_aln_prot_seq
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
                    gene_id,
                    fragment_id,
                    cds_start,
                    cds_end,
                    query_start,
                    query_end,
                    target_start,
                    target_end
                ;
            """
            )
        cursor.execute("""CREATE INDEX IF NOT EXISTS Matches_full_length_idx ON Matches_full_length (fragment_id);""")
        cursor.execute("""ALTER TABLE Matches_full_length ADD COLUMN corrected_target_start INT;""")
        cursor.execute("""ALTER TABLE Matches_full_length ADD COLUMN corrected_target_end INT;""")
        cursor.execute("""ALTER TABLE Matches_full_length ADD COLUMN corrected_dna_identity REAL;""")
        cursor.execute("""ALTER TABLE Matches_full_length ADD COLUMN corrected_prot_identity REAL;""")
        cursor.execute("""ALTER TABLE Matches_full_length ADD COLUMN corrected_target_frame REAL;""")
        cursor.execute("""ALTER TABLE Matches_full_length ADD COLUMN corrected_query_frame REAL;""")
        cursor.execute("""ALTER TABLE Matches_full_length DROP COLUMN target_aln_prot_seq;""")

    def insert_corrected_target_start_end(
            self,
            list_tuples: list
    ) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.executemany(
                """
                UPDATE Matches_full_length
                SET corrected_target_start=?,
                corrected_target_end=?,
                corrected_dna_identity=?,
                corrected_prot_identity=?,
                corrected_target_frame=?,
                corrected_query_frame=?
                WHERE fragment_id=?
                """,
                list_tuples,
            )

    def insert_identity_and_dna_algns_columns(self, list_tuples: list) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute(
                """ ALTER TABLE Matches ADD COLUMN dna_identity REAL;"""
            )
            cursor.execute(
                """ ALTER TABLE Matches ADD COLUMN prot_identity REAL;"""
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
                dna_identity=?,
                prot_identity=?,
                query_dna_seq=?,
                target_dna_seq=?
            WHERE fragment_id=?
            """,
                list_tuples,
            )

    def update_has_duplicate_genes_table(self, list_tuples: list) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.executemany(
                """
                UPDATE Genes
                SET has_duplicated_cds=1
                WHERE gene_id=?
                """,
                list_tuples,
            )

    def insert_event_categ_matches_interdependence_counts(
        self,
        list_tuples: list,
    ) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute(
                """ ALTER TABLE Matches_interdependence_counts ADD COLUMN mode VARCHAR(100);"""
            )
            cursor.executemany(
                """ UPDATE Matches_interdependence_counts SET mode=? WHERE match_trans_id=? """,
                list_tuples,
            )

    def insert_classification_column_to_matches_interdependence_counts_table(
        self, list_tuples: list
    ) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute(
                """ALTER TABLE Matches_interdependence_counts ADD COLUMN classification VARCHAR(100);"""
            )
            cursor.executemany(
                """
                UPDATE Matches_interdependence_counts
                SET classification=?
                WHERE match_trans_id=?
                """,
                list_tuples,
            )

    def insert_percent_query_column_to_fragments(self) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            if not self.check_if_column_in_table_exists(
                table_name="Matches", column_name="percent_query"
            ):
                cursor = db.cursor()
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
                            fragment_id,
                            MAX(f.cds_start, (f.query_start + f.cds_start)) AS intersect_start,
                            MIN(f.cds_end, (f.query_end + f.cds_start)) AS intersect_end,
                            f.cds_end,
                            f.cds_start
                        FROM Matches AS f
                        WHERE f.cds_end >= (f.query_start + f.cds_start)
                        AND f.cds_start <= (f.query_end + f.cds_start)
                    ) AS int
                    WHERE Matches.fragment_id = int.fragment_id;
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
                gene_end
            )
            VALUES (?, ?, ?, ?, ?, ?)
            """
            cursor.execute(insert_gene_table_param, gene_args_tuple)

    def insert_matches(
        self,
        gene_args_tuple: tuple,
        fragments_tuples_list: list,
    ) -> None:
        insert_matches_table_param = """
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
            gene_end
        )
        VALUES (?, ?, ?, ?, ?, ?)
        """
        with contextlib.closing(
            sqlite3.connect(self.results_database_path, timeout=self.timeout_database)
        ) as db:
            with db:
                with contextlib.closing(db.cursor()) as cursor:
                    cursor.execute(insert_gene_table_param, gene_args_tuple)
                    cursor.executemany(
                        insert_matches_table_param, fragments_tuples_list
                    )

    def insert_expansion_table(self, list_tuples: list) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute(
                """ SELECT
                gene_id,
                mode,
                start,
                end,
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
                mode,
                start,
                end,
                degree,
                cluster_id,
                expansion_id
            )
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """
            cursor.executemany(insert_gene_table_param, list_tuples)

    def drop_table(
            self,
            table_name: str
    ):
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute(f"""DROP TABLE IF EXISTS {table_name};""")

    def create_non_reciprocal_fragments_table(self) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute("""
            SELECT sql FROM sqlite_master WHERE type='table' AND name='Matches_full_length';
            """)
            schema = cursor.fetchone()[0]
            new_table_schema = schema.replace(
                "Matches_full_length",
                "Matches_full_length_non_reciprocal"
            )
            cursor.execute(new_table_schema)

    def insert_in_non_reciprocal_fragments_table(
            self,
            fragment_ids_list: list
    ) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            placeholders = ', '.join(['?'] * len(fragment_ids_list))
            query = f"""
                    INSERT INTO Matches_full_length_non_reciprocal
                    SELECT * FROM Matches_full_length
                    WHERE fragment_id IN ({placeholders});
                    """
            cursor.execute(query, fragment_ids_list)

    def insert_matches_interdependence_classification(self, tuples_list: list) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            insert_full_length_event_table_param = """
            INSERT INTO Matches_interdependence_classification (
                match_id,
                gene_id,
                transcript_id,
                cds_start,
                cds_end,
                query_id,
                mode,
                target_id,
                annot_target_start,
                annot_target_end,
                target_start,
                target_end,
                neither,
                query,
                target,
                both
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """
            cursor.executemany(insert_full_length_event_table_param, tuples_list)

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

    def query_concat_categ_pairs(self) -> list:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            if self.check_if_column_in_table_exists(
                table_name="Matches_interdependence_counts", column_name="concat_mode"
            ):
                cursor = db.cursor()
                cursor.execute(
                    """
                SELECT
                    match_trans_id,
                    concat_mode
                FROM Matches_interdependence_counts;
                """
                )
                records = cursor.fetchall()
                cursor.execute("""ALTER TABLE Matches_interdependence_counts DROP COLUMN concat_mode;""")
                return records

    def query_interdependence_counts_matches(
            self
    ) -> list:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute(
                """
                SELECT
                    match_trans_id,
                    gene_id,
                    transcript_count,
                    cum_both,
                    cum_query,
                    cum_target,
                    cum_neither
                    FROM Matches_interdependence_counts
                    ORDER BY gene_id
                """
            )
            return cursor.fetchall()

    def insert_matches_interdependence_expansions_counts(self, tuples_list: list) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.executemany(
                """
                INSERT INTO Matches_interdependence_expansions_counts (
                    gene_id,
                    expansion_id,
                    transcript_count,
                    transcript_id,
                    coding_event_count,
                    total_coding_event_count,
                    missing_coding_event
                )
                VALUES (?, ?, ?, ?, ?, ?, ?)
                """,
                tuples_list,
            )

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
                    f.evalue
                FROM Matches_full_length AS f
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
                    f.evalue
                FROM Matches_full_length AS f
                WHERE f.gene_id='{gene_id}'
                """
            cursor.execute(matches_q)
            return cursor.fetchall()

    def query_expansion_events(
        self,
    ) -> defaultdict:
        expansions_gene_dictionary = defaultdict(lambda: defaultdict(list))
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            fragments_query = """
            SELECT
                mnr.gene_id,
                e.match_id,
                g.gene_start,
                mnr.cds_start,
                mnr.cds_end,
                (COALESCE(mnr.corrected_target_start, mnr.target_start)) AS merged_start,
                (COALESCE(mnr.corrected_target_end, mnr.target_end)) AS merged_end,
                (COALESCE(mnr.corrected_target_start, mnr.target_start) + g.gene_start) AS merged_s,
                (COALESCE(mnr.corrected_target_end, mnr.target_end)+ g.gene_start) AS merged_e,
                e.mode AS mode,
                e.expansion_id
            FROM Matches_full_length_non_reciprocal  as mnr
            INNER JOIN Genes AS g on g.gene_id=mnr.gene_id
            INNER JOIN Expansions AS e
                ON e.gene_id = mnr.gene_id
                AND e.start = merged_start
                and e.end = merged_end
                WHERE mode == "FULL" OR mode == "INSERTION_EXCISION";
            """
            cursor.execute(fragments_query)
            records = cursor.fetchall()
            for record in records:
                (gene_id, match_id, _, cds_start, cds_end,
                 _, _, corr_target_start, corr_target_end,
                 mode, expansion_id) = record
                expansions_gene_dictionary[gene_id][expansion_id].append(
                    (match_id, gene_id, cds_start, cds_end, corr_target_start, corr_target_end)
                )
            return expansions_gene_dictionary

    def query_full_expansion_events(
            self,
    ) -> defaultdict:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute(
                """
            SELECT
                e.gene_id,
                e.start + g.gene_start as cds_start,
                e.end + g.gene_start as cds_end,
                e.expansion_id
            FROM Expansions AS e
            INNER JOIN Genes AS g ON g.gene_id=e.gene_id
            WHERE e.mode="FULL"
            ORDER BY
                e.gene_id, e.expansion_id;
            """
            )
            records = cursor.fetchall()
            expansion_events_dict = defaultdict(lambda: defaultdict(list))
            for record in records:
                gene_id, cds_start, cds_end, cluster_id, expansion_id = record
                expansion_events_dict[gene_id][expansion_id].append(
                    P.open(cds_start, cds_end)
                )
            return expansion_events_dict

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

    def query_genes_with_duplicated_cds(
        self,
    ) -> list:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute("SELECT DISTINCT gene_id FROM Expansions")
            return [record for record in cursor.fetchall()]
