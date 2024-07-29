# ------------------------------------------------------------------------
# This module contains the SqliteHandler class, which is used to handle the
# results database.
# ------------------------------------------------------------------------
import sqlite3
import contextlib
from pathlib import Path
from collections import defaultdict
import portion as P
import pandas as pd


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

    def add_column_to_table(
        self,
        table_name: str,
        column_name: str,
        column_type: str,
    ) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            if self.check_if_table_exists(table_name=table_name):
                if self.check_if_column_in_table_exists(table_name=table_name, column_name=column_name):
                    cursor.execute(f"ALTER TABLE {table_name} DROP COLUMN {column_name};")
                cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN {column_name} {column_type};")

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

    def drop_table(
            self,
            table_name: str
    ):
        with sqlite3.connect(
                self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute(f"""DROP TABLE IF EXISTS {table_name};""")

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
                ON Matches_interdependence_classification (MatchID, GeneID, TranscriptID);
                """
            )
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS Expansions (
                    MatchID INTEGER PRIMARY KEY AUTOINCREMENT,
                    GeneID VARCHAR(100),
                    Mode VARCHAR(100),
                    EventStart INTEGER NOT NULL,
                    EventEnd INTEGER NOT NULL,
                    EventDegree INTEGER NOT NULL,
                    ClusterID INTEGER,
                    ExpansionID INTEGER NOT NULL,
                    FOREIGN KEY (GeneID) REFERENCES Genes(GeneID),
                    UNIQUE (
                        GeneID,
                        EventStart,
                        EventEnd,
                        ExpansionID
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
                            ROW_NUMBER() OVER (ORDER BY fld.GeneID, e.ExpansionID, fld.MatchID) AS ID,
                            fld.GeneID,
                            e.ExpansionID,
                            g.GeneStart,
                            g.TranscriptCount,
                            fld.QueryExonStart,
                            fld.QueryExonEnd,
                            fld.TargetStart,
                            fld.TargetEnd,
                            group_concat(fld.Mode) as ConcatMode,
                            SUM(fld.Both) AS Both,
                            SUM(fld.Query) AS Query,
                            SUM(fld.Target) AS Target,
                            SUM(fld.Neither) AS Neither
                        FROM Matches_interdependence_classification AS fld
                        INNER JOIN Genes AS g ON g.GeneID = fld.GeneID
                        INNER JOIN Expansions AS e ON e.GeneID = fld.GeneID AND e.MatchID = fld.MatchID
                        GROUP BY fld.MatchID, fld.QueryExonStart, fld.QueryExonEnd
                        ORDER BY fld.GeneID, e.ExpansionID
                    );
                """
            )

    def create_filtered_full_length_events_view(
        self,
        query_overlap_threshold: float,
        evalue_threshold: float,
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
                        f.FragmentID,
                        f.GeneID,
                        g.GeneStart,
                        g.GeneEnd,
                        f.QueryExonFrame,
                        f.QueryFrame,
                        f.TargetFrame,
                        g.GeneStrand,
                        f.QueryStrand,
                        f.TargetStrand,
                        f.QueryExonStart,
                        f.QueryExonEnd,
                        f.QueryStart,
                        f.QueryEnd,
                        f.TargetStart,
                        f.TargetEnd,
                        f.Evalue,
                        f.DNAIdentity,
                        f.ProtIdentity,
                        f.QueryAlnProtSeq,
                        f.TargetAlnProtSeq
                    FROM Matches AS f
                    JOIN Genes g ON g.GeneID = f.GeneID
                    WHERE f.AlnQuery >= {query_overlap_threshold}
                    AND f.Evalue <= {evalue_threshold}
                    AND g.GeneStrand = f.TargetStrand
                    AND f.QueryExonFrame = f.QueryFrame
                    ),
                    -- Identify gene_ids+cdss with more than one dupl fragment
                    multi_fragment_genes AS (
                    SELECT
                        cf.GeneID,
                        cf.QueryExonStart,
                        cf.QueryExonEnd
                    FROM in_frame_candidate_fragments AS cf
                    GROUP BY cf.GeneID, cf.QueryExonStart, cf.QueryExonEnd
                    HAVING COUNT(*) > 1
                ),
                    -- Handling multiple fragment genes
                    -- Fragments from genes with more than one fragment are kept
                    overlapping_fragments AS (
                        SELECT
                            cf.*
                        FROM in_frame_candidate_fragments AS cf
                        -- Joining with Genes and filtered gene_ids
                        JOIN multi_fragment_genes AS mfg ON mfg.GeneID = cf.GeneID
                        AND mfg.QueryExonStart = cf.QueryExonStart
                        AND mfg.QueryExonEnd = cf.QueryExonEnd
                    ),
                    filtered_overlapping_fragments AS (
                        SELECT
                            DISTINCT f1.*
                        FROM overlapping_fragments AS f1
                        LEFT JOIN overlapping_fragments AS f2 ON f1.GeneID = f2.GeneID
                        AND f1.QueryExonStart = f2.QueryExonStart
                        AND f1.QueryExonEnd = f2.QueryExonEnd
                        AND f1.FragmentID <> f2.FragmentID
                        AND f1.TargetStart <= f2.TargetEnd
                        AND f1.TargetEnd >= f2.TargetStart
                        -- If step 2 works, introduce step 3, then 4 and so on.
                        WHERE f2.FragmentID IS NULL -- Keep f1 because it lacks an overlapping fragment
                        OR f1.FragmentID = (
                            SELECT
                                FragmentID
                            FROM overlapping_fragments AS ofr
                            WHERE ofr.GeneID = f1.GeneID
                            AND ofr.QueryExonStart = f2.QueryExonStart
                            AND ofr.QueryExonEnd = f2.QueryExonEnd
                            AND ofr.TargetStart <= f2.TargetEnd
                            AND ofr.TargetEnd >= f2.TargetStart
                            ORDER BY
                                CASE WHEN TargetAlnProtSeq NOT LIKE '%*%' THEN 1 ELSE 2 END,
                                Evalue
                                LIMIT 1
                            )
                        ORDER BY f1.FragmentID
                    ),
                    -- Identify gene_ids+cdss with exactly one dupl fragment
                    single_fragment_genes AS (
                        SELECT
                            *
                        FROM in_frame_candidate_fragments AS cf
                        GROUP BY cf.GeneID, cf.QueryExonStart, cf.QueryExonEnd
                        HAVING COUNT(*) = 1
                    ),
                    -- Handling single fragment genes
                    single_gene_fragments AS (
                        SELECT cf.*
                    FROM in_frame_candidate_fragments AS cf
                    JOIN single_fragment_genes sfg ON sfg.GeneID = cf.GeneID
                     AND sfg.FragmentID = cf.FragmentID
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
                    GeneID,
                    FragmentID,
                    QueryExonStart,
                    QueryExonEnd,
                    QueryStart,
                    QueryEnd,
                    TargetStart,
                    TargetEnd
                ;
            """
            )
        cursor.execute("""CREATE INDEX IF NOT EXISTS Matches_full_length_idx ON Matches_full_length (FragmentID);""")
        columns_to_add = [
            ("CorrectedTargetStart", "INTEGER"),
            ("CorrectedTargetEnd", "INTEGER"),
            ("CorrectedDNAIdentity", "REAL"),
            ("CorrectedProtIdentity", "REAL"),
            ("QueryProtSeq", "VARCHAR"),
            ("CorrectedTargetProtSeq", "VARCHAR"),
            ("CorrectedTargetFrame", "INTEGER"),
            ("CorrectedQueryFrame", "INTEGER"),]
        for column_name, column_type in columns_to_add:
            self.add_column_to_table(
                table_name="Matches_full_length",
                column_name=column_name,
                column_type=column_type,
            )

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
                SET CorrectedTargetStart=?,
                CorrectedTargetEnd=?,
                CorrectedDNAIdentity=?,
                CorrectedProtIdentity=?,
                QueryProtSeq=?,
                CorrectedTargetProtSeq=?,
                CorrectedTargetFrame=?,
                CorrectedQueryFrame=?
                WHERE FragmentID=?
                """,
                list_tuples,
            )

    def insert_identity_and_dna_algns_columns(self, list_tuples: list) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute(
                """ ALTER TABLE Matches ADD COLUMN QueryDNASeq VARCHAR;"""
            )
            cursor.execute(
                """ ALTER TABLE Matches ADD COLUMN TargetDNASeq VARCHAR;"""
            )
            cursor.execute(
                """ ALTER TABLE Matches ADD COLUMN DNAIdentity REAL;"""
            )
            cursor.execute(
                """ ALTER TABLE Matches ADD COLUMN ProtIdentity REAL;"""
            )

            cursor.executemany(
                """
            UPDATE Matches
            SET
                QueryDNASeq=?,
                TargetDNASeq=?,
                DNAIdentity=?,
                ProtIdentity=?
            WHERE FragmentID=?
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
                SET Duplication=1
                WHERE GeneID=?
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
                """ ALTER TABLE Matches_interdependence_counts ADD COLUMN Mode VARCHAR(100);"""
            )
            cursor.executemany(
                """ UPDATE Matches_interdependence_counts SET Mode=? WHERE ID=? """,
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
                """ALTER TABLE Matches_interdependence_counts ADD COLUMN Class VARCHAR(100);"""
            )
            cursor.executemany(
                """
                UPDATE Matches_interdependence_counts
                SET Class=?
                WHERE ID=?
                """,
                list_tuples,
            )

    def insert_percent_query_column_to_fragments(self) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            if not self.check_if_column_in_table_exists(
                table_name="Matches", column_name="AlnQuery"
            ):
                cursor = db.cursor()
                cursor.execute(
                    """
                    ALTER TABLE Matches ADD COLUMN AlnQuery DECIMAL(10, 3);
                    """
                )
                cursor.execute(
                    """
                    UPDATE Matches
                    SET AlnQuery =
                     ROUND(
                        CAST(int.intersect_end - int.intersect_start AS REAL) /
                        CAST(int.QueryExonEnd - int.QueryExonStart AS REAL), 3
                    )
                    FROM (
                        SELECT
                            FragmentID,
                            MAX(f.QueryExonStart, (f.QueryStart + f.QueryExonStart)) AS intersect_start,
                            MIN(f.QueryExonEnd, (f.QueryEnd + f.QueryExonStart)) AS intersect_end,
                            f.QueryExonEnd,
                            f.QueryExonStart
                        FROM Matches AS f
                        WHERE f.QueryExonEnd >= (f.QueryStart + f.QueryExonStart)
                        AND f.QueryExonStart <= (f.QueryEnd + f.QueryExonStart)
                    ) AS int
                    WHERE Matches.FragmentID = int.FragmentID;
                """
                )

    def insert_gene_ids_table(self, gene_args_tuple: tuple) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            insert_gene_table_param = """
            INSERT INTO Genes (
                GeneID,
                GeneChrom,
                GeneStrand,
                TranscriptCount,
                GeneStart,
                GeneEnd
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
            GeneID,
            QueryExonStart,
            QueryExonEnd,
            QueryExonFrame,
            QueryFrame,
            QueryStrand,
            TargetFrame,
            TargetStrand,
            Score,
            Bits,
            Evalue,
            AlignmentLength,
            QueryStart,
            QueryEnd,
            TargetStart,
            TargetEnd,
            QueryAlnProtSeq,
            TargetAlnProtSeq,
            Match,
            QueryCountStopCodons,
            TargetCountStopCodons
        )
        VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """
        insert_gene_table_param = """
        INSERT INTO Genes (
            GeneID,
            GeneChrom,
            GeneStrand,
            TranscriptCount,
            GeneStart,
            GeneEnd
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
                GeneID,
                Mode,
                EventStart,
                EventEnd,
                EventDegree,
                ClusterID,
                ExpansionID
                FROM Expansions;"""
            )
            records = cursor.fetchall()
            if records:
                list_tuples = [
                    record for record in list_tuples if record not in records
                ]
            insert_gene_table_param = """
            INSERT INTO Expansions (
                GeneID,
                Mode,
                EventStart,
                EventEnd,
                EventDegree,
                ClusterID,
                ExpansionID
            )
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """
            cursor.executemany(insert_gene_table_param, list_tuples)

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
            SELECT
                FragmentID,
                GeneID,
                GeneStart,
                GeneEnd,
                QueryExonFrame,
                QueryFrame,
                TargetFrame,
                GeneStrand,
                QueryStrand,
                TargetStrand,
                QueryExonStart,
                QueryExonEnd,
                QueryStart,
                QueryEnd,
                TargetStart + GeneStart AS TargetStart,
                TargetEnd + GeneStart AS TargetEnd,
                Evalue,
                DNAIdentity,
                ProtIdentity,
                QueryAlnProtSeq,
                TargetAlnProtSeq,
                COALESCE(CorrectedTargetStart, TargetStart) + GeneStart AS CorrectedTargetStart,
                COALESCE(CorrectedTargetEnd, TargetEnd) + GeneStart AS CorrectedTargetEnd,
                COALESCE(CorrectedDNAIdentity, DNAIdentity) AS CorrectedDNAIdentity,
                COALESCE(CorrectedProtIdentity, ProtIdentity) AS CorrectedProtIdentity,
                COALESCE(QueryProtSeq, QueryAlnProtSeq) AS QueryProtSeq,
                COALESCE(CorrectedTargetProtSeq, TargetAlnProtSeq) AS CorrectedTargetProtSeq,
                COALESCE(CorrectedTargetFrame, TargetFrame) AS CorrectedTargetFrame,
                COALESCE(CorrectedQueryFrame, QueryFrame) AS CorrectedQueryFrame
            FROM Matches_full_length
            WHERE FragmentID IN ({placeholders});
                    """
            cursor.execute(query, fragment_ids_list)

    def insert_matches_interdependence_classification(self, tuples_list: list) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            insert_full_length_event_table_param = """
            INSERT INTO Matches_interdependence_classification (
                MatchID,
                GeneID,
                TranscriptID,
                QueryExonStart,
                QueryExonEnd,
                QueryExonID,
                Mode,
                TargetID,
                TargetIDStart,
                TargetIDEnd,
                TargetStart,
                TargetEnd,
                Neither,
                Query,
                Target,
                Both
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """
            cursor.executemany(insert_full_length_event_table_param, tuples_list)

    def query_concat_categ_pairs(self) -> list:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            if self.check_if_column_in_table_exists(
                table_name="Matches_interdependence_counts", column_name="ConcatMode"
            ):
                cursor = db.cursor()
                cursor.execute(
                    """
                SELECT
                    ID,
                    ConcatMode
                FROM Matches_interdependence_counts;
                """
                )
                records = cursor.fetchall()
                cursor.execute("""ALTER TABLE Matches_interdependence_counts DROP COLUMN ConcatMode;""")
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
                    ID,
                    GeneID,
                    TranscriptCount,
                    Both,
                    Query,
                    Target,
                    Neither
                    FROM Matches_interdependence_counts
                    ORDER BY GeneID
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
                    f.FragmentID,
                    f.GeneID,
                    f.QueryExonStart - f.GeneStart as QueryExonStart,
                    f.QueryExonEnd - f.GeneStart as QueryExonEnd,
                    f.TargetStart,
                    f.TargetEnd,
                    f.Evalue
                FROM Matches_full_length AS f
                ORDER BY
                    f.GeneID;
                """
            else:
                matches_q = f"""
                SELECT
                    f.FragmentID,
                    f.GeneID,
                    f.QueryExonStart - f.GeneStart as QueryExonStart,
                    f.QueryExonEnd - f.GeneStart as QueryExonEnd,
                    f.TargetStart,
                    f.TargetEnd,
                    f.Evalue
                FROM Matches_full_length AS f
                WHERE f.GeneID='{gene_id}'
                """
            cursor.execute(matches_q)
            return cursor.fetchall()

    def query_expansion_coding_events(
        self,
    ) -> defaultdict:
        expansions_gene_dictionary = defaultdict(lambda: defaultdict(list))
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            fragments_query = """
            SELECT
                n.GeneID,
                e.MatchID,
                n.QueryExonStart,
                n.QueryExonEnd,
                n.CorrectedTargetStart,
                n.CorrectedTargetEnd,
                e.Mode,
                e.ExpansionID
             FROM (
                 SELECT
                    mnr.GeneID,
                    g.GeneStart,
                    mnr.QueryExonStart,
                    mnr.QueryExonEnd,
                    mnr.CorrectedTargetStart,
                    mnr.CorrectedTargetEnd
                FROM Matches_full_length_non_reciprocal  as mnr
                INNER JOIN Genes AS g on g.GeneID=mnr.GeneID
                ) AS n
            INNER JOIN Expansions AS e
            ON n.GeneID = e.GeneID
            AND e.EventStart= n.CorrectedTargetStart
            and e.EventEnd = n.CorrectedTargetEnd
            WHERE (e.Mode == "FULL" OR e.Mode == "INSERTION_EXCISION");
            """
            cursor.execute(fragments_query)
            records = cursor.fetchall()
            for record in records:
                (gene_id, match_id, cds_start, cds_end,
                 corr_target_start, corr_target_end,
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
                e.GeneID,
                e.EventStart,
                e.EventEnd,
                e.ExpansionID
            FROM Expansions AS e
            WHERE e.Mode="FULL"
            ORDER BY
                e.GeneID, e.ExpansionID;
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
                f.FragmentID,
                f.GeneID,
                g.GeneStart,
                g.GeneEnd,
                g.GeneChrom,
                f.QueryExonStart,
                f.QueryExonEnd,
                f.QueryStart,
                f.QueryEnd,
                f.TargetStart,
                f.TargetEnd,
                f.QueryStrand,
                f.TargetStrand,
                f.QueryAlnProtSeq,
                f.TargetAlnProtSeq
            FROM Matches as f
            INNER JOIN Genes as g ON g.GeneID=f.GeneID
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
            cursor.execute("SELECT GeneID FROM Genes")
            return [record[0] for record in cursor.fetchall()]

    def query_genes_with_duplicated_cds(
        self,
    ) -> list:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute("SELECT DISTINCT GeneID FROM Expansions")
            return [record for record in cursor.fetchall()]

    def export_all_tables_to_csv(
            self,
            output_dir: Path
    ):
        with sqlite3.connect(
                self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
            tables = [table[0] for table in cursor.fetchall() if "sqlite" not in table[0]]
            for table in tables:
                table_name = table
                df = pd.read_sql_query(f"SELECT * FROM {table_name}", db)
                csv_file_path = output_dir / f"{table_name}.csv"
                df.to_csv(csv_file_path, index=False)
