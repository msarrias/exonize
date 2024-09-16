# ------------------------------------------------------------------------
# This module contains the SqliteHandler class, which is used to handle the
# results database.
# ------------------------------------------------------------------------
import sqlite3
import contextlib
from pathlib import Path
import pandas as pd
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

    @staticmethod
    def batch(iterable, n=1):
        length = len(iterable)
        for ndx in range(0, length, n):
            yield iterable[ndx:min(ndx + n, length)]

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
                CREATE TABLE IF NOT EXISTS Expansions (
                    MatchID INTEGER PRIMARY KEY AUTOINCREMENT,
                    GeneID VARCHAR(100),
                    Mode TEXT CHECK(Mode IN (
                                            'FULL',
                                            'INSERTION',
                                            'EXCISION',
                                            'CANDIDATE',
                                            'TRUNCATION'
                     )),
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
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS Expansions_transcript_interdependence (
                    GeneID VARCHAR(100),
                    ExpansionID INTEGER NOT NULL,
                    NumberTranscripts INTEGER NOT NULL,
                    NumberCodingEvents INTEGER NOT NULL,
                    All_ INTEGER NOT NULL,
                    Present INTEGER NOT NULL,
                    Absent INTEGER NOT NULL,
                    Neither INTEGER NOT NULL,
                    Classification TEXT CHECK(Classification IN (
                                            'OBLIGATE',
                                            'EXCLUSIVE',
                                            'FLEXIBLE',
                                            'OPTIONAL_FLEXIBLE',
                                            'OPTIONAL_EXCLUSIVE',
                                            'OPTIONAL_OBLIGATE'
                     )),
                    ExclusiveEvents TEXT,
                    FOREIGN KEY (GeneID) REFERENCES Genes(GeneID),
                    PRIMARY KEY (ExpansionID, GeneID)
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
        cursor.execute("""
        CREATE INDEX IF NOT EXISTS
            Matches_full_length_idx ON Matches_full_length (FragmentID);
        """)
        columns_to_add = [
            ("CorrectedTargetStart", "INTEGER"),
            ("CorrectedTargetEnd", "INTEGER")
        ]
        for column_name, column_type in columns_to_add:
            self.add_column_to_table(
                table_name="Matches_full_length",
                column_name=column_name,
                column_type=column_type,
            )

    def create_full_expansions_table(
            self,
    ):
        with sqlite3.connect(
                self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute("""
            CREATE TABLE IF NOT EXISTS Expansions_Full AS
            WITH FullExpansionCounts AS (
                SELECT GeneID, ExpansionID, COUNT(*) AS FullCount
            FROM Expansions
            WHERE Mode = 'FULL'
            GROUP BY GeneID, ExpansionID
            HAVING COUNT(*) > 1
            )
            SELECT
                e.*
            FROM Expansions e
            JOIN FullExpansionCounts f ON e.GeneID=f.GeneID AND e.ExpansionID = f.ExpansionID
            WHERE e.Mode='FULL'
            ORDER BY e.GeneID, e.ExpansionID;
            """)

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
                CorrectedTargetEnd=?
                WHERE FragmentID=?
                """,
                list_tuples,
            )

    def insert_identity_and_dna_algns_columns(
            self,
            list_tuples: list
    ) -> None:
        columns_to_add = [
            ("DNAIdentity", "REAL"),
            ("ProtIdentity", "REAL"),
            ("QueryDNASeq", "VARCHAR"),
            ("TargetDNASeq", "VARCHAR"),
        ]
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            for column_name, column_type in columns_to_add:
                self.add_column_to_table(
                    table_name="Matches",
                    column_name=column_name,
                    column_type=column_type,
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

    def update_has_duplicate_genes_table(
            self,
            list_tuples: list
    ) -> None:
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

    def insert_percent_query_column_to_fragments(
            self,
    ) -> None:
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

    def insert_gene_ids_table(
            self,
            gene_args_tuple: tuple,
    ) -> None:
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

    def insert_expansion_table(
            self,
            list_tuples: list,
    ) -> None:
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
        with sqlite3.connect(self.results_database_path, timeout=self.timeout_database) as db:
            cursor = db.cursor()
            cursor.executemany(insert_gene_table_param, list_tuples)

    def insert_expansions_interdependence_classification(
            self,
            list_tuples: list
    ):
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.executemany("""
            INSERT INTO Expansions_transcript_interdependence(
                GeneID,
                ExpansionID,
                NumberTranscripts,
                NumberCodingEvents,
                All_,
                Present,
                Absent,
                Neither,
                Classification,
                ExclusiveEvents
            )
            VALUES (?, ?, ?, ?, ?, ?, ? ,?, ?, ?)
            """, list_tuples)

    def create_non_reciprocal_fragments_table(
            self,
    ) -> None:
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute("""
            SELECT
                sql
            FROM sqlite_master
            WHERE type='table' AND name='Matches_full_length';
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
                COALESCE(CorrectedTargetEnd, TargetEnd) + GeneStart AS CorrectedTargetEnd
            FROM Matches_full_length
            WHERE FragmentID IN ({placeholders});
                    """
            cursor.execute(query, fragment_ids_list)

    def insert_mode_in_non_reciprocal_fragments_table(
            self,
            list_tuples: list
    ) -> None:
        self.add_column_to_table(
            table_name="Matches_full_length_non_reciprocal",
            column_name="Mode",
            column_type="""
            Mode TEXT CHECK(Mode IN (
            'FULL',
            'INSERTION',
            'EXCISION',
            'CANDIDATE',
            'TRUNCATION'
            ))
            """,
        )
        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.executemany(
                """
                UPDATE Matches_full_length_non_reciprocal
                SET Mode=? WHERE FragmentID=?
                """,
                list_tuples,
            )

    def insert_matches_interdependence_classification(
            self,
            tuples_list: list
    ) -> None:
        column_names = [
            'NumberTranscripts',
            'All_', 'Present', 'Absent', 'Neither'
        ]
        for column_name in column_names:
            self.add_column_to_table(
                table_name="Matches_full_length_non_reciprocal",
                column_name=column_name,
                column_type="INTEGER",
            )
        self.add_column_to_table(
            table_name="Matches_full_length_non_reciprocal",
            column_name='Classification',
            column_type="Classification TEXT "
                        "CHECK(Classification IN "
                        "('OBLIGATE','EXCLUSIVE','FLEXIBLE',"
                        "'OPTIONAL_FLEXIBLE','OPTIONAL_EXCLUSIVE','OPTIONAL_OBLIGATE'))"
        )

        with sqlite3.connect(
            self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            insert_full_length_event_table_param = """
            UPDATE Matches_full_length_non_reciprocal
            SET
                NumberTranscripts=?,
                All_=?,
                Present=?,
                Absent=?,
                Neither=?,
                Classification=?
                WHERE FragmentID=?
            """
            cursor.executemany(insert_full_length_event_table_param, tuples_list)

    def query_full_length_events(
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

    def query_raw_matches(
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

    def query_non_reciprocal_coding_matches(
            self,
    ) -> list:
        with sqlite3.connect(
                self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute(
                """
            SELECT
                GeneID,
                FragmentID,
                QueryExonStart,
                QueryExonEnd,
                CorrectedTargetStart,
                CorrectedTargetEnd
            FROM Matches_full_length_non_reciprocal
            WHERE Mode="FULL"
            ORDER BY
                GeneID, FragmentID;
            """)
            return cursor.fetchall()

    def query_coding_expansion_events(
            self,
    ) -> dict:
        with sqlite3.connect(
                self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute(
                """
            SELECT
                GeneID, ExpansionID, EventStart, EventEnd
            FROM Expansions_Full
            ORDER BY
                GeneID, ExpansionID;
            """
            )
            records = cursor.fetchall()
            expansions_dictionary = defaultdict(lambda: defaultdict(list))
            for record in records:
                gene_id, expansion_id, event_start, event_end = record
                expansions_dictionary[gene_id][expansion_id].append(P.open(event_start, event_end))
            return expansions_dictionary

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
    ) -> None:
        with sqlite3.connect(
                self.results_database_path, timeout=self.timeout_database
        ) as db:
            cursor = db.cursor()
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
            tables = [table[0]
                      for table in cursor.fetchall()
                      if ("sqlite" not in table[0] and "Matches" not in table[0])]
            for table in tables:
                table_name = table
                df = pd.read_sql_query(f"SELECT * FROM {table_name}", db)
                csv_file_path = output_dir / f"{table_name}.csv"
                df.to_csv(csv_file_path, index=False)
