import sqlite3


class ExonizeDBHandler:
    def __init__(self, db_path):
        self.db_path = db_path

    def check_if_table_exists(self, table_name):
        with sqlite3.connect(self.db_path) as db:
            cursor = db.cursor()
            cursor.execute(
                f"""
                SELECT name FROM sqlite_master WHERE type='table' AND name='{table_name}';
                """
            )
            return cursor.fetchone() is not None

    def fetch_genes_with_exon_dup(self):
        with sqlite3.connect(self.db_path) as db:
            cursor = db.cursor()
            cursor.execute(
                """
                SELECT
                    GeneID,
                    GeneChrom,
                    GeneStrand,
                    GeneStart,
                    GeneEnd
                FROM Genes
                WHERE Duplication=1;
                """
            )
            return {
                geneid: (chrom, strand, start, end)
                for geneid, chrom, strand, start, end in cursor.fetchall()
            }

    def fetch_all_non_reciprocal_matches(self):
        local_matches = set()
        global_matches = set()
        with sqlite3.connect(self.db_path) as db:
            cursor = db.cursor()
            if self.check_if_table_exists("Local_matches_non_reciprocal"):
                cursor.execute(
                    """
                    SELECT
                        GeneID,
                        QueryExonStart,
                        QueryExonEnd,
                        CorrectedTargetStart,
                        CorrectedTargetEnd,
                        Mode
                    FROM Local_matches_non_reciprocal
                    """
                )
                local_matches = set(cursor.fetchall())
            if self.check_if_table_exists("Global_matches_non_reciprocal"):
                cursor.execute(
                    """
                    SELECT
                        GeneID,
                        QueryExonStart,
                        QueryExonEnd,
                        TargetExonStart,
                        TargetExonEnd
                    FROM Global_matches_non_reciprocal
                    GROUP BY GeneID, QueryExonStart, QueryExonEnd, TargetExonStart, TargetExonEnd;
                    """
                )
                global_matches = set((*res, 'FULL') for res in cursor.fetchall())
        return global_matches.union(local_matches)
