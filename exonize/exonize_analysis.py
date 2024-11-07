import sqlite3
import networkx as nx
import portion as P


class Gene:
    def __init__(
            self,
            gene_id,
            coordinates,
            strand,
            chromosome
    ):
        self.gene_id = gene_id
        self.coordinates = coordinates
        self.strand = strand
        self.chromosome = chromosome
        self.expansions = {}

    def __getitem__(self, expansion_id):
        return self.expansions[expansion_id]


class Expansions:
    def __init__(
            self,
            exonize_db_path: str
    ):
        self.exonize_db_path = exonize_db_path
        self.genes = {}
        self._db_handler = ExonizeDBHandler(self.exonize_db_path)

    def __getitem__(self, gene_id):
        return self.genes[gene_id]

    def add_gene(
            self,
            gene_id: str
    ) -> Gene:
        if gene_id not in self.genes:
            chrom, strand, start, end = self._genes_dict[gene_id]
            return Gene(
                id=gene_id,
                coordinates=(start, end),
                strand=strand,
                chromosome=chrom
            )

    def build_expansion_graph(
            self,
            gene_id: str,
            expansion_id: int
    ):
        expansion_graph = nx.Graph()
        nodes = self._db_handler.gene_expansions_dict[gene_id][expansion_id]['nodes']
        expansion_graph.add_nodes_from(nodes)
        for edge in self._db_handler.gene_expansions_dict[gene_id][expansion_id]['edges']:
            q_start, q_end, t_start, t_end, mode = edge
            expansion_graph.add_edge(
                u_for_edge=(q_start, t_start),
                v_for_edge=(t_start, t_end),
                mode=mode
            )
        return expansion_graph

    def build_expansions(self):
        for gene_id, non_reciprocal_matches in self._db_handler.gene_expansions_dict.items():
            self.genes[gene_id] = self.add_gene(
                gene_id=gene_id
            )
            for expansion_id in non_reciprocal_matches.keys():
                self.genes[gene_id].expansions[expansion_id] = self.build_expansion_graph(
                    gene_id=gene_id,
                    expansion_id=expansion_id
                )


class ExonizeDBHandler:
    def __init__(self, db_path):
        self.db_path = db_path
        self.genes_dict = self.collect_genes(self.fetch_genes_with_exon_dup())
        self.gene_expansions_dict = dict(dict(dict(list)))
        self._expansions_nodes = self.fetch_expansions_nodes()
        self._expansions_edges = self.fetch_expansions_edges()
        self.collect_expansion_nodes(expansions=self._expansions_nodes)
        self.collect_expansions_edges(matches_set=self._expansions_edges)

    def check_if_table_exists(
            self,
            table_name: str
    ) -> bool:
        with sqlite3.connect(self.db_path) as db:
            cursor = db.cursor()
            cursor.execute(
                f"""
                SELECT
                    name
                FROM sqlite_master
                WHERE type='table' AND name='{table_name}';
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
            return cursor.fetchall()

    @staticmethod
    def collect_genes(
            genes_records: list
    ) -> dict:
        return {
            geneid: (chrom, strand, start, end)
            for geneid, chrom, strand, start, end in genes_records
        }

    def collect_expansion_nodes(
            self,
            expansions
    ):
        for gene_id, mode, event_start, event_end, expansion_id in expansions:
            if gene_id not in self.gene_expansions_dict:
                self.gene_expansions_dict[gene_id] = {}
            if expansion_id not in self.gene_expansions_dict[gene_id]:
                self.gene_expansions_dict[gene_id][expansion_id] = dict(nodes=[], edges=[])
            self.gene_expansions_dict[gene_id][expansion_id]['nodes'].append(
                P.open(event_start, event_end)
            )

    def collect_expansions_edges(
            self,
            matches_set: set,
    ):
        for match in matches_set:
            gene_id, q_start, q_end, t_start, t_end, mode = match
            for expansion_id, expansion_atrributes in self.gene_expansions_dict[gene_id].items():
                if P.open(q_start, q_end) and P.open(t_start, t_end) in expansion_atrributes['nodes']:
                    expansion_atrributes['edges'].append(
                        (P.open(q_start, q_end), P.open(t_start, t_end), mode)
                    )

    def fetch_expansions_nodes(
            self
    ):
        with sqlite3.connect(self.db_path) as db:
            cursor = db.cursor()
            cursor.execute("""
                        SELECT
                            GeneID,
                            Mode,
                            EventStart,
                            EventEnd,
                            ExpansionID
                        FROM Expansions
                        """)
        return cursor.fetchall()

    def fetch_expansions_edges(self):
        local_matches = set()
        global_matches = set()
        with sqlite3.connect(self.db_path) as db:
            cursor = db.cursor()
            if self.check_if_table_exists(
                    table_name="Local_matches_non_reciprocal"
            ):
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
            if self.check_if_table_exists(
                    table_name="Global_matches_non_reciprocal"
            ):
                cursor.execute(
                    """
                    SELECT
                        GeneID,
                        QueryExonStart,
                        QueryExonEnd,
                        TargetExonStart,
                        TargetExonEnd
                    FROM Global_matches_non_reciprocal
                    GROUP BY
                        GeneID,
                        QueryExonStart,
                        QueryExonEnd,
                        TargetExonStart,
                        TargetExonEnd;
                    """
                )
                global_matches = set((*res, 'FULL') for res in cursor.fetchall())
        return global_matches.union(local_matches)
