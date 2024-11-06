import sqlite3
import networkx as nx


class ExonizeDBHandler:
    def __init__(self, db_path):
        self.db_path = db_path

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
            return {
                geneid: (chrom, strand, start, end)
                for geneid, chrom, strand, start, end in cursor.fetchall()
            }

    def fetch_expansions_nodes_and_edges(self):
        local_matches = set()
        global_matches = set()
        gene_matches_dict = dict(dict(dict(list)))
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
                    GROUP BY GeneID, QueryExonStart, QueryExonEnd, TargetExonStart, TargetExonEnd;
                    """
                )
                global_matches = set((*res, 'FULL') for res in cursor.fetchall())
            cursor.execute("""
            SELECT
                GeneID,
                Mode,
                EventStart,
                EventEnd,
                ExpansionID
            FROM Expansions
            """)
            expansions = cursor.fetchall()
        for gene_id, mode, event_start, event_end, expansion_id in expansions:
            if gene_id not in gene_matches_dict:
                gene_matches_dict[gene_id] = {}
            if expansion_id not in gene_matches_dict[gene_id]:
                gene_matches_dict[gene_id][expansion_id] = dict(nodes=[], edges=[])
            gene_matches_dict[gene_id][expansion_id]['nodes'].append(P.open(event_start, event_end))
        for match in global_matches.union(local_matches):
            gene_id, q_start, q_end, t_start, t_end, mode = match
            for expansion_id, expansion_atrributes in gene_matches_dict[gene_id].items():
                if P.open(q_start, q_end) and P.open(t_start, t_end) in expansion_atrributes['nodes']:
                    expansion_atrributes['edges'].append(
                        (P.open(q_start, q_end), P.open(t_start, t_end), mode)
                    )
        return gene_matches_dict


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
        self._expansions_dictionary = db_handler.fetch_expansions_nodes_and_edges()
        self._genes_dict = db_handler.fetch_genes_with_exon_dup()

    def __getitem__(self, gene_id):
        return self.genes[gene_id]

    def add_gene(
            self,
            gene_id: str
    ) -> None:
        if gene_id not in self.genes:
            chrom, strand, start, end = self._genes_dict[gene_id]
            self.genes[gene_id]  = Gene(
                id=gene_id,
                coordinates=(start, end),
                strand=strand,
                chromosome=chrom
            )

    def build_expansion_graph(
            self,
            expansion_id: int
    ):
        expansion_graph = nx.Graph()
        nodes = self._expansions_dictionary[gene_id][expansion_id]['nodes']
        expansion_graph.add_nodes_from(nodes)
        for edge in self._expansions_dictionary[gene_id][expansion_id]['edges']:
            q_start, q_end, t_start, t_end, mode = edge
            expansion_graph.add_edge(
                u_for_edge=(q_start, t_start),
                v_for_edge=(t_start, t_end),
                mode=mode
            )
        self.genes[gene_id].expansions[expansion_id] = expansion_graph

    def build_expansions(self):
        for gene_id, non_reciprocal_matches in self._expansions_dictionary.items():
            self.add_gene(gene_id=gene_id)
            for expansion_id in non_reciprocal_matches.keys():
                self.build_expansion_graph(
                    expansion_id=expansion_id
                )
