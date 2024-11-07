import sqlite3

import matplotlib.pyplot as plt
import networkx as nx
import portion as P
from pathlib import Path


class Gene:
    def __init__(
            self,
            gene_id,
            coordinates,
            strand,
            chromosome
    ):
        self.id = gene_id
        self.coordinates = coordinates
        self.strand = strand
        self.chromosome = chromosome
        self.expansions = {}
        self.plot_handler = PlotHandler()

    def __getitem__(self, expansion_id):
        return self.expansions[expansion_id].graph

    def __iter__(self):
        return iter(expansion.graph for expansion in self.expansions.values())

    def __repr__(self):
        return f"<Gene {self.id} with {len(self.expansions)} expansions (iterable of expansion graphs)>"

    def __len__(self):
        return len(self.expansions)

    def draw_expansions_multigraph(
            self,
            figure_path: Path = None
    ):
        self.plot_handler.draw_expansions_multigraph(
            self.expansions.graph,
            figure_path
        )


class Expansion:
    def __init__(
            self,
            expansion_id,
            nodes,
            edges
    ):
        self.expansion_id = expansion_id
        self.graph = nx.Graph()
        self.graph.add_nodes_from([
            (coord, {"type": node_type})
            for coord, node_type in nodes
        ])
        for edge in edges:
            query_coord, target_coord, mode = edge
            self.graph.add_edge(
                u_of_edge=query_coord,
                v_of_edge=target_coord,
                mode=mode
            )


class GenomeExpansions:
    def __init__(
            self,
            exonize_db_path: str
    ):
        self.exonize_db_path = exonize_db_path
        self._genes = {}
        self._db_handler = ExonizeDBHandler(self.exonize_db_path)
        self.build_expansions()

    def __iter__(self):
        return iter(self._genes.values())

    def __contains__(self, n):
        return n in self._genes

    def __getitem__(self, gene_id):
        return self._genes[gene_id]

    def __len__(self):
        return len(self._genes)

    @property
    def genes(self):
        return list(self._genes.keys())

    def add_gene(
            self,
            gene_id: str
    ) -> Gene:
        if gene_id not in self._genes:
            chrom, strand, start, end = self._db_handler.genes_dict[gene_id]
            return Gene(
                gene_id=gene_id,
                coordinates=(start, end),
                strand=strand,
                chromosome=chrom
            )

    def build_expansions(self):
        for gene_id, non_reciprocal_matches in self._db_handler.gene_expansions_dict.items():
            self._genes[gene_id] = self.add_gene(
                gene_id=gene_id
            )
            for expansion_id, data in non_reciprocal_matches.items():
                nodes = data['nodes']
                edges = data['edges']
                expansion = Expansion(
                    expansion_id=expansion_id,
                    nodes=nodes,
                    edges=edges
                )
                self._genes[gene_id].expansions[expansion_id] = expansion


class PlotHandler:
    def __init__(self):
        self._full = 'FULL'
        self._partial_insertion = 'PARTIAL_INSERTION'
        self._partial_excision = 'PARTIAL_EXCISION'
        self._inter_boundary = 'INTER_BOUNDARY'
        self._intronic = 'INTRONIC'
        self.color_map = {
            self._partial_insertion: 'blue',
            self._partial_excision: 'purple',
            self._full: 'green',
            self._intronic: 'red',
            self._inter_boundary: 'orange'
        }

    def draw_expansions_multigraph(
            self,
            gene_graph: nx.MultiGraph,
            figure_path: Path = None,
            draw_overlapping_edges: bool = False,
    ):
        plt.figure(figsize=(16, 8))
        node_colors = [
            self._color_map[node[1]['type']]
            for node in gene_graph.nodes(data=True)
        ]
        components = list(nx.connected_components(gene_graph))
        node_labels = {
            node: f'({node[0]},{node[1]})'
            for node in gene_graph.nodes
        }
        # Create a separate circular layout for each component
        layout_scale = 2
        component_positions = []
        for component in components:
            layout = nx.circular_layout(
                gene_graph.subgraph(component),
                scale=layout_scale
            )
            component_positions.append(layout)
        position_shift = max(layout_scale * 5.5, 15)
        component_position = {}
        for event_idx, layout in enumerate(component_positions):
            for node, position in layout.items():
                shifted_position = (position[0] + event_idx * position_shift, position[1])
                component_position[node] = shifted_position
        if max([len(component) for component in components]) == 2:
            label_positions = component_position
        else:
            label_positions = {
                node: (position[0], position[1] + 0.1)
                for node, position in component_position.items()
            }
        nx.draw_networkx_nodes(
            G=gene_graph,
            node_color=node_colors,
            pos=component_position,
            node_size=350,
        )
        nx.draw_networkx_labels(
            gene_graph,
            label_positions,
            labels=node_labels,
            font_size=8,
            bbox=dict(
                boxstyle="round,pad=0.3",
                edgecolor="white",
                facecolor="white"
            )
        )
        # overlapping_edges:
        for edge in gene_graph.edges(data=True):
            source, target, attributes = edge
            edge_style = attributes.get('style', 'solid')
            edge_color = attributes.get('color', 'black')
            edge_width = attributes.get('width', 1)
            nx.draw_networkx_edges(
                gene_graph,
                component_position,
                edgelist=[(source, target)],
                edge_color=edge_color,
                style=edge_style,
                width=edge_width
            )
        plt.show()
        # overlapping_edges = self.fetch_overlapping_edges(
        #     gene_graph=gene_graph
        # )
        # for cluster in overlapping_edges:
        #     for pair in cluster:
        #         query, targ = pair
        #         nx.draw_networkx_edges(
        #             gene_graph,
        #             component_position,
        #             edgelist=[((query.lower, query.upper),
        #                        (targ.lower, targ.upper))],
        #             edge_color='red',
        #             style='dotted',
        #             width=2
        #         )
        # plt.savefig(figure_path)
        plt.close()


class ExonizeDBHandler:
    def __init__(self, db_path):
        self.db_path = db_path
        self.genes_dict = self.collect_genes(self.fetch_genes_with_exon_dup())
        self.gene_expansions_dict = {}
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
                (P.open(event_start, event_end), mode)
            )

    def collect_expansions_edges(
            self,
            matches_set: set,
    ):
        for match in matches_set:
            gene_id, q_start, q_end, t_start, t_end, mode = match
            for expansion_id, expansion_atrributes in self.gene_expansions_dict[gene_id].items():
                node_coordinates = {coordinate for coordinate, _ in expansion_atrributes['nodes']}
                if P.open(q_start, q_end) and P.open(t_start, t_end) in node_coordinates:
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
                global_matches = set(
                    (*res, 'FULL')
                    for res in cursor.fetchall()
                )
        return global_matches.union(local_matches)
