import sqlite3
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
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

    def build_gene_graph(self):
        gene_graph = nx.Graph(id=self.id)
        for expansion in self.expansions.values():
            for node, data in expansion.graph.nodes(data=True):
                gene_graph.add_node(node, **data)
            for source, target, edge_data in expansion.graph.edges(data=True):
                gene_graph.add_edge(source, target, **edge_data)
        return gene_graph

    def draw_expansions_multigraph(
            self,
            expansion_id: int = None,
            figure_path: Path = None,
            figure_size: tuple[float, float] = (8.0, 8.0),
            legend: bool = True,
            connect_overlapping_nodes: bool = True,
            full_expansion: bool = False
    ):
        if expansion_id is not None:
            G = self.expansions[expansion_id].graph
        else:
            G = self.build_gene_graph()
        self.plot_handler.draw_expansions_multigraph(
            gene_start=self.coordinates.lower,
            gene_graph=G,
            figure_path=figure_path,
            figure_size=figure_size,
            legend=legend,
            connect_overlapping_nodes=connect_overlapping_nodes,
            full_expansion=full_expansion
        )


class Expansion:
    def __init__(
            self,
            expansion_id,
            nodes,
            edges
    ):
        self.graph = nx.Graph()
        self.graph.id = expansion_id
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
                coordinates=P.open(start, end),
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
        self._color_map = {
            self._partial_insertion: 'blue',
            self._partial_excision: 'purple',
            self._full: 'green',
            self._intronic: 'red',
            self._inter_boundary: 'orange'
        }

    @staticmethod
    def component_positions(
            components: list[list],
            gene_graph: nx.MultiGraph
    ):
        layout_scale = 2
        component_positions = []
        large_components = [comp for comp in components if len(comp) > 2]
        small_components = [comp for comp in components if len(comp) <= 2]
        # Generate positions for large components with horizontal shifting
        for component in large_components:
            layout = nx.circular_layout(gene_graph.subgraph(component), scale=layout_scale)
            component_positions.append(layout)
        # Set a reduced horizontal shift between large and small components
        position_shift = layout_scale * 3  # Reduced horizontal shift for a tighter layout
        component_position = {}
        # Special case handling for two-node graph
        if len(gene_graph.nodes) == 2:
            nodex, nodey = list(gene_graph.nodes)
            component_position[nodex] = (0, 0)
            component_position[nodey] = (0.2, 0)  # Small offset to avoid overlap
        else:
            # Place the large components with horizontal shifts
            for event_idx, layout in enumerate(component_positions):
                for node, position in layout.items():
                    x, y = position
                    shifted_position = (x + event_idx * position_shift, y)
                    component_position[node] = shifted_position
            # Calculate the vertical range (y_min and y_max) based on large component positions
            if component_position:
                y_values = [y for x, y in component_position.values()]
                y_min, y_max = min(y_values), max(y_values)
            else:
                y_min, y_max = 0, 0  # Default values if no components are positioned
            # Determine the starting x-position for small components, to the right of large components
            last_large_x = max(x for x, y in component_position.values()) if component_position else 0
            small_component_x = last_large_x + position_shift
            # Calculate even vertical spacing for small components within the y_min and y_max range
            num_small_components = len(small_components)
            if num_small_components > 1:
                small_component_spacing = (y_max - y_min) / (num_small_components - 1)
            else:
                small_component_spacing = 0  # No spacing needed if there's only one component
            # Position small components vertically at the fixed x position,
            # with slight offsets for two-node components
            for idx, component in enumerate(small_components):
                y_position = y_max - idx * small_component_spacing  # even spacing
                if len(component) == 2:
                    # Add a small horizontal offset for two-node components to avoid collapsing
                    x, y = list(component)
                    component_position[x] = (small_component_x - 1, y_position)
                    component_position[y] = (small_component_x + 1, y_position)
                else:
                    for node in component:
                        component_position[node] = (small_component_x, y_position)
        return component_position

    @staticmethod
    def full_expansion(
            graph: nx.Graph,
    ) -> None:
        nodes_to_drop = [
            node
            for node in graph.nodes
            if graph.nodes[node].get('type') != 'FULL'
        ]
        if nodes_to_drop:
            graph.remove_nodes_from(nodes_to_drop)

    def draw_expansions_multigraph(
            self,
            gene_start: int,
            gene_graph: [nx.Graph, nx.MultiGraph],
            figure_size: tuple[float, float] = (8.0, 8.0),
            figure_path: Path = None,
            legend: bool = True,
            connect_overlapping_nodes: bool = True,
            full_expansion: bool = False
    ):
        G = gene_graph.copy()
        if full_expansion:
            self.full_expansion(
                graph=G
            )
        if G.number_of_nodes() > 1:
            plt.figure(figsize=figure_size)
            node_colors = [
                self._color_map[attrib['type']]
                for node, attrib in G.nodes(data=True)
            ]
            components = list(nx.connected_components(G))
            node_labels = {
                node: f'({node.lower - gene_start},{node.upper - gene_start})'
                for node in G.nodes
            }
            component_position = self.component_positions(
                components=components,
                gene_graph=G
            )
            # Adjust label positions
            label_positions = {
                node: (x, y + 0.1)
                for node, (x, y) in component_position.items()
            }
            nx.draw_networkx_nodes(
                G=G,
                node_color=node_colors,
                pos=component_position,
                node_size=350,
            )
            nx.draw_networkx_labels(
                G,
                label_positions,
                labels=node_labels,
                font_size=8,
                bbox=dict(
                    boxstyle="round,pad=0.3",
                    edgecolor="white",
                    facecolor="white"
                )
            )
            for edge in G.edges(data=True):
                source, target, attributes = edge
                edge_style = attributes.get('style', 'solid')
                edge_color = attributes.get('color', 'black')
                edge_width = attributes.get('width', 1)
                nx.draw_networkx_edges(
                    G,
                    component_position,
                    edgelist=[(source, target)],
                    edge_color=edge_color,
                    style=edge_style,
                    width=edge_width
                )
            if connect_overlapping_nodes:
                overlapping_nodes = self.fetch_overlapping_nodes(
                    gene_graph=G
                )
                for cluster in overlapping_nodes:
                    for pair in cluster:
                        nodei, nodej = pair
                        nx.draw_networkx_edges(
                            G,
                            component_position,
                            edgelist=[(nodei, nodej)],
                            edge_color='red',
                            style='dotted',
                            width=2
                        )
            if legend:
                node_attributes = {
                    attrib['type']
                    for _, attrib in G.nodes(data=True)
                }
                legend_elements = [
                    mlines.Line2D(
                        [], [],
                        color=self._color_map[label],
                        marker='o',
                        linestyle='None',
                        markersize=10,
                        label=label
                    )
                    for label in node_attributes
                ]
                plt.legend(
                    handles=legend_elements,
                    loc="upper right",
                    bbox_to_anchor=(1.2, 1),
                    frameon=False
                )

            for spine in plt.gca().spines.values():
                spine.set_visible(False)
            if figure_path:
                plt.savefig(figure_path)
            else:
                plt.show()

    def fetch_overlapping_nodes(
            self,
            gene_graph: nx.MultiGraph
    ):
        overlapping_clusters = self._get_overlapping_clusters(
            target_coordinates_set=set([
                (coordinate, None)
                for coordinate in gene_graph.nodes
            ]),
            threshold=0
        )
        overlapping_clusters = [
            [coordinate for coordinate, _ in cluster]
            for cluster in overlapping_clusters
            if len(cluster) > 1
        ]
        pairs_list = [
            [(I, J)
             for indx_i, I in enumerate(cluster)
             for J in cluster[indx_i + 1:]]
            for cluster in overlapping_clusters
        ]
        return pairs_list

    def _get_overlapping_clusters(
            self,
            target_coordinates_set: set[tuple],
            threshold: float,
    ) -> list[list[tuple]]:
        processed_intervals = set()
        overlapping_clusters = []
        sorted_coordinates = sorted(
            target_coordinates_set,
            key=lambda x: (x[0].lower, x[0].upper)
        )
        for target_coordinate, evalue in sorted_coordinates:
            if target_coordinate not in processed_intervals:
                processed_intervals.add(target_coordinate)
                processed_intervals, cluster = self._find_interval_clusters(
                    sorted_coordinates=sorted_coordinates,
                    processed_intervals=processed_intervals,
                    cluster=[(target_coordinate, evalue)],
                    threshold=threshold
                )
                overlapping_clusters.append(cluster)
        overlapping_clusters.sort(key=len, reverse=True)
        return overlapping_clusters

    def _find_interval_clusters(
            self,
            sorted_coordinates: list,
            processed_intervals: set,
            cluster: list[tuple],
            threshold: float
    ) -> tuple:
        new_cluster = list(cluster)
        for other_coordinate, other_evalue in sorted_coordinates:
            if (other_coordinate not in processed_intervals and all(
                    round(self._min_perc_overlap(
                        intv_i=target_coordinate,
                        intv_j=other_coordinate), 1) >= threshold if threshold > 0 else
                    round(self._min_perc_overlap(
                        intv_i=target_coordinate,
                        intv_j=other_coordinate), 1) > threshold
                    for target_coordinate, evalue in new_cluster
            )):
                new_cluster.append((other_coordinate, other_evalue))
                processed_intervals.add(other_coordinate)
        if new_cluster == cluster:
            return processed_intervals, new_cluster
        else:
            return self._find_interval_clusters(
                sorted_coordinates=sorted_coordinates,
                processed_intervals=processed_intervals,
                cluster=new_cluster,
                threshold=threshold
            )

    @staticmethod
    def _min_perc_overlap(
            intv_i: P.Interval,
            intv_j: P.Interval,
    ) -> float:
        def get_interval_length(
                interval: P.Interval,
        ):
            return sum(intv.upper - intv.lower for intv in interval)
        if intv_i.overlaps(intv_j):
            intersection_span = get_interval_length(intv_i.intersection(intv_j))
            longest_length = max(get_interval_length(intv_i), get_interval_length(intv_j))
            return round(intersection_span / longest_length, 3)
        return 0.0


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
