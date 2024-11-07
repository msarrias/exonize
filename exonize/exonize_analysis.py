import sqlite3
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import networkx as nx
import portion as P
from pathlib import Path


class Gene:
    """
    Gene class is a container for gene expansion graphs.

    Attributes
    ----------
    id : str
        The unique identifier for the gene.
    coordinates : tuple
        The start and end coordinates of the gene on the chromosome.
    strand : str
        The DNA strand ('+' or '-') on which the gene is located.
    chromosome : str
        The chromosome on which the gene is located.
    expansions : dict
        A dictionary where keys are expansion IDs and values are expansion objects, each containing an expansion graph.
    plot_handler : PlotHandler
        An instance of PlotHandler responsible for visualizing the gene's expansions.
    """
    def __init__(
            self,
            gene_id,
            coordinates,
            strand,
            chromosome
    ):
        """
        Initializes a Gene instance.

        Parameters
        ----------
        gene_id : str
            The unique identifier for the gene.
        coordinates : tuple
            The start and end coordinates of the gene on the chromosome.
        strand : str
            The DNA strand ('+' or '-') on which the gene is located.
        chromosome : str
            The chromosome on which the gene is located.
        """
        self.id = gene_id
        self.coordinates = coordinates
        self.strand = strand
        self.chromosome = chromosome
        self.expansions = {}
        self.plot_handler = PlotHandler()

    def __getitem__(self, expansion_id):
        """
        Retrieves the expansion graph for the specified expansion ID.

        Parameters
        ----------
        expansion_id : int
            The ID of the expansion to retrieve.

        Returns
        -------
        networkx.Graph
            The expansion graph associated with the given expansion ID.

        Examples
        --------
        >>> gene[1]  # Retrieves the expansion graph for expansion ID 1
        """
        return self.expansions[expansion_id].graph

    def __iter__(self):
        """
        Returns an iterator over the expansion graphs.

        Returns
        -------
        iterator
            An iterator yielding each expansion graph.

        Examples
        --------
        >>> for graph in gene:
        ...     print(graph)
        """
        return iter(expansion.graph for expansion in self.expansions.values())

    def __repr__(self):
        """
        Returns a string representation of the Gene object.

        Returns
        -------
        str
            A string describing the gene's ID and number of expansions.

        Examples
        --------
        >>> repr(gene)
        '<Gene GENE123 with 0 expansions (iterable of expansion graphs)>'
        """
        return f"<Gene {self.id} with {len(self.expansions)} expansions (iterable of expansion graphs)>"

    def __len__(self):
        """
        Returns the number of expansions associated with the gene.

        Returns
        -------
        int
            The number of expansions.

        Examples
        --------
        >>> len(gene)
        0
        """
        return len(self.expansions)

    def build_gene_graph(self):
        """
        Builds and returns a consolidated gene graph containing nodes and edges from all expansion graphs.

        Returns
        -------
        networkx.Graph
            A combined graph with nodes and edges from all expansions.

        Examples
        --------
        >>> combined_graph = gene.build_gene_graph()
        >>> print(combined_graph.nodes)
        >>> print(combined_graph.edges)
        """
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
        """
        Draws a multi-graph of gene expansions.

        Parameters
        ----------
        expansion_id : int, optional
            The ID of a specific expansion to draw. If None, the combined gene graph is drawn.
        figure_path : Path, optional
            The path to save the figure. If None, the figure is not saved.
        figure_size : tuple of float, optional
            The size of the figure in inches. Default is (8.0, 8.0).
        legend : bool, optional
            Whether to display a legend on the plot. Default is True.
        connect_overlapping_nodes : bool, optional
            Whether to draw edges connecting overlapping nodes in the graph. Default is True.
        full_expansion : bool, optional
            Whether to show the full expansion graph only. Default is False.

        """
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
    """
    Expansion class represents an expansion graph for a specific gene expansion.

    Attributes
    ----------
    graph : networkx.Graph
        A NetworkX graph representing the expansion.
    """
    def __init__(
            self,
            expansion_id,
            nodes,
            edges
    ):
        """
        Initializes an Expansion instance.

        Parameters
        ----------
        expansion_id : int
            The unique identifier for the expansion.
        nodes : list of tuples
            A list of tuples representing the nodes, where each tuple is in the form (coord, node_type).
            - coord: The coordinate of the node (e.g., an integer or other unique identifier).
            - node_type: The type of the node (e.g., 'FULL', 'PARTIAL').
        edges : list of tuples
            A list of tuples representing the edges, where each tuple is in the form (query_coord, target_coord, mode).
            - query_coord: The coordinate of the query node.
            - target_coord: The coordinate of the target node.
            - mode: The match mode.

        """
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
    """
    GenomeExpansions class is a container for managing gene expansions across an entire genome.

     Attributes
    ----------
    exonize_db_path : str
        The file path to the Exonize database.
    _genes : dict
        A dictionary of Gene objects, keyed by gene ID.
    _db_handler : ExonizeDBHandler
        An instance of ExonizeDBHandler for interacting with the database.
    """
    def __init__(
            self,
            exonize_db_path: str
    ):
        """
        Initializes a GenomeExpansions instance and builds expansions from the database.

        Parameters
        ----------
        exonize_db_path : str
            The file path to the Exonize database.
        """
        self.exonize_db_path = exonize_db_path
        self._genes = {}
        self._db_handler = ExonizeDBHandler(self.exonize_db_path)
        self.build_expansions()

    def __iter__(self):
        """
        Returns an iterator over the Gene objects.

        Returns
        -------
        iterator
            An iterator yielding each Gene object.

        Examples
        --------
        >>> for gene in genome_expansions:
        ...     print(gene)
        """
        return iter(self._genes.values())

    def __contains__(self, n):
        """
        Checks if a gene ID exists in the GenomeExpansions.

        Parameters
        ----------
        n : str
            The gene ID to check for existence.

        Returns
        -------
        bool
            True if the gene ID exists, False otherwise.

        Examples
        --------
        >>> "GENE123" in genome_expansions
        True
        """
        return n in self._genes

    def __getitem__(self, gene_id):
        """
        Retrieves a Gene object by gene ID.

        Parameters
        ----------
        gene_id : str
            The ID of the gene to retrieve.

        Returns
        -------
        Gene
            The Gene object associated with the specified gene ID.

        Examples
        --------
        >>> gene = genome_expansions["GENE123"]
        >>> print(gene)
        <Gene GENE123 with 0 expansions (iterable of expansion graphs)>
        """
        return self._genes[gene_id]

    def __len__(self):
        """
        Returns the number of genes in the GenomeExpansions.

        Returns
        -------
        int
            The number of genes in the GenomeExpansions.

        Examples
        --------
        >>> len(genome_expansions)
        18
        """
        return len(self._genes)

    @property
    def genes(self):
        """
        Returns a list of gene IDs.

        Returns
        -------
        list
            A list of gene IDs in the GenomeExpansions.

        Examples
        --------
        >>> genome_expansions.genes
        ['GENE123', 'GENE456', 'GENE789']
        """
        return list(self._genes.keys())

    def add_gene(
            self,
            gene_id: str
    ) -> Gene:
        """
        Adds a Gene object to the genome based on the provided gene ID.

        Parameters
        ----------
        gene_id : str
            The unique identifier for the gene.

        Returns
        -------
        Gene
            The Gene object that was added.

        Examples
        --------
        >>> gene = genome_expansions.add_gene("GENE123")
        >>> print(gene)
        <Gene GENE123 with 0 expansions (iterable of expansion graphs)>
        """
        if gene_id not in self._genes:
            chrom, strand, start, end = self._db_handler.genes_dict[gene_id]
            return Gene(
                gene_id=gene_id,
                coordinates=P.open(start, end),
                strand=strand,
                chromosome=chrom
            )

    def build_expansions(self):
        """
        Constructs the expansion data for each gene by loading nodes and edges from the database.

        This method initializes each Gene object and populates its expansions based on data from
        the Exonize database. Each expansion consists of nodes and edges, forming a graph for each gene.

        Examples
        --------
        >>> genome_expansions.build_expansions()
        >>> print(len(genome_expansions))
        18
        """
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
    """
    PlotHandler is a class for managing and visualizing gene expansion graphs
    with various layout and styling options.

    Attributes
    ----------
    _color_map : dict
        A mapping of node types to specific colors for visualization.
    """
    def __init__(self):
        """
        Initializes the PlotHandler instance with predefined node types and color mappings.
        """
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
        """
        Calculates positions for each component in the graph for plotting.

        Parameters
        ----------
        components : list of lists
            A list of connected components in the gene graph, where each component is a list of nodes.
        gene_graph : nx.MultiGraph
            The graph structure containing nodes and edges of the gene.

        Returns
        -------
        dict
            A dictionary mapping each node to a position (x, y) coordinate for visualization.

        """
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
        """
        Removes nodes from the graph that are not labled as 'FULL' type.

        Parameters
        ----------
        graph : nx.Graph

        """
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
        """
        Draws a multi-graph of gene expansions.

        Parameters
        ----------
        gene_start : int
            The start position of the gene for labeling purposes.
        gene_graph : nx.Graph or nx.MultiGraph
            The graph containing nodes and edges of gene expansions.
        figure_size : tuple of float, optional
            The size of the figure in inches. Default is (8.0, 8.0).
        figure_path : Path, optional
            The path to save the figure. If None, the figure is not saved.
        legend : bool, optional
            Whether to display a legend on the plot. Default is True.
        connect_overlapping_nodes : bool, optional
            Whether to draw edges connecting overlapping nodes in the graph. Default is True.
        full_expansion : bool, optional
            Whether to show the full expansion graph only. Default is False.
        """
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
        """
        Identifies and returns overlapping nodes in the gene graph for visualization.

        Parameters
        ----------
        gene_graph : nx.MultiGraph
            The graph containing nodes and edges of gene expansions.

        Returns
        -------
        list of lists
            A list of overlapping node pairs.
        """
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
    """
        ExonizeDBHandler class manages database interactions for gene and expansion data with exon duplications.

        Attributes
        ----------
        db_path : str
            Path to the SQLite exonize database file.
        genes_dict : dict
            Dictionary of gene information, with gene IDs as keys and tuples of
            chromosome, strand, start, and end as values.
        gene_expansions_dict : dict
            Dictionary of expansions associated with each gene, including nodes and edges for each expansion.
        _expansions_nodes : list
            Expansion nodes fetched from the database.
        _expansions_edges : set
           Expansion edges fetched from the database.
        """
    def __init__(self, db_path):
        """
        Initializes an ExonizeDBHandler instance, fetching genes and expansions from the database.

        Parameters
        ----------
        db_path : str
            Path to the SQLite database file.
        """
        self.db_path = db_path
        self.genes_dict = self.collect_genes(self.fetch_genes_with_exon_dup())
        self.gene_expansions_dict = {}
        self._expansions_nodes = self.fetch_expansions_nodes()
        self._expansions_edges = self.fetch_expansions_edges()
        self.collect_expansion_nodes(expansions=self._expansions_nodes)
        self.collect_expansions_edges(matches_set=self._expansions_edges)

    def _check_if_table_exists(
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
            if self._check_if_table_exists(
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
            if self._check_if_table_exists(
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
