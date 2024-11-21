import sqlite3
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from dna_features_viewer import GraphicFeature, GraphicRecord
import networkx as nx
import portion as P
from pathlib import Path
import numpy as np
import pickle


class ExpansionsContainer:
    """A container for managing gene expansions across an entire genome.

    Attributes:
        exonize_db_path (str): The file path to the Exonize database.
    """
    def __init__(
            self,
            exonize_db_path: str
    ):
        """Initializes a ExpansionsContainer instance and builds ExpansionsContainer from the database.

        Args:
            exonize_db_path (str): The file path to the Exonize database.
        """
        self.exonize_db_path = exonize_db_path
        self._genes = {}
        self._db_handler = _ExonizeDBHandler(self.exonize_db_path)
        self.gene_hierarchy_dictionary = {}
        self.build_expansions()

    def __iter__(
            self
    ) -> iter:
        """Returns an iterator over the Gene objects.

        Returns:
            iter: An iterator yielding each Gene object.

        Examples:
            >>> for gene in expansions_container:
            ...     print(gene)
        """
        return iter(self._genes.values())

    def __contains__(
            self,
            n: str
    ) -> bool:
        """Checks if a gene ID exists in the ExpansionsContainer.

        Args:
            n (str): The gene ID to check for existence.

        Returns:
            bool: True if the gene ID exists, False otherwise.

        Examples:
            >>> "GENE123" in expansions_container
            True
        """
        return n in self._genes

    def __getitem__(
            self,
            gene_id: str
    ) -> Gene:
        """Retrieves a Gene object by gene ID.

        Args:
            gene_id (str): The ID of the gene to retrieve.

        Returns:
            Gene: The Gene object associated with the specified gene ID.

        Examples:
            >>> gene = expansions_container["GENE123"]
            >>> print(gene)
            <Gene GENE123 with 0 expansions (iterable of expansion graphs)>
        """
        return self._genes[gene_id]

    def __len__(
            self
    ) -> int:
        """Returns the number of genes in the ExpansionsContainer.

        Returns:
            int: The number of genes in the ExpansionsContainer.

        Examples:
            >>> len(expansions_container)
            18
        """
        return len(self._genes)

    @property
    def genes(
            self
    ) -> list:
        """Returns a list of gene IDs.

        Returns:
            list: A list of gene IDs in the ExpansionsContainer.

        Examples:
            >>> expansions_container.genes
            ['GENE123', 'GENE456', 'GENE789']
        """
        return list(self._genes.keys())

    def _add_gene(
            self,
            gene_id: str
    ) -> Gene:
        """Adds a Gene object to the genome based on the provided gene ID.

        Args:
            gene_id (str): The unique identifier for the gene.

        Returns:
            Gene: The Gene object that was added.

        Examples:
            >>> gene = expansions_container.add_gene("GENE123")
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

    @staticmethod
    def _read_pkl_file(
            file_path: Path
    ) -> dict:
        """
        read_pkl_file is a function that reads a pickle file and returns
         the object stored in it.
        """
        with open(file_path, 'rb') as handle:
            read_file = pickle.load(handle)
        return read_file

    def _fetch_gene_cdss_set(
            self,
            gene_id: str,
            cds_feature: str = 'CDS'

    ) -> list[tuple]:
        return list(
            set(
                annotation['coordinate']
                for mrna_annotation in self.gene_hierarchy_dictionary[gene_id]['mRNAs'].values()
                for annotation in mrna_annotation['structure']
                if annotation['mode'] == cds_feature
            )
        )

    def _update_genes_cds_coordinates(
            self,
            cds_feature: str = 'CDS'
    ):
        """
        Updates the cds coordinates of the genes in the expansion container.
        """
        for gene in self:
            gene.gene_cds_coordinates = self._fetch_gene_cdss_set(
                gene_id=gene.id,
                cds_feature=cds_feature
            )

    def parse_gene_hierarchy_dictionary(
            self,
            gene_hierarchy_dictionary_path: Path,
            cds_feature: str = 'CDS'
    ):
        if Path.exists(Path(gene_hierarchy_dictionary_path)):
            self.gene_hierarchy_dictionary = self._read_pkl_file(
                file_path=gene_hierarchy_dictionary_path
            )
            self._update_genes_cds_coordinates(
                cds_feature=cds_feature
            )
        else:
            print(f"Path {self.gene_hierarchy_dictionary_path} does not exist.")

    def build_expansions(self):
        """Constructs the gene expansions from the Exonize database.

        This method initializes each Gene object and populates its expansions based on data from
        the Exonize database. Each expansion consists of nodes and edges, forming a graph for each gene.

        Examples:
            >>> expansions_container.build_expansions()
            >>> print(len(expansions_container))
            18
        """
        for gene_id, non_reciprocal_matches in self._db_handler.gene_expansions_dict.items():
            self._genes[gene_id] = self._add_gene(
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


class Gene:
    """Gene class is a container for gene expansion graphs.

    Attributes:
        id (str): The unique identifier for the gene.
        coordinates (portion.Interval): The start and end coordinates of the gene on the chromosome.
        strand (str): The DNA strand ('+' or '-') on which the gene is located.
        chromosome (str): The chromosome on which the gene is located.
        expansions (dict): A dictionary where keys are expansion IDs and values are expansion objects.
    """
    def __init__(
            self,
            gene_id,
            coordinates,
            strand,
            chromosome
    ):
        """Initializes a Gene instance.

        Args:
            gene_id (str): The unique identifier for the gene.
            coordinates (portion.Interval): The start and end coordinates of the gene on the chromosome.
            strand (str): The DNA strand ('+' or '-') on which the gene is located.
            chromosome (str): The chromosome on which the gene is located.
        """
        self.id = gene_id
        self.coordinates = coordinates
        self.strand = strand
        self.chromosome = chromosome
        self.expansions = {}
        self._plot_handler = _PlotHandler()

    def __getitem__(
            self,
            expansion_id: int
    ) -> nx.Graph:
        """Retrieves the expansion graph for the specified expansion ID.

        Args:
            expansion_id (int): The ID of the expansion to retrieve.

        Returns:
            networkx.Graph: The expansion graph associated with the given expansion ID.

        Examples:
            >>> gene[1]  # Retrieves the expansion graph for expansion ID 1
        """
        return self.expansions[expansion_id].graph

    def __iter__(
            self
    ) -> iter:
        """Returns an iterator over the expansion graphs.

        Returns:
            iterator: An iterator yielding each expansion graph.

        Examples:
            >>> for graph in gene:
            ...     print(graph)
        """
        return iter(expansion.graph for expansion in self.expansions.values())

    def __repr__(
            self
    ) -> str:
        """Returns a string representation of the Gene object.

        Returns:
            str: A string describing the gene's ID and number of expansions.

        Examples:
            >>> repr(gene)
            '<Gene GENE123 with 0 expansions (iterable of expansion graphs)>'
        """
        return f"<Gene {self.id} with {len(self.expansions)} expansions (iterable of expansion graphs)>"

    def __len__(
            self
    ) -> int:
        """Returns the number of expansions associated with the gene.

        Returns:
            int: The number of expansions.

        Examples:
            >>> len(gene)
            0
        """
        return len(self.expansions)

    def build_gene_graph(
            self
    ) -> nx.Graph:
        """
        Builds and returns a consolidated gene graph containing nodes and edges from all expansion graphs.

        Returns:
            networkx.Graph: A combined graph with nodes and edges from all expansions.

        Examples:
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
            connect_overlapping_nodes: bool = False,
            color_tandem_pair_edges: bool = False,
            full_expansion: bool = False,
            tandem_edges_color: str = 'blue'
    ) -> None:
        """Draws a multi-graph of gene expansions.

        Args:
            expansion_id (int, optional): The ID of a specific expansion to draw. If None, the  gene graph is drawn.
            figure_path (Path, optional): The path to save the figure. If None, the figure is not saved.
            figure_size (tuple of float, optional): The size of the figure in inches. Default is (8.0, 8.0).
            legend (bool, optional): Whether to display a legend on the plot. Default is True.
            connect_overlapping_nodes (bool, optional): Whether to draw edges connecting overlapping nodes in the graph.
             Default is True.
            color_tandem_pair_edges (bool, optional): Color edges between tandem exon nodes. Default is True.
            full_expansion (bool, optional): Whether to show the full expansion graph only. Default is False.
            tandem_edges_color (str, optional): The color to use for tandem pair edges. Default is 'blue'.
        """
        if expansion_id is not None:
            graph = self.expansions[expansion_id].graph
        else:
            graph = self.build_gene_graph()
        self._plot_handler.draw_expansions_multigraph(
            gene_start=self.coordinates.lower,
            gene_graph=graph,
            figure_path=figure_path,
            figure_size=figure_size,
            legend=legend,
            connect_overlapping_nodes=connect_overlapping_nodes,
            color_tandem_pair_edges=color_tandem_pair_edges,
            full_expansion=full_expansion,
            tandem_edges_color=tandem_edges_color
        )

    def draw_gene_structure(
            self,
            expansion_id: int = None,
            save_path: Path = None
    ):
        """
        Visualize the gene structure, highlighting coding exons and expansion events.

        This method uses the dna_features_viewer library to plot the gene structure, showing the locations
        of coding exons and expansion events within the specified gene. The plot can be saved to a file if
        `save_path` is provided.

        Args:
            expansion_id (int, optional): The ID of the expansion to visualize. Defaults to None, in which
                case the default expansion (if any) will be used.
            save_path (Path, optional): The file path where the plot will be saved. If not provided, the plot
                will be displayed but not saved.

        Returns:
            None
        """
        graph = self.expansions[expansion_id].graph
        self._plot_handler.draw_gene_structure(
            gene_coord=self.coordinates,
            strand=self.strand,
            cds_coords_list=self.gene_cds_coordinates,
            graph=graph,
            save_path=save_path
        )


class Expansion:
    """
    Expansion class represents an expansion graph for a specific gene expansion.

    Attributes:
        graph (networkx.Graph): A NetworkX graph representing the expansion.
    """
    def __init__(
            self,
            expansion_id: int,
            nodes: list[tuple],
            edges: list[tuple]
    ):
        """Initializes an Expansion instance.

        Args:
            expansion_id (int): The unique identifier for the expansion.
            nodes (list of tuples): A list of tuples representing the nodes in the form (coord, node_type).
            edges (list of tuples): A list of tuples representing the edges in the form (q_coord, t_coord, mode).
        """
        self.graph = nx.Graph()
        self.graph.id = expansion_id
        self.graph.add_nodes_from([
            (coord, {"mode": node_type})
            for coord, node_type in nodes
        ])
        for edge in edges:
            query_coord, target_coord, mode, tandem = edge
            self.graph.add_edge(
                u_of_edge=query_coord,
                v_of_edge=target_coord,
                mode=mode,
                tandem=tandem,
                width=1 if not tandem else 2,
            )


class _PlotHandler:
    """_PlotHandler is a class for managing and visualizing gene expansion graphs
    with various layout and styling options.

    Attributes
    ----------
    _color_map : dict
        A mapping of node types to specific colors for visualization.
    """
    def __init__(self):
        """
        Initializes the _PlotHandler instance with predefined node types and color mappings.
        """
        self._full = 'FULL'
        self._partial_insertion = 'PARTIAL_INSERTION'
        self._partial_excision = 'PARTIAL_EXCISION'
        self._inter_boundary = 'INTER_BOUNDARY'
        self._intronic = 'INTRONIC'
        self._color_map = {
            self._partial_insertion: '#0072B2',  # dark blue
            self._partial_excision: '#0072B2',  # burnt orange
            self._full: '#009E73',  # teal green
            self._intronic: '#D55E00',  # reddish pink
            self._inter_boundary: '#E69F00'  # golden yellow
        }

    @staticmethod
    def large_component_position(
            component: list,
            layout_scale: int = 3
    ):
        # Sort nodes within the component by their 'lower' and 'upper' coordinates
        layout = {}
        sorted_nodes = sorted(component, key=lambda coord: (coord.lower, coord.upper))
        n = len(sorted_nodes)
        for i, node in enumerate(sorted_nodes):
            # Circular layout for larger components
            angle = 2 * np.pi * i / n
            layout[node] = (layout_scale * np.cos(angle), layout_scale * np.sin(angle))
        return layout

    @staticmethod
    def _separate_large_small_components(components):
        """Separates large and small components based on their size."""
        large_components = [comp for comp in components if len(comp) > 2]
        small_components = [comp for comp in components if len(comp) <= 2]
        return large_components, small_components

    def _position_large_components(
            self,
            large_components: list[list],
            layout_scale: int = 3
    ) -> list[dict]:
        """Generates positions for large components with layout adjustments."""
        return [self.large_component_position(
            component=component,
            layout_scale=layout_scale
        ) for component in large_components]

    @staticmethod
    def _handle_two_node_case(
            gene_graph: nx.MultiGraph,
    ) -> dict:
        """Handles the special case when the graph has only two nodes."""
        nodex, nodey = list(gene_graph.nodes)
        return {nodex: (0, 0), nodey: (0.2, 0)}

    @staticmethod
    def _place_large_components(
            component_positions,
            position_shift=6
    ):
        """Places large components with horizontal shifts."""
        component_position = {}
        for event_idx, layout in enumerate(component_positions):
            for node, position in layout.items():
                x, y = position
                component_position[node] = (x + event_idx * position_shift, y)
        return component_position

    @staticmethod
    def _get_y_range(
            component_position: dict
    ):
        """Calculates the vertical range (y_min, y_max) of positioned large components."""
        y_values = [y for _, y in component_position.values()]
        return (min(y_values), max(y_values)) if y_values else (0, 0)

    @staticmethod
    def _calculate_small_component_x(
            component_position: dict,
            position_shift: int = 6
    ):
        """Determines the starting x-position for small components."""
        last_large_x = max((x for x, _ in component_position.values()), default=0)
        return last_large_x + position_shift

    @staticmethod
    def _position_small_components(
            small_components: list[list],
            x_position: int,
            y_min: int,
            y_max: int
    ):
        """Positions small components vertically at a fixed x position with even spacing."""
        num_small_components = len(small_components)
        spacing = (y_max - y_min) / (num_small_components - 1) if num_small_components > 1 else 0
        small_component_positions = {}

        for idx, component in enumerate(small_components):
            y_position = y_max - idx * spacing
            if len(component) == 2:
                # Add a small horizontal offset for two-node components to avoid overlap
                x, y = component
                small_component_positions[x] = (x_position - 1, y_position)
                small_component_positions[y] = (x_position + 1, y_position)
            else:
                for node in component:
                    small_component_positions[node] = (x_position, y_position)

        return small_component_positions

    def component_positions(
            self,
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
        large_components, small_components = self._separate_large_small_components(components)
        component_positions = self._position_large_components(
            large_components=large_components,
            layout_scale=layout_scale
        )
        if len(gene_graph.nodes) == 2:
            return self._handle_two_node_case(
                gene_graph=gene_graph
            )
        component_position = self._place_large_components(
            component_positions=component_positions,
            position_shift=layout_scale * 3)
        y_min, y_max = self._get_y_range(component_position)
        small_component_x = self._calculate_small_component_x(component_position)

        small_component_positions = self._position_small_components(
            small_components, small_component_x, y_min, y_max
        )
        component_position.update(small_component_positions)

        return component_position

    @staticmethod
    def full_expansion(
            graph: nx.Graph,
    ) -> None:
        """
        Removes nodes from the graph that are not labled as 'FULL' mode.

        Parameters
        ----------
        graph : nx.Graph

        """
        nodes_to_drop = [
            node
            for node in graph.nodes
            if graph.nodes[node].get('mode') != 'FULL'
        ]
        if nodes_to_drop:
            graph.remove_nodes_from(nodes_to_drop)

    def draw_gene_structure(
            self,
            gene_coord,
            strand,
            cds_coords_list,
            graph,
            save_path=None
    ):
        features = [
            GraphicFeature(
                start=coord.lower - gene_coord.lower,
                end=coord.upper - gene_coord.lower,
                strand=strand
            )
            for coord in cds_coords_list
        ]
        _, ax = plt.subplots(figsize=(9, 2))
        record = GraphicRecord(
            first_index=0,
            sequence_length=gene_coord.upper - gene_coord.lower,
            features=features
        )
        record.plot(ax=ax)
        for coord, attribute in graph.nodes(data=True):
            ax.fill_between(
                (coord.lower - gene_coord.lower, coord.upper - gene_coord.lower),
                +1,
                -10,
                alpha=.4,
                color=self._color_map[attribute['mode']]
            )
        if save_path:
            plt.savefig(save_path, bbox_inches='tight')
        else:
            plt.show()

    def draw_expansions_multigraph(
            self,
            gene_start: int,
            gene_graph: [nx.Graph, nx.MultiGraph],
            figure_size: tuple[float, float] = (8.0, 8.0),
            figure_path: Path = None,
            legend: bool = True,
            color_tandem_pair_edges: bool = False,
            connect_overlapping_nodes: bool = False,
            full_expansion: bool = False,
            tandem_edges_color: str = 'blue'
    ):
        """
        Draws a multi-graph of gene expansions.
        """
        graph = gene_graph.copy()
        if full_expansion:
            self.full_expansion(graph=graph)

        if graph.number_of_nodes() <= 1:
            return

        plt.figure(figsize=figure_size)

        node_colors = self._get_node_colors(graph)
        components = list(nx.connected_components(graph))
        node_labels = self._generate_node_labels(graph, gene_start)
        component_position = self.component_positions(components=components, gene_graph=graph)
        label_positions = self._adjust_label_positions(component_position)

        self._draw_nodes(graph, component_position, node_colors)
        self._draw_labels(graph, label_positions, node_labels)
        self._draw_edges(graph, component_position, color_tandem_pair_edges, tandem_edges_color)

        if connect_overlapping_nodes:
            self._draw_overlapping_edges(graph, component_position)

        if legend:
            self._draw_legend(graph, color_tandem_pair_edges, tandem_edges_color)

        self._finalize_figure(figure_path)

    def _get_node_colors(
            self,
            graph: nx.Graph
    ):
        """Get colors for each node based on their mode."""
        return [self._color_map[attrib['mode']] for _, attrib in graph.nodes(data=True)]

    @staticmethod
    def _generate_node_labels(
            graph: nx.Graph,
            gene_start: int
    ):
        """Generate labels for each node based on their position relative to gene start."""
        return {node: f'({node.lower - gene_start},{node.upper - gene_start})' for node in graph.nodes}

    @staticmethod
    def _adjust_label_positions(
            component_position: dict
    ):
        """Adjust label positions slightly above each node."""
        return {node: (x, y + 0.1) for node, (x, y) in component_position.items()}

    @staticmethod
    def _draw_nodes(
            graph: nx.Graph,
            component_position: dict,
            node_colors: list
    ) -> None:
        """Draws the nodes of the graph."""
        nx.draw_networkx_nodes(G=graph, pos=component_position, node_color=node_colors, node_size=350)

    @staticmethod
    def _draw_labels(
            graph: nx.Graph,
            label_positions: dict,
            node_labels: dict
    ):
        """Draws labels for the nodes in the graph."""
        nx.draw_networkx_labels(
            graph, label_positions, labels=node_labels, font_size=8,
            bbox=dict(boxstyle="round,pad=0.3", edgecolor="white", facecolor="white")
        )

    @staticmethod
    def _draw_edges(
            graph: nx.Graph,
            component_position: dict,
            color_tandem_pair_edges: bool,
            tandem_edges_color: str
    ):
        """Draws the edges of the graph with specified styles and colors."""
        for source, target, attributes in graph.edges(data=True):
            edge_color = tandem_edges_color if color_tandem_pair_edges and attributes.get('tandem') else 'black'
            edge_style = attributes.get('style', 'solid')
            edge_width = attributes.get('width', 1)
            nx.draw_networkx_edges(
                graph, component_position, edgelist=[(source, target)],
                edge_color=edge_color, style=edge_style, width=edge_width
            )

    def _draw_overlapping_edges(
            self,
            graph: nx.Graph,
            component_position: dict
    ):
        """Draws edges between overlapping nodes with a dotted red style."""
        overlapping_nodes = self.fetch_overlapping_nodes(
            gene_graph=graph
        )
        for cluster in overlapping_nodes:
            for nodei, nodej in cluster:
                nx.draw_networkx_edges(
                    graph, component_position, edgelist=[(nodei, nodej)],
                    edge_color='red', style='dotted', width=2
                )

    def _draw_legend(
            self,
            graph: nx.Graph,
            color_tandem_pair_edges: bool,
            tandem_edges_color: str
    ):
        """Draws the legend for the graph, showing node types and optional tandem edges."""
        node_attributes = {attrib['mode'] for _, attrib in graph.nodes(data=True)}
        legend_elements = [
            mlines.Line2D([], [], color=self._color_map[label], marker='o', linestyle='None', markersize=10,
                          label=label)
            for label in node_attributes
        ]
        if color_tandem_pair_edges:
            legend_elements.append(
                mlines.Line2D([], [], color=tandem_edges_color, linewidth=2, label='TANDEM')
            )
        plt.legend(handles=legend_elements, loc="upper right", bbox_to_anchor=(1.2, 1), frameon=False)

    @staticmethod
    def _finalize_figure(
            figure_path: Path = None
    ):
        """Finalizes the figure by saving or displaying it, and hiding spines."""
        for spine in plt.gca().spines.values():
            spine.set_visible(False)
        if figure_path:
            plt.savefig(figure_path, bbox_inches='tight')
        else:
            plt.show()

    def fetch_overlapping_nodes(
            self,
            gene_graph: nx.Graph
    ):
        """
        Identifies and returns overlapping nodes in the gene graph for visualization.

        Parameters
        ----------
        gene_graph : nx.Graph
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


class _ExonizeDBHandler:
    """
        _ExonizeDBHandler class manages database interactions for gene and expansion data with exon duplications.

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
        Initializes an _ExonizeDBHandler instance, fetching genes and expansions from the database.

        Parameters
        ----------
        db_path : str
            Path to the SQLite database file.
        """
        self.db_path = db_path
        self.genes_dict = self.collect_genes(self.fetch_genes_with_exon_dup())
        self.tandem_pairs_dict = self.fetch_tandem_pairs()
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
                WHERE mode='table' AND name='{table_name}';
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

    def fetch_tandem_pairs(
            self,
    ) -> dict:
        exonize_tandem_pairs = {}
        with sqlite3.connect(self.db_path) as db:
            cursor = db.cursor()
            cursor.execute(
                """
                SELECT
                    GeneID,
                    PredecessorStart,
                    PredecessorEnd,
                    SuccessorStart,
                    SuccessorEnd
                FROM Expansions_full_tandem
                WHERE TandemPair=1;
                """
            )
            tandem_pairs = cursor.fetchall()
        for record in tandem_pairs:
            gene_id, pred_s, pred_e, succ_s, succ_e = record
            if gene_id not in exonize_tandem_pairs:
                exonize_tandem_pairs[gene_id] = set()
            pair = tuple(sorted((P.open(pred_s, pred_e), P.open(succ_s, succ_e)), key=lambda x: (x.lower, x.upper)))
            exonize_tandem_pairs[gene_id].add(pair)
        return exonize_tandem_pairs

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

    def check_for_tandemness(
            self,
            sorted_pair: tuple,
            gene_id: str
    ) -> bool:
        if gene_id in self.tandem_pairs_dict:
            if sorted_pair in self.tandem_pairs_dict[gene_id]:
                return True
        return False

    def collect_expansions_edges(
            self,
            matches_set: set,
    ):
        for match in matches_set:
            gene_id, q_start, q_end, t_start, t_end, mode = match
            for expansion_id, expansion_atrributes in self.gene_expansions_dict[gene_id].items():
                node_coordinates = {coordinate for coordinate, _ in expansion_atrributes['nodes']}
                sorted_pair = tuple(
                    sorted((P.open(q_start, q_end), P.open(t_start, t_end)), key=lambda x: (x.lower, x.upper))
                )
                if all(coord in node_coordinates for coord in sorted_pair):
                    is_tandem = self.check_for_tandemness(sorted_pair=sorted_pair, gene_id=gene_id)
                    expansion_atrributes['edges'].append(
                        (P.open(q_start, q_end), P.open(t_start, t_end), mode, is_tandem)
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
