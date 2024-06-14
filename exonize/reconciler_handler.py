# ------------------------------------------------------------------------
# This module contains the CounterHandler class, which is used to handle
# the counter object in the BlastEngine class.
# ------------------------------------------------------------------------
import networkx as nx
from collections import defaultdict
import portion as P
import matplotlib.pyplot as plt
from pathlib import Path


class ReconcilerHandler(object):
    def __init__(
            self,
            blast_engine: object,
            cds_overlapping_threshold: float,
    ):
        self.environment = blast_engine.environment
        self.data_container = blast_engine.data_container
        self.database_interface = blast_engine.database_interface
        self.blast_engine = blast_engine
        self.cds_overlapping_threshold = cds_overlapping_threshold

    @staticmethod
    def compute_average(
            a_list: list
    ) -> float:
        return sum(a_list) / len(a_list)

    def get_candidate_cds_reference(
            self,
            cds_coordinates_list: list[P.Interval],
            overlapping_coordinates_list: list[tuple[P.Interval, float]]
    ) -> P.Interval:
        """
        get_candidate_cds_reference is a function that given a list
        of overlapping target coordinates, returns the best candidate CDS reference

        Example:
         >>> cds_list = [P.open(0, 100), P.open(200, 300)]
         >>> overlapping_list = [(P.open(5, 100), 0.9), (P.open(10, 100), 0.9)]
         >>> get_candidate_cds_reference(cds_list, overlapping_list)
         >>> P.open(0, 100)

        :param cds_coordinates_list: list of CDS coordinates across all transcripts,
         sorted by the lower bound. The intervals can overlap.
        :param overlapping_coordinates_list: list of overlapping target coordinates

        """
        # We want to find the CDS that overlaps with the overlapping targets by
        # the highest percentage. We will compute the average of the two-way overlapping
        # percentage between the CDS and the overlapping targets and take the
        # CDS with the highest average.
        cand_reference_list = [
            (cds_coordinate,
             self.compute_average(
                 a_list=[
                     self.blast_engine.get_average_overlap_percentage(
                         intv_i=cds_coordinate,
                         intv_j=target_coordinate
                     )
                     for target_coordinate, _ in overlapping_coordinates_list]
             )
             ) for cds_coordinate in cds_coordinates_list
            # the two-way overlapping percentage should be higher than the threshold
            # maybe this is a too strict condition, consider using the average instead
            if all(
                [self.blast_engine.min_perc_overlap(
                    intv_i=cds_coordinate,
                    intv_j=target_coordinate) > self.cds_overlapping_threshold
                 for target_coordinate, _ in overlapping_coordinates_list]
            )
        ]
        if cand_reference_list:
            # Since we are considering all CDSs across all transcripts,
            # we might end up with more than one CDS candidate.
            return max(cand_reference_list, key=lambda x: x[1])[0]
        return P.open(0, 0)

    def get_overlapping_clusters(
            self,
            target_coordinates_set: set[tuple[P.Interval, float]],
            threshold: float,
    ) -> list[list[tuple]]:
        """
        Get overlapping clusters of target coordinates.
        """
        processed_intervals = set()
        overlapping_clusters = []

        # Sort the coordinates by their lower bound for efficient processing
        sorted_coordinates = sorted(target_coordinates_set, key=lambda x: x[0].lower)
        for target_coordinate, evalue in sorted_coordinates:
            if target_coordinate not in processed_intervals:
                # Create a cluster for the current coordinate
                cluster = [(target_coordinate, evalue)]
                # Find overlapping coordinates and add to the cluster
                for other_coordinate, other_evalue in sorted_coordinates:
                    if (
                            target_coordinate != other_coordinate and
                            (self.blast_engine.min_perc_overlap(
                                intv_i=target_coordinate,
                                intv_j=other_coordinate) > threshold
                             or target_coordinate.contains(other_coordinate))
                    ):
                        cluster.append((other_coordinate, other_evalue))
                        processed_intervals.add(other_coordinate)
                overlapping_clusters.append(cluster)
                processed_intervals.add(target_coordinate)
        overlapping_clusters.sort(key=len, reverse=True)
        return overlapping_clusters

    def build_reference_dictionary(
            self,
            cds_candidates_dictionary: dict,
            clusters_list: list[list]
    ) -> dict:
        # this dictionary should be for finding CDS reference and intron reference.
        # This should be applied to the clusters.
        reference_dictionary = dict()
        for coordinates_cluster in clusters_list:
            # First: let's look for targets with a high two-ways overlap percentace with a CDS
            candidate_reference = self.get_candidate_cds_reference(
                cds_coordinates_list=cds_candidates_dictionary['candidates_cds_coordinates'],
                overlapping_coordinates_list=coordinates_cluster
            )
            if candidate_reference:
                ref_type = 'FULL'
            # Second: let's look for targets contained in introns if all
            # the targets in the cluster are contained in introns
            # we take the hit with the lowest evalue as reference
            else:
                candidate_reference = min(coordinates_cluster, key=lambda x: x[1])[0] if all(
                    not target_coordinate.overlaps(cds_coordinate)  # we don't want any overlap with CDS
                    for cds_coordinate in cds_candidates_dictionary['candidates_cds_coordinates']
                    for target_coordinate, _ in coordinates_cluster
                ) else None
                ref_type = 'INACTIVE_UNANNOTATED' if candidate_reference else None
            if candidate_reference:
                for target_coordinate, _ in coordinates_cluster:
                    reference_dictionary[target_coordinate] = {
                        'reference': candidate_reference,
                        'mode': ref_type
                    }
            else:
                for target_coordinate, _ in coordinates_cluster:
                    # Process separately if no shared reference
                    individual_reference = self.get_candidate_cds_reference(
                        cds_coordinates_list=cds_candidates_dictionary['candidates_cds_coordinates'],
                        overlapping_coordinates_list=[(target_coordinate, _)]
                    )
                    ref_type = 'FULL' if individual_reference else None

                    if not ref_type:
                        contained_match_in_cds = [
                            cds_coordinate
                            for cds_coordinate in
                            cds_candidates_dictionary['candidates_cds_coordinates']
                            if cds_coordinate.contains(target_coordinate)
                        ]
                        if contained_match_in_cds:
                            ref_type = 'INSERTION_EXCISION'

                        elif all([
                            not cds_coordinate.overlaps(target_coordinate)
                            for cds_coordinate in cds_candidates_dictionary['candidates_cds_coordinates']]
                        ):
                            ref_type = 'INACTIVE_UNANNOTATED'
                        else:
                            ref_type = 'TRUNCATION_ACQUISITION'
                        individual_reference = target_coordinate

                    reference_dictionary[target_coordinate] = {
                        'reference': individual_reference,
                        'mode': ref_type
                    }
        return reference_dictionary

    @staticmethod
    def create_events_multigraph(
            reference_coordinates_dictionary: dict,
            query_coordinates_set: set,
            tblastx_records_set: set,
    ) -> nx.MultiGraph:
        gene_graph = nx.MultiGraph()
        target_coordinates_set = set([
            (reference['reference'], reference['mode'])
            for reference in reference_coordinates_dictionary.values()
        ])

        set_of_nodes = set([
            ((node_coordinate.lower, node_coordinate.upper), coordinate_type)
            for node_coordinate, coordinate_type in [
                *[(coordinate, 'FULL') for coordinate in query_coordinates_set],
                *target_coordinates_set
            ]
        ])
        gene_graph.add_nodes_from(
            [node_coordinate for node_coordinate, _ in set_of_nodes]
        )
        for node in set_of_nodes:
            node_coordinate, coordinate_type = node
            gene_graph.nodes[node_coordinate]['type'] = coordinate_type

        for event in tblastx_records_set:
            (fragment_id, _, cds_start, cds_end, target_start, target_end, evalue) = event
            target_coordinate = P.open(target_start, target_end)  # exact target coordinates
            # we take the "reference target coordinates"
            reference_coordinate = reference_coordinates_dictionary[target_coordinate]['reference']
            mode = reference_coordinates_dictionary[target_coordinate]['mode']
            gene_graph.add_edge(
                u_for_edge=(cds_start, cds_end),
                v_for_edge=(reference_coordinate.lower, reference_coordinate.upper),
                fragment_id=fragment_id,
                target=(target_start, target_end),
                evalue=evalue,
                mode=mode,
                color='black',
                width=2
            )
        return gene_graph

    @staticmethod
    def draw_event_multigraph(
            gene_graph: nx.MultiGraph,
            figure_path: Path,
    ):
        color_map = {
            'INSERTION_EXCISION': 'blue',
            'FULL': 'green',
            'INACTIVE_UNANNOTATED': 'red',
            'TRUNCATION_ACQUISITION': 'orange'
        }

        plt.figure(figsize=(16, 8))
        node_colors = [
            color_map[node[1]['type']]
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

        # Draw the graph with edge attributes
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
        # Draw edges with different styles and colors
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
        plt.savefig(figure_path)
        plt.close()

    @staticmethod
    def build_mode_dictionary(
            reference_coordinates_dictionary: dict
    ) -> dict:
        return {
            reference_coordinate['reference']: reference_coordinate['mode']
            for reference_coordinate in reference_coordinates_dictionary.values()
        }

    @staticmethod
    def build_event_coordinates_dictionary(
            component: set[tuple],
            mode_dict: dict,
            gene_graph: nx.MultiGraph
    ) -> dict:
        event_coord_dict = {}
        for start, end in component:
            node_coord = P.open(int(start), int(end))
            mode = mode_dict.get(node_coord, "FULL")  # queries are not included in the mode dictionary
            degree = gene_graph.degree((start, end))
            event_coord_dict[node_coord] = [mode, degree, None]  # cluster_id is None initially
        return event_coord_dict

    def assign_cluster_ids_to_event_coordinates(
            self,
            event_coordinates_dictionary: dict
    ) -> None:
        # Convert event coordinates to a set of tuples for clustering
        node_coordinates_set = set((node_coordinate, 0)
                                   for node_coordinate in event_coordinates_dictionary.keys()
                                   )
        # get clusters of overlapping nodes
        node_coordinates_clusters = [
            [node_coordinate for node_coordinate, _ in cluster]
            for cluster in self.get_overlapping_clusters(
                target_coordinates_set=node_coordinates_set,
                threshold=0  # we want to get all overlaps regardless of the percentage
            )
            if len(cluster) > 1
        ]
        cluster_id = 0
        # Assign cluster id to each cluster within the component
        for cluster in node_coordinates_clusters:
            for node_coordinate in cluster:
                event_coordinates_dictionary[node_coordinate][2] = cluster_id
            cluster_id += 1

    @staticmethod
    def map_edges_to_records(
            graph: nx.MultiGraph,
            expansion_id_counter: int,
            component: set[tuple]
    ) -> list[tuple]:
        comp_event_list = []
        subgraph = graph.subgraph(component)
        # in here we are mapping the event id to each fragment id (event)
        # where each fragment represents a BLAST hit between the query and the target
        # hence, in the graph each fragment is represented by an edge
        for edge in subgraph.edges(data=True):
            # edge is a tuple (node1, node2, attributes)
            comp_event_list.append(
                (expansion_id_counter, edge[-1]['fragment_id'])
            )
        return comp_event_list

    @staticmethod
    def build_events_list(
            gene_id: str,
            event_coordinates_dictionary: dict,
            expansion_id_counter: int,
    ) -> list[tuple]:
        return [
            (gene_id,
             mode,
             node_coordinate.lower,
             node_coordinate.upper,
             degree,
             cluster_id,
             expansion_id_counter
             )
            for node_coordinate, (mode, degree, cluster_id) in event_coordinates_dictionary.items()
        ]

    def get_reconciled_graph_and_expansion_events_tuples(
            self,
            reference_coordinates_dictionary: dict,
            gene_id: str,
            gene_graph: nx.MultiGraph
    ):
        mode_dictionary = self.build_mode_dictionary(
            reference_coordinates_dictionary=reference_coordinates_dictionary
        )
        # each disconnected component will represent a duplication event within a gene
        # and the number of nodes in the component is the number of duplicated exons
        # i.e., each event is described by a number of duplicated exons
        disconnected_components = list(nx.connected_components(gene_graph))
        expansion_events_list = []
        expansion_id_counter = 0
        for component in disconnected_components:
            # First: Assign event id to each component
            event_coordinates_dictionary = self.build_event_coordinates_dictionary(
                component=component,
                mode_dict=mode_dictionary,
                gene_graph=gene_graph
            )
            # Second : Assign cluster ids within events
            self.assign_cluster_ids_to_event_coordinates(
                event_coordinates_dictionary=event_coordinates_dictionary
            )
            expansion_events_list.extend(
                self.build_events_list(
                    gene_id=gene_id,
                    event_coordinates_dictionary=event_coordinates_dictionary,
                    expansion_id_counter=expansion_id_counter
                )
            )
            expansion_id_counter += 1
        return gene_graph, expansion_events_list

    @staticmethod
    def get_gene_events_dictionary(
            tblastx_full_matches_list: list[tuple],
    ):
        full_matches_dictionary = defaultdict(set)
        # group full matches by gene id
        for match in tblastx_full_matches_list:
            gene_id = match[1]
            full_matches_dictionary[gene_id].add(match)
        return full_matches_dictionary

    @staticmethod
    def get_hits_query_and_target_coordinates(
            tblastx_records_set: set,
    ):
        query_coordinates = set()
        target_coordinates = set()
        for record in tblastx_records_set:
            fragment_id, gene_id, cds_start, cds_end, target_start, target_end, evalue = record
            query_coordinates.add(P.open(cds_start, cds_end))
            target_coordinates.add((P.open(target_start, target_end), evalue))
        return query_coordinates, target_coordinates

    @staticmethod
    def center_and_sort_cds_coordinates(
            cds_coordinates: list,
            gene_start: int
    ) -> list[P.Interval]:
        return sorted(
                list(set([P.open(i.lower - gene_start, i.upper - gene_start) for i in cds_coordinates])),
                key=lambda x: (x.lower, x.upper)
        )

    def align_target_coordinates(
            self,
            gene_id: str,
            tblastx_records_set: set[tuple],
    ) -> tuple[set[P.Interval], dict]:
        gene_start = self.data_container.gene_hierarchy_dictionary[gene_id]['coordinate'].lower
        cds_candidates_dictionary = self.blast_engine.get_candidate_cds_coordinates(
            gene_id=gene_id
        )
        # center cds coordinates to gene start
        cds_candidates_dictionary['candidates_cds_coordinates'] = self.center_and_sort_cds_coordinates(
            cds_coordinates=cds_candidates_dictionary['candidates_cds_coordinates'],
            gene_start=gene_start
        )
        query_coordinates, target_coordinates = self.get_hits_query_and_target_coordinates(
            tblastx_records_set=tblastx_records_set
        )
        overlapping_targets = self.get_overlapping_clusters(
            target_coordinates_set=target_coordinates,
            threshold=self.cds_overlapping_threshold
        )
        reference_coordinates_dictionary = self.build_reference_dictionary(
            cds_candidates_dictionary=cds_candidates_dictionary,
            clusters_list=overlapping_targets
        )
        return query_coordinates, reference_coordinates_dictionary
