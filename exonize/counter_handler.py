# ------------------------------------------------------------------------
# This module contains the CounterHandler class, which is used to handle
# the counter object in the BlastEngine class.
# ------------------------------------------------------------------------
import networkx as nx
from collections import defaultdict
import portion as P
import random
import re
import tempfile


class CounterHandler(object):
    def __init__(
            self,
            blast_engine: object,
            cds_overlapping_threshold: float,
            ):
        self.environment = blast_engine.environment
        self.data_container = blast_engine.data_container
        self.blast_engine = blast_engine
        self.cds_overlapping_threshold = cds_overlapping_threshold

    @staticmethod
    def compute_average(
            a_list: list
    ) -> float:
        return sum(a_list) / len(a_list)

    def two_way_overlapping(
            self,
            intv_i: P.Interval,
            intv_j: P.Interval,
            threshold: float
    ) -> bool:
        return (self.blast_engine.get_overlap_percentage(
                    intv_i=intv_i,
                    intv_j=intv_j) > threshold
                and
                self.blast_engine.get_overlap_percentage(
                    intv_i=intv_j,
                    intv_j=intv_i) > threshold
                )

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
                [self.two_way_overlapping(
                    intv_i=cds_coordinate,
                    intv_j=target_coordinate,
                    threshold=self.cds_overlapping_threshold
                )
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
                            self.two_way_overlapping(
                                intv_i=target_coordinate,
                                intv_j=other_coordinate,
                                threshold=threshold
                            )
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
    ) -> dict[dict]:
        # this dictionary should be for finding CDS reference and intron reference.
        # This should be applied to the clusters.
        reference_dictionary = dict()
        for coordinates_cluster in clusters_list:
            # First: let's look for targets overlapping with a CDS
            candidate_reference = self.get_candidate_cds_reference(
                cds_coordinates_list=cds_candidates_dictionary['candidates_cds_coordinates'],
                overlapping_coordinates_list=coordinates_cluster
            )
            if candidate_reference:
                ref_type = 'coding'
            # Second: let's look for targets overlapping with introns
            else:
                candidate_reference = min(coordinates_cluster, key=lambda x: x[1])[0] if all(
                    not target_coordinate.overlaps(cds_coordinate)  # we don't want any overlap with CDS
                    for cds_coordinate in cds_candidates_dictionary['candidates_cds_coordinates']
                    for target_coordinate, _ in coordinates_cluster
                ) else None
                ref_type = 'non_coding' if candidate_reference else None
            for target_coordinate, _ in coordinates_cluster:
                if candidate_reference:
                    reference_dictionary[target_coordinate] = {
                        'reference_coordinate': candidate_reference,
                        'reference_type': ref_type
                    }
                else:
                    # Process separately if no shared reference
                    individual_reference = self.get_candidate_cds_reference(
                        cds_coordinates_list=cds_candidates_dictionary['candidates_cds_coordinates'],
                        overlapping_coordinates_list=[(target_coordinate, _)]
                    )
                    reference_dictionary[target_coordinate] = {
                        'reference_coordinate': individual_reference if individual_reference else target_coordinate,
                        'reference_type': 'coding' if individual_reference else 'non_coding'
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
            reference['reference_coordinate']
            for reference in reference_coordinates_dictionary.values()
        ])
        gene_graph.add_nodes_from(
            set([(node_coordinate.lower, node_coordinate.upper)
                 for node_coordinate in [*query_coordinates_set, *target_coordinates_set]
                 ])
        )
        for event in tblastx_records_set:
            (fragment_id, _, source_start, source_end,
             target_start, target_end, evalue, event_type) = event[:8]
            target_coordinate = P.open(target_start, target_end)  # exact target coordinates
            reference_type = reference_coordinates_dictionary[target_coordinate]['reference_type']
            gene_graph.add_edge(
                u_for_edge=(source_start, source_end),
                v_for_edge=(target_start, target_end),
                fragment_id=fragment_id,
                query_CDS=(source_start, source_end),
                target=(target_start, target_end),
                evalue=evalue,
                event_type=event_type,
                reference_type=reference_type,
                color='black',
                width=2
            )
        return gene_graph

    @staticmethod
    def build_reference_type_dictionary(
            reference_coordinates_dictionary: dict
    ) -> dict:
        return {
            reference_coordinate['reference_coordinate']: reference_coordinate['reference_type']
            for reference_coordinate in reference_coordinates_dictionary.values()
        }

    @staticmethod
    def build_event_coordinates_dictionary(
            component: set[tuple],
            reference_type_dict: dict,
            gene_graph: nx.MultiGraph
    ) -> dict:
        event_coord_dict = {}
        for start, end in component:
            node_coord = P.open(int(start), int(end))
            ref_type = reference_type_dict.get(node_coord, 'coding')
            degree = gene_graph.degree((start, end))
            event_coord_dict[node_coord] = [ref_type, degree, None]  # cluster_id is None initially
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
            event_id_counter: int,
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
                (event_id_counter, edge[-1]['fragment_id'])
            )
        return comp_event_list

    @staticmethod
    def build_events_list(
            gene_id: str,
            event_coordinates_dictionary: dict,
            event_id_counter: int,
    ) -> list[tuple]:
        return [
            (gene_id,
             node_attributes[0],
             node_coordinate.lower,
             node_coordinate.upper,
             *node_attributes[1:],
             event_id_counter
             )
            for node_coordinate, node_attributes in event_coordinates_dictionary.items()
        ]

    def get_events_tuples_from_multigraph(
            self,
            reference_coordinates_dictionary: dict,
            gene_id: str,
            gene_graph: nx.MultiGraph
    ):
        reference_type_dictionary = self.build_reference_type_dictionary(
            reference_coordinates_dictionary=reference_coordinates_dictionary
        )
        # each disconnected component will represent a duplication event within a gene
        # and the number of nodes in the component is the number of duplicated exons
        # i.e., each event is described by a number of duplicated exons
        disconnected_components = list(nx.connected_components(gene_graph))
        gene_fragments_with_event_ids_list = []
        gene_events_list = []
        event_id_counter = 0
        for component in disconnected_components:
            # First: Assign event id to each component
            event_coordinates_dictionary = self.build_event_coordinates_dictionary(
                component=component,
                reference_type_dict=reference_type_dictionary,
                gene_graph=gene_graph
            )
            # Second : Assign cluster ids within events
            self.assign_cluster_ids_to_event_coordinates(
                event_coordinates_dictionary=event_coordinates_dictionary
            )
            # Third: map each BLAST hit to an event id
            component_fragments_with_event_ids_list = self.map_edges_to_records(
                graph=gene_graph,
                event_id_counter=event_id_counter,
                component=component
            )
            gene_fragments_with_event_ids_list.extend(
                component_fragments_with_event_ids_list
            )
            gene_events_list.extend(
                self.build_events_list(
                    gene_id=gene_id,
                    event_coordinates_dictionary=event_coordinates_dictionary,
                    event_id_counter=event_id_counter
                )
            )
            event_id_counter += 1
        return gene_fragments_with_event_ids_list, gene_events_list

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
        query_coordinates = set(
            [P.open(record[2], record[3]) for record in tblastx_records_set]
        )
        target_coordinates = set(
            [(P.open(record[4], record[5]), record[-2]) for record in tblastx_records_set]
        )
        return query_coordinates, target_coordinates

    def assign_event_ids(
            self,
            tblastx_full_matches_list: list[tuple],
    ) -> (list[tuple], set[tuple]):
        genes_events_tuples = list()
        genes_events_set = set()
        # group full matches by gene id
        full_matches_dictionary = self.get_gene_events_dictionary(
            tblastx_full_matches_list=tblastx_full_matches_list
        )
        for gene_id, tblastx_records_set in full_matches_dictionary.items():
            gene_start = self.data_container.gene_hierarchy_dictionary[gene_id]['coordinate'].lower
            cds_candidates_dictionary = self.blast_engine.get_candidate_cds_coordinates(
                gene_id=gene_id
            )
            # center cds coordinates to gene start
            cds_candidates_dictionary['candidates_cds_coordinates'] = sorted(
                list(set([
                    P.open(i.lower - gene_start, i.upper - gene_start)
                    for i in cds_candidates_dictionary['candidates_cds_coordinates']])),
                key=lambda x: (x.lower, x.upper)
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
            gene_graph = self.create_events_multigraph(
                reference_coordinates_dictionary=reference_coordinates_dictionary,
                query_coordinates_set=query_coordinates,
                tblastx_records_set=tblastx_records_set
            )
            gene_fragments_with_event_ids_list, gene_events_set = self.get_events_tuples_from_multigraph(
                reference_coordinates_dictionary=reference_coordinates_dictionary,
                gene_id=gene_id,
                gene_graph=gene_graph
            )
            if len(gene_fragments_with_event_ids_list) != len(tblastx_records_set):
                self.environment.logger.exception(
                    f'{gene_id}: {len(gene_fragments_with_event_ids_list)} events found,'
                    f' {len(tblastx_records_set)} expected.'
                )
            genes_events_tuples.extend(gene_fragments_with_event_ids_list)
            genes_events_set.update(gene_events_set)
        if len(genes_events_tuples) != len(tblastx_full_matches_list):
            self.environment.logger.exception(
                f'{len(genes_events_tuples)} events found,'
                f' {len(tblastx_full_matches_list)} expected.'
            )
        return genes_events_tuples, genes_events_set
