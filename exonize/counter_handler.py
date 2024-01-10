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

    def get_shorter_intv_overlapping_percentage(
            self,
            a: P.Interval,
            b: P.Interval
    ) -> float:
        """
        get_shorter_intv_overlapping_percentage is a function
        that given two intervals, returns the percentage of
        overlap of the shorter interval with the longer interval.
        """
        shorter, longer = self.blast_engine.get_shorter_longer_interv(a, b)
        return self.blast_engine.get_overlap_percentage(longer, shorter)

    @staticmethod
    def compute_average(
            a_list: list
    ) -> float:
        return sum(a_list) / len(a_list)

    def two_way_overlapping(
            self,
            intv_i: P.Interval,
            intv_j: P.Interval,
            threshold: float = self.cds_overlapping_threshold
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
         - cds_coordinates_list = [P.open(0, 100), P.open(200, 300)]
         - overlapping_coordinates_list = [(P.open(5, 100), 0.9), (P.open(10, 100), 0.9)]
         - get_candidate_cds_reference
         - returns P.open(0, 100)

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
                [self.two_way_overlapping(cds_coordinate, target_coordinate)
                 for target_coordinate, _ in overlapping_coordinates_list]
            )
        ]
        if cand_reference_list:
            # Since we are considering all CDSs across all transcripts,
            # we might end up with more than one CDS candidate.
            return max(cand_reference_list, key=lambda x: x[1])[0]
        return P.open(0, 0)

    @staticmethod
    def get_overlapping_clusters(
            target_coordinates_set: set[tuple[P.Interval, float]]
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
                    if (target_coordinate != other_coordinate and
                            self.two_way_overlapping(
                                intv_i=target_coordinate,
                                intv_j=other_coordinate
                            )):
                        cluster.append((other_coordinate, other_evalue))
                        processed_intervals.add(other_coordinate)

                overlapping_clusters.append(cluster)
                processed_intervals.add(target_coordinate)

        # Optionally, sort clusters by their size in descending order
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
                    for cds_coordinate in sorted_cds_coordinates_list
                    for target_coordinate, _ in coordinates_cluster
                ) else None
                ref_type = 'non_coding' if candidate_reference else None
            for target_coordinate, _ in coordinates_cluster:
                if candidate_reference:
                    reference_dictionary[target_coordinate] = {
                        'intv_ref': candidate_reference,
                        'ref': ref_type
                    }
                else:
                    # Process separately if no shared reference
                    individual_reference = self.get_candidate_cds_reference(
                        cds_coordinates_list=sorted_cds_coordinates_list,
                        overlapping_coordinates_list=[(target_coordinate, _)]
                    )
                    reference_dictionary[target_coordinate] = {
                        'intv_ref': individual_reference if individual_reference else target_coordinate,
                        'ref': 'coding' if individual_reference else 'non_coding'
                    }
        return reference_dictionary

    @staticmethod
    def create_events_multigraph(
            reference_coordinates_dictionary: dict,
            query_coords_set: set,
            records_set: set,
    ) -> nx.MultiGraph:
        gene_G = nx.MultiGraph()
        target_coords_set = set([i['intv_ref'] for i in reference_coordinates_dictionary.values()])
        gene_G.add_nodes_from(
            set([(i.lower, i.upper)
                 for i in [*query_coords_set, *target_coords_set]
                 ])
        )
        for event in records_set:
            source = (event[2], event[3])  # exact CDS coordinates
            target = reference_coordinates_dictionary[P.open(event[4], event[5])]['intv_ref']  # we need a reference
            ref_des = reference_coordinates_dictionary[P.open(event[4], event[5])]['ref']
            gene_G.add_edge(
                source,
                (target.lower, target.upper),
                fragment_id=event[0],
                query_CDS=(event[2], event[3]),
                target=(event[4], event[5]),
                evalue=event[6],
                event_type=event[7],
                ref=ref_des,
                color='black',
                width=2
            )
        return gene_G

    def get_events_tuples_from_multigraph(
            self,
            reference_coordinates_dictionary: dict,
            gene_id: str,
            gene_graph: nx.MultiGraph
    ) -> (list[tuple], set[tuple]):
        reference = {i['intv_ref']: i['ref'] for i in reference_coordinates_dictionary.values()}
        disconnected_components = list(nx.connected_components(gene_graph))
        events_tuples, events_list, event_id_counter, cluster_counter = list(), list(), 0, 0
        for component in disconnected_components:
            temp, comp_event_list = dict(), list()
            for node_x, node_y in component:
                ref = 'coding'
                node_intv = P.open(int(node_x), int(node_y))
                if node_intv in reference:
                    ref = reference[node_intv]
                temp[node_intv] = [ref,  # either coding/non_coding/coding_non_coding
                                   sum(1 for _ in gene_graph.neighbors((node_x, node_y))),  # degree
                                   None,  # cluster_id
                                   event_id_counter]
            temp_list_tuples = set((node, 0) for node in temp.keys())
            node_clusters = [[i[0] for i in cluster] for cluster in
                             self.get_overlapping_clusters(
                                 target_coordinates_set=temp_list_tuples
                             )
                             if len(cluster) > 1]
            for cluster in node_clusters:
                for node in cluster:
                    temp[node][2] = cluster_counter
                cluster_counter += 1
            subgraph = gene_graph.subgraph(component)
            for edge in subgraph.edges(data=True):
                # edge is a tuple (node1, node2, attributes)
                node1, node2, attributes = edge
                comp_event_list.append(
                    (event_id_counter, attributes['fragment_id'])
                )
            event_id_counter += 1
            events_tuples.extend(comp_event_list)
            events_list.extend([
                (gene_id, value[0], coord.lower, coord.upper, *value[1:])
                for coord, value in temp.items()])
        return events_tuples, events_list

    def assign_event_ids(
            self,
            full_matches_list: list[tuple],
    ) -> (list[tuple], set[tuple]):
        genes_events_tuples, genes_events_set = list(), set()
        full_matches_dict = defaultdict(set)

        for match in full_matches_list:
            full_matches_dict[match[1]].add(match)

        for gene_id, records in full_matches_dict.items():
            gene_start = self.data_container.gene_hierarchy_dictionary[gene_id]['coordinate'].lower
            cds_candidates = self.blast_engine.get_candidate_cds_coordinates(gene_id=gene_id)
            # center cds coordinates to gene start
            cds_candidates['candidates_cds_coordinates'] = sorted(
                list(set([
                    P.open(i.lower - gene_start, i.upper - gene_start)
                    for i in cds_candidates['candidates_cds_coordinates']])),
                key=lambda x: (x.lower, x.upper)
            )
            query_coordinates = set(
                [P.open(record[2], record[3]) for record in records]
            )
            target_coordinates = set(
                [(P.open(record[4], record[5]), record[-2]) for record in records]
            )

            overlapping_targets = self.get_overlapping_clusters(
                target_coordinates_set=target_coordinates,
                threshold=self.cds_overlapping_threshold
            )
            reference_coordinates_dictionary = self.build_reference_dictionary(
                cds_candidates_dictionary=cds_candidates,
                clusters_list=overlapping_targets
            )
            G = self.create_events_multigraph(
                reference_coordinates_dictionary=reference_coordinates_dictionary,
                query_coords_set=query_coordinates,
                records_set=records
            )
            gene_events_tuples, gene_events_set = self.get_events_tuples_from_multigraph(
                reference_coordinates_dictionary=reference_coordinates_dictionary,
                gene_id=gene_id,
                gene_graph=G
            )
            if len(gene_events_tuples) != len(records):
                self.environment.logger.exception(
                    f'{gene_id}: {len(gene_events_tuples)} events found,'
                    f' {len(records)} expected.'
                )
            genes_events_tuples.extend(gene_events_tuples)
            genes_events_set.update(gene_events_set)
        if len(genes_events_tuples) != len(full_matches_list):
            self.environment.logger.exception(
                f'{len(genes_events_tuples)} events found,'
                f' {len(full_matches_list)} expected.'
            )
        return genes_events_tuples, genes_events_set
