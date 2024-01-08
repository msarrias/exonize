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
            blast_engine: BLASTsearcher,
            cds_overlapping_threshold: float,
            ):
        self.environment = blast_engine.environment
        self.data_container = blast_engine.data_container
        self.blast_engine = blast_engine
        self.cds_overlapping_threshold = cds_overlapping_threshold

    def get_average_overlap_percentage(
            self,
            intv_a: P.Interval,
            intv_b: P.Interval
    ) -> float:
        return sum([self.blast_engine.get_overlap_percentage(intv_a, intv_b),
                    self.blast_engine.get_overlap_percentage(intv_b, intv_a)]) / 2

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

    def get_candidate_reference_dictionary(
            self,
            coordinates: list[P.Interval],
            intv_list: list[tuple[P.Interval, float]]
    ) -> P.Interval:
        cand_ref = [query_intv for query_intv in coordinates
                    if all(
                        self.get_average_overlap_percentage(query_intv, intv_i) >= self.cds_overlapping_threshold
                        for intv_i, _ in intv_list
                    )]
        if len(cand_ref) == 1:
            return cand_ref[0]
        elif cand_ref:
            cand_ref = [(cand_ref_intv,
                         sum([self.get_average_overlap_percentage(cand_ref_intv, intv_i)
                              for intv_i, _ in intv_list]) / len(intv_list))
                        for cand_ref_intv in cand_ref]
            return max(cand_ref, key=lambda x: x[1])[0]
        return P.open(0, 0)

    def overlap_condition(
            self,
            coordinate,
            x,
            threshold
    ):
        perc = self.get_shorter_intv_overlapping_percentage(a=coordinate, b=x)
        return perc >= threshold and x != coordinate

    def get_overlapping_clusters(
            self,
            coordinates_set: set[tuple[P.Interval, float]],
            threshold: float,
    ) -> list[list[tuple]]:
        overlapping_targets_list, skip_events = list(), list()
        for coord, evalue in coordinates_set:
            if coord not in skip_events:
                temp = [(i, i_evalue) for i, i_evalue in coordinates_set
                        if self.overlap_condition(coordinate=coord,
                                                  x=i,
                                                  threshold=threshold
                                                  )
                        ]
                if temp:
                    skip_events.extend([coord, *[i for i, _ in temp]])
                    flat_temp = [*temp, (coord, evalue)]
                    flat_temp.sort(key=lambda x: (x[0].lower, x[0].upper))
                    overlapping_targets_list.append(flat_temp)
                else:
                    overlapping_targets_list.append([(coord, evalue)])
        overlapping_targets_list.sort(key=len, reverse=True)
        return overlapping_targets_list

    def build_reference_dictionary(
            self,
            cds_candidates_dictionary: dict,
            overlapping_targets_list: list[list]
    ) -> dict[dict]:
        # this dictionary should be for finding CDS reference and intron reference.
        # This should be applied to the clusters.
        ref_dict = dict()
        for intv_list in overlapping_targets_list:
            # First: let's look for targets overlapping with a CDS
            sorted_CDS_coords_list = sorted(cds_candidates_dictionary['set_coords'],
                                            key=lambda x: (x.lower, x.upper))
            cand_ref = self.get_candidate_reference_dictionary(
                coordinates=sorted_CDS_coords_list,
                intv_list=intv_list
            )
            if cand_ref:
                for intv_i, _ in intv_list:
                    ref_dict[intv_i] = dict(intv_ref=cand_ref, ref='coding')
            # Second: let's look for targets overlapping with introns
            if all(not target_intv.overlaps(cds_intv)
                   for cds_intv in sorted_CDS_coords_list
                   for target_intv, _ in intv_list):
                cand_ref = min(intv_list, key=lambda x: x[1])[0]
                for intv_i, _ in intv_list:
                    ref_dict[intv_i] = dict(intv_ref=cand_ref, ref='non_coding')
            # if there is no shared reference, we take it separetly
            elif not cand_ref:
                for intv_i, i_evalue in intv_list:
                    cand_ref = self.get_candidate_reference_dictionary(
                        coordinates=sorted_CDS_coords_list,
                        intv_list=[(intv_i, i_evalue)]
                    )
                    if cand_ref:
                        ref_dict[intv_i] = dict(intv_ref=cand_ref, ref='coding')
                    else:
                        ref_dict[intv_i] = dict(intv_ref=intv_i, ref='non_coding')
        return ref_dict

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
                                 coordinates_set=temp_list_tuples,
                                 threshold=0.01
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
            cds_candidates['set_coords'] = set([
                P.open(i.lower - gene_start, i.upper - gene_start)
                for i in cds_candidates['set_coords']
            ])
            query_coordinates = set([
                P.open(record[2], record[3])
                for record in records])
            target_coordinates = set([
                (P.open(record[4], record[5]), record[-2])
                for record in records])

            overlapping_targets = self.get_overlapping_clusters(
                coordinates_set=target_coordinates,
                threshold=self.cds_overlapping_threshold
            )
            reference_coordinates_dictionary = self.build_reference_dictionary(
                cds_candidates_dictionary=cds_candidates,
                overlapping_targets_list=overlapping_targets
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
