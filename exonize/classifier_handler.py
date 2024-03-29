# ------------------------------------------------------------------------
# This module contains the ClassifierHandler class.
# The ClassifierHandler class is a class that contains the methods used to
# classify tblastx hits as full-length duplications.
# ------------------------------------------------------------------------
import sys
import random
import re
import tempfile
import portion as P
from typing import Union


class ClassifierHandler(object):
    def __init__(
            self,
            blast_engine: object,
            cds_overlapping_threshold: float,
    ):
        self.data_container = blast_engine.data_container
        self.blast_engine = blast_engine
        self.database_interface = blast_engine.database_interface
        self.cds_overlapping_threshold = cds_overlapping_threshold

        # Initialize other internal attributes
        self.__neither, self.__query, self.__target = 0, 0, 0
        self.__both, self.__target_full, self.__target_insertion = 0, 0, 0
        self.__annot_target_start, self.__annot_target_end, self.__target_type = None, None, None
        self.__query_cds, self.__target_cds = "-", "-"
        self.__query_cds_frame, self.__target_cds_frame = " ", " "
        self.__found = False
        self.tuples_full_length_duplications = list()
        self.tuples_truncation_events = list()
        self.tuples_obligatory_events = list()

    def initialize_variables(
            self,
    ) -> None:
        """
        initializes variables used in the identify_full_length_duplications function
        """
        self.__neither, self.__query, self.__target = 0, 0, 0
        self.__both, self.__target_full, self.__target_insertion = 0,  0, 0
        self.__annot_target_start, self.__annot_target_end, self.__target_type = None, None, None
        self.__query_cds, self.__target_cds = "-", "-"
        self.__query_cds_frame, self.__target_cds_frame = " ", " "
        self.__found = False

    def initialize_list_of_tuples(
            self,
    ) -> None:
        """
        initializes the list of tuples used to store the identified events in the
         identify_full_length_duplications function
        """
        self.tuples_full_length_duplications = list()
        self.tuples_obligatory_events = list()
        self.tuples_truncation_events = list()

    @staticmethod
    def recover_query_cds(
            transcript_dictionary: dict,
            query_coordinate: P.Interval,
    ) -> Union[tuple, None]:
        """
        Find the CDS id in the transcript dictionary that matches the given coordinate.
        :param transcript_dictionary: Dictionary containing CDS information.
        :param query_coordinate: The coordinate to search for.
        :return: The ID of the matching CDS or None if no match is found.
        """
        for cds in transcript_dictionary["structure"]:
            if cds["coordinate"] == query_coordinate and cds["type"] == "CDS":
                return cds["id"], cds["coordinate"], cds["frame"]
        return None

    def identify_query(
            self,
            transcript_dictionary: dict,
            cds_coordinate: P.Interval,
    ) -> None:
        """
        identify_query is a function that identifies the tblastx
        query CDS in the gene transcript.
        A transcript cannot have overlapping CDSs. If the query CDS overlaps
        with more than one CDS, the program exits.
        """
        query_cds = self.recover_query_cds(
            transcript_dictionary=transcript_dictionary,
            query_coordinate=cds_coordinate
        )
        # if the cds belongs to the transcript
        if query_cds:
            self.__query_cds, _, self.__query_cds_frame = query_cds
            self.__query = 1
            self.__target_type = 'QUERY_ONLY'

    def target_out_of_mrna(
            self,
            transcript_coordinate: P.Interval,
            mrna_id: str,
            row_tuple: list[tuple],
    ) -> bool:
        """
        target_out_of_mrna is a function that identifies tblastx hits
        that are outside the mRNA transcript.
        """
        (fragment_id, gene_id, gene_start,
         _, cds_start, cds_end, query_start, query_end,
         target_start, target_end, evalue
         ) = row_tuple
        target_coordinate = P.open(target_start + gene_start, target_end + gene_start)
        if (target_coordinate.upper < transcript_coordinate.lower
                or transcript_coordinate.upper < target_coordinate.lower):
            if (self.__query + self.__target) == 0:
                self.__neither = 1
                self.tuples_full_length_duplications.append((
                    fragment_id, gene_id, mrna_id,
                    cds_start, cds_end, self.__query_cds, query_start, query_end,
                    "OUT_OF_MRNA",
                    self.__target_cds, self.__annot_target_start,
                    self.__annot_target_end,
                    target_start, target_end,
                    self.__neither, self.__query, self.__target, self.__both,
                    evalue
                ))
                return True
        return False

    def find_overlapping_annotation(
            self,
            transcript_dictionary: dict,
            cds_coordinate: P.Interval,
    ) -> Union[tuple, None]:
        """
        find_overlapping_annotation is a function that given a transcript dictionary
        and a CDS interval, returns a list of tuples with the following structure:
        (CDS_id, CDS_coord, CDS_frame) for all CDSs that overlap with the query CDS
        interval.
        Note:
        - The set of representative CDS only contains information about
        the coordinates of the CDSs, not their IDs.Hence,
        if we find hits for a CDS, we need to identify the CDS in the gene transcript.
        Recall that not all CDS are part of all transcripts (e.g. alternative splicing).
        - Following the classification criteria described in
        'get_candidate_cds_coordinates', we identify the query CDS described in
        case a (get_candidate_cds_coordinates), i.e, when they are considerd to overlap.
        """
        overlaps = [
            (
                self.blast_engine.get_average_overlap_percentage(
                    intv_i=annotation['coordinate'],
                    intv_j=cds_coordinate
                ),
                (annotation['id'], annotation['coordinate'], int(annotation['frame']))
            )
            for annotation in transcript_dictionary['structure']
            if annotation['type'] == 'CDS'
            and self.blast_engine.min_perc_overlap(
                intv_i=annotation['coordinate'],
                intv_j=cds_coordinate
            ) > self.cds_overlapping_threshold
        ]

        if not overlaps:
            return None

        # Find the tuple with the maximum average overlap
        _, best_matching_cds = max(overlaps, key=lambda x: x[0])

        return best_matching_cds

    def indetify_full_target(
            self,
            transcript_dictionary: dict[dict],
            target_coordinate: P.Interval,
    ) -> None:
        """
        indetify_full_target is a function that identifies tblastx
        hits that are full-length duplications as described
        in self.get_candidate_query_cds
        """
        target_only = self.find_overlapping_annotation(
            transcript_dictionary=transcript_dictionary,
            cds_coordinate=target_coordinate
        )
        if target_only:
            self.__target_cds, target_cds_coordinate, self.__target_cds_frame = target_only
            self.__target_full = 1
            self.__found = True
            self.__target_type = "FULL"
            self.__annot_target_start = target_cds_coordinate.lower
            self.__annot_target_end = target_cds_coordinate.upper

    def indentify_insertion_target(
            self,
            transcript_dictionary: dict,
            target_coordinate: P.Interval
    ) -> None:
        """
        Identifies tblastx hits that are insertion duplications these can be:
        - INS_CDS: if the insertion is in a CDS
        - INS_UTR: if the insertion is in an untranslated region
        - DEACTIVATED: if the insertion is in an intron
        """
        overlapping_annotations = [
            (annotation['id'], annotation['coordinate'], annotation['type'])
            for annotation in transcript_dictionary['structure']
            if annotation['coordinate'].contains(target_coordinate)
        ]
        # There will only be one overlapping annotation
        # since the feature annotations that we are interested in
        # do not overlap in the same transcript
        # Determine the type of the hit based on the filtered annotations
        if overlapping_annotations:
            self.__target_cds, target_cds_coordinate, annot_type = overlapping_annotations[0]
            self.__found = True
            self.__annot_target_start = target_cds_coordinate.lower
            self.__annot_target_end = target_cds_coordinate.upper

            if annot_type == 'CDS':
                self.__target_insertion = 1
                self.__target_type = "INS_CDS"
            elif annot_type == 'UTR':
                self.__target_type = "INS_UTR"
            else:  # assuming the only other type is 'intron'
                self.__target_type = "DEACTIVATED"

    def indentify_truncation_target(
            self,
            transcript_dictionary: dict[str],
            mrna_id: str,
            row_tuple: list
    ) -> None:
        """
        Identifies tblastx hits that are truncation duplications.
        These are hits that span across more than one annotation
        in the transcript architecture. We record a line per annotation
        that is truncated.
        """
        (fragment_id, gene_id, gene_start,
         _, cds_start, cds_end, query_start, query_end,
         target_start, target_end, evalue) = row_tuple
        target_coordinate = P.open(target_start + gene_start, target_end + gene_start)
        coordinate_dictionary = self.get_interval_dictionary(
            transcript_dictionary=transcript_dictionary['structure'],
            target_coordinate=target_coordinate,
            transcript_coordinate=transcript_dictionary['coordinate']
        )
        if coordinate_dictionary:
            self.__target_type = "TRUNC"
            self.__target_cds = None
            self.__annot_target_start = None
            self.__annot_target_end = None
            for seg_b, annotation_dictionary in coordinate_dictionary.items():
                coord_b = annotation_dictionary['coordinate']
                self.tuples_truncation_events.append(
                    (
                        fragment_id, gene_id, mrna_id,
                        transcript_dictionary['coordinate'].lower,
                        transcript_dictionary['coordinate'].upper,
                        self.__query_cds, cds_start, cds_end, query_start, query_end,
                        target_start, target_end,
                        annotation_dictionary['id'],
                        annotation_dictionary['type'],
                        coord_b.lower, coord_b.upper, seg_b.lower, seg_b.upper
                    )
                )

    def identify_obligate_pair(
            self,
            transcript_coordinate: P.Interval,
            mrna_id: str,
            row_tuple: list[tuple]
    ) -> None:
        """
        Identifies tblastx hits that are obligate pairs.
        These are hits where the query and target show as CDSs
        in the transcript in question.
        """
        (fragment_id, gene_id, _, _,
         cds_start, cds_end, query_start, query_end,
         target_start, target_end, evalue) = row_tuple
        self.__both = 1
        self.__query, self.__target = 0, 0
        self.tuples_obligatory_events.append((
            fragment_id, gene_id, mrna_id,
            transcript_coordinate.lower, transcript_coordinate.upper,
            cds_start, cds_end, self.__query_cds, self.__query_cds_frame,
            query_start, query_end, self.__target_cds, self.__target_cds_frame,
            self.__annot_target_start, self.__annot_target_end,
            target_start, target_end, self.__target_type
        ))

    def identify_neither_pair(self) -> None:
        """
        Identifies tblastx hits that are neither pairs.
        These are hits where the query and target do not show as CDSs
        in the transcript in question.
        """
        self.__neither = 1
        self.__target_type = 'NEITHER'

    def insert_full_length_duplication_tuple(
            self,
            mrna_id: str,
            row_tuple: list[tuple]
    ):
        """
        insert_full_length_duplication_tuple is a function that appends
         the "row" event to the list of tuples.
        """
        (fragment_id, gene_id, _, _,
         cds_start, cds_end, query_start, query_end,
         target_start, target_end, evalue) = row_tuple
        self.tuples_full_length_duplications.append((
            fragment_id, gene_id, mrna_id,
            cds_start, cds_end, self.__query_cds, query_start, query_end,
            self.__target_type, self.__target_cds, self.__annot_target_start,
            self.__annot_target_end,
            target_start, target_end, self.__neither, self.__query,
            self.__target,
            self.__both, evalue
        ))

    def insert_classified_tuples_in_results_database(self) -> None:
        """
        insert_tuples_in_results_database is a function that inserts the list
        of tuples collected by the identify_full_length_duplications
        function into the results database.
        """
        self.database_interface.instert_full_length_event(
            tuples_list=self.tuples_full_length_duplications
        )
        self.database_interface.instert_obligate_event(
            tuples_list=self.tuples_obligatory_events
        )
        self.database_interface.insert_truncation_event(
            tuples_list=self.tuples_truncation_events
        )

    def identify_full_length_duplications(
            self,
    ) -> None:
        """
        identify_full_length_duplications is a function that identifies full-length
        duplications following our classification model.
        The function iterates over all representative tblastx hits and for each transcript
        associated with the gene harboring the event it identifies the following events:
        - I. Full exon duplication: the match fully (following our coverage criteria)
        overlaps with a CDS annotation.
        - II. Insertion: the match is found within a larger CDS.
        - III. Deactivation or unnanotated: the match is found in an intron or UTR.
        - IV. Trunctation: the match spans more than one annotation (e.g., CDS, intron, UTR).
        The identify_full_length_duplications function also looks for:
         - I. Obligate events: defined as events where the query CDS and the target CDS are
          included within the same transcript.
        - II. Neither events: defined as events where the query CDS is not found in the
         transcript and the target figures a trunctation or deactivation.
        The identified events are stored in the results database (self.results_database) in the tables:
         - ObligatoryEvents: identifies both, query and target CDSs. A single record is stored per event.
         - TruncationEvents: identifies all covered annotations. The number of records per event
          will correspond to the number of annotations covered by the event.
         - FullLengthDuplications tables: identifies full-length duplications.
          The number of records per event will correspond to the number of transcripts associated
          with the gene harboring the event.
        """
        rows = self.database_interface.query_filtered_full_duplication_events()
        self.initialize_list_of_tuples()
        for row in rows:
            (_, gene_id, gene_start, _,
             cds_start, cds_end, _, _,
             target_start, target_end, _) = row
            cds_coordinate = P.open(cds_start, cds_end)
            target_coordinate = P.open(target_start + gene_start, target_end + gene_start)
            for mrna_id, transcript_dictionary \
                    in self.data_container.gene_hierarchy_dictionary[gene_id]['mRNAs'].items():
                # we want to classify each transcript independently
                self.initialize_variables()
                transcript_coordinate = transcript_dictionary['coordinate']
                # ####### QUERY ONLY - FULL LENGTH #######
                self.identify_query(
                    transcript_dictionary=transcript_dictionary,
                    cds_coordinate=cds_coordinate
                )
                # ###### CHECK: TARGET REGION NOT IN mRNA #######
                if self.target_out_of_mrna(
                        transcript_coordinate=transcript_coordinate,
                        mrna_id=mrna_id,
                        row_tuple=row
                ):
                    continue
                # ####### TARGET ONLY - FULL LENGTH #######
                self.indetify_full_target(
                    transcript_dictionary=transcript_dictionary,
                    target_coordinate=target_coordinate
                )
                if self.__target_full == 0:
                    # ####### INSERTION #######
                    self.indentify_insertion_target(
                        transcript_dictionary=transcript_dictionary,
                        target_coordinate=target_coordinate
                    )
                    if not self.__found:
                        # ####### TRUNCATION #######
                        self.indentify_truncation_target(
                            transcript_dictionary=transcript_dictionary,
                            mrna_id=mrna_id,
                            row_tuple=row
                        )
                self.__target = self.__target_full + self.__target_insertion
                # ####### OBLIGATE PAIR #######
                if self.__query + self.__target == 2:
                    self.identify_obligate_pair(
                        transcript_coordinate=transcript_dictionary['coordinate'],
                        mrna_id=mrna_id,
                        row_tuple=row
                    )
                # ####### NEITHER PAIR #######
                elif self.__query + self.__target == 0:
                    self.identify_neither_pair()
                self.insert_full_length_duplication_tuple(
                    mrna_id=mrna_id,
                    row_tuple=row
                )

    @staticmethod
    def get_interval_dictionary(
            transcript_dictionary: dict,
            target_coordinate: P.Interval,
            transcript_coordinate: P.Interval,
    ) -> dict:

        def sort_key_intervals_dictionary(
                intv_dict: dict
        ) -> dict:
            sorted_intervals = sorted(
                list(intv_dict.keys()),
                key=lambda item: (item.lower, item.upper)
            )
            return {coord: intv_dict[coord] for coord in sorted_intervals}

        utr_features = ['five_prime_UTR', 'three_prime_UTR']
        coordinates_dictionary = {}
        out_of_utr_coordinate = None
        if target_coordinate.lower < transcript_coordinate.lower:
            out_of_utr_coordinate = P.open(target_coordinate.lower, transcript_coordinate.lower)
        elif transcript_coordinate.upper < target_coordinate.upper:
            out_of_utr_coordinate = P.open(transcript_coordinate.upper, target_coordinate.upper)
        if out_of_utr_coordinate:
            coordinates_dictionary[out_of_utr_coordinate] = dict(
                id='out_of_UTR',
                type='out_of_UTR',
                coordinate=out_of_utr_coordinate
            )
        for annotation in transcript_dictionary:
            if (annotation['coordinate'].overlaps(target_coordinate)
                    and not annotation['coordinate'].contains(target_coordinate)):
                feature_coordinate = annotation['coordinate']
                annotation_dictionary = dict(
                    id=annotation['id'],
                    type=annotation['type'],
                    coordinate=annotation['coordinate']
                    )
                intersection_coordinate = feature_coordinate & target_coordinate
                if intersection_coordinate not in coordinates_dictionary:
                    coordinates_dictionary[intersection_coordinate] = annotation_dictionary
                elif coordinates_dictionary[intersection_coordinate]['type'] in utr_features:
                    continue
                elif (coordinates_dictionary[intersection_coordinate]['type'] == 'CDS'
                      and annotation['type'] == 'exon'):
                    continue
                elif (coordinates_dictionary[intersection_coordinate]['type'] == 'intron'
                      and annotation['type'] in utr_features):
                    coordinates_dictionary[intersection_coordinate] = annotation_dictionary
                elif (coordinates_dictionary[intersection_coordinate]['type'] == 'exon'
                      and annotation['type'] in utr_features):
                    coordinates_dictionary[intersection_coordinate] = annotation_dictionary
                elif (coordinates_dictionary[intersection_coordinate]['type'] == 'exon'
                      and annotation['type'] == 'CDS'):
                    coordinates_dictionary[intersection_coordinate] = annotation_dictionary
                else:
                    print('check here')
        return sort_key_intervals_dictionary(coordinates_dictionary)
