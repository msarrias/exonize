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
            blast_engine,
            cds_overlapping_threshold
    ):
        self.data_container = blast_engine.data_container
        self.database_interface = blast_engine.database_interface
        self.cds_overlapping_threshold = cds_overlapping_threshold

        # Initialize other internal attributes
        self.__neither, self.__query, self.__target = 0, 0, 0
        self.__both, self.__target_full, self.__target_insertion = 0, 0, 0
        self.__annot_target_start, self.__annot_target_end, self.__target_type = None, None, None
        self.__query_cds, self.__target_cds = "-", "-"
        self.__query_cds_frame, self.__target_cds_frame = " ", " "
        self.__found = False
        self.__tuples_full_length_duplications = list()
        self.__tuples_insertion_duplications = list()
        self.__tuples_truncation_events = list()
        self.__tuples_obligatory_events = list()

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
        self.__tuples_full_length_duplications = list()
        self.__tuples_obligatory_events = list()
        self.__tuples_truncation_events = list()

    @staticmethod
    def get_overlap_percentage(
            intv_i: P.Interval,
            intv_j: P.Interval,
    ) -> float:
        """
        Given two intervals, the function
        get_overlap_percentage returns the percentage
        of the overlapping region relative to an interval b.
        """
        intersection = intv_i & intv_j
        if intersection:
            return (intersection.upper - intersection.lower) / (intv_j.upper - intv_j.lower)
        return 0

    def find_overlapping_annot(
            self,
            transcript_dictionary: dict,
            cds_coordinate: P.Interval,
    ) -> list:
        """
        find_overlapping_annot is a function that given a transcript dictionary
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
        return [(annotation['id'], annotation['coordinate'], int(annotation['frame']))
                for annotation in transcript_dictionary['structure']
                if annotation['type'] == 'CDS'
                and all(
                [
                    self.get_overlap_percentage(
                        annotation['coordinate'],
                        cds_coordinate
                    ) >= self.cds_overlapping_threshold,
                    self.get_overlap_percentage(
                        cds_coordinate,
                        annotation['coordinate']
                    ) >= self.cds_overlapping_threshold
                ])
                ]

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
        query_only = self.find_overlapping_annot(
            transcript_dictionary=transcript_dictionary,
            cds_coordinate=cds_coordinate
        )
        if len(query_only) > 1:
            self.data_container.log.logger.error(
                f'Overlapping query CDSs: {query_only}, '
                f'please review your GFF3 file'
            )
            sys.exit()
        elif query_only:
            self.__query_cds, _, self.__query_cds_frame = query_only[0]
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
                self.__tuples_full_length_duplications.append((
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

    def indetify_full_target(
            self,
            transcript_dictionary: dict[dict],
            target_coordinate: P.Interval,
    ) -> None:
        """
        indetify_full_target is a function that identifies tblastx
        hits that are full-length duplications as described
        in self.find_overlapping_annot
        """
        target_only = self.find_overlapping_annot(
            transcript_dictionary=transcript_dictionary,
            cds_coordinate=target_coordinate
        )
        if len(target_only) > 1:
            self.data_container.log.logger.error(
                f'overlapping target CDSs: {target_only}'
            )
            sys.exit()
        elif target_only:
            self.__target_cds, target_cds_coordinate, self.__target_cds_frame = target_only[0]
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
        - INS_UTR: if the insertion is in an untralated region
        - DEACTIVATED: if the insertion is in an intron
        """
        def filter_structure_by_interval_and_type(
                structure: dict,
                t_intv_: P.Interval,
                annot_type: str
        ) -> list:
            return [(i['id'], i['coordinate']) for i in structure
                    if (i['coordinate'].contains(t_intv_)
                    and annot_type in i['type'])
                    ]

        insertion_CDS_ = filter_structure_by_interval_and_type(
            structure=transcript_dictionary['structure'],
            t_intv_=target_coordinate,
            annot_type='CDS'
        )
        if insertion_CDS_:
            self.__target_cds, t_CDS_coord_ = insertion_CDS_[0]
            self.__target_insertion = 1
            self.__target_type = "INS_CDS"
            self.__found = True
            self.__annot_target_start, self.__annot_target_end = t_CDS_coord_.lower, t_CDS_coord_.upper
        else:
            insertion_UTR_ = filter_structure_by_interval_and_type(
                structure=transcript_dictionary['structure'],
                t_intv_=target_coordinate,
                annot_type='UTR'
            )
            if insertion_UTR_:
                self.__target_cds, t_CDS_coord_ = insertion_UTR_[0]
                self.__target_type = "INS_UTR"
                self.__found = True
                self.__annot_target_start, self.__annot_target_end = t_CDS_coord_.lower, t_CDS_coord_.upper
            else:
                insertion_intron_ = filter_structure_by_interval_and_type(
                    structure=transcript_dictionary['structure'],
                    t_intv_=target_coordinate,
                    annot_type='intron'
                )
                if insertion_intron_:
                    self.__target_cds, t_CDS_coord_ = insertion_intron_[0]
                    self.__target_type = "DEACTIVATED"
                    self.__found = True
                    self.__annot_target_start, self.__annot_target_end = t_CDS_coord_.lower, t_CDS_coord_.upper

    def indentify_truncation_target(
            self,
            transcript_dictionary: dict[str],
            mrna_id: str,
            row_tuple: list
    ) -> None:
        """
        Identifies tblastx hits that are truncation duplications.
        These are hits that span across more than one annotation
        in the transcript architecture. We record a line per annotation that is truncated.
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
            for seg_b, value in coordinate_dictionary.items():
                coord_b = value['coordinate']
                self.__tuples_truncation_events.append((
                    fragment_id, gene_id, mrna_id,
                    transcript_dictionary['coordinate'].lower, transcript_dictionary['coordinate'].upper,
                    self.__query_cds, cds_start, cds_end, query_start, query_end,
                    target_start, target_end, value['id'], value['type'],
                    coord_b.lower, coord_b.upper, seg_b.lower, seg_b.upper)
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
        self.__tuples_obligatory_events.append((
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
        self.__tuples_full_length_duplications.append((
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
            tuples_list=self.__tuples_full_length_duplications
        )
        self.database_interface.instert_obligatory_event(
            tuples_list=self.__tuples_obligatory_events
        )
        self.database_interface.instert_truncation_event(
            tuples_list=self.__tuples_truncation_events
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
        self.insert_classified_tuples_in_results_database()

    @staticmethod
    def get_interval_dictionary(
            transcript_dictionary: dict,
            target_coordinate,
            transcript_coordinate
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
