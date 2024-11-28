# ------------------------------------------------------------------------
# This module contains the Exonize class, which is the main class of the package.
# The Exonize class contains all the methods necessary to run the exon duplication
# search pipeline.
# The pipeline is completed with the aid of the following classes:
# - EnvironmentSetup: This class sets up the environment and creates the working
#  directory.
# - DataPreprocessor: This class contains all the methods necessary to preprocess
#  the data.
# - SqliteHandler: This class contains all the methods necessary to interact with
#  the SQLite database.
# - BLASTsearcher: This class contains all the methods necessary to perform the
# tblastx search.
# - ClassifierHandler: This class contains all the methods necessary to classify
# the tblastx hits.
# - CounterHandler: This class contains all the methods necessary to count the
# duplication events.
# The pipeline is executed by the run_exonize_pipeline method.
# ------------------------------------------------------------------------

import os
from itertools import permutations
from typing import Any, Sequence, Iterator
from datetime import datetime
import sys
from pathlib import Path
import time
import random
import portion as P
from exonize.environment_setup import EnvironmentSetup
from exonize.data_preprocessor import DataPreprocessor
from exonize.sqlite_handler import SqliteHandler
from exonize.searcher import Searcher
from exonize.classifier_handler import ClassifierHandler
from exonize.reconciler_handler import ReconcilerHandler


class Exonize(object):
    def __init__(
            self,
            gff_file_path: Path,
            genome_file_path: Path,
            pdb_ids_mapping_db_path: Path,
            pdb_structures_path: Path,
            output_prefix: str,
            output_directory_path: Path,
            cpus_number: int = os.cpu_count(),
            gene_annot_feature: str = 'gene',
            cds_annot_feature: str = 'CDS',
            transcript_annot_feature: str = 'mRNA',
            sequence_base: int = 1,
            frame_base: int = 0,
            min_exon_length: int = 30,
            evalue_threshold: float = 1e-3,
            self_hit_threshold: float = 0.5,
            query_coverage_threshold: float = 0.9,
            exon_clustering_overlap_threshold: float = 0.9,
            targets_clustering_overlap_threshold: float = 0.9,
            pair_coverage_threshold: float = 0.9,
            fraction_of_aligned_positions: float = 0.9,
            peptide_identity_threshold: float = 0.4,
            plddt_threshold: float = 0.7,
            RMSD_threshold: float = 2,
            TM_score_threshold: float = 0.5,
            csv: bool = False,
            enable_debug: bool = False,
            soft_force: bool = False,
            hard_force: bool = False,
            global_search: bool = False,
            local_search: bool = False,
            structural_search: bool = False,
            sleep_max_seconds: int = 60,
            timeout_database: int = 160
    ):
        self.tic = datetime.now()
        self.local_full_matches_dictionary = {}
        self.global_full_matches_dictionary = {}
        self.structural_full_matches_dictionary = {}

        self.environment = EnvironmentSetup(
            genome_file_path=genome_file_path,
            gff_file_path=gff_file_path,
            pdb_ids_mapping_db_path=pdb_ids_mapping_db_path,
            pdb_structures_path=pdb_structures_path,
            output_prefix=output_prefix,
            output_directory_path=output_directory_path,
            gene_annot_feature=gene_annot_feature,
            cds_annot_feature=cds_annot_feature,
            transcript_annot_feature=transcript_annot_feature,
            sequence_base=sequence_base,
            frame_base=frame_base,
            min_exon_length=min_exon_length,
            evalue_threshold=evalue_threshold,
            peptide_identity_threshold=peptide_identity_threshold,
            fraction_of_aligned_positions=fraction_of_aligned_positions,
            pair_coverage_threshold=pair_coverage_threshold,
            exon_clustering_overlap_threshold=exon_clustering_overlap_threshold,
            targets_clustering_overlap_threshold=targets_clustering_overlap_threshold,
            query_coverage_threshold=query_coverage_threshold,
            self_hit_threshold=self_hit_threshold,
            plddt_threshold=plddt_threshold,
            RMSD_threshold=RMSD_threshold,
            TM_score_threshold=TM_score_threshold,
            global_search=global_search,
            local_search=local_search,
            structural_search=structural_search,
            hard_force=hard_force,
            soft_force=soft_force,
            debug_mode=enable_debug,
            csv=csv,
            sleep_max_seconds=sleep_max_seconds,
            timeout_database=timeout_database,
            cpus_number=cpus_number
        )
        self.database_interface = SqliteHandler(
            environment=self.environment
        )
        self.data_container = DataPreprocessor(
            database_interface=self.database_interface
        )
        self.search_engine = Searcher(
            data_container=self.data_container
        )
        self.event_classifier = ClassifierHandler(
            search_engine=self.search_engine
        )
        self.event_reconciler = ReconcilerHandler(
            search_engine=self.search_engine
            )

    def generate_unique_events_list(
            self,
            events_list: list,
            event_type_idx: int
    ) -> list:
        new_events_list = []
        events_ids = []
        for event in events_list:
            mrna_concat_event_types = list(set(event[event_type_idx].rsplit(',')))
            mrna_events_perm = self.generate_combinations(strings=mrna_concat_event_types)
            keys = [i for i in mrna_events_perm if i in events_ids]
            if not keys:
                events_ids.append(mrna_events_perm[0])
                event_n = mrna_events_perm[0]
            else:
                event_n = keys[0]
            new_events_list.append((event_n, event[0]))
        return new_events_list

    @staticmethod
    def generate_combinations(
            strings: list[str]
    ) -> list[str]:
        result = set()
        for perm in permutations(strings):
            result.add('-'.join(perm))
        return list(result)

    def log_search_progress(
            self,
            unprocessed_gene_ids_list,
            local_search: bool = False,
            global_search: bool = False,
            structural_search: bool = False
    ):
        search = ''
        if local_search:
            search = 'local'
        elif global_search:
            search = 'global'
        elif structural_search:
            search = 'structural'
        gene_ids_list = list(self.data_container.gene_hierarchy_dictionary.keys())
        gene_count = len(gene_ids_list)
        out_message = (
            f'Resuming {search} search'
            if len(unprocessed_gene_ids_list) != gene_count
            else f'Starting {search} search'
        )
        self.environment.logger.info(
            f'{out_message} for {len(unprocessed_gene_ids_list)}/{gene_count} genes'
        )
        self.environment.logger.info('Exonizing: this may take a while...')

    def log_completed_search(
            self,
            local_search: bool = False,
            global_search: bool = False,
            structural_search: bool = False):
        search = ''
        if local_search:
            search = 'Local'
        elif global_search:
            search = 'Global'
        elif structural_search:
            search = 'Structural'
        self.environment.logger.info(
            f'{search} search has been completed. '
            'If you wish to re-run the analysis, '
            'consider using the hard-force/soft-force flag'
        )

    def cleanup_local_search(self):
        if self.environment.LOCAL_SEARCH:
            self.database_interface.clear_results_database(
                except_tables=[
                    'Genes',
                    'Search_monitor',
                    'Local_matches',
                    'Structural_matches_non_reciprocal'
                ]
            )
            self.data_container.initialize_database()

    def local_search_complete_identity_and_coverage(self):
        self.database_interface.insert_percent_query_column_to_fragments()
        matches_list = self.database_interface.query_raw_matches()
        identity_and_sequence_tuples = self.search_engine.get_identity_and_dna_seq_tuples(
            matches_list=matches_list
        )
        self.database_interface.insert_identity_and_dna_algns_columns(
            list_tuples=identity_and_sequence_tuples
        )

    def local_search(
            self
    ):
        unprocessed_gene_ids_list = self.database_interface.query_to_process_gene_ids(local_search=True)
        if unprocessed_gene_ids_list:
            self.log_search_progress(
                unprocessed_gene_ids_list=unprocessed_gene_ids_list,
                local_search=True
            )
            status: int
            forks: int = 0
            for balanced_batch in self.even_batches(
                    data=list(unprocessed_gene_ids_list),
                    number_of_batches=self.environment.FORKS_NUMBER,
            ):
                if os.fork():
                    forks += 1
                    if forks >= self.environment.FORKS_NUMBER:
                        _, status = os.wait()
                        code = os.waitstatus_to_exitcode(status)
                        assert code in (os.EX_OK, os.EX_TEMPFAIL, os.EX_SOFTWARE)
                        assert code != os.EX_SOFTWARE
                        forks -= 1
                else:
                    status = os.EX_OK
                    try:
                        self.search_engine.cds_local_search(gene_id_list=list(balanced_batch))
                    except Exception as exception:
                        self.environment.logger.exception(str(exception))
                        status = os.EX_SOFTWARE
                    finally:
                        # This prevents the child process forked above to keep along the for loop upon completion
                        # of the try/except block. If this was not present, it would resume were it left off, and
                        # fork in turn its own children, duplicating the work done, and creating a huge mess.
                        # We do not want that, so we gracefully exit the process when it is done.
                        os._exit(status)  # https://docs.python.org/3/library/os.html#os._exit
            # guarantees that all forked processes will be terminated before proceeding with the rest
            while forks > 0:
                _, status = os.wait()
                code = os.waitstatus_to_exitcode(status)
                assert code in (os.EX_OK, os.EX_TEMPFAIL, os.EX_SOFTWARE)
                assert code != os.EX_SOFTWARE
                forks -= 1
            self.local_search_complete_identity_and_coverage()
        else:
            self.log_completed_search(local_search=True)
        self.cleanup_local_search()
        self.database_interface.create_filtered_full_length_events_view(
            query_overlap_threshold=self.environment.query_coverage_threshold,
            evalue_threshold=self.environment.evalue_threshold,
        )

    def cleanup_global_search(self):
        except_tables = [
            'Genes', 'Search_monitor',
            'Global_matches_non_reciprocal',
            'Structural_matches_non_reciprocal'
        ]
        if self.environment.SEARCH_ALL:
            except_tables.append('Local_matches')
        self.database_interface.clear_results_database(except_tables=except_tables)
        self.data_container.initialize_database()

    def global_search(
            self
    ):
        genes_to_update = []
        unprocessed_gene_ids_list = self.database_interface.query_to_process_gene_ids(global_search=True)
        if unprocessed_gene_ids_list:
            self.log_search_progress(
                unprocessed_gene_ids_list=unprocessed_gene_ids_list,
                global_search=True
            )
            status: int
            forks: int = 0
            for balanced_batch in self.even_batches(
                    data=list(unprocessed_gene_ids_list),
                    number_of_batches=self.environment.FORKS_NUMBER,
            ):
                if os.fork():
                    forks += 1
                    if forks >= self.environment.FORKS_NUMBER:
                        _, status = os.wait()
                        code = os.waitstatus_to_exitcode(status)
                        assert code in (os.EX_OK, os.EX_TEMPFAIL, os.EX_SOFTWARE)
                        assert code != os.EX_SOFTWARE
                        forks -= 1
                else:
                    status = os.EX_OK
                    try:
                        self.search_engine.cds_global_search(genes_list=list(balanced_batch))
                    except Exception as exception:
                        self.environment.logger.exception(str(exception))
                        status = os.EX_SOFTWARE
                    finally:
                        os._exit(status)  # https://docs.python.org/3/library/os.html#os._exit
            while forks > 0:
                _, status = os.wait()
                code = os.waitstatus_to_exitcode(status)
                assert code in (os.EX_OK, os.EX_TEMPFAIL, os.EX_SOFTWARE)
                assert code != os.EX_SOFTWARE
                forks -= 1
            genes_to_update = self.database_interface.query_gene_ids_global_search()
        else:
            self.log_completed_search(global_search=True)
        self.cleanup_global_search()
        if genes_to_update:
            self.database_interface.update_has_duplicate_genes_table(
                list_tuples=[(gene,) for gene in genes_to_update]
            )

    def structural_search(
            self
    ):
        unprocessed_gene_ids_list = self.data_container.gene_hierarchy_dictionary.keys()
        if unprocessed_gene_ids_list:
            self.log_search_progress(
                unprocessed_gene_ids_list=unprocessed_gene_ids_list,
                structural_search=True
            )
            status: int
            forks: int = 0
            for balanced_batch in self.even_batches(
                    data=list(unprocessed_gene_ids_list),
                    number_of_batches=self.environment.FORKS_NUMBER,
            ):
                if os.fork():
                    forks += 1
                    if forks >= self.environment.FORKS_NUMBER:
                        _, status = os.wait()
                        code = os.waitstatus_to_exitcode(status)
                        assert code in (os.EX_OK, os.EX_TEMPFAIL, os.EX_SOFTWARE)
                        assert code != os.EX_SOFTWARE
                        forks -= 1
                else:
                    status = os.EX_OK
                    try:
                        self.search_engine.cds_structural_search(gene_id_list=list(balanced_batch))
                    except Exception as exception:
                        self.environment.logger.exception(str(exception))
                        status = os.EX_SOFTWARE
                    finally:
                        os._exit(status)  # https://docs.python.org/3/library/os.html#os._exit
            while forks > 0:
                _, status = os.wait()
                code = os.waitstatus_to_exitcode(status)
                assert code in (os.EX_OK, os.EX_TEMPFAIL, os.EX_SOFTWARE)
                assert code != os.EX_SOFTWARE
                forks -= 1
        else:
            self.log_completed_search(structural_search=True)

    def classify_matches_transcript_interdependence(
            self,
            non_reciprocal_coding_matches_list: list
    ) -> list[tuple]:
        match_interdependence_tuples = []
        for match in non_reciprocal_coding_matches_list:
            (gene_id,
             match_id,
             query_start,
             query_end,
             corrected_target_start,
             corrected_target_end) = match
            match_interdependence_classification = self.event_classifier.classify_coding_match_interdependence(
                gene_id=gene_id,
                match_id=match_id,
                query_coordinates=P.open(query_start, query_end),
                target_coordinates=P.open(corrected_target_start, corrected_target_end)
            )
            match_interdependence_tuples.append(match_interdependence_classification)
        return match_interdependence_tuples

    def fetch_tandem_tuples(
            self,
            full_events_list: list
    ) -> list:
        tandemness_tuples = []
        if full_events_list:
            expansions_dictionary = self.event_reconciler.build_expansion_dictionary(
                records=full_events_list
            )
            tandemness_tuples = self.event_reconciler.get_gene_full_events_tandemness_tuples(
                expansions_dictionary=expansions_dictionary
            )
        return tandemness_tuples

    def reconcile(
            self,
            genes_list: list,
    ) -> None:
        for gene_id in genes_list:
            global_records_set = set()
            local_records_set = set()
            structural_records_set = set()
            query_coordinates = set()
            targets_reference_coordinates_dictionary = {}
            if self.environment.SEARCH_ALL or (self.environment.LOCAL_SEARCH and not self.environment.SEARCH_ALL):
                if gene_id in self.local_full_matches_dictionary:
                    local_records_set = self.local_full_matches_dictionary[gene_id]
                    cds_candidates_dictionary = self.search_engine.get_candidate_cds_coordinates(
                        gene_id=gene_id
                    )
                    (query_coordinates,
                     targets_reference_coordinates_dictionary
                     ) = self.event_reconciler.align_target_coordinates(
                        gene_id=gene_id,
                        local_records_set=local_records_set,
                        cds_candidates_dictionary=cds_candidates_dictionary
                    )
                    corrected_coordinates_tuples = self.event_reconciler.get_matches_corrected_coordinates_and_identity(
                        gene_id=gene_id,
                        local_records_set=local_records_set,
                        targets_reference_coordinates_dictionary=targets_reference_coordinates_dictionary
                    )
                    attempt = False
                    while not attempt:
                        try:
                            self.database_interface.insert_corrected_target_start_end(
                                list_tuples=corrected_coordinates_tuples
                            )
                            attempt = True
                        except Exception as e:
                            if "locked" in str(e):
                                time.sleep(random.randrange(start=0, stop=self.environment.sleep_max_seconds))
                            else:
                                self.environment.logger.exception(e)
                                sys.exit()
            if self.environment.SEARCH_ALL or (self.environment.GLOBAL_SEARCH and not self.environment.SEARCH_ALL):
                global_records_set = self.global_full_matches_dictionary[gene_id]
            if self.environment.STRUCTURAL_SEARCH:
                structural_records_set = self.structural_full_matches_dictionary[gene_id]
            gene_graph = self.event_reconciler.create_events_multigraph(
                targets_reference_coordinates_dictionary=targets_reference_coordinates_dictionary,
                query_local_coordinates_set=query_coordinates,
                local_records_set=local_records_set,
                global_records_set=global_records_set,
                structural_records_set=structural_records_set
            )
            (gene_events_list,
             non_reciprocal_fragment_ids_list,
             full_events_list
             ) = self.event_reconciler.get_reconciled_graph_and_expansion_events_tuples(
                targets_reference_coordinates_dictionary=targets_reference_coordinates_dictionary,
                gene_id=gene_id,
                gene_graph=gene_graph
            )
            tandemness_tuples = self.fetch_tandem_tuples(full_events_list=full_events_list)
            attempt = False
            while not attempt:
                try:
                    self.database_interface.insert_expansion_table(
                        list_tuples=gene_events_list,
                        list_tuples_full=full_events_list,
                        list_tuples_tandemness=tandemness_tuples,
                    )
                    attempt = True
                except Exception as e:
                    if "locked" in str(e):
                        time.sleep(random.randrange(start=0, stop=self.environment.sleep_max_seconds))
                    else:
                        self.environment.logger.exception(e)
                        sys.exit()
            if self.environment.SEARCH_ALL or (self.environment.LOCAL_SEARCH and not self.environment.SEARCH_ALL):
                attempt = False
                while not attempt:
                    try:
                        self.database_interface.insert_in_non_reciprocal_fragments_table(
                            fragment_ids_list=non_reciprocal_fragment_ids_list,
                            gene_id=gene_id
                        )
                        attempt = True
                    except Exception as e:
                        if "locked" in str(e):
                            time.sleep(random.randrange(start=0, stop=self.environment.sleep_max_seconds))
                        else:
                            self.environment.logger.exception(e)
                            sys.exit()

    def events_reconciliation(
            self,
    ):
        genes_to_process = set()
        if self.environment.SEARCH_ALL or (self.environment.LOCAL_SEARCH and not self.environment.SEARCH_ALL):
            self.database_interface.create_non_reciprocal_fragments_table()
            local_full_matches_list = self.database_interface.query_full_length_events()
            self.local_full_matches_dictionary = self.event_reconciler.get_gene_events_dictionary(
                local_full_matches_list=local_full_matches_list
            )
            genes_to_process = genes_to_process.union(set(self.local_full_matches_dictionary.keys()))
        if self.environment.SEARCH_ALL or (self.environment.GLOBAL_SEARCH and not self.environment.SEARCH_ALL):
            self.global_full_matches_dictionary = self.database_interface.query_global_cds_events()
            genes_to_process = genes_to_process.union(set(self.global_full_matches_dictionary.keys()))
        if self.environment.STRUCTURAL_SEARCH:
            self.structural_full_matches_dictionary = self.database_interface.query_structural_cds_events()
        status: int
        forks: int = 0
        for balanced_batch in self.even_batches(
                data=list(genes_to_process),
                number_of_batches=self.environment.FORKS_NUMBER,
        ):
            if os.fork():
                forks += 1
                if forks >= self.environment.FORKS_NUMBER:
                    _, status = os.wait()
                    code = os.waitstatus_to_exitcode(status)
                    assert code in (os.EX_OK, os.EX_TEMPFAIL, os.EX_SOFTWARE)
                    assert code != os.EX_SOFTWARE
                    forks -= 1
            else:
                status = os.EX_OK
                try:
                    self.reconcile(genes_list=list(balanced_batch))
                except Exception as exception:
                    self.environment.logger.exception(str(exception))
                    status = os.EX_SOFTWARE
                finally:
                    os._exit(status)
        while forks > 0:
            _, status = os.wait()
            code = os.waitstatus_to_exitcode(status)
            assert code in (os.EX_OK, os.EX_TEMPFAIL, os.EX_SOFTWARE)
            assert code != os.EX_SOFTWARE
            forks -= 1
        self.database_interface.drop_table(table_name='Matches_full_length')
        genes_with_duplicates = self.database_interface.query_genes_with_duplicated_cds()
        self.database_interface.update_has_duplicate_genes_table(list_tuples=genes_with_duplicates)

    def transcript_interdependence_classification(
            self
    ):
        # MATCHES INTERDEPENDENCE CLASSIFICATION
        if self.environment.SEARCH_ALL or (self.environment.GLOBAL_SEARCH and not self.environment.SEARCH_ALL):
            cds_global_matches_list = self.database_interface.query_cds_global_matches()
            transcripts_iterdependence_global_matches_tuples = self.classify_matches_transcript_interdependence(
                non_reciprocal_coding_matches_list=cds_global_matches_list
            )
            self.database_interface.insert_matches_interdependence_classification(
                tuples_list=[
                    (n_mrnas, all_, present, absent, neither, category, frag_id)
                    for _, frag_id, n_mrnas, _, all_, present, absent, neither, category, _
                    in transcripts_iterdependence_global_matches_tuples
                ],
                table_name='Global_matches_non_reciprocal',
                table_identifier_column='GlobalFragmentID'
            )

        if self.environment.SEARCH_ALL or (self.environment.LOCAL_SEARCH or not self.environment.SEARCH_ALL):
            non_reciprocal_coding_matches_list = self.database_interface.query_non_reciprocal_coding_matches()
            transcripts_iterdependence_tuples = self.classify_matches_transcript_interdependence(
                non_reciprocal_coding_matches_list=non_reciprocal_coding_matches_list
            )
            self.database_interface.insert_matches_interdependence_classification(
                tuples_list=[
                    (n_mrnas, all_, present, absent, neither, category, frag_id)
                    for _, frag_id, n_mrnas, _, all_, present, absent, neither, category, _
                    in transcripts_iterdependence_tuples
                ],
                table_name='Local_matches_non_reciprocal',
                table_identifier_column='LocalFragmentID'
            )
        # EXPANSION INTERDEPENDENCE CLASSIFICATION
        expansions_dictionary = self.database_interface.query_coding_expansion_events()
        expansion_interdependence_tuples = self.event_classifier.classify_expansion_interdependence(
            expansions_dictionary=expansions_dictionary
        )
        self.database_interface.insert_expansions_interdependence_classification(
            list_tuples=expansion_interdependence_tuples
            )
        if self.environment.GLOBAL_SEARCH:
            self.database_interface.drop_table('Expansions')

    def runtime_logger(
            self,
    ):
        gene_ids_list = list(self.data_container.gene_hierarchy_dictionary.keys())
        runtime_hours = round((datetime.now() - self.tic).total_seconds() / 3600, 2)
        with open(self.environment.log_file_name, 'w') as f:
            f.write(self.environment.exonize_pipeline_settings)
            f.write(
                f'\nRuntime (hours):              {runtime_hours}'
                f'\nNumber of processed genes:    {len(gene_ids_list)}'
            )

    @staticmethod
    def even_batches(
            data: Sequence[Any],
            number_of_batches: int = 1,
    ) -> Iterator[Sequence[Any]]:
        """
        Given a list and a number of batches, returns 'number_of_batches'
         consecutive subsets elements of 'data' of even size each, except
         for the last one whose size is the remainder of the division of
        the length of 'data' by 'number_of_batches'.
        """
        even_batch_size = (len(data) // number_of_batches) + 1
        for batch_number in range(number_of_batches):
            batch_start_index = batch_number * even_batch_size
            batch_end_index = min((batch_number + 1) * even_batch_size, len(data))
            yield data[batch_start_index:batch_end_index]

    def run_exonize_pipeline(
            self,
    ) -> None:
        """
        run_exonize_pipeline iterates over all genes in the gene_hierarchy_dictionary
        attribute and performs a tblastx search for each representative
        CDS (see get_candidate_CDS_coords).
        The steps are the following:
        - 1. The function checks if the gene has already been processed.
         If so, the gene is skipped.
        - 2. For each gene, the function iterates over all representative
         CDSs and performs a tblastx search.
        - 3. The percent_query (hit query coverage) column is inserted
        in the Fragments table.
        - 4. The Filtered_full_length_events View is created. This view
        contains all tblastx hits that have passed the filtering step.
        - 5. The Mrna_counts View is created. This view contains the number
        of transcripts associated with each gene.
        - 6. The function creates the Cumulative_counts table. This table
        contains the cumulative counts of the different event types across
        transcripts.
        - 7. The function collects all the raw concatenated event types
        (see query_concat_categ_pairs).
        - 8. The function generates the unique event types
         (see generate_unique_events_list) so that no event type is repeated.
        - 9. The function inserts the unique event types in the
        Event_categ_full_length_events_cumulative_counts table.
        - 10. The function collects the identity and DNA sequence tuples
         (see get_identity_and_dna_seq_tuples) and inserts them in the Fragments table.
        - 11. The function collects all events in the
        Full_length_events_cumulative_counts table.
        - 12. The function reconciles the events by assigning a "pair ID"
         to each event (see assign_pair_ids).
        - 13. The function creates the Exclusive_pairs view. This view contains
         all the events that follow the mutually exclusive category.
        """
        self.environment.logger.info(f'Running Exonize for: {self.environment.output_prefix}')
        self.data_container.prepare_data()
        if self.environment.STRUCTURAL_SEARCH:
            self.structural_search()
            self.data_container.structural_search_cdss_to_query = self.data_container.fetch_structural_queried_cdss()
        if self.environment.SEARCH_ALL:
            self.local_search()
            self.global_search()
        elif not self.environment.SEARCH_ALL and self.environment.GLOBAL_SEARCH:
            self.global_search()
        elif not self.environment.SEARCH_ALL and self.environment.LOCAL_SEARCH:
            self.local_search()
        self.environment.logger.info('Starting reconciliation and classification...')
        self.environment.logger.info('Reconciling local matches')
        self.events_reconciliation()
        self.environment.logger.info('Classifying events')
        self.transcript_interdependence_classification()
        self.runtime_logger()
        self.environment.logger.info('Process completed successfully')
        if self.environment.CSV:
            self.database_interface.export_all_tables_to_csv(
                output_dir=self.environment.csv_path
            )
        self.data_container.clear_working_directory()
