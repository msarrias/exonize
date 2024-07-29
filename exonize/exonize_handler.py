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
from datetime import date, datetime
import sys
from pathlib import Path
import cProfile
import gc
from exonize.profiling import get_run_performance_profile
from exonize.environment_setup import EnvironmentSetup
from exonize.data_preprocessor import DataPreprocessor
from exonize.sqlite_handler import SqliteHandler
from exonize.blast_searcher import BLASTsearcher
from exonize.classifier_handler import ClassifierHandler
from exonize.reconciler_handler import ReconcilerHandler


class Exonize(object):
    def __init__(
            self,
            gff_file_path: Path,
            genome_file_path: Path,
            output_prefix: str,
            draw_event_multigraphs: bool,
            csv: bool,
            enable_debug: bool,
            soft_force: bool,
            hard_force: bool,
            evalue_threshold: float,
            sleep_max_seconds: int,
            min_exon_length: int,
            cds_overlapping_threshold: float,
            query_overlapping_threshold: float,
            self_hit_threshold: float,
            cpus_number: int,
            timeout_database: int,
            output_directory_path: Path
    ):
        self._DEBUG_MODE = enable_debug
        self.SOFT_FORCE = soft_force
        self.HARD_FORCE = hard_force
        self.FORKS_NUMBER = cpus_number

        self.gff_file_path = gff_file_path
        self.genome_file_path = genome_file_path
        self.output_prefix = output_prefix
        self.evalue_threshold = evalue_threshold
        self.cds_overlapping_threshold = cds_overlapping_threshold
        self.query_overlapping_threshold = query_overlapping_threshold
        self.self_hit_threshold = self_hit_threshold
        self.min_exon_length = min_exon_length
        self.sleep_max_seconds = sleep_max_seconds
        self.timeout_database = timeout_database
        self.draw_event_multigraphs = draw_event_multigraphs
        self.csv = csv
        self.tic = datetime.now()
        self.full_matches_dictionary = {}

        if not self.output_prefix:
            self.output_prefix = gff_file_path.stem

        if output_directory_path:
            self.working_directory = output_directory_path / f'{self.output_prefix}_exonize'
        else:
            self.working_directory = Path(f'{self.output_prefix}_exonize')
        self.results_database_path = self.working_directory / f'{self.output_prefix}_results.db'
        self.log_file_name = self.working_directory / f"exonize_settings_{datetime.now():%Y%m%d_%H%M%S}.log"
        self.PROFILE_PATH = self.working_directory / 'cProfile_dump_stats.dmp'

        # Initialize logger and set up environment
        self.environment = EnvironmentSetup(
            hard_force=self.HARD_FORCE,
            soft_force=self.SOFT_FORCE,
            working_directory=self.working_directory,
            results_database_path=self.results_database_path,
        )
        self.database_interface = SqliteHandler(
            results_database_path=self.results_database_path,
            timeout_database=self.timeout_database,
        )
        self.data_container = DataPreprocessor(
            logger_obj=self.environment,
            database_interface=self.database_interface,
            working_directory=self.working_directory,
            gff_file_path=self.gff_file_path,
            output_prefix=self.output_prefix,
            genome_file_path=self.genome_file_path,
            debug_mode=self._DEBUG_MODE,
            evalue_threshold=evalue_threshold,
            self_hit_threshold=self.self_hit_threshold,
            cds_overlapping_threshold=self.cds_overlapping_threshold,
            query_overlapping_threshold=self.query_overlapping_threshold,
            min_exon_length=self.min_exon_length,
            draw_event_multigraphs=self.draw_event_multigraphs,
            csv=self.csv,
        )

        self.blast_engine = BLASTsearcher(
            data_container=self.data_container,
            sleep_max_seconds=self.sleep_max_seconds,
            self_hit_threshold=self.self_hit_threshold,
            min_exon_length=self.min_exon_length,
            cds_overlapping_threshold=self.cds_overlapping_threshold,
            evalue_threshold=self.evalue_threshold,
            debug_mode=self._DEBUG_MODE,
        )
        self.event_classifier = ClassifierHandler(
            blast_engine=self.blast_engine,
            cds_overlapping_threshold=self.cds_overlapping_threshold,
        )
        self.event_reconciler = ReconcilerHandler(
            blast_engine=self.blast_engine,
            cds_overlapping_threshold=self.cds_overlapping_threshold,
        )
        self.exonize_pipeline_settings = f"""
Exonize - settings
--------------------------------
Date:                        {date.today()}
python version:              {sys.version}
cpu count:                   {self.FORKS_NUMBER}
--------------------------------
Indentifier:                 {self.output_prefix}
GFF file:                    {gff_file_path}
Genome file:                 {genome_file_path}
--------------------------------
tblastx e-value threshold:   {evalue_threshold}
CDS overlapping threshold:   {cds_overlapping_threshold}
Query overlapping threshold: {query_overlapping_threshold}
Self-hit threshold:          {self_hit_threshold}
Min exon length (bps):       {min_exon_length}
--------------------------------
Exonize results database:   {self.results_database_path.name}
        """

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

    def search(
            self
    ):
        gene_ids_list = list(self.data_container.gene_hierarchy_dictionary.keys())
        processed_gene_ids_list = set(
            self.database_interface.query_gene_ids_in_results_database()
        )
        unprocessed_gene_ids_list = list(set(gene_ids_list) - processed_gene_ids_list)
        if unprocessed_gene_ids_list:
            gene_count = len(gene_ids_list)
            out_message = 'Starting exon duplication search'
            if len(unprocessed_gene_ids_list) != gene_count:
                out_message = 'Resuming search'
            self.environment.logger.info(
                f'{out_message} for'
                f' {len(unprocessed_gene_ids_list)}/{gene_count} genes'
            )

            self.environment.logger.info(
                'Exonizing: this may take a while ...'
            )
            for balanced_batch in self.even_batches(
                    data=unprocessed_gene_ids_list,
                    number_of_batches=self.FORKS_NUMBER,
            ):
                # This part effectively forks a child process, independent of the parent process, and that
                # will be responsible for processing the genes in the batch, parallel to the other children forked
                # in the same way during the rest of the loop.
                # A note to understand why this works: os.fork() returns 0 in the child process, and the PID
                # of the child process in the parent process. So the parent process always goes to the part
                # evaluated to True (if os.fork()), and the child process always goes to the part evaluated
                # to false (0 is evaluated to False).
                # The parallel part happens because the main parent process will keep along the for loop, and will
                # fork more children, until the number of children reaches the maximum number of children allowed,
                # doing nothing else but forking until 'FORKS_NUMBER' is reached.
                # # Benchmark without any parallel computation:
                # pr = cProfile.Profile()
                # pr.enable()
                # for gene_id in unprocessed_gene_ids_list:
                #     self.blast_engine.find_coding_exon_duplicates(gene_id)
                # pr.disable()
                # get_run_performance_profile(self.PROFILE_PATH, pr)
                # Benchmark with parallel computation using os.fork:
                pr = cProfile.Profile()
                pr.enable()
                gc.collect()
                gc.freeze()
                # transactions_pks: set[int]
                status: int
                code: int
                forks: int = 0
                if os.fork():
                    forks += 1
                    if forks >= self.FORKS_NUMBER:
                        _, status = os.wait()
                        code = os.waitstatus_to_exitcode(status)
                        assert code in (os.EX_OK, os.EX_TEMPFAIL, os.EX_SOFTWARE)
                        assert code != os.EX_SOFTWARE
                        forks -= 1
                else:
                    status = os.EX_OK
                    try:
                        self.blast_engine.find_coding_exon_duplicates(
                            gene_id_list=list(balanced_batch)
                        )
                    except Exception as exception:
                        print(exception)
                        self.environment.logger.exception(
                            str(exception)
                        )
                        status = os.EX_SOFTWARE
                    finally:
                        # This prevents the child process forked above to keep along the for loop upon completion
                        # of the try/except block. If this was not present, it would resume were it left off, and
                        # fork in turn its own children, duplicating the work done, and creating a huge mess.
                        # We do not want that, so we gracefully exit the process when it is done.
                        os._exit(status)  # https://docs.python.org/3/library/os.html#os._exit
                # This blocks guarantees that all forked processes will be terminated before proceeding with the rest
                while forks > 0:
                    _, status = os.wait()
                    code = os.waitstatus_to_exitcode(status)
                    assert code in (os.EX_OK, os.EX_TEMPFAIL, os.EX_SOFTWARE)
                    assert code != os.EX_SOFTWARE
                    forks -= 1
                    gc.unfreeze()
                    pr.disable()
                    get_run_performance_profile(self.PROFILE_PATH, pr)
            self.database_interface.insert_percent_query_column_to_fragments()
            matches_list = self.database_interface.query_fragments()
            identity_and_sequence_tuples = self.blast_engine.get_identity_and_dna_seq_tuples(
                matches_list=matches_list
            )
            self.database_interface.insert_identity_and_dna_algns_columns(
                list_tuples=identity_and_sequence_tuples
            )
        else:
            self.environment.logger.info(
                'All genes have been processed. If you want to re-run the analysis, '
                'consider using the hard-force/soft-force flag'
            )
            self.database_interface.clear_results_database()
            self.database_interface.connect_create_results_database()
            # self.database_interface.create_matches_interdependence_expansions_counts_table()
        self.database_interface.create_filtered_full_length_events_view(
            query_overlap_threshold=self.query_overlapping_threshold,
            evalue_threshold=self.evalue_threshold,
        )

    def classify_expansion_events_interdependence(
            self,
    ) -> list[tuple]:
        expansion_interdependence_tuples = []
        expansion_events_dict = self.database_interface.query_full_expansion_events()
        for gene_id, expansions_dict in expansion_events_dict.items():
            for expansion_id, coding_events_coordinates_list in expansions_dict.items():
                if len(coding_events_coordinates_list) > 1:
                    expansion_interdependence_tuples.extend(
                        self.event_classifier.classify_expansion_transcript_interdependence(
                            gene_id=gene_id,
                            expansion_id=expansion_id,
                            coding_coordinates_list=coding_events_coordinates_list
                        )
                    )
        return expansion_interdependence_tuples

    def classify_matches_transcript_interdependence(
            self
    ) -> list[tuple]:
        # Classify matches based on the mode and interdependence
        expansions_gene_dictionary = self.database_interface.query_expansion_coding_events()
        self.event_classifier.initialize_list_of_tuples()
        for gene_id, expansions_dict in expansions_gene_dictionary.items():
            for expansion_id, records in expansions_dict.items():
                for record in records:
                    self.event_classifier.classify_match_interdependence(
                        row_tuple=record
                    )
        return self.event_classifier.tuples_match_transcript_interdependence

    def update_mode_cumulative_counts_table(
            self,
    ):
        query_concat_categ_pair_list = self.database_interface.query_concat_categ_pairs()
        if query_concat_categ_pair_list:
            reduced_event_types_tuples = self.generate_unique_events_list(
                events_list=query_concat_categ_pair_list,
                event_type_idx=-1
            )
            self.database_interface.insert_event_categ_matches_interdependence_counts(
                list_tuples=reduced_event_types_tuples
            )

    def reconcile(
            self,
            gene_id: str,
    ) -> None:
        tblastx_records_set = self.full_matches_dictionary[gene_id]
        cds_candidates_dictionary = self.blast_engine.get_candidate_cds_coordinates(
            gene_id=gene_id
        )
        (query_coordinates,
         reference_coordinates_dictionary) = self.event_reconciler.align_target_coordinates(
            gene_id=gene_id,
            tblastx_records_set=tblastx_records_set
        )
        corrected_coordinates_tuples = self.event_reconciler.get_matches_corrected_coordinates_and_identity(
            gene_id=gene_id,
            tblastx_records_set=tblastx_records_set,
            reference_coordinates_dictionary=reference_coordinates_dictionary,
            cds_candidates_dictionary=cds_candidates_dictionary
        )
        self.database_interface.insert_corrected_target_start_end(
            list_tuples=corrected_coordinates_tuples
        )
        gene_graph = self.event_reconciler.create_events_multigraph(
            tblastx_records_set=tblastx_records_set,
            query_coordinates_set=query_coordinates,
            reference_coordinates_dictionary=reference_coordinates_dictionary
        )
        if self.draw_event_multigraphs:
            self.event_reconciler.draw_event_multigraph(
                gene_graph=gene_graph,
                figure_path=self.data_container.multigraphs_path / f'{gene_id}.png'
            )
        (gene_events_list,
         non_reciprocal_fragment_ids_list
         ) = self.event_reconciler.get_reconciled_graph_and_expansion_events_tuples(
            reference_coordinates_dictionary=reference_coordinates_dictionary,
            gene_id=gene_id,
            gene_graph=gene_graph
        )
        self.database_interface.insert_expansion_table(
            list_tuples=gene_events_list
        )
        self.database_interface.insert_in_non_reciprocal_fragments_table(
            fragment_ids_list=non_reciprocal_fragment_ids_list
        )

    def events_reconciliation(
            self,
    ):
        tblastx_full_matches_list = self.database_interface.query_full_events()
        self.database_interface.create_non_reciprocal_fragments_table()
        # group full matches by gene id
        self.full_matches_dictionary = self.event_reconciler.get_gene_events_dictionary(
            tblastx_full_matches_list=tblastx_full_matches_list
        )
        genes_to_process = list(self.full_matches_dictionary.keys())
        for gene_id in genes_to_process:
            self.reconcile(gene_id)
        self.database_interface.drop_table('Matches_full_length')
        genes_with_duplicates = self.database_interface.query_genes_with_duplicated_cds()
        self.database_interface.update_has_duplicate_genes_table(
            list_tuples=genes_with_duplicates
        )

    def events_classification(
            self
    ):
        transcripts_iterdependence_tuples = self.classify_matches_transcript_interdependence()
        self.database_interface.insert_matches_interdependence_classification(
            tuples_list=transcripts_iterdependence_tuples
        )
        self.database_interface.create_matches_interdependence_counts_table()
        self.update_mode_cumulative_counts_table()
        records = self.database_interface.query_interdependence_counts_matches()
        records_to_insert = self.event_classifier.classify_transcript_interdependence_counts(
            records
        )
        self.database_interface.insert_classification_column_to_matches_interdependence_counts_table(
            list_tuples=records_to_insert
        )
        # # CLASSIFY MATCHES TRANSCRIPT INTERDEPENDENCE
        # expansion_interdependence_tuples = self.classify_expansion_events_interdependence()
        # self.database_interface.insert_matches_interdependence_expansions_counts(
        #     tuples_list=expansion_interdependence_tuples
        # )

    def runtime_logger(self):
        gene_ids_list = list(self.data_container.gene_hierarchy_dictionary.keys())
        runtime_hours = round((datetime.now() - self.tic).total_seconds() / 3600, 2)
        if runtime_hours < 1:
            runtime_hours = ' < 1'
        with open(self.log_file_name, 'w') as f:
            f.write(self.exonize_pipeline_settings)
            f.write(
                f'\nRuntime (hours):             {runtime_hours}'
                f'\nNumber of processed genes:   {len(gene_ids_list)}'
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
        # We round up to the upper integer value to guarantee that there
        # will be 'number_of_batches' batches
        even_batch_size = (len(data) // number_of_batches) + 1
        for batch_number in range(number_of_batches):
            batch_start_index = batch_number * even_batch_size
            batch_end_index = min((batch_number + 1) * even_batch_size, len(data))
            yield data[batch_start_index:batch_end_index]

    def run_exonize_pipeline(self) -> None:
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
        self.environment.logger.info(f'Running Exonize for: {self.output_prefix}')
        self.data_container.prepare_data()
        self.search()
        self.environment.logger.info('Reconciling matches')
        self.events_reconciliation()
        self.environment.logger.info('Classifying events')
        self.events_classification()
        self.runtime_logger()
        self.environment.logger.info('Process completed successfully')
        if self.csv:
            self.database_interface.export_all_tables_to_csv(
                output_dir=self.data_container.csv_path
            )
        self.data_container.clear_working_directory()
