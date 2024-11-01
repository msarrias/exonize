import logging
import os
import shutil
from pathlib import Path
from datetime import datetime


class EnvironmentSetup(object):
    logger: logging.Logger

    def __init__(
            self,
            genome_file_path: Path,
            gff_file_path: Path,
            output_prefix: str,
            output_directory_path: Path,
            gene_annot_feature: str,
            cds_annot_feature: str,
            transcript_annot_feature: str,
            sequence_base: int,
            frame_base: int,
            min_exon_length: int,
            evalue_threshold: float,
            peptide_identity_threshold: float,
            fraction_of_aligned_positions: float,
            exon_clustering_overlap_threshold: float,
            targets_clustering_overlap_threshold: float,
            query_coverage_threshold: float,
            self_hit_threshold: float,
            global_search: bool,
            local_search: bool,
            hard_force: bool,
            soft_force: bool,
            debug_mode: bool,
            csv: bool,
            sleep_max_seconds: int,
            timeout_database: int,
            cpus_number: int,
    ):
        self.__FILE_ONLY_INFO = 9
        self.DEBUG_MODE = debug_mode
        self.SOFT_FORCE = soft_force
        self.HARD_FORCE = hard_force
        self.FORKS_NUMBER = cpus_number
        self.GLOBAL_SEARCH = global_search
        self.LOCAL_SEARCH = local_search
        self.SEARCH_ALL = not self.GLOBAL_SEARCH and not self.LOCAL_SEARCH
        self.CSV = csv

        self.gff_file_path = gff_file_path,
        self.genome_file_path = genome_file_path
        self.output_directory_path = output_directory_path
        self.output_prefix = output_prefix

        # Search criteria parameters
        self.evalue_threshold = evalue_threshold
        self.min_exon_length = min_exon_length
        self.self_hit_threshold = self_hit_threshold
        self.query_coverage_threshold = query_coverage_threshold
        self.exon_clustering_overlap_threshold = exon_clustering_overlap_threshold
        self.targets_clustering_overlap_threshold = targets_clustering_overlap_threshold
        self.fraction_of_aligned_positions = fraction_of_aligned_positions
        self.peptide_identity_threshold = peptide_identity_threshold
        self.sequence_base = sequence_base
        self.frame_base = frame_base

        # other
        self.sleep_max_seconds = sleep_max_seconds
        self.timeout_database = timeout_database

        # Annotation features
        self.gene_annot_feature = gene_annot_feature
        self.cds_annot_feature = cds_annot_feature
        self.transcript_annot_feature = transcript_annot_feature

        # Constants - mode classification
        self.full = 'FULL'
        self.partial_insertion = 'PARTIAL_INSERTION'
        self.partial_excision = 'PARTIAL_EXCISION'
        self.inter_boundary = 'INTER_BOUNDARY'
        self.intronic = 'INTRONIC'

        self.color_map = {
            self.partial_insertion: 'blue',
            self.partial_excision: 'purple',
            self.full: 'green',
            self.intronic: 'red',
            self.inter_boundary: 'orange'

        }
        if not self.output_prefix:
            self.output_prefix = self.gff_file_path.stem
        if self.output_directory_path:
            self.working_directory = self.output_directory_path / f'{self.output_prefix}_exonize'
        else:
            self.working_directory = Path(f'{self.output_prefix}_exonize')
        self.results_database_path = self.working_directory / f'{self.output_prefix}_results.db'
        self.log_file_name = self.working_directory / f"{self.output_prefix}.log"
        self.PROFILE_PATH = self.working_directory / 'cProfile_dump_stats.dmp'

        # Derived attributes that depend on initial parameters
        self.genome_database_path = self.working_directory / f'{self.output_prefix}_genome_annotations.db'
        self.protein_database_path = self.working_directory / f'{self.output_prefix}_protein.db'
        self.gene_hierarchy_path = self.working_directory / f"{self.output_prefix}_gene_hierarchy.pkl"
        if self.CSV:
            self.csv_path = self.working_directory / "csvs"
        if self.DEBUG_MODE:
            os.makedirs(self.working_directory / 'input', exist_ok=True)
            os.makedirs(self.working_directory / 'output', exist_ok=True)
        if self.CSV:
            os.makedirs(self.csv_path, exist_ok=True)
        self.setup_environment()

    def configure_logger(self):
        """
        configure_logger is a function that configures the logger.
        INFO level is used for the log file and WARNING and ERROR
        level for the console.
        """
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        # Define console handler for the "INFO" level and above ("DEBUG", "WARNING", "EXCEPTION").
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(logging.Formatter('[%(levelname)s]: %(message)s'))

        logging.basicConfig(
            level=logging.INFO, handlers=[console_handler]
        )

    def setup_environment(self):
        if self.HARD_FORCE:
            if self.working_directory.exists():
                shutil.rmtree(self.working_directory)
        elif self.SOFT_FORCE:
            if self.results_database_path.exists():
                os.remove(self.results_database_path)
        os.makedirs(self.working_directory, exist_ok=True)
        self.configure_logger()
