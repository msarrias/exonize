import gffutils
import gzip
import os
import pickle
import subprocess
import sys
import shutil
from Bio import SeqIO
import portion as P
from pathlib import Path
import tarfile


class DataPreprocessor(object):
    utr_features = ['five_prime_UTR', 'three_prime_UTR']
    features_of_interest = ['CDS', 'exon', 'intron'] + utr_features

    def __init__(
            self,
            logger_obj: object,
            database_interface: object,
            working_directory: Path,
            gff_file_path: Path,
            output_prefix: str,
            genome_file_path: Path,
            debug_mode: bool,
            evalue_threshold: float,
            self_hit_threshold: float,
            cds_overlapping_threshold: float,
            query_overlapping_threshold: float,
            min_exon_length: int,
            draw_event_multigraphs: bool,
            csv: bool
    ):
        self.environment = logger_obj
        self.database_interface = database_interface
        self.working_directory = working_directory
        self.gff_file_path = gff_file_path
        self.output_prefix = output_prefix
        self.genome_file_path = genome_file_path
        self.evalue_threshold = evalue_threshold
        self.self_hit_threshold = self_hit_threshold
        self.cds_overlapping_threshold = cds_overlapping_threshold
        self.query_overlapping_threshold = query_overlapping_threshold
        self.min_exon_length = min_exon_length
        self.timeout_database = database_interface.timeout_database
        self.results_database = database_interface.results_database_path
        self.draw_event_multigraphs = draw_event_multigraphs
        self.csv = csv
        self._DEBUG_MODE = debug_mode

        self.database_features = None
        self.old_filename = None
        self.genome_database = None
        self.genome_dictionary = dict()
        self.gene_hierarchy_dictionary = dict()

        # Derived attributes that depend on initial parameters
        self.genome_database_path = self.working_directory / f'{self.output_prefix}_genome_annotations.db'
        self.protein_database_path = self.working_directory / f'{self.output_prefix}_protein.db'
        self.gene_hierarchy_path = self.working_directory / f"{self.output_prefix}_gene_hierarchy.pkl"
        if self.draw_event_multigraphs:
            self.multigraphs_path = self.working_directory / 'multigraphs'
        if self.csv:
            self.csv_path = self.working_directory / "csvs"

    @staticmethod
    def dump_pkl_file(
            out_file_path: Path,
            records_dictionary: dict
    ) -> None:
        """
        dump_pkl_file is a function that dumps an object into a pickle file.
        """
        with open(out_file_path, 'wb') as handle:
            pickle.dump(records_dictionary, handle)

    @staticmethod
    def read_pkl_file(
            file_path: Path
    ) -> dict:
        """
        read_pkl_file is a function that reads a pickle file and returns
         the object stored in it.
        """
        with open(file_path, 'rb') as handle:
            read_file = pickle.load(handle)
        return read_file

    @staticmethod
    def sort_list_intervals_dict(
            list_dictionaries: list,
            reverse=False
    ) -> list:
        """
        sort_list_intervals_dict is a function that sorts a list
        of dictionaries based on the coordinates of the intervals
        present in the dictionaries.
        The list is sorted in ascending order by default.
        """
        return sorted(
            list_dictionaries,
            key=lambda x: (x['coordinate'].lower, x['coordinate']),
            reverse=reverse
        )

    @staticmethod
    def min_perc_overlap(
            intv_i: P.Interval,
            intv_j: P.Interval,
    ) -> float:
        def get_interval_length(
                interval: P.Interval,
        ):
            return sum(intv.upper - intv.lower for intv in interval)

        """
        Given two intervals, the function returns the percentage of the overlapping
        region relative to the longest interval. The percentage overlap of the shortest
        interval will always be greater or equal than that of the longest interval.
        :param intv_i:
        :param intv_j:
        :return:
        """
        if intv_i.overlaps(intv_j):
            intersection_span = get_interval_length(intv_i.intersection(intv_j))
            longest_length = max(get_interval_length(intv_i), get_interval_length(intv_j))
            return round(intersection_span / longest_length, 3)
        return 0.0

    @staticmethod
    def reverse_sequence_bool(
            gene_strand: str,
    ) -> bool:
        """
        reverse_sequence_bool checks if the gene is in the negative
        strand and returns True if it is.
        :param gene_strand: strand
        """
        return gene_strand == '-'

    def convert_gtf_to_gff(
            self,
    ) -> None:
        """
        Convert a GTF file to GFF format using gffread. Flags description:
        -'O': This flag is used to enable the output of the file in GFF3 format.
        -'o': This flag is used to specify the output file name.
        """
        gffread_command = ["gffread", self.old_filename, "-O", "-o", self.gff_file_path]
        subprocess.call(gffread_command)

    def create_genome_database(
            self,
    ) -> None:
        """
        create_genome_database is a function that creates a gffutils
        database from a GFF3 file.
        Args:
        - dbfn: path to the database file
        - force: if True, the database will be overwritten if it
        already exists
        - keep_order: if True, the order of the features in the GFF
        file will be preserved
        - merge_strategy: if 'create_unique', the database will be
         created with unique IDs
        - sort_attribute_values: if True, the attribute values will
         be sorted
        - disable_infer_genes: if True, the function will not attempt
         to automatically infer gene features
        - disable_infer_transcripts: if True, the function will not
        attempt to automatically infer transcript features
        """
        try:
            self.environment.logger.info(
                "Parsing annotations - This may take a while..."
            )
            self.genome_database = gffutils.create_db(
                data=str(self.gff_file_path),
                dbfn=str(self.genome_database_path),
                force=True,
                keep_order=True,
                merge_strategy='create_unique',
                sort_attribute_values=True,
                disable_infer_genes=True,
                disable_infer_transcripts=True
            )
        except ValueError as e:
            self.environment.logger.critical(
                f"Incorrect genome annotations file {e}"
            )
            sys.exit()

    def search_create_intron_annotations(
            self,
    ) -> None:
        """
        search_create_intron_annotations is a function that verifies
        that the gffutils database contains intron annotations, if not,
        it attempts to write them.
        Some of the db.update() parameters description:
        - make_backup: if True, a backup of the database will be
        created before updating it.
        """
        if 'intron' not in self.database_features:
            self.environment.logger.info(
                "The GFF file does not contain intron annotations - "
                "attempting to write intron annotations in database",
            )
            try:
                self.genome_database.update(
                    list(self.genome_database.create_introns()),
                    make_backup=False
                )
            except ValueError as e:
                self.environment.logger.critical(
                    f"failed to write intron annotations in database. "
                    f"Please provide a GFF3 file with intron annotations {e}"
                )
                sys.exit()

    def create_parse_or_update_database(
            self,
    ) -> None:
        """
        create_parse_or_update_database is a function that in the
        absence of a database it:
        (i)   Provided a GTF file, it converts the file to a gff format.
        (ii)  Creates DB with the gff file by means of the gffutils library.
        Once the db exists it is loaded and the following step is performed:
        (i) Verifies that the database contains intron annotations, if not,
        it attempts to write them.
        """
        if not self.gene_hierarchy_path.exists():
            if not self.genome_database_path.exists():
                if '.gtf' in self.gff_file_path.suffix:
                    self.old_filename = self.gff_file_path.stem
                    self.gff_file_path = Path(f"{self.old_filename}.gff")
                    self.convert_gtf_to_gff()
                    self.environment.logger.info(
                        'the GTF file has been converted into a GFF3 file'
                    )
                    self.environment.logger.info(
                        f'with filename: {self.gff_file_path}'
                    )
                self.create_genome_database()
            if not self.genome_database:
                self.environment.logger.info(
                    "Reading annotations database"
                )
                self.load_genome_database()
            self.database_features = list(self.genome_database.featuretypes())
            self.search_create_intron_annotations()

    def load_genome_database(
            self,
    ) -> None:
        """
        load_genome_database is a function that loads a gffutils database.
        - dbfn: path to the database file
        - keep_order: This is a parameter that is passed when creating
        the FeatureDB instance. When keep_order is set to True, the order
        of attributes in the GFF/GTF file will be preserved when they are
        retrieved from the database.
        """
        try:
            self.genome_database = gffutils.FeatureDB(
                str(self.genome_database_path),
                keep_order=True
            )
        except ValueError as e:
            self.environment.logger.critical(
                f"Incorrect data base path {e}"
            )
            sys.exit()

    def read_genome(
            self,
    ) -> None:
        """
        read_genome is a function that reads a FASTA file and stores
        the masked/unmasked genome sequence in a dictionary.
        The dictionary has the following structure: {chromosome: sequence}
        """
        self.environment.logger.info("Reading genome")
        try:
            self.genome_dictionary = {}
            if self.genome_file_path.suffix == '.gz':
                with gzip.open(self.genome_file_path, mode='rt') as genome_file:  # 'rt' for textmode
                    parsed_genome = SeqIO.parse(genome_file, 'fasta')
                    for record in parsed_genome:
                        self.genome_dictionary[record.id] = str(record.seq)
            else:
                with open(self.genome_file_path, mode='r') as genome_file:
                    parsed_genome = SeqIO.parse(genome_file, 'fasta')
                    for record in parsed_genome:
                        self.genome_dictionary[record.id] = str(record.seq)
        except (ValueError, FileNotFoundError) as e:
            self.environment.logger.critical(
                f"Incorrect genome file path: {e}"
            )
            sys.exit()
        except OSError as e:
            self.environment.logger.critical(
                f"There was an error reading the genome file: {e}"
            )
            sys.exit()

    def create_gene_hierarchy_dictionary(
            self,
    ) -> None:
        """
        Constructs a nested dictionary to represent the hierarchical structure
        and attributes of genes and their related mRNA transcripts based on genomic
        feature data. The created hierarchy is stored in the attribute
        `self.gene_hierarchy_dictionary` and is also saved as a pickle file.
        Note:
        - GFF coordinates are 1-based. Thus, 1 is subtracted from the start position
        to convert them to 0-based coordinates.
        - If the gene is in the negative strand the direction of transcription and
        translation is opposite to the direction the DNA sequence is represented
        meaning that translation starts from the last CDS
        Structure of `self.gene_hierarchy_dictionary`:
        {
        gene_id_1: {
            coordinate: gene_coord_1,
            'chrom': chromosome_1,
            'strand': strand_1,
            'mRNAs': {
                mRNA_id_1: {
                    coordinate: mRNA_coord_1,
                    'strand': strand_1,
                    'structure': [
                        {
                            'id': feature_id_1,
                            'coordinate': feature_coord_1,
                            'frame': frame_1,
                            'type': feature_type_1,
                            'attributes': attribute_dict_1
                        },
                        ...
                    ]
                },
                ...
            }
        },
        ...
        }
        """
        self.environment.logger.info(
            "Fetching gene-hierarchy data from genome annotations"
        )
        for gene in self.genome_database.features_of_type('gene'):
            mrna_transcripts = [
                mrna_transcript for mrna_transcript
                in self.genome_database.children(
                    gene.id,
                    featuretype='mRNA',
                    order_by='start'
                )
            ]
            if mrna_transcripts:
                gene_coordinate = P.open(gene.start - 1, gene.end)
                mrna_dictionary = dict(
                    coordinate=gene_coordinate,
                    chrom=gene.chrom,
                    strand=gene.strand,
                    attributes=dict(gene.attributes),
                    mRNAs=dict()
                )
                for mrna_annot in mrna_transcripts:
                    mrna_coordinate = P.open(mrna_annot.start - 1, mrna_annot.end)
                    mrna_dictionary['mRNAs'][mrna_annot.id] = dict(
                        coordinate=mrna_coordinate,
                        strand=gene.strand,
                        structure=list()
                    )
                    mrna_transcripts_list = list()
                    for child in self.genome_database.children(
                            mrna_annot.id,
                            featuretype=self.features_of_interest,
                            order_by='start'
                    ):
                        child_coordinate = P.open(child.start - 1, child.end)
                        if child_coordinate:
                            mrna_transcripts_list.append(
                                dict(
                                    id=child.id,  # ID attribute
                                    coordinate=child_coordinate,  # ID coordinate starting at 0
                                    frame=child.frame,  # One of '0', '1' or '2'.
                                    type=child.featuretype,   # feature type name
                                    attributes=dict(child.attributes)
                                )   # feature attributes
                            )
                    # if the gene is in the negative strand the direction of
                    # transcription and translation is opposite to the direction the
                    # DNA sequence is represented meaning that translation starts
                    # from the last CDS
                    reverse = self.reverse_sequence_bool(gene_strand=gene.strand)
                    mrna_dictionary['mRNAs'][mrna_annot.id]['structure'] = self.sort_list_intervals_dict(
                        list_dictionaries=mrna_transcripts_list,
                        reverse=reverse,
                    )
                self.gene_hierarchy_dictionary[gene.id] = mrna_dictionary

    def fetch_gene_cdss_set(
            self,
            gene_id: str
    ) -> list[tuple]:
        return list(
            set(
                (coordinate, annotation_structure['frame'])
                for mrna_annotation in self.gene_hierarchy_dictionary[gene_id]['mRNAs'].values()
                for annotation_structure in mrna_annotation['structure']
                for coordinate in (annotation_structure['coordinate'],)
                if (annotation_structure['type'] == 'CDS'
                    and (coordinate.upper - coordinate.lower) >= self.min_exon_length)
            )
        )

    @staticmethod
    def flatten_clusters_representative_exons(
            cluster_list: list
    ) -> list:
        return [
            cluster[0][0] if len(cluster) == 1 else min(cluster, key=lambda x: x[0].upper - x[0].lower)[0]
            for cluster in cluster_list
        ]

    def get_overlapping_clusters(
            self,
            target_coordinates_set: set[tuple],
            threshold: float,
    ) -> list[list[tuple]]:
        processed_intervals = set()
        overlapping_clusters = []
        sorted_coordinates = sorted(target_coordinates_set, key=lambda x: x[0].lower)
        for target_coordinate, evalue in sorted_coordinates:
            if target_coordinate not in processed_intervals:
                processed_intervals.add(target_coordinate)
                processed_intervals, cluster = self.find_interval_clusters(
                    sorted_coordinates=sorted_coordinates,
                    processed_intervals=processed_intervals,
                    cluster=[(target_coordinate, evalue)],
                    threshold=threshold
                )
                overlapping_clusters.append(cluster)
        overlapping_clusters.sort(key=len, reverse=True)
        return overlapping_clusters

    def find_interval_clusters(
            self,
            sorted_coordinates: list,
            processed_intervals: set,
            cluster: list[tuple],
            threshold: float
    ) -> tuple:
        new_cluster = list(cluster)
        for other_coordinate, other_evalue in sorted_coordinates:
            if (other_coordinate not in processed_intervals and all(
                    round(self.min_perc_overlap(
                        intv_i=target_coordinate,
                        intv_j=other_coordinate), 1) >= threshold if threshold > 0 else
                    round(self.min_perc_overlap(
                        intv_i=target_coordinate,
                        intv_j=other_coordinate), 1) > threshold
                    for target_coordinate, evalue in new_cluster
            )):
                new_cluster.append((other_coordinate, other_evalue))
                processed_intervals.add(other_coordinate)
        if new_cluster == cluster:
            return processed_intervals, new_cluster
        else:
            return self.find_interval_clusters(
                sorted_coordinates=sorted_coordinates,
                processed_intervals=processed_intervals,
                cluster=new_cluster,
                threshold=threshold
            )

    @staticmethod
    def compress_directory(
            source_dir: Path
    ) -> None:
        output_filename = source_dir.with_suffix('.tar.gz')
        with tarfile.open(output_filename, "w:gz") as tar:
            base_dir = source_dir.name
            tar.add(source_dir, arcname=base_dir)

    def clear_working_directory(
            self,
    ) -> None:
        if self.gene_hierarchy_path.exists() and self.genome_database_path.exists():
            os.remove(self.genome_database_path)
        if self.draw_event_multigraphs:
            self.compress_directory(source_dir=self.multigraphs_path)
            shutil.rmtree(self.multigraphs_path)
        if self.csv:
            self.compress_directory(source_dir=self.csv_path)
            shutil.rmtree(self.csv_path)

    def prepare_data(
            self,
    ) -> None:
        """
        prepare_data is a wrapper function that:
        (i)   creates the database with the genomic annotations (if it does not exist)
        (ii)  reads or creates the gene hierarchy dictionary
        (iii) reads the genome sequence
        (iv)  connects or creates the results database
        """
        if self._DEBUG_MODE:
            os.makedirs(self.working_directory / 'input', exist_ok=True)
            os.makedirs(self.working_directory / 'output', exist_ok=True)
        if self.draw_event_multigraphs:
            os.makedirs(self.multigraphs_path, exist_ok=True)
        if self.csv:
            os.makedirs(self.csv_path, exist_ok=True)
        self.create_parse_or_update_database()
        self.read_genome()
        if self.gene_hierarchy_path.exists():
            self.gene_hierarchy_dictionary = self.read_pkl_file(
                file_path=self.gene_hierarchy_path
            )
        else:
            self.create_gene_hierarchy_dictionary()
            self.dump_pkl_file(
                out_file_path=self.gene_hierarchy_path,
                records_dictionary=self.gene_hierarchy_dictionary
            )
            os.remove(self.genome_database_path)
        self.database_interface.connect_create_results_database()
        # self.database_interface.create_matches_interdependence_expansions_counts_table()
        if self._DEBUG_MODE:
            self.environment.logger.warning(
                "All tblastx io files will be saved."
                " This may take a large amount of disk space."
            )
