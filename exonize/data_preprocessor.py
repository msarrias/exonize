import gffutils
import os
import pickle
import re
import subprocess
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import portion as P
from datetime import date


class DataPreprocessor(object):
    utr_features = ['five_prime_UTR', 'three_prime_UTR']
    features_of_interest = ['CDS', 'exon', 'intron'] + utr_features

    def __init__(
            self,
            logger_obj: object,
            database_interface: object,
            working_directory: str,
            gff_file_path: str,
            specie_identifier: str,
            genome_file_path: str,
            genome_pickled_file_path: str,
            debug_mode: bool,
            hard_masking: bool,
            evalue_threshold: float,
            self_hit_threshold: float,
            cds_overlapping_threshold: float,
            query_overlapping_threshold: float,
            min_exon_length: int,
            masking_percentage_threshold: float,
    ):
        self.environment = logger_obj
        self.database_interface = database_interface
        self.working_directory = working_directory
        self.gff_file_path = gff_file_path
        self.specie_identifier = specie_identifier
        self.genome_file_path = genome_file_path
        self.genome_pickled_file_path = genome_pickled_file_path
        self.evalue_threshold = evalue_threshold
        self.self_hit_threshold = self_hit_threshold
        self.cds_overlapping_threshold = cds_overlapping_threshold
        self.query_overlapping_threshold = query_overlapping_threshold
        self.min_exon_length = min_exon_length
        self.masking_percentage_threshold = masking_percentage_threshold
        self.timeout_database = database_interface.timeout_database
        self.results_database = database_interface.results_database_path
        self._DEBUG_MODE = debug_mode
        self._HARD_MASKING = hard_masking

        self.database_features = None
        self.old_filename = None
        self.genome_database = None
        self.genome_dictionary = dict()
        self.gene_hierarchy_dictionary = dict()

        # Derived attributes that depend on initial parameters
        self.genome_database_path = os.path.join(
            self.working_directory,
            f'{self.specie_identifier}_genome_annotations.db'
        )
        self.protein_database_path = os.path.join(
            self.working_directory,
            f'{self.specie_identifier}_protein.db'
        )
        self.gene_hierarchy_path = os.path.join(
            self.working_directory,
            f"{self.specie_identifier}_gene_hierarchy.pkl"
        )

    @staticmethod
    def dump_pkl_file(
            out_file_path: str,
            records_dictionary: dict
    ) -> None:
        """
        dump_pkl_file is a function that dumps an object into a pickle file.
        """
        with open(out_file_path, 'wb') as handle:
            pickle.dump(records_dictionary, handle)

    @staticmethod
    def read_pkl_file(file_path: str) -> dict:
        """
        read_pkl_file is a function that reads a pickle file and returns
         the object stored in it.
        """
        with open(file_path, 'rb') as handle:
            read_file = pickle.load(handle)
        return read_file

    @staticmethod
    def sort_list_intervals_dict(list_dictionaries: list, reverse=False) -> list:
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
    def reverse_sequence_bool(gene_strand: str):
        """
        reverse_sequence_bool checks if the gene is in the negative
        strand and returns True if it is.
        :param gene_strand: strand
        """
        return gene_strand == '-'

    @staticmethod
    def get_tblastx_version():
        result = subprocess.run(
            ["tblastx", "-version"],
            capture_output=True,
            text=True
        )
        if result.returncode == 0:
            return result.stdout.strip()
        else:
            raise Exception("Error executing tblastx: " + result.stderr)

    def convert_gtf_to_gff(self,) -> None:
        """
        Convert a GTF file to GFF format using gffread. Flags description:
        -'O': This flag is used to enable the output of the file in GFF3 format.
        -'o': This flag is used to specify the output file name.
        """
        gffread_command = ["gffread", self.old_filename, "-O", "-o", self.gff_file_path]
        subprocess.call(gffread_command)

    def create_genome_database(self,) -> None:
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
            self.environment.logger.info("Creating annotations database")
            self.genome_database = gffutils.create_db(
                self.gff_file_path,
                dbfn=self.genome_database_path,
                force=True,
                keep_order=True,
                merge_strategy='create_unique',
                sort_attribute_values=True,
                disable_infer_genes=True,
                disable_infer_transcripts=True
            )
        except ValueError as e:
            self.environment.logger.exception(f"Incorrect genome annotations file {e}")
            sys.exit()

    def search_create_intron_annotations(self,) -> None:
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
                self.genome_database.update(list(self.genome_database.create_introns()), make_backup=False)
            except ValueError as e:
                self.environment.logger.exception(
                    f"failed to write intron annotations in database. "
                    f"Please provide a GFF3 file with intron annotations {e}"
                )
                sys.exit()

    def create_parse_or_update_database(self,) -> None:
        """
        create_parse_or_update_database is a function that in the
        absence of a database it:
        (i)   Provided a GTF file, it converts the file to a gff format.
        (ii)  Creates DB with the gff file by means of the gffutils library.
        Once the db exists it is loaded and the following step is performed:
        (i) Verifies that the database contains intron annotations, if not,
        it attempts to write them.
        """
        if not os.path.exists(self.genome_database_path):
            if 'gtf' in self.gff_file_path:
                self.old_filename = self.gff_file_path
                self.gff_file_path = f"{self.old_filename.rsplit('.gtf')[0]}.gff"
                self.convert_gtf_to_gff()
                self.environment.logger.info('the GTF file has been converted into a GFF3 file')
                self.environment.logger.info(f'with filename: {self.gff_file_path}')
            self.create_genome_database()
        if not self.genome_database:
            self.environment.logger.info("Reading annotations database")
            self.load_genome_database()
        self.database_features = list(self.genome_database.featuretypes())
        self.search_create_intron_annotations()

    def load_genome_database(self,) -> None:
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
                self.genome_database_path,
                keep_order=True
            )
        except ValueError as e:
            self.environment.logger.exception(
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
        hard_masking_regex = re.compile(pattern='[a-z]')
        self.environment.logger.info("Reading genome")
        if (self.genome_pickled_file_path is not None
                and os.path.exists(path=self.genome_pickled_file_path)):
            self.genome_dictionary = self.read_pkl_file(file_path=self.genome_pickled_file_path)
        else:
            try:
                with open(self.genome_file_path) as genome_file:
                    parsed_genome = SeqIO.parse(genome_file, format='fasta')
                    if self._HARD_MASKING:
                        self.genome_dictionary = {
                            fasta.id: hard_masking_regex.sub(repl='N', string=str(fasta.seq))
                            for fasta in parsed_genome
                        }
                    else:
                        self.genome_dictionary = {
                            fasta.id: str(fasta.seq)
                            for fasta in parsed_genome
                        }
            except (ValueError, FileNotFoundError) as e:
                self.environment.logger.exception(
                    f"Incorrect genome file path {e}"
                )
                sys.exit()
            if self.genome_pickled_file_path is not None:
                self.dump_pkl_file(
                    out_file_path=self.genome_pickled_file_path,
                    records_dictionary=self.genome_dictionary,
                )

    def construct_mrna_sequence(
            self,
            chromosome: str,
            gene_strand: str,
            cds_coordinates_list: list[dict],
    ) -> str:
        """
        construct_mrna_sequence is a function that constructs the mRNA sequence
        based on genomic coordinates of the CDSs and the strand of the gene.
        """
        mrna_sequence = ''
        for cds_coordinate in cds_coordinates_list:
            start = cds_coordinate['coordinate'].lower
            end = cds_coordinate['coordinate'].upper
            cds_sequence = self.genome_dictionary[chromosome][start:end]
            if gene_strand == '-':
                cds_sequence = str(Seq(cds_sequence).reverse_complement())
            mrna_sequence += cds_sequence
        return mrna_sequence

    @staticmethod
    def trim_sequence_to_codon_length(
            sequence: str,
            is_final_cds: bool,
            gene_id: str,
            transcript_id: str
    ) -> str:
        """
        Trim the sequence to ensure its length is a multiple of 3.
        If it's not the final CDS and an overhang is present, an exception is raised.

        :param sequence: The DNA sequence to be checked.
        :param is_final_cds: Boolean indicating if this is the final CDS in the sequence.
        :param gene_id: The gene ID.
        :param transcript_id: The transcript ID.
        :return: Trimmed sequence or the original sequence if no trimming is needed.
        """
        overhang = len(sequence) % 3
        if overhang and is_final_cds:
            return sequence[:-overhang]
        elif overhang:
            raise ValueError(
                f' {gene_id}, {transcript_id}, non-final CDS has an overhang:'
                f' {len(sequence)} is not divisible by 3'
            )
        return sequence

    def construct_peptide_sequences(
            self,
            gene_id: str,
            transcript_id: str,
            mrna_sequence: str,
            cds_coordinates_list: list[dict],
    ) -> tuple[str, list[tuple]]:
        """
        Construct a protein sequence from transcriptomic coordinates and collect
        corresponding CDSs in both DNA and protein formats.
        Given the transcriptomic coordinates and considering the reading frames
        specified in a GFF file, this function constructs a protein sequence while
        managing the intricacies of exon stitching and reading frame maintenance across
        possibly intron-interrupted CDS regions. In the context of a GFF file, feature
        coordinates are generally relative to the genomic sequence, which implies that
        reading frames may be disrupted by introns.
        ----------------
        Example:
            Consider the following genomic coordinates of CDSs and their reading frames:
            cds1: (0,127)      reading frame: 0
            cds2: (4545,4682)  reading frame: 2
            cds3: (6460,6589)  reading frame: 0
            cds4: (7311,7442)  reading frame: 0

            When translated to transcriptomic coordinates (while still referring to
            genomic positions), considering reading frames, the coordinates become:
            cds1: (0, 127), (4545, 4547)
            cds2: (4547, 4682)
            cds3: (6460, 6589)
            cds4: (7311, 7442)
        ----------------
        Note:
            - The first two nucleotides of cds2 complete the last codon of cds1, and thus,
             it is in line with the specified reading frame of 2 for cds2.
            - It is not uncommon that the length of the last CDS is not a multiple of 3.
             This is because the reading frames of different CDSs across transcripts are
             not necessarily aligned. So that the extra nucleotides are necesary for satisfying
             completition of all transcripts.
        """
        mrna_peptide_sequence = ''
        start_coord = 0
        n_coords = len(cds_coordinates_list) - 1
        cds_list_tuples = list()

        for coord_idx, cds_dictionary in enumerate(cds_coordinates_list):
            frame_cds = int(cds_dictionary['frame'])
            start = cds_dictionary['coordinate'].lower
            end = cds_dictionary['coordinate'].upper
            len_coord = end - start
            frame_next_cds = 0
            end_coord = len_coord + start_coord
            if coord_idx != n_coords:
                frame_next_cds = int(cds_coordinates_list[coord_idx + 1]['frame'])
            cds_dna_sequence = self.trim_sequence_to_codon_length(
                sequence=mrna_sequence[start_coord + frame_cds:  end_coord + frame_next_cds],
                is_final_cds=coord_idx == n_coords,
                gene_id=gene_id,
                transcript_id=transcript_id
            )
            cds_peptide_sequence = str(Seq(cds_dna_sequence).translate())
            mrna_peptide_sequence += cds_peptide_sequence
            start_coord = end_coord
            frame_cds = frame_next_cds
            cds_list_tuples.append(
                (gene_id,
                 transcript_id,
                 coord_idx,
                 cds_dictionary['id'],
                 frame_cds,
                 start, end,
                 cds_dna_sequence,
                 cds_peptide_sequence)
            )
        return mrna_peptide_sequence, cds_list_tuples

    def create_gene_hierarchy_dictionary(self) -> None:
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
            "Fetching gene-hierarchy data and writing protein database"
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
                    mRNAs=dict()
                )
                gene_tuples_list_peptide_transcripts = list()
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
                    cds_annotations_list = [cds_annotation for cds_annotation in
                                            mrna_dictionary['mRNAs'][mrna_annot.id]['structure']
                                            if cds_annotation['type'] == 'CDS'
                                            ]
                    mrna_dna_sequence = self.construct_mrna_sequence(
                        chromosome=gene.chrom,
                        gene_strand=gene.strand,
                        cds_coordinates_list=cds_annotations_list,
                    )
                    peptide_sequence, cds_list_tuples = self.construct_peptide_sequences(
                        gene_id=gene.id,
                        transcript_id=mrna_annot.id,
                        mrna_sequence=mrna_dna_sequence,
                        cds_coordinates_list=cds_annotations_list,
                    )
                    self.database_interface.insert_into_cdss_table(
                        database_path=self.protein_database_path,
                        gene_args_tuple_list=cds_list_tuples
                    )
                    gene_tuples_list_peptide_transcripts.append(
                        (gene.id, gene.chrom, gene.strand,
                         gene_coordinate.lower, gene_coordinate.upper,
                         mrna_annot.id,
                         mrna_coordinate.lower, mrna_coordinate.upper,
                         peptide_sequence)
                    )
                self.gene_hierarchy_dictionary[gene.id] = mrna_dictionary

                self.database_interface.insert_into_proteins_table(
                    database_path=self.protein_database_path,
                    gene_args_list_tuple=gene_tuples_list_peptide_transcripts
                )

    def prepare_data(self) -> None:
        """
        prepare_data is a wrapper function that:
        (i)   creates the database with the genomic annotations (if it does not exist)
        (ii)  reads or creates the gene hierarchy dictionary
        (iii) reads the genome sequence
        (iv)  connects or creates the results database
        """
        if self._DEBUG_MODE:
            os.makedirs(os.path.join(self.working_directory, 'input'), exist_ok=True)
            os.makedirs(os.path.join(self.working_directory, 'output'), exist_ok=True)
        self.create_parse_or_update_database()
        self.database_interface.create_protein_table(
            database_path=self.protein_database_path
        )
        self.read_genome()
        if os.path.exists(self.gene_hierarchy_path):
            self.gene_hierarchy_dictionary = self.read_pkl_file(file_path=self.gene_hierarchy_path)
        else:
            self.create_gene_hierarchy_dictionary()
            self.dump_pkl_file(
                out_file_path=self.gene_hierarchy_path,
                records_dictionary=self.gene_hierarchy_dictionary
            )
        self.database_interface.connect_create_results_database()
        self.database_interface.insert_pipeline_settings(
            settings_list=[
                self.specie_identifier,
                date.today(),
                self.get_tblastx_version(),
                self.genome_file_path,
                self.gff_file_path,
                self.results_database,
                self.evalue_threshold,
                self.self_hit_threshold,
                self.cds_overlapping_threshold,
                self.query_overlapping_threshold,
                self.min_exon_length,
                self.masking_percentage_threshold
            ]
        )
        if self._DEBUG_MODE:
            self.environment.logger.warning(
                "All tblastx io files will be saved."
                " This may take a large amount of disk space."
            )
