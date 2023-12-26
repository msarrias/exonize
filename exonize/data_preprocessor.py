from exonize.sqlite_utils import *
import gffutils
import os
import pickle
import re
import subprocess
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import portion as P


class DataPreprocessor(object):
    UTR_features = ['five_prime_UTR', 'three_prime_UTR']
    feat_of_interest = ['CDS', 'exon', 'intron'] + UTR_features

    def __init__(self, logger_, **attributes):
        self.log = logger_
        for key, value in attributes.items():
            setattr(self, key, value)
        self.db_features = None
        self.old_filename = None
        self.db = None
        self.genome = None
        self.gene_hierarchy_dict = None

        # Derived attributes that depend on initial parameters
        self.db_path = os.path.join(
            self.working_dir,
            f'{self.specie_identifier}_genome_annotations.db'
        )
        self.protein_db_path = os.path.join(
            self.working_dir,
            f'{self.specie_identifier}_protein.db'
        )
        self.gene_hierarchy_path = os.path.join(
            self.working_dir,
            f"{self.specie_identifier}_gene_hierarchy.pkl"
        )

    @staticmethod
    def dump_pkl_file(out_filepath: str, obj: dict) -> None:
        """
        dump_pkl_file is a function that dumps an object into a pickle file.
        """
        with open(out_filepath, 'wb') as handle:
            pickle.dump(obj, handle)

    @staticmethod
    def read_pkl_file(filepath: str) -> dict:
        """
        read_pkl_file is a function that reads a pickle file and returns
         the object stored in it.
        """
        with open(filepath, 'rb') as handle:
            read_file = pickle.load(handle)
        return read_file

    @staticmethod
    def sort_list_intervals_dict(list_dicts: list, reverse=False) -> list:
        """
        sort_list_intervals_dict is a function that sorts a list
        of dictionaries based on the coordinates of the intervals
        present in the dictionaries.
        The list is sorted in ascending order by default.
        """
        return sorted(
            list_dicts,
            key=lambda x: (x['coord'].lower, x['coord']),
            reverse=reverse
        )

    @staticmethod
    def reverse_sequence_bool(gene_strand_: str):
        """
        reverse_sequence_bool checks if the gene is in the negative
        strand and returns True if it is.
        :param gene_strand_: strand
        """
        return gene_strand_ == '-'

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
            self.log.logger.info("Creating annotations database")
            self.db = gffutils.create_db(
                self.gff_file_path,
                dbfn=self.db_path,
                force=True,
                keep_order=True,
                merge_strategy='create_unique',
                sort_attribute_values=True,
                disable_infer_genes=True,
                disable_infer_transcripts=True
            )
        except ValueError as e:
            self.log.logger.exception(f"Incorrect genome annotations file {e}")
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
        if 'intron' not in self.db_features:
            self.log.logger.info(
                "The GFF file does not contain intron annotations - "
                "attempting to write intron annotations in database",
            )
            try:
                self.db.update(list(self.db.create_introns()), make_backup=False)
            except ValueError as e:
                self.log.logger.exception(
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
        if not os.path.exists(self.db_path):
            if 'gtf' in self.gff_file_path:
                self.old_filename = self.gff_file_path
                self.gff_file_path = f"{self.old_filename.rsplit('.gtf')[0]}.gff"
                self.convert_gtf_to_gff()
                self.log.logger.info('the GTF file has been converted into a GFF3 file')
                self.log.logger.info(f'with filename: {self.gff_file_path}')
            self.create_genome_database()
        if not self.db:
            self.log.logger.info("Reading annotations database")
            self.load_db()
        self.db_features = list(self.db.featuretypes())
        self.search_create_intron_annotations()

    def load_db(self,) -> None:
        """
        load_db is a function that loads a gffutils database.
        - dbfn: path to the database file
        - keep_order: This is a parameter that is passed when creating
        the FeatureDB instance. When keep_order is set to True, the order
        of attributes in the GFF/GTF file will be preserved when they are
        retrieved from the database.
        """
        try:
            self.db = gffutils.FeatureDB(self.db_path, keep_order=True)
        except ValueError as e:
            self.log.logger.exception(f"Incorrect data base path {e}")
            sys.exit()

    def read_genome(self,) -> None:
        """
        read_genome is a function that reads a FASTA file and stores
        the masked/unmasked genome sequence in a dictionary.
        The dictionary has the following structure: {chromosome: sequence}
        """
        hard_masking_regex = re.compile('[a-z]')
        self.log.logger.info("Reading genome file")
        if self.genome_pickled_file_path is not None and os.path.exists(self.genome_pickled_file_path):
            self.genome = self.read_pkl_file(self.genome_pickled_file_path)
        else:
            try:
                with open(self.genome_path) as genome_file:
                    parsed_genome = SeqIO.parse(genome_file, 'fasta')
                    if self._HARD_MASKING:
                        self.genome = {
                            fasta.id: hard_masking_regex.sub('N', str(fasta.seq))
                            for fasta in parsed_genome
                        }
                    else:
                        self.genome = {
                            fasta.id: str(fasta.seq)
                            for fasta in parsed_genome
                        }
            except (ValueError, FileNotFoundError) as e:
                self.log.logger.exception(f"Incorrect genome file path {e}")
                sys.exit()
            if self.genome_pickled_file_path is not None:
                self.dump_pkl_file(
                    out_filepath=self.genome_pickled_file_path,
                    obj=self.genome,
                )

    def construct_mRNA_sequence(
            self,
            chrom_,
            strand_,
            coords_
    ) -> str:
        """
        construct_mRNA_sequence is a function that constructs the mRNA sequence
        based on genomic coordinates of the CDSs and the strand of the gene.
        """
        seq = ''
        for coord_ in coords_:
            start, end = coord_['coord'].lower, coord_['coord'].upper
            exon_ = self.genome[chrom_][start:end]
            if strand_ == '-':
                exon_ = str(Seq(exon_).reverse_complement())
            seq += exon_
        return seq

    def check_for_overhangs(
            self,
            seq: str,
            cds_idx: int,
            n_CDSs: int,
            gene_id: str,
            trans_id: str,
    ) -> str:
        overhang = len(seq) % 3
        if overhang:
            if cds_idx == (n_CDSs - 1):
                seq = seq[:-overhang]
            else:
                self.log.logger.error(
                    f'check here: {gene_id}, {trans_id}'
                )
                sys.exit()
        return seq

    def construct_protein_sequence(
            self,
            gene_id,
            trans_id_,
            mRNA_seq_,
            coords_,
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
        prot_seq_, temp, start_coord = '', list(), 0
        n_coords = len(coords_)
        cds_list_tuples = list()
        for coord_idx, coord_ in enumerate(coords_):
            frame = int(coord_['frame'])
            s, e = coord_['coord'].lower, coord_['coord'].upper
            len_coord, frame_next = e - s, 0
            end_coord = len_coord + start_coord
            if coord_idx != len(coords_) - 1:
                frame_next = int(coords_[coord_idx + 1]['frame'])
            exon_seq = self.check_for_overhangs(
                seq=mRNA_seq_[start_coord + frame:  end_coord + frame_next],
                cds_idx=coord_idx,
                n_CDSs=n_coords,
                gene_id=gene_id,
                trans_id=trans_id_,
            )
            exon_prot = str(Seq(exon_seq).translate())
            prot_seq_ += exon_prot
            start_coord, frame = end_coord, frame_next
            cds_list_tuples.append(
                (gene_id,
                 trans_id_,
                 coord_idx,
                 coord_['id'],
                 frame,
                 s, e,
                 exon_seq,
                 exon_prot)
            )
        return prot_seq_, cds_list_tuples

    def create_gene_hierarchy_dict(self) -> None:
        """
        Constructs a nested dictionary to represent the hierarchical structure
        and attributes of genes and their related mRNA transcripts based on genomic
        feature data. The created hierarchy is stored in the attribute
        `self.gene_hierarchy_dict` and is also saved as a pickle file.
        Note:
        - GFF coordinates are 1-based. Thus, 1 is subtracted from the start position
        to convert them to 0-based coordinates.
        - If the gene is in the negative strand the direction of transcription and
        translation is opposite to the direction the DNA sequence is represented
        meaning that translation starts from the last CDS
        Structure of `self.gene_hierarchy_dict`:
        {
        gene_id_1: {
            'coord': gene_coord_1,
            'chrom': chromosome_1,
            'strand': strand_1,
            'mRNAs': {
                mRNA_id_1: {
                    'coord': mRNA_coord_1,
                    'strand': strand_1,
                    'structure': [
                        {
                            'id': feature_id_1,
                            'coord': feature_coord_1,
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
        self.log.logger.info(
            "Fetching gene-hierarchy data and writing protein database"
        )
        self.gene_hierarchy_dict = dict()
        for gene in self.db.features_of_type('gene'):
            mrna_transcripts = [mRNA_t for mRNA_t
                                in self.db.children(gene.id, featuretype='mRNA', order_by='start')
                                ]
            if mrna_transcripts:
                gene_coord = P.open(gene.start - 1, gene.end)
                mrna_dict = dict(
                    coord=gene_coord,
                    chrom=gene.chrom,
                    strand=gene.strand,
                    mRNAs=dict()
                )
                protein_arg_list_tuples = list()
                for mrna_annot in mrna_transcripts:
                    mrna_coord = P.open(mrna_annot.start - 1, mrna_annot.end)
                    mrna_dict['mRNAs'][mrna_annot.id] = dict(
                        coord=mrna_coord,
                        strand=gene.strand,
                        structure=list()
                    )
                    temp_mrna_transcript = list()
                    for child in self.db.children(
                            mrna_annot.id,
                            featuretype=self.feat_of_interest,
                            order_by='start'
                    ):
                        coord = P.open(child.start - 1, child.end)
                        if coord:
                            temp_mrna_transcript.append(dict(
                                id=child.id,  # ID attribute
                                coord=coord,  # ID coordinate starting at 0
                                frame=child.frame,  # One of '0', '1' or '2'.
                                type=child.featuretype,   # feature type name
                                attributes=dict(child.attributes))   # feature attributes
                            )
                    # if the gene is in the negative strand the direction of
                    # transcription and translation is opposite to the direction the
                    # DNA sequence is represented meaning that translation starts
                    # from the last CDS
                    reverse = self.reverse_sequence_bool(gene_strand_=gene.strand)
                    mrna_dict['mRNAs'][mrna_annot.id]['structure'] = self.sort_list_intervals_dict(
                        list_dicts=temp_mrna_transcript, reverse=reverse,
                    )
                    list_cds_annot = [i for i in
                                      mrna_dict['mRNAs'][mrna_annot.id]['structure']
                                      if i['type'] == 'CDS'
                                      ]
                    mRNA_seq = self.construct_mRNA_sequence(
                        chrom_=gene.chrom,
                        strand_=gene.strand,
                        coords_=list_cds_annot,
                    )
                    prot_seq, CDS_list_tuples = self.construct_protein_sequence(
                        gene_id=gene.id,
                        trans_id_=gene.strand,
                        mRNA_seq_=mRNA_seq,
                        coords_=list_cds_annot,
                    )
                    insert_into_CDSs(db_path=self.protein_db_path,
                                     timeout_db=self.timeout_db,
                                     gene_args_tuple_list=CDS_list_tuples
                                     )
                    protein_arg_list_tuples.append((
                        gene.id, gene.chrom, gene.strand,
                        gene_coord.lower, gene_coord.upper,
                        mrna_annot.id, mrna_coord.lower, mrna_coord.upper,
                        prot_seq)
                    )
                self.gene_hierarchy_dict[gene.id] = mrna_dict
                insert_into_proteins(
                    db_path=self.protein_db_path,
                    timeout_db=self.timeout_db,
                    gene_args_list_tuple=protein_arg_list_tuples
                )
        self.dump_pkl_file(out_filepath=self.gene_hierarchy_path, obj=self.gene_hierarchy_dict)

    def prepare_data(self) -> None:
        """
        prepare_data is a wrapper function that:
        (i)   creates the database with the genomic annotations (if it does not exist)
        (ii)  reads or creates the gene hierarchy dictionary
        (iii) reads the genome sequence
        (iv)  connects or creates the results database
        """
        if self._DEBUG_MODE:
            os.makedirs(os.path.join(self.working_dir, 'input'), exist_ok=True)
            os.makedirs(os.path.join(self.working_dir, 'output'), exist_ok=True)
        self.create_parse_or_update_database()
        create_protein_table(db_path=self.protein_db_path, timeout_db=self.timeout_db)
        self.read_genome()
        if os.path.exists(self.gene_hierarchy_path):
            self.gene_hierarchy_dict = self.read_pkl_file(self.gene_hierarchy_path)
        else:
            self.create_gene_hierarchy_dict()
        connect_create_results_db(db_path=self.results_db, timeout_db=self.timeout_db)
        if self._DEBUG_MODE:
            self.log.logger.warning(
                "All tblastx io files will be saved."
                " This may take a large amount of disk space."
            )
