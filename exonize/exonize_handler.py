import cProfile
import gc
import gffutils
import os
import pickle
import random
import re
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from collections.abc import Sequence, Iterator
from typing import Any, Union

import networkx as nx
import portion as P
import sqlite3
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from exonize.profiling import get_run_performance_profile, PROFILE_PATH
from exonize.sqlite_utils import *
from exonize.utils import *
from exonize.environment_setup import EnvironmentSetup


class Exonize(object):
    stop_codons = ["TAG", "TGA", "TAA"]
    UTR_features = ['five_prime_UTR', 'three_prime_UTR']
    feat_of_interest = ['CDS', 'exon', 'intron'] + UTR_features

    def __init__(
            self,
            gff_file_path,
            genome_path,
            specie_identifier,
            enable_debug,
            hard_masking,
            soft_force,
            hard_force,
            evalue_threshold,
            sleep_max_seconds,
            min_exon_length,
            cds_overlapping_threshold,
            masking_perc_threshold,
            self_hit_threshold,
            timeout_db,
            genome_pickled_file_path,
    ):
        self._DEBUG_MODE = enable_debug
        self.SOFT_FORCE = soft_force
        self.HARD_FORCE = hard_force
        self.HARD_MASKING = hard_masking

        self.gff_file_path = gff_file_path
        self.genome_path = genome_path
        self.genome_pickled_file_path = genome_pickled_file_path
        self.specie_identifier = specie_identifier
        self.evalue = evalue_threshold
        self.cds_overlapping_threshold = cds_overlapping_threshold
        self.masking_perc_threshold = masking_perc_threshold
        self.self_hit_threshold = self_hit_threshold
        self.min_exon_len = min_exon_length
        self.secs = sleep_max_seconds
        self.timeout_db = timeout_db

        self.genome = None
        self.gene_hierarchy_dict = None
        self.db_features = None
        self.old_filename = None
        self.db = None

        # Derived attributes that depend on initial parameters
        self.working_dir = f'{self.specie_identifier}_exonize'
        self.db_path = os.path.join(self.working_dir, f'{self.specie_identifier}_genome_annotations.db')
        self.protein_db_path = os.path.join(self.working_dir, f'{self.specie_identifier}_protein.db')
        self.results_db = os.path.join(self.working_dir, f'{self.specie_identifier}_results.db')
        self.gene_hierarchy_path = os.path.join(self.working_dir, f"{self.specie_identifier}_gene_hierarchy.pkl")

        # Initialize other internal attributes
        self.__neither, self.__query, self.__target, self.__target_full, self.__target_insertion = 0, 0, 0, 0, 0
        self.__both = 0
        self.__annot_target_start, self.__annot_target_end, self.__target_t = None, None, None
        self.__query_CDS, self.__target_CDS, self.__query_CDS_frame, self.__target_CDS_frame = "-", "-", " ", " "
        self.__found = False
        self.__tuples_full_length_duplications = list()
        self.__tuples_insertion_duplications = list()
        self.__tuples_truncation_events = list()

        # Set up environment
        self.log = EnvironmentSetup(**self.__dict__)
        self.log.setup_environment()

    @staticmethod
    def dump_pkl_file(out_filepath: str, obj: dict) -> None:
        """
        dump_pkl_file is a function that dumps an object into a pickle file.
        """
        with open(out_filepath, 'wb') as handle:
            pickle.dump(obj, handle)

    @staticmethod
    def dump_fasta_file(out_filepath: str, seq_dict: dict) -> None:
        """
        dump_fasta_file is a function that dumps a dictionary with sequences into a FASTA file.
        """
        with open(out_filepath, "w") as handle:
            for annot_id, annot_seq in seq_dict.items():
                record = SeqRecord(Seq(annot_seq), id=str(annot_id), description='')
                SeqIO.write(record, handle, "fasta")

    @staticmethod
    def read_pkl_file(filepath: str) -> dict:
        """
        read_pkl_file is a function that reads a pickle file and returns the object stored in it.
        """
        with open(filepath, 'rb') as handle:
            read_file = pickle.load(handle)
        return read_file

    @staticmethod
    def reverse_sequence_bool(gene_strand_: str):
        """
        reverse_sequence_bool checks if the gene is in the negative strand and returns True if it is.
        :param gene_strand_: strand
        """
        return gene_strand_ == '-'

    @staticmethod
    def get_shorter_longer_interv(a: P.Interval, b: P.Interval) -> tuple:
        """
        Given two intervals, the function
        get_shorter_longer_interv returns the smaller
        and the larger interval in length.
        """
        return (a, b) if a.upper - a.lower < b.upper - b.lower else (b, a)

    @staticmethod
    def get_overlap_percentage(a: P.Interval, b: P.Interval) -> float:
        """
        Given two intervals, the function
        get_overlap_percentage returns the percentage
        of the overlapping region relative to an interval b.
        """
        intersection = a & b
        return (intersection.upper - intersection.lower) / (b.upper - b.lower) if intersection else 0

    @staticmethod
    def sort_list_intervals_dict(list_dicts: list, reverse=False) -> list:
        """
        sort_list_intervals_dict is a function that sorts a list of dictionaries based on the coordinates of the
        intervals present in the dictionaries. The list is sorted in ascending order by default.
        """
        return sorted(list_dicts, key=lambda x: (x['coord'].lower, x['coord']), reverse=reverse)

    def create_parse_or_update_database(self) -> None:
        """
        create_parse_or_update_database is a function that in the absence of a database it:
        (i)   Provided a GTF file, it converts the file to a gff format.
        (ii)  Creates DB with the gff file by means of the gffutils library.
        Once the db exists it is loaded and the following step is performed:
        (i) Verifies that the database contains intron annotations, if not, it attempts to write them.
        """
        def convert_gtf_to_gff() -> None:
            """
            Convert a GTF file to GFF format using gffread. Flags description:
            -'O': This flag is used to enable the output of the file in GFF3 format.
            -'o': This flag is used to specify the output file name.
            """
            gffread_command = ["gffread", self.old_filename, "-O", "-o", self.gff_file_path]
            subprocess.call(gffread_command)

        def create_genome_database() -> None:
            """
            create_genome_database is a function that creates a gffutils database from a GFF3 file.
            Args:
            - dbfn: path to the database file
            - force: if True, the database will be overwritten if it already exists
            - keep_order: if True, the order of the features in the GFF file will be preserved
            - merge_strategy: if 'create_unique', the database will be created with unique IDs
            - sort_attribute_values: if True, the attribute values will be sorted
            - disable_infer_genes: if True, the function will not attempt to automatically infer gene features
            - disable_infer_transcripts: if True, the function will not attempt to automatically infer transcript
              features
            """
            try:
                self.log.logger.info("Creating annotations database")
                self.db = gffutils.create_db(self.gff_file_path,
                                             dbfn=self.db_path,
                                             force=True,
                                             keep_order=True,
                                             merge_strategy='create_unique',
                                             sort_attribute_values=True,
                                             disable_infer_genes=True,
                                             disable_infer_transcripts=True)
            except ValueError as e:
                self.log.logger.exception(f"Incorrect genome annotations file {e}")
                sys.exit()

        def search_create_intron_annotations() -> None:
            """
            search_create_intron_annotations is a function that verifies that the gffutils database contains intron
            annotations, if not, it attempts to write them. Some of the db.update() parameters description:
            - make_backup: if True, a backup of the database will be created before updating it.
            - id_spec: dictionary with the following structure: {feature_type: [intron_id]} where the intron_id function
              returns the intron identifier based on the feature attribute (self.id_spec_attribute) present in all
              annotations in the gff file. Common choices are: "ID" or "Parent".
            """
            if 'intron' not in self.db_features:
                self.log.logger.info(
                    "The GFF file does not contain intron annotations - "
                    "attempting to write intron annotations in database",
                )
                try:
                    self.db.update(list(self.db.create_introns()), make_backup=False)
                except ValueError as e:
                    self.log.logger.exception(f"failed to write intron annotations in database. "
                                          f"Please provide a GFF3 file with intron annotations {e}")
                    sys.exit()
        if not os.path.exists(self.db_path):
            if 'gtf' in self.gff_file_path:
                self.old_filename = self.gff_file_path
                self.gff_file_path = f"{self.old_filename.rsplit('.gtf')[0]}.gff"
                convert_gtf_to_gff()
                self.log.logger.info('the GTF file has been converted into a GFF3 file')
                self.log.logger.info(f'with filename: {self.gff_file_path}')
            create_genome_database()
        if not self.db:
            self.log.logger.info("Reading annotations database")
            self.load_db()
        self.db_features = list(self.db.featuretypes())
        search_create_intron_annotations()

    def load_db(self) -> None:
        """
        load_db is a function that loads a gffutils database.
        - dbfn: path to the database file
        - keep_order: This is a parameter that is passed when creating the FeatureDB instance. When keep_order is set to
        True, the order of attributes in the GFF/GTF file will be preserved when they are retrieved from the database.
        """
        try:
            self.db = gffutils.FeatureDB(self.db_path, keep_order=True)
        except ValueError as e:
            self.log.logger.exception(f"Incorrect data base path {e}")
            sys.exit()

    def read_genome(self) -> None:
        """
        read_genome is a function that reads a FASTA file and stores the masked/unmasked genome sequence in a
        dictionary.
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
                self.dump_pkl_file(self.genome_pickled_file_path, self.genome)

    def create_gene_hierarchy_dict(self) -> None:
        """
        Constructs a nested dictionary to represent the hierarchical structure and attributes of genes and
        their related mRNA transcripts based on genomic feature data. The created hierarchy is stored in the attribute
        `self.gene_hierarchy_dict` and is also saved as a pickle file.
        Note:
        - GFF coordinates are 1-based. Thus, 1 is subtracted from the start position to convert them to 0-based
          coordinates.
        - If the gene is in the negative strand the direction of transcription and translation is opposite to the
        direction the DNA sequence is represented meaning that translation starts from the last CDS
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
        def construct_mRNA_sequence(chrom_, strand_, coords_):
            """
            construct_mRNA_sequence is a function that constructs the mRNA sequence based on genomic coordinates of the
            CDSs and the strand of the gene.
            """
            seq = ''
            for coord_ in coords_:
                start, end = coord_['coord'].lower, coord_['coord'].upper
                exon_ = self.genome[chrom_][start:end]
                if strand_ == '-':
                    exon_ = str(Seq(exon_).reverse_complement())
                seq += exon_
            return seq

        def check_for_overhangs(seq: str, cds_idx: int, n_CDSs: int, gene_id: str, trans_id: str) -> str:
            overhang = len(seq) % 3
            if overhang:
                if cds_idx == (n_CDSs - 1):
                    seq = seq[:-overhang]
                else:
                    self.log.logger.error(f'check here: {gene_id}, {trans_id}')
                    sys.exit()
            return seq

        def construct_protein_sequence(gene_id, trans_id_, mRNA_seq_, coords_):
            """
            Construct a protein sequence from transcriptomic coordinates and collect corresponding CDSs in both DNA and
            protein formats.
            Given the transcriptomic coordinates and considering the reading frames specified in a GFF file, this
            function constructs a protein sequence while managing the intricacies of exon stitching and reading frame
            maintenance across possibly intron-interrupted CDS regions.
            In the context of a GFF file, feature coordinates are generally relative to the genomic sequence, which
            implies that reading frames may be disrupted by introns.
            ----------------
            Example:
                Consider the following genomic coordinates of CDSs and their reading frames:
                cds1: (0,127)      reading frame: 0
                cds2: (4545,4682)  reading frame: 2
                cds3: (6460,6589)  reading frame: 0
                cds4: (7311,7442)  reading frame: 0

                When translated to transcriptomic coordinates (while still referring to genomic positions), considering
                reading frames, the coordinates become:
                cds1: (0, 127), (4545, 4547)
                cds2: (4547, 4682)
                cds3: (6460, 6589)
                cds4: (7311, 7442)
            ----------------
            Note:
                - The first two nucleotides of cds2 complete the last codon of cds1, and thus, it is in line with the
                specified reading frame of 2 for cds2.
                - It is not uncommon that the length of the last CDS is not a multiple of 3. This is because the reading
                frames of different CDSs across transcripts are not necessarily aligned. So that the extra nucleotides
                are necesary for satisfying completition of all transcripts.
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
                exon_seq = check_for_overhangs(mRNA_seq_[start_coord + frame:  end_coord + frame_next],
                                               coord_idx, n_coords, gene_id, trans_id_)
                exon_prot = str(Seq(exon_seq).translate())
                prot_seq_ += exon_prot
                start_coord, frame = end_coord, frame_next
                cds_list_tuples.append((gene_id, trans_id_, coord_idx, coord_['id'], frame, s, e, exon_seq, exon_prot))
            return prot_seq_, cds_list_tuples

        self.log.logger.info("Fetching gene-hierarchy data and writing protein database")
        self.gene_hierarchy_dict = dict()
        for gene in self.db.features_of_type('gene'):
            mrna_transcripts = [mRNA_t for mRNA_t in self.db.children(gene.id, featuretype='mRNA', order_by='start')]
            if mrna_transcripts:
                gene_coord = P.open(gene.start - 1, gene.end)
                mrna_dict = dict(coord=gene_coord, chrom=gene.chrom, strand=gene.strand, mRNAs=dict())
                protein_arg_list_tuples = list()
                for mrna_annot in mrna_transcripts:
                    mrna_coord = P.open(mrna_annot.start - 1, mrna_annot.end)
                    mrna_dict['mRNAs'][mrna_annot.id] = dict(coord=mrna_coord, strand=gene.strand, structure=list())
                    temp_mrna_transcript = list()
                    for child in self.db.children(mrna_annot.id, featuretype=self.feat_of_interest, order_by='start'):
                        coord = P.open(child.start - 1, child.end)
                        if coord:
                            temp_mrna_transcript.append(dict(id=child.id,  # ID attribute
                                                             coord=coord,  # ID coordinate starting at 0
                                                             frame=child.frame,  # One of '0', '1' or '2'.
                                                             type=child.featuretype,   # feature type name
                                                             attributes=dict(child.attributes)))  # feature attributes
                    # if the gene is in the negative strand the direction of transcription and translation
                    # is opposite to the direction the DNA sequence is represented meaning that translation starts
                    # from the last CDS
                    reverse = self.reverse_sequence_bool(gene.strand)
                    mrna_dict['mRNAs'][mrna_annot.id]['structure'] = self.sort_list_intervals_dict(
                        temp_mrna_transcript, reverse,
                    )
                    list_cds_annot = [i for i in mrna_dict['mRNAs'][mrna_annot.id]['structure'] if i['type'] == 'CDS']
                    mRNA_seq = construct_mRNA_sequence(gene.chrom, gene.strand, list_cds_annot)
                    prot_seq, CDS_list_tuples = construct_protein_sequence(
                        gene.id, gene.strand, mRNA_seq, list_cds_annot,
                    )
                    insert_into_CDSs(self.protein_db_path, self.timeout_db, CDS_list_tuples)
                    protein_arg_list_tuples.append((gene.id, gene.chrom, gene.strand,
                                                    gene_coord.lower, gene_coord.upper,
                                                    mrna_annot.id, mrna_coord.lower, mrna_coord.upper,
                                                    prot_seq))
                self.gene_hierarchy_dict[gene.id] = mrna_dict
                insert_into_proteins(self.protein_db_path, self.timeout_db, protein_arg_list_tuples)
        self.dump_pkl_file(self.gene_hierarchy_path, self.gene_hierarchy_dict)

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
        create_protein_table(self.protein_db_path, self.timeout_db)
        self.read_genome()
        if os.path.exists(self.gene_hierarchy_path):
            self.gene_hierarchy_dict = self.read_pkl_file(self.gene_hierarchy_path)
        else:
            self.create_gene_hierarchy_dict()
        connect_create_results_db(self.results_db, self.timeout_db)
        if self._DEBUG_MODE:
            self.log.logger.warning("All tblastx io files will be saved. This may take a large amount of disk space.")

    def execute_tblastx(self, query_filename: str, target_filename: str, output_file: str):
        """
        execute_tblastx is a function that executes a tblastx search with the following parameters:
        - tblastx: A BLAST tool that compares the six-frame translations of a nucleotide query sequence
        against the six-frame translations of a nucleotide sequence database.
        - query: query file name
        - subject: subject file name
        - evalue: Expectation value (E) threshold for reporting matches against database sequences.
        - qcov_hsp_perc: Percent query coverage per hsp (high-scoring pair).
        - outfmt: alignment view options - 5: XML output format
        - out: output file name
        """
        tblastx_command = ['tblastx',
                           '-query', query_filename,
                           '-subject', target_filename,
                           '-evalue', str(self.evalue),
                           '-qcov_hsp_perc', str(self.cds_overlapping_threshold * 100),
                           '-outfmt', '5',  # XML output format
                           '-out', output_file]
        subprocess.run(tblastx_command)

    def align_CDS(self, gene_id: str, query_seq: str, hit_seq: str, query_coord: P.Interval, cds_frame: str) -> dict:
        """
        align_CDS is a function that performs a tblastx search between a query sequence (CDS) and a target sequence (gene).
        :param gene_id: gene identifier
        :param query_seq: query sequence (CDS)
        :param hit_seq: target sequence (gene)
        :param query_coord: query coordinates (CDS) interval
        :param cds_frame: frame of the CDS
        :return: dict with the following structure: {hsp_id: {'score': '', 'bits': '','evalue': '',...}}
        """
        def tblastx_with_saved_io(
                ident: str,
                gene_id_: str,
                hit_seq_: str,
                query_seq_: str,
                query_coord_: P.Interval,
                gene_coord_: P.Interval,
                cds_frame_: str,
        ) -> dict:
            """
            tblastx_with_saved_io is a function that executes a tblastx search saving input and output files. This
            function is used for debugging purposes. The input and output files are saved in the following paths:
            - input: input/{ident}_query.fa and input/{gene_id_}_target.fa where ident is the identifier of the query
            sequence (CDS) and gene_id_ is the identifier of the target sequence (gene).
            - output: output/{ident}_output.xml where ident is the identifier of the query sequence (CDS).
            """
            output_file = os.path.join(self.working_dir, f'output/{ident}_output.xml')
            if not os.path.exists(output_file):
                query_filename = os.path.join(self.working_dir, f'input/{ident}_query.fa')
                target_filename = os.path.join(self.working_dir, f'input/{gene_id_}_target.fa')
                if not os.path.exists(target_filename):
                    self.dump_fasta_file(target_filename, {f"{gene_id_}": hit_seq_})
                self.dump_fasta_file(query_filename, {ident: query_seq_})
                self.execute_tblastx(query_filename, target_filename, output_file)
            with open(output_file, "r") as result_handle:
                blast_records = NCBIXML.parse(result_handle)
                try:
                    temp_ = self.parse_tblastx_output(blast_records, query_coord_, gene_coord_, cds_frame_)
                except Exception as e:
                    self.log.logger.exception(e)
                    sys.exit()
            return temp_

        def execute_tblastx_using_tempfiles(
                hit_seq_: str,
                query_seq_: str,
                query_coord_: P.Interval,
                gene_coord_: P.Interval,
                cds_frame_: str,
        ) -> dict:
            """
            execute_tblastx_using_tempfiles is a function that executes a tblastx search using temporary files.
            """
            with tempfile.TemporaryDirectory(dir=self.working_dir) as tmpdirname:
                query_filename = f'{tmpdirname}/query.fa'
                target_filename = f'{tmpdirname}/target.fa'
                self.dump_fasta_file(query_filename, {'query': query_seq_})
                self.dump_fasta_file(target_filename, {'target': hit_seq_})
                output_file = f'{tmpdirname}/output.xml'
                self.execute_tblastx(query_filename, target_filename, output_file)
                with open(output_file, 'r') as result_handle:
                    blast_records = NCBIXML.parse(result_handle)
                    try:
                        temp_ = self.parse_tblastx_output(blast_records, query_coord_, gene_coord_, cds_frame_)
                    except Exception as e:
                        self.log.logger.exception(e)
                        sys.exit()
                return temp_

        chrom = self.gene_hierarchy_dict[gene_id]['chrom']
        gene_coord = self.gene_hierarchy_dict[gene_id]['coord']
        identifier = f'{gene_id}_{chrom}_{str(query_coord.lower)}_{query_coord.upper}'.replace(':', '_')
        if self._DEBUG_MODE:
            temp = tblastx_with_saved_io(identifier, gene_id, hit_seq, query_seq, query_coord, gene_coord, cds_frame)
        else:
            temp = execute_tblastx_using_tempfiles(hit_seq, query_seq, query_coord, gene_coord, cds_frame)
        return temp

    def parse_tblastx_output(
            self,
            blast_records: dict,
            q_coord: P.Interval,
            hit_coord: P.Interval,
            cds_frame_: str,
    ) -> dict:
        """
        the parse_tblastx_output function parses the output of a tblastx search, where a single sequence (CDS)
        has been queried against a single target (gene). Meaning that we only expect to find one BLAST record.
        We only  consider hits that:
            (i)   have an e-value lower than the threshold,
            (ii)  have a minimum alignment length percentage of the query sequence and
            (iii) that do not overlap with the query sequence (self-hit), with a maximum overlap (self.self_hit_threshold)
                  of the query sequence.
        :param blast_records: 'Record' object with blast records
        :param q_coord: query coordinates (CDS) interval
        :param hit_coord: hit coordinates
        :param cds_frame_: frame of the CDS
        :return: dict with the following structure: {target_id {hsp_id: {'score': '', 'bits': '','evalue': '',...}}}
        """
        def get_hsp_dict(
                hsp,
                cds_frame: str,
        ) -> dict:
            return dict(
                cds_frame=cds_frame,
                score=hsp.score,
                bits=hsp.bits,
                evalue=hsp.expect,
                alignment_len=hsp.align_length * 3,
                hit_frame=hsp.frame,
                query_start=hsp.query_start - 1,
                query_end=hsp.query_end,
                target_start=hsp.sbjct_start - 1,
                target_end=hsp.sbjct_end,
                query_aln_prot_seq=hsp.query,
                target_aln_prot_seq=hsp.sbjct,
                query_num_stop_codons=hsp.query.count('*'),
                target_num_stop_codons=hsp.sbjct.count('*'),
                match=hsp.match)

        res_tblastx = dict()
        # since we are performing a single query against a single subject, there's only one blast_record
        for blast_record in blast_records:
            if len(blast_record.alignments) == 0:
                continue
            alignment = blast_record.alignments[0]  # Assuming only one alignment per blast_record
            if len([aln for aln in blast_record.alignments]) > 1:
                self.log.logger.error("More than one alignment per blast_record")
                sys.exit()
            for hsp_idx, hsp_rec in enumerate(alignment.hsps):
                blast_target_coord = P.open((hsp_rec.sbjct_start - 1) + hit_coord.lower, hsp_rec.sbjct_end + hit_coord.lower)
                if (self.get_overlap_percentage(q_coord, blast_target_coord) < self.self_hit_threshold
                        and self.get_overlap_percentage(blast_target_coord, q_coord) < self.self_hit_threshold):
                    res_tblastx[hsp_idx] = get_hsp_dict(hsp_rec, cds_frame_)
        return res_tblastx

    def get_gene_tuple(
            self,
            gene_id: str,
            has_dup_bin: int,
    ) -> tuple:
        """
        get_gene_tuple is a function that given a gene_id, returns a tuple with the following structure:
        (gene_id, chromosome, strand, start_coord, end_coord, 1 if it has a duplication event 0 otherwise)
        """
        gene_coord = self.gene_hierarchy_dict[gene_id]['coord']
        return (gene_id,
                self.gene_hierarchy_dict[gene_id]['chrom'],
                self.gene_hierarchy_dict[gene_id]['strand'],
                gene_coord.lower,
                gene_coord.upper,
                has_dup_bin)

    def get_candidate_CDS_coords(
            self,
            gene_id: str,
    ) -> dict:
        """
        get_candidate_CDS_coords is a function that given a gene_id, collects all the CDS coordinates with a length
        greater than the minimum exon length (self.min_exon_len) across all transcript.
        If there are overlapping CDS coordinates, they are resolved according to the following criteria:
            - a: if they overlap by more than the overlapping threshold (self.cds_overlapping_threshold) of both CDS
             lengths the shortest CDS is selected.
            - b: both CDSs are selected otherwise
        The rationale behind choosing the shortest CDS is that it will reduce exclusion of short duplications due
        to coverage thresholds.
        :param gene_id: gene identifier
        :return: list of representative CDS coordinates for the gene across all transcripts
        """

        def get_first_overlapping_intervals(
            sorted_intervals: list[P.Interval],
        ) -> Union[tuple[P.Interval, P.Interval], tuple[None, None]]:
            """
            Given a list of intervals, returns the first consecutive interval tuples with the overlapping pairs
            described in case a (get_candidate_CDS_coords).
            :param sorted_intervals: list of intervals
            :return: list of tuples with pairs of overlapping intervals
            """
            first_overlap_index = 0
            while first_overlap_index < len(sorted_intervals) - 1:
                current_interval = sorted_intervals[first_overlap_index]
                next_interval = sorted_intervals[first_overlap_index + 1]
                if (
                    self.get_overlap_percentage(current_interval, next_interval) >= self.cds_overlapping_threshold
                    and self.get_overlap_percentage(next_interval, current_interval) >= self.cds_overlapping_threshold
                ):
                    return current_interval, next_interval
                first_overlap_index += 1
            return None, None

        def resolve_overlaps_coords_list(
                sorted_cds_coordinates: list[P.Interval],
        ) -> list[P.Interval]:
            """
            resolve_overlaps_coords_list is a function that given a list of coordinates, resolves overlaps according
            to the criteria described above. Since the pairs described in case b (get_candidate_CDS_coords) are not
            considered to be overlapping, they are both included in the search.
            :param sorted_cds_coordinates: list of coordinates sorted according to the following criteria:
                (i) lower coordinate, (ii) upper coordinate
            :return: list of coordinates without overlaps
            """
            new_list = list()
            # We only process the first pair of overlapping intervals since the resolved overlap could also overlap
            # with the next interval in the list.
            intv_a, intv_b = get_first_overlapping_intervals(sorted_cds_coordinates)
            if all((intv_a, intv_b)):
                shorter, _ = self.get_shorter_longer_interv(intv_a, intv_b)
                # Note: new_list is populated through enumeration of 'sorted_cds_coordinates', so it will also be sorted
                # according to the same criteria used for 'sorted_cds_coordinates'.
                for idx, coord in enumerate(sorted_cds_coordinates):
                    if coord != intv_a:
                        new_list.append(coord)
                    else:
                        new_list.append(shorter)
                        new_list.extend(sorted_cds_coordinates[idx + 2:])
                        return resolve_overlaps_coords_list(new_list)
            return sorted_cds_coordinates

        cds_coords_and_frames = list(set(
            (coordinate, annotation_structure['frame'])
            for mrna_annotation in self.gene_hierarchy_dict[gene_id]['mRNAs'].values()
            for annotation_structure in mrna_annotation['structure']
            for coordinate in (annotation_structure['coord'],)
            if (
                annotation_structure['type'] == 'CDS' and
                (coordinate.upper - coordinate.lower) >= self.min_exon_len
            )
        ))
        if cds_coords_and_frames:
            repr_cds_frame_dict = dict()
            for cds_coord, frame in cds_coords_and_frames:
                if cds_coord in repr_cds_frame_dict:
                    repr_cds_frame_dict[cds_coord]['frame'] += f'_{str(frame)}'
                else:
                    repr_cds_frame_dict[cds_coord] = {'frame': frame}
            sorted_cds_coords_list: list[P.Interval] = sorted(
                [coordinate for coordinate, _ in cds_coords_and_frames],
                key=lambda coordinate: (coordinate.lower, coordinate.upper),
            )
            return {
                'set_coords': resolve_overlaps_coords_list(sorted_cds_coords_list),
                'cds_frame_dict': repr_cds_frame_dict,
            }
        return dict()

    def find_coding_exon_duplicates(
            self,
            gene_id: str,
    ) -> None:
        """
        find_coding_exon_duplicates is a function that given a gene_id, performs a tblastx for each representative CDS
        (see get_candidate_CDS_coords). If the tblastx search returns hits, they are stored in the "results" database,
        otherwise the gene is recorded as having no duplication event.
        :param gene_id: gene identifier
        """
        def check_for_masking(
                chrom_: str,
                gene_id_: str,
                seq_: str,
                coord_: P.Interval,
                type_='gene',
        ):
            """
            check_for_masking is a function that checks if it is hardmasked. If the sequence is a gene and the percentage
            of hardmasking is greater than the threshold (self.masking_perc_threshold) the CDS duplication search is aborted
            and the gene is recorded as having no duplication event. If the sequence is a CDS and the percentage of
            hardmasking is greater than the threshold, the CDS is not queried. Logs for hardmasked genes and CDS are stored
            in the logs attribute and later dumped into a file.
            :param chrom_: chromosome identifier
            :param gene_id_: gene identifier
            :param seq_: sequence
            :param coord_: coordinates
            :param type_: type of sequence (gene or CDS)
            """

            try:
                masking_perc = round(seq_.count('N') / len(seq_), 3)
                if masking_perc > self.masking_perc_threshold:
                    seq_ = ''
                    if type_ == 'gene':
                        self.log.logger.file_only_info(f'Gene {gene_id_} in chromosome {chrom_} and coordinates {str(coord_.lower)}, '
                                                   f'{str(coord_.upper)} is hardmasked.')
                        insert_gene_ids_table(self.results_db, self.timeout_db, self.get_gene_tuple(gene_id_, 0))
                    if type_ == 'CDS':
                        self.log.logger.file_only_info(f'Gene {gene_id_} - {masking_perc * 100} of CDS {str(coord_.lower)},'
                                                   f' {str(coord_.upper)} located in chromosome {chrom_} is hardmasked.')
                return seq_

            except KeyError as e_:
                self.log.logger.exception(f'Either there is missing a chromosome in the genome file '
                                      f'or the chromosome identifiers in the GFF3 and FASTA files do not match {e_}')
                sys.exit()

        CDS_blast_dict = dict()
        time.sleep(random.randrange(0, self.secs))
        chrom, gene_coord, gene_strand = (self.gene_hierarchy_dict[gene_id]['chrom'],
                                          self.gene_hierarchy_dict[gene_id]['coord'],
                                          self.gene_hierarchy_dict[gene_id]['strand'])
        temp_gs = str(Seq(self.genome[chrom][gene_coord.lower:gene_coord.upper]))
        gene_seq = check_for_masking(chrom, gene_id, temp_gs, gene_coord, type_='gene')
        if gene_seq:
            CDS_coords_dict = self.get_candidate_CDS_coords(gene_id)
            if CDS_coords_dict:
                for cds_coord in CDS_coords_dict['set_coords']:
                    temp_cds = str(Seq(self.genome[chrom][cds_coord.lower:cds_coord.upper]))
                    cds_seq = check_for_masking(chrom, gene_id, temp_cds, cds_coord, type_='CDS')
                    if cds_seq:
                        cds_frame_ = CDS_coords_dict['cds_frame_dict'][cds_coord]['frame']
                        tblastx_o = self.align_CDS(gene_id, cds_seq, gene_seq, cds_coord, cds_frame_)
                        if tblastx_o:
                            CDS_blast_dict[cds_coord] = tblastx_o
                if CDS_blast_dict:
                    self.insert_fragments_table(gene_id, CDS_blast_dict)
                else:
                    insert_gene_ids_table(self.results_db, self.timeout_db, self.get_gene_tuple(gene_id, 0))

    def insert_fragments_table(
            self,
            gene_id: str,
            blast_cds_dict: dict,
    ) -> None:
        """
        insert_fragments_table is a function that given a gene_id and a dictionary with the tblastx results for each
        CDS, inserts in the (i) Genes, and (ii) Fragments table of the results (self.results_db) database.
        These two steps are done in paralell to avoid incomplete data in case of a crash.
        A single entry is recorded in the Genes table, while multiple entries can be recorded in the Fragments table.
        The fragments refer to the tblastx HSPs (high-scoring pairs) that have passed the filtering criteria.
        The fragments tuple has the following structure:
        (gene_id, cds_coord_start, cds_coord_end, hit_query_frame, hit_query_strand, hit_target_frame, hit_target_strand,
         score, bits, evalue, alignment_length, query_start, query_end, target_start, target_end, query_aln_prot_seq,
         target_aln_prot_seq, match (alignment sequence), query_num_stop_codons (count of "*" in target alignment),
         target_num_stop_codons (count of "*" in target alignment))
        :param gene_id: gene identifier
        :param blast_cds_dict: dictionary with the tblastx results. The dictionary has the following structure:
        {cds_coord: {hsp_idx: {'score': '', 'bits': '','evalue': '',...}}}
        """

        def reformat_frame_strand(
                frame: int,
        ) -> tuple:
            """
            reformat_frame_strand is a function that converts the frame to a 0-based index and defines a strand variable
            based on the frame sign.
            """
            n_frame = abs(frame) - 1
            n_strand = '-' if frame < 0 else '+'
            return n_frame, n_strand

        def get_fragment_tuple(
                gene_id_: str,
                cds_coord: P.Interval,
                blast_hits: dict,
                hsp_idx: int,
        ) -> tuple:
            hsp_dict = blast_hits[hsp_idx]
            hit_q_frame, hit_t_frame = hsp_dict['hit_frame']
            hit_q_f, hit_q_s = reformat_frame_strand(hit_q_frame)
            hit_t_f, hit_t_s = reformat_frame_strand(hit_t_frame)
            return (gene_id_, cds_coord.lower, cds_coord.upper, '_'.join(list(hsp_dict['cds_frame'])),
                    hit_q_f, hit_q_s, hit_t_f, hit_t_s,
                    hsp_dict['score'], hsp_dict['bits'], hsp_dict['evalue'],
                    hsp_dict['alignment_len'], hsp_dict['query_start'], hsp_dict['query_end'],
                    hsp_dict['target_start'], hsp_dict['target_end'], hsp_dict['query_aln_prot_seq'],
                    hsp_dict['target_aln_prot_seq'], hsp_dict['match'],
                    hsp_dict['query_num_stop_codons'], hsp_dict['target_num_stop_codons'])

        tuple_list = [get_fragment_tuple(gene_id, cds_coord, blast_hits, hsp_idx)
                      for cds_coord, blast_hits in blast_cds_dict.items() for hsp_idx, hsp_dict in blast_hits.items()]
        insert_fragments_table_param, insert_gene_table_param = insert_fragments_calls()
        with sqlite3.connect(self.results_db, timeout=self.timeout_db) as db:
            cursor = db.cursor()
            cursor.execute(insert_gene_table_param, self.get_gene_tuple(gene_id, 1))
            cursor.executemany(insert_fragments_table_param, tuple_list)
            db.commit()

    def find_overlapping_annot(
            self,
            trans_dict: dict,
            cds_intv: P.Interval,
    ) -> list:
        """
        find_overlapping_annot is a function that given a transcript dictionary and a CDS interval, returns a list of
        tuples with the following structure: (CDS_id, CDS_coord, CDS_frame) for all CDSs that overlap with the query CDS
        interval.
        Note:
        - The set of representative CDS only contains information about the coordinates of the CDSs, not their IDs.Hence,
        if we find hits for a CDS, we need to identify the CDS in the gene transcript. Recall that not all CDS are part
        of all transcripts (e.g. alternative splicing).
        - Following the classification criteria described in 'get_candidate_CDS_coords', we identify the query CDS described in
        case a (get_candidate_CDS_coords), i.e, when they are considerd to overlap.
        """
        return [(i['id'], i['coord'], int(i['frame'])) for i in trans_dict['structure']
                if i['type'] == 'CDS'
                and all([self.get_overlap_percentage(i['coord'], cds_intv) >= self.cds_overlapping_threshold,
                         self.get_overlap_percentage(cds_intv, i['coord']) >= self.cds_overlapping_threshold])]

    def identify_full_length_duplications(self) -> None:
        """
        identify_full_length_duplications is a function that identifies full-length duplications following our classification
        model. The function iterates over all representative tblastx hits and for each transcript associated with the gene harboring
        the event it identifies the following events:
        - I. Full exon duplication: the match fully (following our coverage criteria) overlaps with a CDS annotation.
        - II. Insertion: the match is found within a larger CDS.
        - III. Deactivation or unnanotated: the match is found in an intron or UTR.
        - IV. Trunctation: the match spans more than one annotation (e.g., CDS, intron, UTR).
        The identify_full_length_duplications function also looks for:
         - I. Obligate events: defined as events where the query CDS and the target CDS are included within the same
          transcript.
        - II. Neither events: defined as events where the query CDS is not found in the transcript and the target
        figures a trunctation or deactivation.
        The identified events are stored in the results database (self.results_db) in the tables:
         - ObligatoryEvents: identifies both, query and target CDSs. A single record is stored per event.
         - TruncationEvents: identifies all covered annotations. The number of records per event will correspond to the
         number of annotations covered by the event.
         - FullLengthDuplications tables: identifies full-length duplications. The number of records per event will correspond
         to the number of transcripts associated with the gene harboring the event.
        """
        def initialize_variables() -> None:
            """
            initializes variables used in the identify_full_length_duplications function
            """
            self.__neither, self.__query, self.__target, self.__target_full, self.__target_insertion = 0, 0, 0, 0, 0
            self.__both = 0
            self.__annot_target_start, self.__annot_target_end, self.__target_t = None, None, None
            self.__query_CDS, self.__target_CDS, self.__query_CDS_frame, self.__target_CDS_frame = "-", "-", " ", " "
            self.__found = False

        def initialize_list_of_tuples() -> None:
            """
            initializes the list of tuples used to store the identified events in the identify_full_length_duplications function
            """
            self.__tuples_full_length_duplications, self.__tuples_obligatory_events, self.__tuples_truncation_events = list(), list(), list()

        def identify_query(
                trans_dict_: dict,
                cds_intv_: P.Interval,
        ) -> None:
            """
            identify_query is a function that identifies the tblastx query CDS in the gene transcript. A transcript
            cannot have overlapping CDSs. If the query CDS overlaps with more than one CDS, the program exits.
            """
            query_only_ = self.find_overlapping_annot(trans_dict_, cds_intv_)
            if len(query_only_) > 1:
                self.log.logger.error(f'Overlapping query CDSs: {query_only_}, please review your GFF3 file')
                sys.exit()
            elif query_only_:
                self.__query_CDS, _, self.__query_CDS_frame = query_only_[0]
                self.__query = 1
                self.__target_t = 'QUERY_ONLY'

        def target_out_of_mRNA(
                trans_coord_: P.Interval,
                mrna_: str,
                row_: list[tuple],
        ) -> bool:
            """
            target_out_of_mRNA is a function that identifies tblastx hits that are outside the mRNA transcript.
            """
            fragment_id_, gene_id_, gene_s_, _, cds_s_, cds_e_, query_s_, query_e_, target_s_, target_e_, evalue_ = row_
            target_intv_ = P.open(target_s_ + gene_s_, target_e_ + gene_s_)
            if target_intv_.upper < trans_coord_.lower or trans_coord_.upper < target_intv_.lower:
                if (self.__query + self.__target) == 0:
                    self.__neither = 1
                    self.__tuples_full_length_duplications.append((fragment_id_, gene_id_, mrna_,
                                                                   cds_s_, cds_e_, self.__query_CDS, query_s_, query_e_, "OUT_OF_MRNA",
                                                                   self.__target_CDS, self.__annot_target_start, self.__annot_target_end,
                                                                   target_s_, target_e_,
                                                                   self.__neither, self.__query, self.__target, self.__both, evalue_))
                    return True
            return False

        def indetify_full_target(trans_dict_: dict[dict], target_intv_: P.Interval) -> None:
            """
            indetify_full_target is a function that identifies tblastx hits that are full-length duplications as described
            in self.find_overlapping_annot
            """
            target_only_ = self.find_overlapping_annot(trans_dict_, target_intv_)
            if len(target_only_) > 1:
                self.log.logger.error(f'overlapping target CDSs: {target_only_}')
                sys.exit()
            elif target_only_:
                self.__target_CDS, t_CDS_coord_, self.__target_CDS_frame = target_only_[0]
                self.__target_full = 1
                self.__found = True
                self.__target_t = "FULL"
                self.__annot_target_start, self.__annot_target_end = t_CDS_coord_.lower, t_CDS_coord_.upper

        def indentify_insertion_target(trans_dict_: dict, target_intv_: P.Interval) -> None:
            """
            Identifies tblastx hits that are insertion duplications these can be:
            - INS_CDS: if the insertion is in a CDS
            - INS_UTR: if the insertion is in an untralated region
            - DEACTIVATED: if the insertion is in an intron
            """
            def filter_structure_by_interval_and_type(structure: dict, t_intv_, annot_type: str) -> list:
                return [(i['id'], i['coord']) for i in structure
                        if (i['coord'].contains(t_intv_) and annot_type in i['type'])]

            insertion_CDS_ = filter_structure_by_interval_and_type(trans_dict_['structure'], target_intv_, 'CDS')
            if insertion_CDS_:
                self.__target_CDS, t_CDS_coord_ = insertion_CDS_[0]
                self.__target_insertion = 1
                self.__target_t = "INS_CDS"
                self.__found = True
                self.__annot_target_start, self.__annot_target_end = t_CDS_coord_.lower, t_CDS_coord_.upper
            else:
                insertion_UTR_ = filter_structure_by_interval_and_type(trans_dict_['structure'], target_intv_, 'UTR')
                if insertion_UTR_:
                    self.__target_CDS, t_CDS_coord_ = insertion_UTR_[0]
                    self.__target_t = "INS_UTR"
                    self.__found = True
                    self.__annot_target_start, self.__annot_target_end = t_CDS_coord_.lower, t_CDS_coord_.upper
                else:
                    insertion_intron_ = filter_structure_by_interval_and_type(trans_dict_['structure'], target_intv_, 'intron')
                    if insertion_intron_:
                        self.__target_CDS, t_CDS_coord_ = insertion_intron_[0]
                        self.__target_t = "DEACTIVATED"
                        self.__found = True
                        self.__annot_target_start, self.__annot_target_end = t_CDS_coord_.lower, t_CDS_coord_.upper

        def indentify_truncation_target(trans_dict_: dict[str], mrna_: str, row_: list) -> None:
            """
            Identifies tblastx hits that are truncation duplications. These are hits that span across more than one annotation
            in the transcript architecture. We record a line per annotation that is truncated.
            """
            fragment_id_, gene_id_, gene_s_, _, cds_s_, cds_e_, query_s_, query_e_, target_s_, target_e_, evalue_ = row_
            target_intv_ = P.open(target_s_ + gene_s_, target_e_ + gene_s_)
            intv_dict_ = get_interval_dictionary(trans_dict_['structure'], target_intv_, trans_dict_['coord'])
            if intv_dict_:
                self.__target_t = "TRUNC"
                self.__target_CDS, self.__annot_target_start, self.__annot_target_end = None, None, None
                for seg_b, value in intv_dict_.items():
                    coord_b = value['coord']
                    self.__tuples_truncation_events.append(
                        (fragment_id_, gene_id_, mrna_, trans_dict_['coord'].lower, trans_dict_['coord'].upper,
                         self.__query_CDS, cds_s_, cds_e_, query_s_, query_e_,
                         target_s_, target_e_, value['id'], value['type'],
                         coord_b.lower, coord_b.upper, seg_b.lower, seg_b.upper))

        def identify_obligate_pair(trans_coord_: P.Interval, mrna_: str, row_: list[tuple]) -> None:
            """
            Identifies tblastx hits that are obligate pairs. These are hits where the query and target show as CDSs in the
            transcript in question.
            """
            fragment_id_, gene_id_, _, _, cds_s_, cds_e_, query_s_, query_e_, target_s_, target_e_, evalue_ = row_
            self.__both = 1
            self.__query, self.__target = 0, 0
            self.__tuples_obligatory_events.append((fragment_id_, gene_id_, mrna_, trans_coord_.lower, trans_coord_.upper,
                                                    cds_s_, cds_e_, self.__query_CDS, self.__query_CDS_frame,
                                                    query_s_, query_e_, self.__target_CDS, self.__target_CDS_frame,
                                                    self.__annot_target_start, self.__annot_target_end,
                                                    target_s_, target_e_, self.__target_t))

        def identify_neither_pair() -> None:
            """
            Identifies tblastx hits that are neither pairs. These are hits where the query and target do not show as CDSs
            in the transcript in question.
            """
            self.__neither = 1
            self.__target_t = 'NEITHER'

        def insert_full_length_duplication_tuple(mrna_: str, row_: list[tuple]):
            """
            insert_full_length_duplication_tuple is a function that appends the "row" event to the list of tuples.
            """
            fragment_id_, gene_id_, _, _, cds_s_, cds_e_, query_s_, query_e_, target_s_, target_e_, evalue_ = row_
            self.__tuples_full_length_duplications.append((fragment_id_, gene_id_, mrna_,
                                                           cds_s_, cds_e_, self.__query_CDS, query_s_, query_e_,
                                                           self.__target_t, self.__target_CDS, self.__annot_target_start,
                                                           self.__annot_target_end,
                                                           target_s_, target_e_, self.__neither, self.__query, self.__target,
                                                           self.__both, evalue_))

        def insert_tuples_in_results_db() -> None:
            """
            insert_tuples_in_results_db is a function that inserts the list of tuples collected by the identify_full_length_duplications
            function into the results database.
            """
            instert_full_length_event(self.results_db, self.timeout_db, self.__tuples_full_length_duplications)
            instert_obligatory_event(self.results_db, self.timeout_db, self.__tuples_obligatory_events)
            instert_truncation_event(self.results_db, self.timeout_db, self.__tuples_truncation_events)

        rows = query_filtered_full_duplication_events(self.results_db, self.timeout_db)
        initialize_list_of_tuples()
        for row in rows:
            _, gene_id, gene_s, _, cds_s, cds_e, _, _, target_s, target_e, _ = row
            cds_intv, target_intv = P.open(cds_s, cds_e), P.open(target_s + gene_s, target_e + gene_s)
            for mrna, trans_dict in self.gene_hierarchy_dict[gene_id]['mRNAs'].items():
                initialize_variables()  # we want to classify each transcript independently
                trans_coord = trans_dict['coord']
                # ####### QUERY ONLY - FULL LENGTH #######
                identify_query(trans_dict, cds_intv)
                # ###### CHECK: TARGET REGION NOT IN mRNA #######
                if target_out_of_mRNA(trans_coord, mrna, row):
                    continue
                # ####### TARGET ONLY - FULL LENGTH #######
                indetify_full_target(trans_dict, target_intv)
                if self.__target_full == 0:
                    # ####### INSERTION #######
                    indentify_insertion_target(trans_dict, target_intv)
                    if not self.__found:
                        # ####### TRUNCATION #######
                        indentify_truncation_target(trans_dict, mrna, row)
                self.__target = self.__target_full + self.__target_insertion
                # ####### OBLIGATE PAIR #######
                if self.__query + self.__target == 2:
                    identify_obligate_pair(trans_dict['coord'], mrna, row)
                # ####### NEITHER PAIR #######
                elif self.__query + self.__target == 0:
                    identify_neither_pair()
                insert_full_length_duplication_tuple(mrna, row)
        insert_tuples_in_results_db()

    def assign_event_ids(self, full_matches_list: list[tuple]) -> (list[tuple], set[tuple]):
        def get_average_overlapping_percentage(intv_a: P.Interval, intv_b: P.Interval) -> float:
            return sum([self.get_overlap_percentage(intv_a, intv_b), self.get_overlap_percentage(intv_b, intv_a)]) / 2

        def get_shorter_intv_overlapping_percentage(a: P.Interval, b: P.Interval) -> float:
            """
            get_shorter_intv_overlapping_percentage is a function that given two intervals, returns the percentage of
            overlap of the shorter interval with the longer interval.
            """
            shorter, longer = self.get_shorter_longer_interv(a, b)
            return self.get_overlap_percentage(longer, shorter)

        def get_candidate_reference_dict(coordinates: list[P.Interval],
                                         intv_list: list[tuple[P.Interval, float]]) -> P.Interval:
            cand_ref = [query_intv for query_intv in coordinates
                        if all(get_average_overlapping_percentage(query_intv, intv_i) >= self.cds_overlapping_threshold
                               for intv_i, _ in intv_list)]
            if len(cand_ref) == 1:
                return cand_ref[0]
            elif cand_ref:
                cand_ref = [(cand_ref_intv,
                             sum([get_average_overlapping_percentage(cand_ref_intv, intv_i)
                                  for intv_i, _ in intv_list]) / len(intv_list)) for cand_ref_intv in cand_ref]
                return max(cand_ref, key=lambda x: x[1])[0]
            return P.open(0, 0)

        def get_overlapping_clusters(
            coordinates_set: set[tuple[P.Interval, float]],
            threshold: float,
        ) -> list[list[tuple]]:
            def overlap_condition(coordinate, x):
                perc = get_shorter_intv_overlapping_percentage(coordinate, x)
                return perc >= threshold and x != coordinate

            overlapping_targets_list, skip_events = list(), list()
            for coord, evalue in coordinates_set:
                if coord not in skip_events:
                    temp = [(i, i_evalue) for i, i_evalue in coordinates_set if overlap_condition(coord, i)]
                    if temp:
                        skip_events.extend([coord, *[i for i, _ in temp]])
                        tr = [*temp, (coord, evalue)]
                        tr.sort(key=lambda x: (x[0].lower, x[0].upper))
                        overlapping_targets_list.append(tr)
                    else:
                        overlapping_targets_list.append([(coord, evalue)])
            overlapping_targets_list.sort(key=len, reverse=True)
            return overlapping_targets_list

        def build_reference_dictionary(cds_candidates_dict: dict, overlapping_targets_list: list[list]) -> dict[dict]:
            # this dictionary should be for finding CDS reference and intron reference.
            # This should be applied to the clusters.
            ref_dict = dict()
            for intv_list in overlapping_targets_list:
                # First: let's look for targets overlapping with a CDS
                sorted_CDS_coords_list = sorted(cds_candidates_dict['set_coords'],
                                                key=lambda x: (x.lower, x.upper))
                cand_ref = get_candidate_reference_dict(sorted_CDS_coords_list, intv_list)
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
                        cand_ref = get_candidate_reference_dict(sorted_CDS_coords_list, [(intv_i, i_evalue)])
                        if cand_ref:
                            ref_dict[intv_i] = dict(intv_ref=cand_ref, ref='coding')
                        else:
                            ref_dict[intv_i] = dict(intv_ref=intv_i, ref='non_coding')
            return ref_dict

        def create_events_multigraph(ref_coord_dic: dict,
                                     query_coords_set: set,
                                     records_set: set) -> nx.MultiGraph:
            gene_G = nx.MultiGraph()
            target_coords_set = set([i['intv_ref'] for i in ref_coord_dic.values()])
            gene_G.add_nodes_from(set([(i.lower, i.upper) for i in [*query_coords_set, *target_coords_set]]))
            for event in records_set:
                source = (event[2], event[3])  # exact CDS coordinates
                target = ref_coord_dic[P.open(event[4], event[5])]['intv_ref']  # we need a reference
                ref_des = ref_coord_dic[P.open(event[4], event[5])]['ref']
                gene_G.add_edge(source,
                                (target.lower, target.upper),
                                fragment_id=event[0],
                                query_CDS=(event[2], event[3]),
                                target=(event[4], event[5]),
                                evalue=event[6],
                                event_type=event[7],
                                ref=ref_des,
                                color='black', width=2)
            return gene_G

        def get_events_tuples_from_multigraph(reference_dict: dict,
                                              gene_identif: str,
                                              gene_G: nx.MultiGraph) -> (list[tuple], set[tuple]):
            reference = {i['intv_ref']: i['ref'] for i in reference_dict.values()}
            disconnected_components = list(nx.connected_components(gene_G))
            events_tuples, events_list, event_id_counter, cluster_counter = list(), list(), 0, 0
            for component in disconnected_components:
                temp, comp_event_list = dict(), list()
                for node_x, node_y in component:
                    ref = 'coding'
                    node_intv = P.open(int(node_x), int(node_y))
                    if node_intv in reference:
                        ref = reference[node_intv]
                    temp[node_intv] = [ref,  # either coding/non_coding/coding_non_coding
                                       sum(1 for _ in gene_G.neighbors((node_x, node_y))),  # degree
                                       None,  # cluster_id
                                       event_id_counter]
                temp_list_tuples = set((node, 0) for node in temp.keys())
                node_clusters = [[i[0] for i in cluster] for cluster
                                 in get_overlapping_clusters(temp_list_tuples, 0.01)
                                 if len(cluster) > 1]
                for cluster in node_clusters:
                    for node in cluster:
                        temp[node][2] = cluster_counter
                    cluster_counter += 1
                subgraph = gene_G.subgraph(component)
                for edge in subgraph.edges(data=True):
                    # edge is a tuple (node1, node2, attributes)
                    node1, node2, attributes = edge
                    comp_event_list.append((event_id_counter, attributes['fragment_id']))
                event_id_counter += 1
                events_tuples.extend(comp_event_list)
                events_list.extend([(gene_identif,
                                     value[0],
                                     coord.lower, coord.upper,
                                     *value[1:])
                                    for coord, value in temp.items()])
            return events_tuples, events_list

        genes_events_tuples, genes_events_set = list(), set()
        full_matches_dict = defaultdict(set)
        for match in full_matches_list:
            full_matches_dict[match[1]].add(match)
        for gene_id, records in full_matches_dict.items():
            gene_start = self.gene_hierarchy_dict[gene_id]['coord'].lower
            cds_candidates = self.get_candidate_CDS_coords(gene_id)
            cds_candidates['set_coords'] = set([P.open(i.lower - gene_start, i.upper - gene_start)
                                                for i in cds_candidates['set_coords']])
            query_coordinates, target_coordinates = (set([P.open(i[2], i[3]) for i in records]),
                                                     set([(P.open(i[4], i[5]), i[-2]) for i in records]))
            overlapping_targets = get_overlapping_clusters(target_coordinates, self.cds_overlapping_threshold)
            ref_coord_dict = build_reference_dictionary(cds_candidates, overlapping_targets)
            G = create_events_multigraph(ref_coord_dict, query_coordinates, records)
            gene_events_tuples, gene_events_set = get_events_tuples_from_multigraph(ref_coord_dict, gene_id, G)
            if len(gene_events_tuples) != len(records):
                self.log.logger.exception(f'{gene_id}: {len(gene_events_tuples)} events found, {len(records)} expected.')
            genes_events_tuples.extend(gene_events_tuples)
            genes_events_set.update(gene_events_set)
        if len(genes_events_tuples) != len(full_matches_list):
            self.log.logger.exception(f'{len(genes_events_tuples)} events found, {len(full_matches_list)} expected.')
        return genes_events_tuples, genes_events_set

    def get_identity_and_dna_seq_tuples(self) -> list[tuple]:
        """
        Retrieves DNA sequences from tblastx query and target, computes their DNA and amino acid identity,
        and returns a list of tuples with the structure:
        (DNA_identity, AA_identity, query_dna_seq, target_dna_seq, fragment_id)
        """
        def fetch_dna_sequence(chrom: str,
                               annot_start: int,
                               annot_end: int,
                               pos_annot_s: int,
                               pos_annot_e: int,
                               strand: str) -> str:
            """
            Retrieve a subsequence from a genomic region, reverse complementing it if on the negative strand.
            """
            seq = Seq(self.genome[chrom][annot_start:annot_end][pos_annot_s:pos_annot_e])
            return str(seq.reverse_complement()) if self.reverse_sequence_bool(strand) else str(seq)

        def compute_identity(seq_1: str, seq_2: str) -> float:
            """
            Compute the identity between two sequences (seq_1 and seq_2) using the Hamming distance method.
            """
            # Calculate the Hamming distance and return it
            return round(sum(i == j for i, j in zip(seq_1, seq_2)) / len(seq_2), 3)

        def process_fragment(fragment: list) -> tuple[float, float, str, str, int]:
            """
            process fragment recovers the query/target DNA sequences (since the tblastx alignment is done in amino acids)
            and computes the DNA and amino acid identity. This information is returned in a tuple and later used to update the
            Fragments table.
            """
            (fragment_id, gene_id, gene_start, gene_end, gene_chrom,
             CDS_start, CDS_end, query_start, query_end, target_start, target_end,
             query_strand, target_strand, query_aln_prot_seq, target_aln_prot_seq) = fragment
            query_dna_seq = fetch_dna_sequence(gene_chrom, CDS_start, CDS_end, query_start, query_end, query_strand)
            target_dna_seq = fetch_dna_sequence(gene_chrom, gene_start, gene_end, target_start, target_end,
                                                target_strand)

            if len(query_dna_seq) != len(target_dna_seq):
                raise ValueError(f'{gene_id}: CDS {(CDS_start, CDS_end)} search - sequences must have the same length.')

            return (compute_identity(query_dna_seq, target_dna_seq),
                    compute_identity(query_aln_prot_seq, target_aln_prot_seq),
                    query_dna_seq, target_dna_seq, fragment_id)

        return [process_fragment(fragment) for fragment in query_fragments(self.results_db, self.timeout_db)]

    def run_exonize_pipeline(self) -> None:
        """
        run_exonize_pipeline iterates over all genes in the gene_hierarchy_dict attribute and performs a tblastx search
        for each representative CDS (see get_candidate_CDS_coords).
        The steps are the following:
        - 1. The function checks if the gene has already been processed. If so, the gene is skipped.
        - 2. For each gene, the function iterates over all representative CDSs and performs a tblastx search.
        - 3. The percent_query (hit query coverage) column is inserted in the Fragments table.
        - 4. The Filtered_full_length_events View is created. This view contains all tblastx hits that have passed the
        filtering step.
        - 5. The Mrna_counts View is created. This view contains the number of transcripts associated with each gene.
        - 6. The function creates the Cumulative_counts table. This table contains the cumulative counts of the different
        event types across transcripts.
        - 7. The function collects all the raw concatenated event types (see query_concat_categ_pairs).
        - 8. The function generates the unique event types (see generate_unique_events_list) so that no event type is repeated.
        - 9. The function inserts the unique event types in the Event_categ_full_length_events_cumulative_counts table.
        - 10. The function collects the identity and DNA sequence tuples (see get_identity_and_dna_seq_tuples) and inserts
        them in the Fragments table.
        - 11. The function collects all events in the Full_length_events_cumulative_counts table.
        - 12. The function reconciles the events by assigning a "pair ID" to each event (see assign_pair_ids).
        - 13. The function creates the Exclusive_pairs view. This view contains all the events that follow the mutually exclusive
        category.
        """

        def even_batches(
            data: Sequence[Any],
            number_of_batches: int = 1,
        ) -> Iterator[Sequence[Any]]:
            """
            Given a list and a number of batches, returns 'number_of_batches' consecutive subsets elements of 'data' of
            even size each, except for the last one whose size is the remainder of the division of the length of 'data'
            by 'number_of_batches'.
            """
            # We round up to the upper integer value to guarantee that there will be 'number_of_batches' batches
            even_batch_size = (len(data) // number_of_batches) + 1
            for batch_number in range(number_of_batches):
                batch_start_index = batch_number * even_batch_size
                batch_end_index = min((batch_number + 1) * even_batch_size, len(data))
                yield data[batch_start_index:batch_end_index]

        self.prepare_data()
        gene_ids = list(self.gene_hierarchy_dict.keys())
        processed_gene_ids = set(query_gene_ids_in_res_db(self.results_db, self.timeout_db))
        unprocessed_gene_ids = list(set(gene_ids) - processed_gene_ids)
        if unprocessed_gene_ids:
            gene_n = len(gene_ids)
            self.log.logger.info(f'Starting exon duplication search for {len(unprocessed_gene_ids)}/{gene_n} genes.')
            with tqdm(total=len(unprocessed_gene_ids), position=0, leave=True, ncols=50) as progress_bar:
                # Benchmark without any parallel computation:
                # pr = cProfile.Profile()
                # pr.enable()
                # for gene_id in unprocessed_gene_ids:
                #     self.find_coding_exon_duplicates(gene_id)
                #     progress_bar.update(1)
                # pr.disable()
                # pr.dump_stats(PROFILE_PATH)
                # get_run_performance_profile(PROFILE_PATH)

                # Benchmark with parallel computation using os.fork:
                pr = cProfile.Profile()
                pr.enable()
                gc.collect()
                gc.freeze()
                transactions_pks: set[int]
                status: int
                code: int
                forks: int = 0
                FORKS_NUMBER = os.cpu_count()  # This is pretty greedy, could be changed and put in a config file
                for balanced_genes_batch in even_batches(
                    data=unprocessed_gene_ids,
                    number_of_batches=FORKS_NUMBER,
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
                    if os.fork():
                        forks += 1
                        if forks >= FORKS_NUMBER:
                            _, status = os.wait()
                            code = os.waitstatus_to_exitcode(status)
                            assert code in (os.EX_OK, os.EX_TEMPFAIL, os.EX_SOFTWARE)
                            assert code != os.EX_SOFTWARE
                            forks -= 1
                    else:
                        status = os.EX_OK
                        try:
                            for gene_id in balanced_genes_batch:
                                self.find_coding_exon_duplicates(gene_id)
                            progress_bar.update(len(balanced_genes_batch))
                        except Exception as exception:
                            sys.stdout.write(exception)
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
                pr.dump_stats(PROFILE_PATH)
                get_run_performance_profile(PROFILE_PATH)
        else:
            self.log.logger.info('All genes have been processed. If you want to re-run the analysis, '
                             'consider using the hard-force/soft-force flag')
        insert_percent_query_column_to_fragments(self.results_db, self.timeout_db)
        create_filtered_full_length_events_view(self.results_db, self.timeout_db)
        create_mrna_counts_view(self.results_db, self.timeout_db)
        self.log.logger.info('Classifying events')
        self.identify_full_length_duplications()
        create_cumulative_counts_table(self.results_db, self.timeout_db)
        query_concat_categ_pair_list = query_concat_categ_pairs(self.results_db, self.timeout_db)
        reduced_event_types_tuples = generate_unique_events_list(query_concat_categ_pair_list, -1)
        insert_event_categ_full_length_events_cumulative_counts(self.results_db, self.timeout_db, reduced_event_types_tuples)
        identity_and_sequence_tuples = self.get_identity_and_dna_seq_tuples()
        insert_identity_and_dna_algns_columns(self.results_db, self.timeout_db, identity_and_sequence_tuples)
        full_matches = query_full_events(self.results_db, self.timeout_db)
        fragments_tuples, events_set = self.assign_event_ids(full_matches)
        insert_event_id_column_to_full_length_events_cumulative_counts(self.results_db, self.timeout_db, fragments_tuples)
        insert_events_table(self.results_db, self.timeout_db, events_set)
        create_exclusive_pairs_view(self.results_db, self.timeout_db)
        self.log.logger.info('Process completed successfully.')
