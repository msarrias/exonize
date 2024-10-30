# ------------------------------------------------------------------------
# Purpose: This module contains the BLASTsearcher class, which is used to
# perform tblastx searches between CDSs and genes.
# ------------------------------------------------------------------------
import os
import random
import subprocess
import sys
import tempfile
from pathlib import Path
import time
import portion as P
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio import AlignIO


class Searcher(object):
    def __init__(
            self,
            data_container: object,
            sleep_max_seconds: int,
            self_hit_threshold: float,
            min_exon_length: int,
            evalue_threshold: float,
            query_coverage_threshold: float,
            exon_clustering_overlap_threshold: float,
            peptide_identity_threshold: float = 0.4,
            fraction_of_aligned_positions: float = 0.9,
            debug_mode: bool = False
    ):
        self.data_container = data_container
        self.database_interface = data_container.database_interface
        self.environment = data_container.environment
        self.evalue_threshold = evalue_threshold
        self.query_coverage_threshold = query_coverage_threshold
        self.fraction_of_aligned_positions = fraction_of_aligned_positions
        self.peptide_identity_threshold = peptide_identity_threshold
        self.exon_clustering_overlap_threshold = exon_clustering_overlap_threshold
        self.sleep_max_seconds = sleep_max_seconds
        self.self_hit_threshold = self_hit_threshold
        self.min_exon_length = min_exon_length

        self._DEBUG_MODE = debug_mode

    @staticmethod
    def dump_fasta_file(
            out_file_path: Path,
            seq_dictionary: dict,
    ) -> None:
        """
        dump_fasta_file is a function that dumps a  dictionary with
        sequences into a FASTA file.
        :param out_file_path: output file path
        :param seq_dictionary: dictionary with sequences with the following structure:
        {sequence_id: sequence}
        """
        with open(out_file_path, "w") as handle:
            for annotation_id, annotation_sequence in seq_dictionary.items():
                record = SeqRecord(
                    Seq(annotation_sequence),
                    id=str(annotation_id),
                    description=''
                )
                SeqIO.write(record, handle, "fasta")

    @staticmethod
    def get_overlap_percentage(
            intv_i: P.Interval,
            intv_j: P.Interval,
    ) -> float:
        """
        Given two intervals, the function get_overlap_percentage returns the percentage
        of the overlapping region relative to an interval j.
        """
        intersection = intv_i & intv_j
        if intersection:
            return (intersection.upper - intersection.lower) / (intv_j.upper - intv_j.lower)
        return 0

    @staticmethod
    def compute_identity(
            sequence_i: str,
            sequence_j: str
    ) -> float:
        """
        Compute the identity between two sequences
        (seq_1 and seq_2) using the Hamming distance method.
        """
        # Calculate the Hamming distance and return it
        if len(sequence_i) != len(sequence_j):
            raise ValueError('Undefined for sequences of unequal length')
        return round(sum(i == j for i, j in zip(sequence_i, sequence_j)) / len(sequence_j), 3)

    @staticmethod
    def reformat_tblastx_frame_strand(
            frame: int,
    ) -> tuple:
        """
        reformat_tblastx_frame_strand is a function that converts the frame to
        a 0-based index and defines a strand variable based on the frame sign.
        :param frame: 6 different frames are possible: +1, +2, +3, -1, -2, -3
        """
        n_frame = abs(frame) - 1
        n_strand = '-' if frame < 0 else '+'
        return n_frame, n_strand

    @staticmethod
    def reverse_sequence_bool(
            strand: str
    ):
        """
        reverse_sequence_bool checks if the gene is in the negative
        strand and returns True if it is.
        :param strand: + or -
        """
        return strand == '-'

    @staticmethod
    def execute_muscle(
            seq_file_path: Path,
            output_file_path: Path
    ):
        muscle_command = [
            "muscle",
            "-in",
            seq_file_path,
            "-out",
            output_file_path
        ]
        subprocess.run(
            muscle_command,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
        )

    def perform_msa(
            self,
            query: str,
            target: str
    ):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = Path(temp_dir, 'input.fasta')
            output_file = Path(temp_dir, 'align.fasta')
            self.dump_fasta_file(
                out_file_path=input_file,
                seq_dictionary={'query': query, 'target': target}
            )
            self.execute_muscle(
                seq_file_path=input_file,
                output_file_path=output_file
            )
            # Read the alignment result from the output file
            if os.path.exists(output_file):
                alignment = AlignIO.read(output_file, "fasta")
                if alignment:
                    return [str(i.seq) for i in alignment]
            else:
                return []

    @staticmethod
    def get_hsp_dictionary(
            hsp,
            cds_frame: str,
    ) -> dict:
        return dict(
            cds_frame=cds_frame,
            score=hsp.score,
            bits=hsp.bits,
            evalue=hsp.expect,
            alignment_lenth=hsp.align_length * 3,
            hit_frame=hsp.frame,
            query_start=hsp.query_start - 1,
            query_end=hsp.query_end,
            target_start=hsp.sbjct_start - 1,
            target_end=hsp.sbjct_end,
            query_aln_prot_seq=hsp.query,
            target_aln_prot_seq=hsp.sbjct,
            query_num_stop_codons=hsp.query.count('*'),
            target_num_stop_codons=hsp.sbjct.count('*'),
            match=hsp.match
        )

    def execute_tblastx(
            self,
            query_file_path: Path,
            target_file_path: Path,
            output_file_path: Path,
            strand: str,
    ):
        """
        execute_tblastx is a function that executes a tblastx
        search with the following parameters:
        - tblastx: A BLAST tool that compares the six-frame translations
         of a nucleotide query sequence
        against the six-frame translations of a nucleotide sequence database.
        - query: query file name
        - subject: subject file name
        - evalue: Expectation value (E) threshold for reporting
         matches against database sequences.
        - qcov_hsp_perc: Percent query coverage per hsp (high-scoring pair).
        - outfmt: alignment view options - 5: XML output format
        - out: output file name
        """
        tblastx_command = [
            'tblastx',
            '-query',
            query_file_path,
            '-subject',
            target_file_path,
            '-evalue',
            str(self.evalue_threshold),
            '-strand',
            strand,
            '-qcov_hsp_perc',
            str(self.query_coverage_threshold * 100),
            '-outfmt',
            '5',  # XML output format
            '-out',
            output_file_path
        ]
        subprocess.run(tblastx_command)

    def parse_tblastx_output(
            self,
            blast_records: dict,
            q_coord: P.Interval,
            hit_coord: P.Interval,
            cds_frame: str,
    ) -> dict:
        """
        the parse_tblastx_output function parses the output of a tblastx search,
         where a single sequence (CDS) has been queried against a single target (gene).
         Meaning that we only expect to find one BLAST record.
        We only  consider hits that:
            (i)   have an e-value lower than the threshold,
            (ii)  have a minimum alignment length percentage of the query sequence and
            (iii) that do not overlap with the query sequence (self-hit),
            with a maximum overlap (self.self_hit_threshold) of the query sequence.
        :param blast_records: 'Record' object with blast records
        :param q_coord: query coordinates (CDS) interval
        :param hit_coord: hit coordinates
        :param cds_frame: frame of the CDS
        :return: dict with the following structure:
         {target_id {hsp_id: {'score': '', 'bits': '','evalue': '',...}}}
        """
        res_tblastx = dict()
        # since we are performing a single query against a single subject,
        # there's only one blast_record
        for blast_record in blast_records:
            if len(blast_record.alignments) == 0:
                continue
            # Assuming only one alignment per blast_record
            alignment = blast_record.alignments[0]
            if len([aln for aln in blast_record.alignments]) > 1:
                self.environment.logger.error(
                    "More than one alignment per blast_record"
                )
                sys.exit()
            for hsp_idx, hsp_record in enumerate(alignment.hsps):
                blast_target_coord = P.open(
                    (hsp_record.sbjct_start - 1) + hit_coord.lower,
                    hsp_record.sbjct_end + hit_coord.lower
                )
                if self.data_container.min_perc_overlap(
                        intv_i=q_coord,
                        intv_j=blast_target_coord
                ) <= self.self_hit_threshold:
                    res_tblastx[hsp_idx] = self.get_hsp_dictionary(
                        hsp=hsp_record,
                        cds_frame=cds_frame
                    )
        return res_tblastx

    def tblastx_with_saved_io(
            self,
            identifier: str,
            gene_id: str,
            strand: str,
            cds_frame: str,
            hit_sequence: str,
            query_sequence: str,
            query_coordinate: P.Interval,
            gene_coordinate: P.Interval
    ) -> dict:
        """
        tblastx_with_saved_io is a function that executes a tblastx
        search saving input and output files. This
        function is used for debugging purposes. The input and output files
         are saved in the following paths:
        - input: input/{ident}_query.fa and input/{gene_id_}_target.fa where
         ident is the identifier of the query
        sequence (CDS) and gene_id_ is the identifier of the target sequence (gene).
        - output: output/{ident}_output.xml where ident is
         the identifier of the query sequence (CDS).
        """
        output_file_path = self.data_container.working_directory / f'{identifier}_output.xml'
        if not os.path.exists(output_file_path):
            query_file_path = self.data_container.working_directory / f'input/{identifier}_query.fa'
            target_file_path = self.data_container.working_directory / f'input/{gene_id}_target.fa'
            if not target_file_path.exists():
                self.dump_fasta_file(
                    out_file_path=target_file_path,
                    seq_dictionary={f"{gene_id}": hit_sequence}
                )
            self.dump_fasta_file(
                out_file_path=query_file_path,
                seq_dictionary={identifier: query_sequence}
            )
            self.execute_tblastx(
                query_file_path=query_file_path,
                target_file_path=target_file_path,
                output_file_path=output_file_path,
                strand=strand
            )
        with open(output_file_path, "r") as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            try:
                tblastx_output_dictionary = self.parse_tblastx_output(
                    blast_records=blast_records,
                    q_coord=query_coordinate,
                    hit_coord=gene_coordinate,
                    cds_frame=cds_frame
                )
            except Exception as e:
                self.environment.logger.exception(e)
                sys.exit()
        return tblastx_output_dictionary

    def execute_tblastx_using_tempfiles(
            self,
            hit_sequence: str,
            query_sequence: str,
            query_coordinate: P.Interval,
            gene_coordinate: P.Interval,
            cds_frame: str,
            strand: str
    ) -> dict:
        """
        execute_tblastx_using_tempfiles is a function that executes
        a tblastx search using temporary files.
        """
        with tempfile.TemporaryDirectory(dir=self.data_container.working_directory) as temporary_directory:
            query_file_path = Path(temporary_directory, 'query.fa')
            target_file_path = Path(temporary_directory, 'target.fa')
            self.dump_fasta_file(
                out_file_path=query_file_path,
                seq_dictionary={'query': query_sequence}
            )
            self.dump_fasta_file(
                out_file_path=target_file_path,
                seq_dictionary={'target': hit_sequence}
            )
            output_file_path = Path(temporary_directory, 'output.xml')
            self.execute_tblastx(
                query_file_path=query_file_path,
                target_file_path=target_file_path,
                output_file_path=output_file_path,
                strand=strand
            )
            with open(output_file_path, 'r') as result_handle:
                blast_records = NCBIXML.parse(result_handle)
                try:
                    tblastx_output_dictionary = self.parse_tblastx_output(
                        blast_records=blast_records,
                        q_coord=query_coordinate,
                        hit_coord=gene_coordinate,
                        cds_frame=cds_frame
                    )
                except Exception as e:
                    self.environment.logger.exception(e)
                    sys.exit()
            return tblastx_output_dictionary

    def get_candidate_cds_coordinates(
            self,
            gene_id: str,
    ) -> dict:
        """
        get_candidate_cds_coordinates is a function that given a gene_id,
        collects all the CDS coordinates with a length greater than the
        minimum exon length across all transcript.
        If there are overlapping CDS coordinates, they are resolved
        according to the following criteria:

        * If they overlap by more than the overlapping threshold
          of both CDS lengths the one with the highest overlapping
          percentage is selected.
        * Both CDSs are selected otherwise

        The rationale behind choosing the interval with the higher percentage
        of overlap is that it will favor the selection of shorter intervals
        reducing the exclusion of short duplications due to coverage thresholds.

        :param gene_id: gene identifier
        :return: list of representative CDS coordinates across all transcripts
        """
        # collect all CDS coordinates and frames across all transcripts
        # we are interested in the frames to account for the unlikely event
        # that two CDS with same coordinates in different transcripts
        # have different frames
        cds_coordinates_and_frames = self.data_container.fetch_gene_cdss_set(
            gene_id=gene_id
        )
        if cds_coordinates_and_frames:
            clusters = self.data_container.get_overlapping_clusters(
                target_coordinates_set=set(
                    (coordinate, None) for coordinate, frame in cds_coordinates_and_frames
                    if coordinate.upper - coordinate.lower >= self.min_exon_length),
                threshold=self.exon_clustering_overlap_threshold
            )
            if clusters:
                representative_cdss = self.data_container.flatten_clusters_representative_exons(
                    cluster_list=clusters,
                )
                representative_cds_frame_dictionary = dict()
                for cds_coordinate, frame in cds_coordinates_and_frames:
                    if cds_coordinate in representative_cdss:
                        if cds_coordinate in representative_cds_frame_dictionary:
                            representative_cds_frame_dictionary[cds_coordinate] += f'_{str(frame)}'
                        else:
                            representative_cds_frame_dictionary[cds_coordinate] = str(frame)
                return dict(
                    candidates_cds_coordinates=representative_cdss,
                    cds_frame_dict=representative_cds_frame_dictionary
                )
        return dict()

    def align_cds(
            self,
            gene_id: str,
            query_sequence: str,
            hit_sequence: str,
            query_coordinate: P.Interval,
            cds_frame: str,
    ) -> dict:
        """
        align_cds is a function that performs a tblastx search between a query
        sequence (CDS) and a target sequence (gene).
        :param gene_id: gene identifier
        :param query_sequence: query sequence (CDS)
        :param hit_sequence: target sequence (gene)
        :param query_coordinate: query coordinates (CDS) interval
        :param cds_frame: frame of the CDS
        :return: dict with the following structure:
         {hsp_id: {'score': '', 'bits': '','evalue': '',...}}
        """
        chromosome = self.data_container.gene_hierarchy_dictionary[gene_id]['chrom']
        gene_coordinate = self.data_container.gene_hierarchy_dictionary[gene_id]['coordinate']
        gene_strand = self.data_container.gene_hierarchy_dictionary[gene_id]['strand']
        tblastx_strand_input = 'plus' if gene_strand == '+' else 'minus'
        identifier = (
            f'{gene_id}_{chromosome}_'
            f'{str(query_coordinate.lower)}_'
            f'{query_coordinate.upper}'
        ).replace(':', '_')
        if self._DEBUG_MODE:
            tblastx_o = self.tblastx_with_saved_io(
                identifier=identifier,
                gene_id=gene_id,
                hit_sequence=hit_sequence,
                query_sequence=query_sequence,
                query_coordinate=query_coordinate,
                gene_coordinate=gene_coordinate,
                cds_frame=cds_frame,
                strand=tblastx_strand_input
            )
        else:
            tblastx_o = self.execute_tblastx_using_tempfiles(
                hit_sequence=hit_sequence,
                query_sequence=query_sequence,
                query_coordinate=query_coordinate,
                gene_coordinate=gene_coordinate,
                cds_frame=cds_frame,
                strand=tblastx_strand_input
            )
        return tblastx_o

    def get_fragment_tuple(
            self,
            gene_id: str,
            cds_coordinate: P.Interval,
            blast_hits: dict,
            hsp_idx: int,
    ) -> tuple:
        hsp_dictionary = blast_hits[hsp_idx]
        hit_query_frame, hit_target_frame = hsp_dictionary['hit_frame']
        hit_query_frame, hit_query_strand = self.reformat_tblastx_frame_strand(
            frame=hit_query_frame
        )
        hit_target_frame, hit_target_strand = self.reformat_tblastx_frame_strand(
            frame=hit_target_frame
        )
        return (
            gene_id,
            cds_coordinate.lower, cds_coordinate.upper,
            '_'.join(list(hsp_dictionary['cds_frame'])),
            hit_query_frame, hit_query_strand,
            hit_target_frame, hit_target_strand,
            hsp_dictionary['score'],
            hsp_dictionary['bits'],
            hsp_dictionary['evalue'],
            hsp_dictionary['alignment_lenth'],
            hsp_dictionary['query_start'],
            hsp_dictionary['query_end'],
            hsp_dictionary['target_start'],
            hsp_dictionary['target_end'],
            hsp_dictionary['query_aln_prot_seq'],
            hsp_dictionary['target_aln_prot_seq'],
            hsp_dictionary['match'],
            hsp_dictionary['query_num_stop_codons'],
            hsp_dictionary['target_num_stop_codons']
        )

    def get_gene_tuple(
            self,
            gene_id: str
    ) -> tuple:
        """
        get_gene_tuple is a function that given a gene_id,
         returns a tuple with the following structure:
        (gene_id, chromosome, strand, start_coord,
        end_coord, 1 if it has a duplication event 0 otherwise)
        """
        gene_coordinate = self.data_container.gene_hierarchy_dictionary[gene_id]['coordinate']
        return (
            gene_id,
            self.data_container.gene_hierarchy_dictionary[gene_id]['chrom'],
            self.data_container.gene_hierarchy_dictionary[gene_id]['strand'],
            len(self.data_container.gene_hierarchy_dictionary[gene_id]['mRNAs']),
            gene_coordinate.lower,
            gene_coordinate.upper
        )

    def local_search(
            self,
            gene_id_list: list[str],
    ) -> None:
        """
        local_search is a function that given a gene_id,
        performs a tblastx for each representative CDS (see get_candidate_cds_coordinates).
        If the tblastx search returns hits, they are stored in the "results" database,
        otherwise the gene is recorded as having no duplication event.
        :param gene_id_list: gene identifier
        """
        for gene_id in gene_id_list:
            blast_hits_dictionary = dict()
            chromosome, gene_coordinate = (
                self.data_container.gene_hierarchy_dictionary[gene_id]['chrom'],
                self.data_container.gene_hierarchy_dictionary[gene_id]['coordinate']
            )
            gene_dna_sequence = str(
                Seq(self.data_container.genome_dictionary[chromosome][gene_coordinate.lower:gene_coordinate.upper])
            )
            cds_coordinates_dictionary = self.get_candidate_cds_coordinates(gene_id=gene_id)
            no_dup = True
            if cds_coordinates_dictionary:
                for cds_coordinate in cds_coordinates_dictionary['candidates_cds_coordinates']:
                    # note that we are not accounting for the frame at this stage, that will be part of
                    # the filtering step (since tblastx alignments account for 3 frames)
                    cds_dna_sequence = str(
                        Seq(self.data_container.genome_dictionary[chromosome][
                            cds_coordinate.lower:cds_coordinate.upper])
                    )
                    cds_frame = cds_coordinates_dictionary['cds_frame_dict'][cds_coordinate]
                    tblastx_o = self.align_cds(
                        gene_id=gene_id,
                        query_sequence=cds_dna_sequence,
                        hit_sequence=gene_dna_sequence,
                        query_coordinate=cds_coordinate,
                        cds_frame=cds_frame
                    )
                    if tblastx_o:
                        blast_hits_dictionary[cds_coordinate] = tblastx_o
                attempt = False
                if blast_hits_dictionary:
                    no_dup = False
                    while not attempt:
                        try:
                            self.populate_fragments_table(
                                gene_id=gene_id,
                                blast_results_dictionary=blast_hits_dictionary
                            )
                            attempt = True
                        except Exception as e:
                            if "locked" in str(e):
                                time.sleep(random.randrange(start=0, stop=self.sleep_max_seconds))
                            else:
                                self.environment.logger.exception(e)
                                sys.exit()
            if no_dup:
                attempt = False
                while not attempt:
                    try:
                        self.database_interface.insert_gene_ids_table(
                            gene_args_tuple=self.get_gene_tuple(
                                gene_id=gene_id
                            )
                        )
                        attempt = True
                    except Exception as e:
                        if "locked" in str(e):
                            time.sleep(random.randrange(start=0, stop=self.sleep_max_seconds))
                        else:
                            self.environment.logger.exception(e)
                            sys.exit()

    def populate_fragments_table(
            self,
            gene_id: str,
            blast_results_dictionary: dict,
    ) -> None:
        """
        populate_fragments_table is a function that given a gene_id and a dictionary
        with the tblastx results for each CDS, inserts in the
        (i) Genes, and (ii) Fragments table of the results (self.results_database_path) database.
        These two steps are done in paralell to avoid incomplete data in case of a crash.
        A single entry is recorded in the Genes table, while multiple entries
        can be recorded in the Fragments table.
        The fragments refer to the tblastx HSPs (high-scoring pairs) that have passed
        the filtering criteria. The fragments tuple has the following structure:
        (gene_id, cds_coord_start, cds_coord_end, hit_query_frame,
         hit_query_strand, hit_target_frame, hit_target_strand,
         score, bits, evalue, alignment_length, query_start,
         query_end, target_start, target_end, query_aln_prot_seq,
         target_aln_prot_seq, match (alignment sequence),
          query_num_stop_codons (count of "*" in target alignment),
         target_num_stop_codons (count of "*" in target alignment))
        :param gene_id: gene identifier
        :param blast_results_dictionary: dictionary with the tblastx results.
        The dictionary has the following structure:
        {cds_coord: {hsp_idx: {'score': '', 'bits': '','evalue': '',...}}}
        """
        tuple_list = [
            self.get_fragment_tuple(
                gene_id=gene_id,
                cds_coordinate=cds_coord,
                blast_hits=blast_hits,
                hsp_idx=hsp_idx
            )
            for cds_coord, blast_hits in blast_results_dictionary.items()
            for hsp_idx, hsp_dictionary in blast_hits.items()
        ]
        self.database_interface.insert_matches(
            gene_args_tuple=self.get_gene_tuple(
                gene_id=gene_id
            ),
            fragments_tuples_list=tuple_list
        )

    def fetch_dna_sequence(
            self,
            chromosome: str,
            annotation_start: int,
            annotation_end: int,
            trim_start: int,
            trim_end: int,
            strand: str
    ) -> str:
        """
        Retrieve a subsequence from a genomic region,
         reverse complementing it if on the negative strand.
        """
        sequence = Seq(
            self.data_container.genome_dictionary[chromosome][annotation_start:annotation_end][trim_start:trim_end]
        )
        if self.reverse_sequence_bool(strand=strand):
            return str(sequence.reverse_complement())
        return str(sequence)

    def process_fragment(
            self,
            fragment: list
    ) -> tuple:
        """
        process fragment recovers the query/target DNA sequences
        (since the tblastx alignment is done in amino acids)
        and computes the DNA and amino acid identity.
        This information is returned in a tuple and later used to update the
        Fragments table.
        """
        (fragment_id, gene_id, gene_start,
         gene_end, gene_chrom, cds_start,
         cds_end, query_start, query_end,
         target_start, target_end,
         query_strand, target_strand,
         query_aln_prot_seq, target_aln_prot_seq) = fragment

        query_dna_seq = self.fetch_dna_sequence(
            chromosome=gene_chrom,
            annotation_start=cds_start,
            annotation_end=cds_end,
            trim_start=query_start,
            trim_end=query_end,
            strand=query_strand
        )
        target_dna_seq = self.fetch_dna_sequence(
            chromosome=gene_chrom,
            annotation_start=gene_start,
            annotation_end=gene_end,
            trim_start=target_start,
            trim_end=target_end,
            strand=target_strand
        )

        if len(query_dna_seq) != len(target_dna_seq):
            self.environment.logger.exception(
                f'{gene_id}: CDS {(cds_start, cds_end)} search '
                f'- sequences must have the same length.'
            )

        return (
            query_dna_seq,
            target_dna_seq,
            self.compute_identity(
                sequence_i=query_dna_seq,
                sequence_j=target_dna_seq
            ),
            self.compute_identity(
                sequence_i=query_aln_prot_seq,
                sequence_j=target_aln_prot_seq
            ),
            fragment_id
        )

    def get_identity_and_dna_seq_tuples(
            self,
            matches_list: list
    ) -> list[tuple]:
        """
        Retrieves DNA sequences from tblastx query and target,
         computes their DNA and amino acid identity,
        and returns a list of tuples with the structure:
        (DNA_identity, AA_identity, query_dna_seq, target_dna_seq, fragment_id)
        """
        return [self.process_fragment(fragment=fragment)
                for fragment in matches_list]

    @staticmethod
    def get_lengths_ratio(
            intvi: P.Interval,
            intvj: P.Interval
    ) -> float:
        short, long = sorted([intvi.upper - intvi.lower, intvj.upper - intvj.lower])
        return short / long

    def fetch_pairs_for_global_alignments(
            self,
            cds_list: list[tuple]
    ) -> set:
        pairs = set()
        for idxi, (coordi, _) in enumerate(cds_list):
            for coordj, _ in cds_list[idxi:]:
                self_overlap = self.data_container.min_perc_overlap(coordi, coordj)
                pair_length_coverage = self.get_lengths_ratio(coordi, coordj)
                if self_overlap == 0 and pair_length_coverage >= self.query_coverage_threshold:
                    pair = tuple(sorted((coordi, coordj), key=lambda x: x.lower - x.upper))
                    pairs.add(pair)
        return pairs

    def cds_global_search(
            self,
            genes_list: list[str]
    ):
        for gene_id in genes_list:
            retain_pairs = set()
            gene_chrom = self.data_container.gene_hierarchy_dictionary[gene_id]['chrom']
            gene_strand = self.data_container.gene_hierarchy_dictionary[gene_id]['strand']
            candidate_cdss = self.get_candidate_cds_coordinates(gene_id=gene_id)
            if candidate_cdss:
                cds_list = [
                    (coord, frame)
                    for coord, frame in self.data_container.fetch_gene_cdss_set(gene_id=gene_id)
                    if coord in candidate_cdss['candidates_cds_coordinates']
                ]
                if len(cds_list) > 1:
                    gene_pairs = self.fetch_pairs_for_global_alignments(
                        cds_list=cds_list
                    )
                    try:
                        transcripts_dictionary = self.data_container.get_transcript_seqs_dict(gene_id=gene_id)
                    except ValueError as e:
                        self.environment.logger.warning(f"{gene_id}: {str(e)}")
                        continue
                    for pair in gene_pairs:
                        coord_i, coord_j = pair
                        seqs_i = self.data_container.recover_prot_dna_seq(
                            cds_coordinate=coord_i,
                            transcript_dict=transcripts_dictionary
                        )
                        seqs_j = self.data_container.recover_prot_dna_seq(
                            cds_coordinate=coord_j,
                            transcript_dict=transcripts_dictionary
                        )
                        if not all([
                            len(seqs) == len(set([frames for *_, frames in seqs]))
                            for seqs in [seqs_i, seqs_j]
                        ]):
                            print('check here', pair)
                        else:
                            for seq_i in seqs_i:
                                dna_i, prot_i, frame_i_tuple = seq_i
                                prev_frame_i, frame_i = frame_i_tuple
                                for seq_j in seqs_j:
                                    dna_j, prot_j, frame_j_tuple = seq_j
                                    prev_frame_j, frame_j = frame_j_tuple
                                    align_dna = self.perform_msa(dna_i, dna_j)
                                    if align_dna:
                                        align_di, align_dj = align_dna
                                        identd = self.compute_identity(
                                            sequence_i=align_di,
                                            sequence_j=align_dj
                                        )
                                        align_prot = self.perform_msa(prot_i, prot_j)
                                        if align_prot:
                                            align_pi, align_pj = align_prot
                                            identp = self.compute_identity(
                                                sequence_i=align_pi,
                                                sequence_j=align_pj
                                            )
                                            align_pos_fract = align_pi.count('-')/len(align_pi)
                                            if (identp > self.peptide_identity_threshold and
                                                    align_pos_fract < (1 - self.fraction_of_aligned_positions)):
                                                retain_pairs.add((
                                                    gene_id, gene_chrom, gene_strand,
                                                    coord_i.lower, coord_i.upper,
                                                    coord_j.lower, coord_j.upper,
                                                    prev_frame_i, frame_i,
                                                    prev_frame_j, frame_j,
                                                    identd, identp,
                                                    align_di, align_dj,
                                                    align_pi, align_pj
                                                ))
                    if retain_pairs:
                        attempt = False
                        while not attempt:
                            try:
                                self.database_interface.insert_global_cds_alignments(
                                    list_tuples=list(retain_pairs)
                                )
                                attempt = True
                            except Exception as e:
                                if "locked" in str(e):
                                    time.sleep(random.randrange(start=0, stop=self.sleep_max_seconds))
                                else:
                                    self.environment.logger.exception(e)
                                    sys.exit()
