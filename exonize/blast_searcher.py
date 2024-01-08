# ------------------------------------------------------------------------
# Purpose: This module contains the BLASTsearcher class, which is used to
# perform tblastx searches between CDSs and genes.
# ------------------------------------------------------------------------
import os
import random
import subprocess
import sys
import tempfile
import time
import portion as P
from typing import Union
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML


class BLASTsearcher(object):
    def __init__(
            self,
            data_container: object,
            masking_percentage_threshold: float,
            sleep_max_seconds: int,
            self_hit_threshold: float,
            min_exon_length: int,
            cds_overlapping_threshold: float,
            evalue_threshold: float,
            debug_mode: bool,
    ):
        self.data_container = data_container
        self.database_interface = data_container.database_interface
        self.environment = data_container.environment
        self.masking_percentage_threshold = masking_percentage_threshold
        self.sleep_max_seconds = sleep_max_seconds
        self.self_hit_threshold = self_hit_threshold
        self.min_exon_length = min_exon_length
        self.cds_overlapping_threshold = cds_overlapping_threshold
        self.evalue_threshold = evalue_threshold
        self._DEBUG_MODE = debug_mode

    @staticmethod
    def dump_fasta_file(
            out_file_path: str,
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
        Given two intervals, the functionget_overlap_percentage returns the percentage
        of the overlapping region relative to an interval j.
        """
        intersection = intv_i & intv_j
        if intersection:
            return (intersection.upper - intersection.lower) / (intv_j.upper - intv_j.lower)
        return 0

    @staticmethod
    def get_candidate_cds(
            intv_i: P.Interval,
            intv_j: P.Interval,
    ) -> P.Interval:
        """
        Given two sorted intervals (sorted by lower bounds), this function returns
        the interval with the smaller length if they have different lengths. If both
        intervals are of equal length, it returns the first interval (intv_i).
        Assumption: intv_i and intv_j are sorted such that intv_i.lower < intv_j.lower.
        param intv_i: the first interval
        param intv_j: the second interval
        Returns: P.interval
        """
        if (intv_j.upper - intv_j.lower) <= (intv_i.upper - intv_i.lower):
            return intv_j
        return intv_i

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

    def check_for_masking(
            self,
            chromosome: str,
            gene_id: str,
            sequence: str,
            coordinate: P.Interval,
            annotation_type='gene',
    ):
        """
        check_for_masking is a function that checks if it is hardmasked.
        If the sequence is a gene and the percentage of hardmasking is greater than
        the threshold (self.masking_perc_threshold) the CDS duplication search is aborted
        and the gene is recorded as having no duplication event.
        If the sequence is a CDS and the percentage of hardmasking is greater than the
        threshold, the CDS is not queried. Logs for hardmasked genes and CDS are stored
        in the logs attribute and later dumped into a file.
        :param chromosome: chromosome identifier
        :param gene_id: gene identifier
        :param sequence: sequence
        :param coordinate: coordinates
        :param annotation_type: type of sequence (gene or CDS)
        """
        try:
            masking_percentage = round(sequence.count('N') / len(sequence), 3)
            if masking_percentage > self.masking_percentage_threshold:
                sequence = ''
                if annotation_type == 'gene':
                    self.environment.logger.file_only_info(
                        f'Gene {gene_id} in chromosome {chromosome} '
                        f'and coordinates {str(coordinate.lower)}, '
                        f'{str(coordinate.upper)} is hardmasked.'
                    )
                    self.database_interface.insert_gene_ids_table(
                        gene_args_tuple=self.get_gene_tuple(
                            gene_id=gene_id,
                            has_duplication_binary=0
                        )
                    )
                if annotation_type == 'CDS':
                    self.environment.logger.file_only_info(
                        f'Gene {gene_id} - {masking_percentage * 100} '
                        f'of CDS {str(coordinate.lower)},'
                        f' {str(coordinate.upper)} located in chromosome {chromosome} '
                        f'is hardmasked.'
                    )
            return sequence

        except KeyError as e:
            self.environment.logger.exception(
                f'Either there is missing a chromosome in the genome file '
                f'or the chromosome identifiers in the GFF3 and FASTA'
                f' files do not match {e}'
            )
            sys.exit()

    def execute_tblastx(
            self,
            query_file_path: str,
            target_file_path: str,
            output_file_path: str,
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
            '-qcov_hsp_perc',
            str(self.cds_overlapping_threshold * 100),
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
                if (
                        self.get_overlap_percentage(
                            intv_i=q_coord,
                            intv_j=blast_target_coord
                        ) < self.self_hit_threshold
                        and
                        self.get_overlap_percentage(
                            intv_i=blast_target_coord,
                            intv_j=q_coord
                        ) < self.self_hit_threshold
                ):
                    res_tblastx[hsp_idx] = self.get_hsp_dictionary(
                        hsp=hsp_record,
                        cds_frame=cds_frame
                    )
        return res_tblastx

    def tblastx_with_saved_io(
            self,
            identifier: str,
            gene_id: str,
            hit_sequence: str,
            query_sequence: str,
            query_coordinate: P.Interval,
            gene_coordinate: P.Interval,
            cds_frame: str,
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
        output_file_path = os.path.join(
            self.data_container.working_directory,
            f'output/{identifier}_output.xml'
        )
        if not os.path.exists(output_file_path):
            query_file_path = os.path.join(
                self.data_container.working_directory,
                f'input/{identifier}_query.fa'
            )
            target_file_path = os.path.join(
                self.data_container.working_directory,
                f'input/{gene_id}_target.fa'
            )
            if not os.path.exists(target_file_path):
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
                output_file_path=output_file_path
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
    ) -> dict:
        """
        execute_tblastx_using_tempfiles is a function that executes
        a tblastx search using temporary files.
        """
        with tempfile.TemporaryDirectory(dir=self.data_container.working_directory) as temporary_directory:
            query_file_path = f'{temporary_directory}/query.fa'
            target_file_path = f'{temporary_directory}/target.fa'
            self.dump_fasta_file(
                out_file_path=query_file_path,
                seq_dictionary={'query': query_sequence}
            )
            self.dump_fasta_file(
                out_file_path=target_file_path,
                seq_dictionary={'target': hit_sequence}
            )
            output_file_path = f'{temporary_directory}/output.xml'
            self.execute_tblastx(
                query_file_path=query_file_path,
                target_file_path=target_file_path,
                output_file_path=output_file_path
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

    def get_first_overlapping_intervals(
            self,
            sorted_intervals: list[P.Interval],
    ) -> Union[tuple[P.Interval, P.Interval], tuple[None, None]]:
        """
        Given a list of intervals, returns the first consecutive
        interval tuples with the overlapping pairs described in
        case a (get_candidate_cds_coordinates).
        :param sorted_intervals: list of intervals
        :return: list of tuples with pairs of overlapping intervals
        """
        first_overlap_index = 0
        while first_overlap_index < len(sorted_intervals) - 1:
            current_interval = sorted_intervals[first_overlap_index]
            next_interval = sorted_intervals[first_overlap_index + 1]
            if (
                    self.get_overlap_percentage(
                        intv_i=current_interval,
                        intv_j=next_interval
                    ) >= self.cds_overlapping_threshold
                    and
                    self.get_overlap_percentage(
                        intv_i=next_interval,
                        intv_j=current_interval
                    ) >= self.cds_overlapping_threshold
            ):
                return current_interval, next_interval
            first_overlap_index += 1
        return None, None

    def resolve_overlaps_between_coordinates(
            self,
            sorted_cds_coordinates: list[P.Interval],
    ) -> list[P.Interval]:
        """
        resolve_overlaps_between_coordinates is a recursive function that given
        a list of coordinates, resolves overlaps according to the following criteria:

        * a: if they have the distinct lengths, and  they overlap by more than
          the overlapping threshold of both CDS lengths the shortest CDS is selected.
        * b: if they have the same length, and they overlap by more than the
          overlapping threshold of both CDS lengths the first CDS is selected.
        * c: both CDSs are selected otherwise

        Since the pairs described in case c are not considered to be overlapping,
        they are both included in the search.
        :param sorted_cds_coordinates: list of unique and sorted coordinates.
        :return: list of coordinates without "overlaps"
        """
        new_list = list()
        # We only process the first pair of overlapping intervals since
        # the resolved overlap could also overlap with the next interval in the list.
        intv_i, intv_j = self.get_first_overlapping_intervals(
            sorted_intervals=sorted_cds_coordinates
        )
        if all((intv_i, intv_j)):
            candidate, _ = self.get_candidate_cds(intv_i=intv_i, intv_j=intv_j)
            # Note: new_list is populated through enumeration
            # of 'sorted_cds_coordinates', so it will also be sorted
            # according to the same criteria used for 'sorted_cds_coordinates'.
            for idx, cds_coordinate in enumerate(sorted_cds_coordinates):
                if cds_coordinate != intv_i:
                    new_list.append(cds_coordinate)
                else:
                    new_list.append(candidate)
                    new_list.extend(sorted_cds_coordinates[idx + 2:])
                    return self.resolve_overlaps_between_coordinates(
                        sorted_cds_coordinates=new_list
                    )
        return sorted_cds_coordinates

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

        * a: if they overlap by more than the overlapping threshold
          of both CDS lengths the shortest CDS is selected.
        * b: both CDSs are selected otherwise

        The rationale behind choosing the shortest CDS is that it will
        reduce exclusion of short duplications due to coverage thresholds.

        :param gene_id: gene identifier
        :return: list of representative CDS coordinates across all transcripts
        """
        # collect all CDS coordinates and frames across all transcripts
        # we are interested in the frames to account for the unlikely event
        # that two CDS with same coordinates in different transcripts
        # have different frames
        cds_coordinates_and_frames: list[tuple[P.Interval, str]] = list(
            set(
                (coordinate, annotation_structure['frame'])
                for mrna_annotation in self.data_container.gene_hierarchy_dictionary[gene_id]['mRNAs'].values()
                for annotation_structure in mrna_annotation['structure']
                for coordinate in (annotation_structure['coordinate'],)
                if (
                        annotation_structure['type'] == 'CDS' and
                        (coordinate.upper - coordinate.lower) >= self.min_exon_length
                )
            )
        )
        if cds_coordinates_and_frames:
            representative_cds_frame_dictionary = dict()
            for cds_coordinate, frame in cds_coordinates_and_frames:
                if cds_coordinate in representative_cds_frame_dictionary:
                    representative_cds_frame_dictionary[cds_coordinate]['frame'] += f'_{str(frame)}'
                else:
                    representative_cds_frame_dictionary[cds_coordinate] = {'frame': frame}
            sorted_cds_coordinates_list: list[P.Interval] = sorted(
                [coordinate for coordinate, _ in cds_coordinates_and_frames],
                key=lambda coordinate: (coordinate.lower, coordinate.upper),
            )
            return dict(
                set_coords=self.resolve_overlaps_between_coordinates(
                    sorted_cds_coordinates=sorted_cds_coordinates_list
                ),
                cds_frame_dict=representative_cds_frame_dictionary,
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
                cds_frame=cds_frame
            )
        else:
            tblastx_o = self.execute_tblastx_using_tempfiles(
                hit_sequence=hit_sequence,
                query_sequence=query_sequence,
                query_coordinate=query_coordinate,
                gene_coordinate=gene_coordinate,
                cds_frame=cds_frame
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
            gene_id: str,
            has_duplication_binary: int,
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
            gene_coordinate.lower,
            gene_coordinate.upper,
            has_duplication_binary
        )

    def find_coding_exon_duplicates(
            self,
            gene_id: str,
    ) -> None:
        """
        find_coding_exon_duplicates is a function that given a gene_id,
        performs a tblastx for each representative CDS (see get_candidate_cds_coordinates).
        If the tblastx search returns hits, they are stored in the "results" database,
        otherwise the gene is recorded as having no duplication event.
        :param gene_id: gene identifier
        """
        blast_hits_dictionary = dict()
        # time.sleep(random.randrange(start=0, stop=self.sleep_max_seconds))
        chromosome, gene_coordinate, gene_strand = (
            self.data_container.gene_hierarchy_dictionary[gene_id]['chrom'],
            self.data_container.gene_hierarchy_dictionary[gene_id]['coordinate'],
            self.data_container.gene_hierarchy_dictionary[gene_id]['strand']
        )
        temp_gene_sequence = str(
            Seq(self.data_container.genome_dictionary[chromosome][gene_coordinate.lower:gene_coordinate.upper])
        )
        gene_dna_sequence = self.check_for_masking(
            chromosome=chromosome,
            gene_id=gene_id,
            sequence=temp_gene_sequence,
            coordinate=gene_coordinate,
            annotation_type='gene'
        )
        if gene_dna_sequence:
            cds_coordinates_dictionary = self.get_candidate_cds_coordinates(gene_id=gene_id)
            if cds_coordinates_dictionary:
                for cds_coordinate in cds_coordinates_dictionary['set_coords']:
                    temp_dna_cds_sequence = str(
                        Seq(self.data_container.genome_dictionary[chromosome][cds_coordinate.lower:cds_coordinate.upper])
                    )
                    cds_dna_sequence = self.check_for_masking(
                        chromosome=chromosome,
                        gene_id=gene_id,
                        sequence=temp_dna_cds_sequence,
                        coordinate=cds_coordinate,
                        annotation_type='CDS'
                    )
                    if cds_dna_sequence:
                        cds_frame = cds_coordinates_dictionary['cds_frame_dict'][cds_coordinate]['frame']
                        tblastx_o = self.align_cds(
                            gene_id=gene_id,
                            query_sequence=cds_dna_sequence,
                            hit_sequence=gene_dna_sequence,
                            query_coordinate=cds_coordinate,
                            cds_frame=cds_frame
                        )
                        if tblastx_o:
                            blast_hits_dictionary[cds_coordinate] = tblastx_o
                if blast_hits_dictionary:
                    self.populate_fragments_table(
                        gene_id=gene_id,
                        blast_results_dictionary=blast_hits_dictionary
                    )
                else:
                    self.database_interface.insert_gene_ids_table(
                        gene_args_tuple=self.get_gene_tuple(
                            gene_id=gene_id,
                            has_duplication_binary=0
                        )
                    )

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
        self.database_interface.insert_fragments(
            gene_args_tuple=self.get_gene_tuple(gene_id=gene_id, has_duplication_binary=1),
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
    ) -> tuple[float, float, str, str, int]:
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

        return (self.compute_identity(sequence_i=query_dna_seq, sequence_j=target_dna_seq),
                self.compute_identity(sequence_i=query_aln_prot_seq, sequence_j=target_aln_prot_seq),
                query_dna_seq,
                target_dna_seq,
                fragment_id)

    def get_identity_and_dna_seq_tuples(
            self,
    ) -> list[tuple]:
        """
        Retrieves DNA sequences from tblastx query and target,
         computes their DNA and amino acid identity,
        and returns a list of tuples with the structure:
        (DNA_identity, AA_identity, query_dna_seq, target_dna_seq, fragment_id)
        """
        return [self.process_fragment(fragment=fragment)
                for fragment in self.database_interface.query_fragments()]
