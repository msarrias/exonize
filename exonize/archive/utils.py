
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
import tempfile
import re
import random


def get_dna_seq(seq: str, frame: int, start: int, end: int) -> str:
    _nucleotide_correspondence_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

    def reverse_complement(seq_: str) -> str:
        return ''.join([_nucleotide_correspondence_dict[nucleotide] for nucleotide in seq_][::-1])

    if frame < 0:
        rev_comp = reverse_complement(seq)
        return rev_comp[start:end]
    else:
        return seq[start:end]


def codon_alignment(dna_seq_a: str, dna_seq_b: str, peptide_seq_a: str, peptide_seq_b: str) -> tuple:
    dna_align_a = gaps_from_peptide(peptide_seq_a, dna_seq_a)
    dna_align_b = gaps_from_peptide(peptide_seq_b, dna_seq_b)
    return dna_align_a, dna_align_b


def gaps_from_peptide(peptide_seq: str, nucleotide_seq: str) -> str:
    """ Transfers gaps from aligned peptide seq into codon partitioned nucleotide seq (codon alignment)
          - peptide_seq is an aligned peptide sequence with gaps that need to be transferred to nucleotide seq
          - nucleotide_seq is an un-aligned dna sequence whose codons translate to peptide seq"""

    def chunks(seq: str, n: int) -> list:
        """ Yield successive n-sized chunks from l."""
        for i in range(0, len(seq), n):
            yield seq[i:i + n]
    codons = [codon for codon in chunks(nucleotide_seq, 3)]  # splits nucleotides into codons (triplets)
    gaped_codons = []
    codon_count = 0
    for aa in peptide_seq:  # adds '---' gaps to nucleotide seq corresponding to peptide
        if aa != '*':
            gaped_codons.append(codons[codon_count])
            codon_count += 1
        else:
            gaped_codons.append('---')
    return ''.join(gaped_codons)


def muscle_pairwise_alignment(seq1: str, seq2: str, muscle_exe="muscle"):
    with tempfile.TemporaryDirectory() as temp_dir_name:
        input_file = f'{temp_dir_name}/input.fa'
        out_file = f'{temp_dir_name}/out.fa'
        with open(input_file, "w") as f:
            f.write(seq1 + "\n")
            f.write(seq2 + "\n")
        # Run MUSCLE for pairwise alignment
        muscle_cline = MuscleCommandline(muscle_exe, input=input_file, out=out_file)
        muscle_cline()
        # Parse the alignment from the output file
        alignment = AlignIO.read(out_file, "fasta")
    return alignment


def exclude_terminal_gaps_from_pairwise_alignment(seq1: str, seq2: str) -> tuple:
    if len(seq1) == len(seq2):
        p = re.compile(r'^-*([^-\s].*?[^-])-*$')
        seq1_match, seq2_match = p.search(seq1), p.search(seq2)
        s1, e1 = seq1_match.start(1), seq1_match.end(1)
        s2, e2 = seq2_match.start(1), seq2_match.end(1)
        if (e2-s2) == (e1-s1):
            return seq1, seq2
        elif (e2-s2) < (e1-s1):
            return seq1[s2:e2], seq2[s2:e2]
        else:
            return seq1[s1:e1], seq2[s1:e1]
    else:
        print('The input sequences are not aligned')

#
# The following two functions are not used and depends on an undefined
# function:get_overlap_percentage
#
# def get_average_overlapping_percentage(intv_a, intv_b) -> float:
#     return sum([get_overlap_percentage(intv_a, intv_b), get_overlap_percentage(intv_b, intv_a)]) / 2
#
#
# def get_overlapping_set_of_coordinates(list_coords: list, overlapping_threshold=0.9) -> dict:
#     overlapping_coords = {}
#     skip_coords = []
#     for idx, i in enumerate(list_coords):
#         if i not in skip_coords:
#             if i not in overlapping_coords:
#                 overlapping_coords[i] = []
#             for j in list_coords[idx+1:]:
#                 if j not in skip_coords:
#                     if all([get_overlap_percentage(i, j) > overlapping_threshold,
#                             get_overlap_percentage(j, i) > overlapping_threshold]):
#                         overlapping_coords[i].append(j)
#                         skip_coords.append(j)
#             skip_coords.append(i)
#     return overlapping_coords


def check_if_non_overlapping(intv_list: list) -> bool:
    for idx, interval in enumerate(intv_list[:-1]):
        if interval.overlaps(intv_list[idx + 1]):
            return False
    return True


def get_non_overlapping_coords_set(overlapping_coords_dict: dict) -> list:
    return [max([interval, *overlp_intv], key=lambda x: x.upper - x.lower)
            for interval, overlp_intv in overlapping_coords_dict.items()]


def sample_color(n=1) -> list:
    return ["#"+''.join([random.choice('0123456789ABCDEF') for _ in range(6)]) for _ in range(n)]


def strand_string_to_integer(strand: str) -> int:
    if strand == '-':
        return -1
    return 1


def check_for_consecutive_intervals(list_intervals: list, target_intv) -> bool:
    first, last = list_intervals[0], list_intervals[-1]
    if first.lower <= target_intv.lower and target_intv.upper <= last.upper:
        if len(list_intervals) <= 1:
            return True
        for i in range(1, len(list_intervals)):
            if abs(list_intervals[i].lower - list_intervals[i - 1].upper) > 2:
                return False
        return True
    else:
        return False


def get_unmatched_events(string1: str, string2: str) -> list:
    a = string1.rsplit('-')
    b = string2.rsplit('-')
    a_copy = list(a)
    for i in a:
        if i in b:
            b.remove(i)
            a_copy.remove(i)
    a_copy.extend(b)
    return a_copy


def get_overlapping_dict(interval_dict: dict) -> dict:
    """
    gets the next overlap to the right
    for each interval in the interval
    dictionary keys `interval_dict`
    """
    overlapping_dict = {interval: [] for interval in interval_dict}
    intervals_list = list(interval_dict.keys())
    for idx, (feat_interv, feature_annot) in enumerate(interval_dict.items()):
        if idx != (len(interval_dict) - 1) and feat_interv.overlaps(intervals_list[idx + 1]):
            overlapping_dict[feat_interv] = intervals_list[idx+1]
    return overlapping_dict
