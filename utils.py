import pickle  # saving and loading dictionaries
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO


def read_pkl_file(filepath: str) -> dict:
    with open(filepath, 'rb') as handle:
        read_file = pickle.load(handle)
    return read_file

# write a function display_time that given a time of start "tic" and end "tac" prints out the time in seconds
# if (tac - tic) is less than a minute, minutes if (tac - tic) is less than an hour


def dump_pkl_file(out_filepath: str, obj: dict) -> None:
    with open(out_filepath, 'wb') as handle:
        pickle.dump(obj, handle)


def dump_fasta_file(out_filepath: str, seq_dict: dict) -> None:
    with open(out_filepath, "w") as handle:
        for annot_id, annot_seq in seq_dict.items():
            record = SeqRecord(Seq(annot_seq), id=str(annot_id), description='')
            SeqIO.write(record, handle, "fasta")


def sort_intervals_dict(interval_dict: dict) -> dict:
    sorted_intervals = sorted(list(interval_dict.keys()), key=lambda item: (item.lower, item.upper))
    return {coord: interval_dict[coord] for coord in sorted_intervals}


def hamming_distance(seq_a: str, seq_b: str) -> float:
    return sum([i == j for i, j in zip(seq_a, seq_b)]) / len(seq_b)


def reverse_complement(seq):
    return ''.join([{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}[B] for B in seq][::-1])


def get_dna_seq(seq, frame, start, end):
    if frame < 0:
        rev_comp = reverse_complement(seq)
        return rev_comp[start:end]
    else:
        return seq[start:end]


def seq_record(seq, seq_id):
    return SeqRecord(Seq(seq), id=seq_id)


def batch(iterable, n=1):
    it_length = len(iterable)
    for ndx in range(0, it_length, n):
        yield iterable[ndx:min(ndx + n, it_length)]


def codon_alignment(dna_seq_a, dna_seq_b, peptide_seq_a, peptide_seq_b):
    dna_align_a = gaps_from_peptide(peptide_seq_a, dna_seq_a)
    dna_align_b = gaps_from_peptide(peptide_seq_b, dna_seq_b)
    return dna_align_a, dna_align_b


def get_small_large_interval(a, b):
    """
    Given two intervals, the function
    get_small_large_interval returns the smaller
    and the larger interval in length.
    """
    len_a = a.upper - a.lower
    len_b = b.upper - b.lower
    if len_a < len_b:
        return a, b
    else:
        return b, a


def get_overlapping_percentage(a, b):
    """
    Given two intervals, the function
    get_overlapping_percentage returns the percentage
    of the overlapping region relative to an interval b.
    """
    intersection = a & b
    if intersection:
        return (intersection.upper - intersection.lower) / (b.upper - b.lower)
    else:
        return 0


def strand_to_int(strand):
    if strand == '+':
        return 1
    else:
        return -1


def check_int_type(fragment_dict):
    return all(isinstance(x, int) for x in [fragment_dict['query_coord'].lower,
                                            fragment_dict['query_coord'].upper,
                                            fragment_dict['hit_coord'].lower,
                                            fragment_dict['hit_coord'].upper])


def gaps_from_peptide(peptide_seq, nucleotide_seq):
    """ Transfers gaps from aligned peptide seq into codon partitioned nucleotide seq (codon alignment)
          - peptide_seq is an aligned peptide sequence with gaps that need to be transferred to nucleotide seq
          - nucleotide_seq is an un-aligned dna sequence whose codons translate to peptide seq"""
    def chunks(l, n):
        """ Yield successive n-sized chunks from l."""
        for i in range(0, len(l), n):
            yield l[i:i+n]
    codons = [codon for codon in chunks(nucleotide_seq, 3)]  # splits nucleotides into codons (triplets)
    gapped_codons = []
    codon_count = 0
    for aa in peptide_seq:  # adds '---' gaps to nucleotide seq corresponding to peptide
        if aa != '*':
            gapped_codons.append(codons[codon_count])
            codon_count += 1
        else:
            gapped_codons.append('---')
    return ''.join(gapped_codons)

