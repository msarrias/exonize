import pickle  # saving and loading dictionaries
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
import tempfile
import re
import random


def read_pkl_file(filepath: str) -> dict:
    with open(filepath, 'rb') as handle:
        read_file = pickle.load(handle)
    return read_file


def dump_pkl_file(out_filepath: str, obj: dict) -> None:
    with open(out_filepath, 'wb') as handle:
        pickle.dump(obj, handle)


def dump_fasta_file(out_filepath: str, seq_dict: dict) -> None:
    with open(out_filepath, "w") as handle:
        for annot_id, annot_seq in seq_dict.items():
            record = SeqRecord(Seq(annot_seq), id=str(annot_id), description='')
            SeqIO.write(record, handle, "fasta")


def sort_list_intervals_dict(list_dicts: list) -> list:
    return sorted(list_dicts, key=lambda x: (x['coord'].lower, x['coord']))


def sort_key_intervals_dict(interval_dict: dict) -> dict:
    sorted_intervals = sorted(list(interval_dict.keys()), key=lambda item: (item.lower, item.upper))
    return {coord: interval_dict[coord] for coord in sorted_intervals}


def hamming_distance(seq_a: str, seq_b: str) -> float:
    return sum([i == j for i, j in zip(seq_a, seq_b)]) / len(seq_b)


def reverse_complement(seq: str) -> str:
    return ''.join([{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}[nucleotide] for nucleotide in seq][::-1])


def get_dna_seq(seq: str, frame: int, start: int, end: int) -> str:
    if frame < 0:
        rev_comp = reverse_complement(seq)
        return rev_comp[start:end]
    else:
        return seq[start:end]


def batch(iterable: list, n=1) -> list:
    it_length = len(iterable)
    for ndx in range(0, it_length, n):
        yield iterable[ndx:min(ndx + n, it_length)]


def get_overlap_percentage(a, b) -> float:
    """
    Given two intervals, the function
    get_overlap_percentage returns the percentage
    of the overlapping region relative to an interval b.
    """
    intersection = a & b
    if intersection:
        return (intersection.upper - intersection.lower) / (b.upper - b.lower)
    else:
        return 0


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


def get_fragment_tuple(gene_id: str, cds_coord, blast_hits: dict, hsp_idx: int) -> tuple:
    hsp_dict = blast_hits[hsp_idx]
    hit_q_frame, hit_t_frame = hsp_dict['hit_frame']
    hit_q_f, hit_q_s = reformat_frame_strand(hit_q_frame)
    hit_t_f, hit_t_s = reformat_frame_strand(hit_t_frame)
    return (gene_id,
            cds_coord.lower,
            cds_coord.upper,
            hit_q_f,
            hit_q_s,
            hit_t_f,
            hit_t_s,
            hsp_dict['score'],
            hsp_dict['bits'],
            hsp_dict['evalue'],
            hsp_dict['alignment_len'],
            hsp_dict['query_start'],
            hsp_dict['query_end'],
            hsp_dict['target_start'],
            hsp_dict['target_end'],
            hsp_dict['query_dna_seq'],
            hsp_dict['target_dna_seq'],
            hsp_dict['query_aln_prot_seq'],
            hsp_dict['target_aln_prot_seq'],
            hsp_dict['match'],
            hsp_dict['query_num_stop_codons'],
            hsp_dict['target_num_stop_codons'],
            hsp_dict['dna_perc_identity'],
            hsp_dict['prot_perc_identity'])


def get_hsp_dict(hsp, query_seq: str, hit_seq: str) -> dict:
    q_frame, t_frame = hsp.frame
    q_dna_seq = get_dna_seq(query_seq, q_frame, (hsp.query_start - 1), hsp.query_end)
    t_dna_seq = get_dna_seq(hit_seq, t_frame, (hsp.sbjct_start - 1), hsp.sbjct_end)
    prot_perc_identity = round(hamming_distance(hsp.query, hsp.sbjct), 3)
    dna_perc_identity = round(hamming_distance(q_dna_seq, t_dna_seq), 3)
    return dict(
            score=hsp.score,
            bits=hsp.bits,
            evalue=hsp.expect,
            alignment_len=hsp.align_length * 3,
            hit_frame=hsp.frame,
            query_start=hsp.query_start - 1,
            query_end=hsp.query_end,
            target_start=hsp.sbjct_start - 1,
            target_end=hsp.sbjct_end,
            query_dna_seq=q_dna_seq,
            target_dna_seq=t_dna_seq,
            query_aln_prot_seq=hsp.query,
            target_aln_prot_seq=hsp.sbjct,
            query_dna_align=q_dna_seq,
            target_dna_align=t_dna_seq,
            query_num_stop_codons=hsp.query.count('*'),
            target_num_stop_codons=hsp.sbjct.count('*'),
            match=hsp.match,
            prot_perc_identity=prot_perc_identity,
            dna_perc_identity=dna_perc_identity
            )


def reformat_frame_strand(frame: int) -> tuple:
    n_frame = abs(frame) - 1
    n_strand = '+'
    if frame < 0:
        n_strand = '-'
    return n_frame, n_strand


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


def get_overlapping_set_of_coordinates(list_coords: list) -> dict:
    overlapping_coords = {}
    skip_coords = []
    for idx, i in enumerate(list_coords):
        if i not in skip_coords:
            if i not in overlapping_coords:
                overlapping_coords[i] = []
            for j in list_coords[idx+1:]:
                if j not in skip_coords:
                    if all([get_overlap_percentage(i, j) > 0.9,
                            get_overlap_percentage(j, i) > 0.9]):
                        overlapping_coords[i].append(j)
                        skip_coords.append(j)
            skip_coords.append(i)
    return overlapping_coords


def get_non_overlapping_coords_set(overlapping_coords_dict: dict) -> list:
    return [max([intv, *overlp_intv], key=lambda x: x.upper - x.lower)
            for intv, overlp_intv in overlapping_coords_dict.items()]


def sample_color(n=1) -> list:
    return ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
            for i in range(n)]


def strand_string_to_integer(strand: str) -> int:
    if strand == '-':
        return -1
    return 1


def identify_obligatory_events():
    pass


def identify_truncation_events():
    pass


def identify_trunctation():
    pass


def identify_deactivated_duplication():
    pass


def identify_MEXs_transcripts():
    pass



def exonize_asci_art():
    exonize_ansi_regular = """
    
    ███████╗██╗  ██╗ ██████╗ ███╗   ██╗██╗███████╗███████╗
    ██╔════╝╚██╗██╔╝██╔═══██╗████╗  ██║██║╚══███╔╝██╔════╝
    █████╗   ╚███╔╝ ██║   ██║██╔██╗ ██║██║  ███╔╝ █████╗
    ██╔══╝   ██╔██╗ ██║   ██║██║╚██╗██║██║ ███╔╝  ██╔══╝
    ███████╗██╔╝ ██╗╚██████╔╝██║ ╚████║██║███████╗███████╗
    ╚══════╝╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚═╝╚══════╝╚══════╝
    """
    print(exonize_ansi_regular)
