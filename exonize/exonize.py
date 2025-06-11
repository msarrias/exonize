import argparse
from pathlib import Path
import os
import sys
from exonize.exonize_handler import Exonize
from exonize import __version__


def exonize_ascii_art_logo() -> None:
    exonize_ansi_regular = """
    ███████╗██╗  ██╗ ██████╗ ███╗   ██╗██╗███████╗███████╗
    ██╔════╝╚██╗██╔╝██╔═══██╗████╗  ██║██║╚══███╔╝██╔════╝
    █████╗   ╚███╔╝ ██║   ██║██╔██╗ ██║██║  ███╔╝ █████╗
    ██╔══╝   ██╔██╗ ██║   ██║██║╚██╗██║██║ ███╔╝  ██╔══╝
    ███████╗██╔╝ ██╗╚██████╔╝██║ ╚████║██║███████╗███████╗
    ╚══════╝╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚═╝╚══════╝╚══════╝
        """
    print(exonize_ansi_regular, file=sys.stderr)
    print(f"Exonize {__version__}\n"
          "    authors: Marina Herrera Sarrias, Department of Mathematics, Stockholm University,\n"
          "             Christopher Wheat, Department of Zoology, Stockholm University\n"
          "             Liam M. Longo, Earth-Life Science Institute (ELSI), Institute of Science Tokyo\n"
          "             Lars Arvestad, Department of Mathematics, Stockholm University\n"

          "maintainers: Marina Herrera Sarrias, Department of Mathematics, Stockholm University,\n"
          "             Lars Arvestad, Department of Mathematics, Stockholm University\n"
          "    Contact: arvestad@math.su.se\n"
          "     GitHub: https://github.com/msarrias/exonize\n"
          "\n",
          file=sys.stderr)


def argument_parser():
    parser = argparse.ArgumentParser(
        description='exonize: A tool for discovering exon duplications.'
    )
    # Required Arguments
    parser.add_argument(
        'gff_file_path',
        type=Path,
        help='Path to GFF file.'
    )
    parser.add_argument(
        'genome_file_path',
        type=Path,
        help='Path to genome file.'
    )
    # Optional Arguments for GFF annotations
    parser.add_argument(
        '--pdb_ids_mapping_db_path',
        default=None,
        type=Path,
        help='Path to the PDB IDs mapping database.'
    )
    parser.add_argument(
        '--pdb_structures_path',
        default=None,
        type=Path,
        help='Path to the PDB structures.'
    )
    parser.add_argument(
        '-gfeat',
        '--gene-annot-feature',
        default='gene',
        help='Gene feature in annotation. Default is gene.'
    )
    parser.add_argument(
        '-cdsfeat',
        '--cds-annot-feature',
        default='CDS',
        help='CDS feature in annotation. Default is CDS.'
    )
    parser.add_argument(
        '-transfeat',
        '--transcript-annot-feature',
        default='transcript',
        help='Transcript feature in annotation. Default is transcript.'
    )
    parser.add_argument(
        '-sb',
        '--sequence-base',
        type=int,
        default=1,
        help='Annotation sequence numbering base (0/1). Default is 1.'
    )
    parser.add_argument(
        '-fb',
        '--frame-base',
        type=int,
        default=0,
        help='Frame base (0/1). Default is 0.'
    )
    parser.add_argument(
        '-l',
        '--min-exon-length',
        default=30,
        type=int,
        help='Minimum exon length. Default is 30.'
    )
    parser.add_argument(
        '-e',
        '--evalue-threshold',
        default=1e-3,
        type=float,
        help='E-value threshold. Default is 1e-3.'
    )
    parser.add_argument(
        '-ts',
        '--self-hit-threshold',
        default=0.5,
        type=float,
        help='Self-hit threshold. Default is 0.5.'
    )
    parser.add_argument(
        '-te',
        '--query-coverage-threshold',
        default=0.8,
        type=float,
        help='query coverage threshold. Default is 0.8.'
    )
    parser.add_argument(
        '-ce',
        '--exon-clustering-overlap-threshold',
        default=0.9,
        type=float,
        help='Exon clustering overlap threshold. Default is 0.9.'
    )
    parser.add_argument(
        '-ct',
        '--targets-clustering-overlap-threshold',
        default=0.9,
        type=float,
        help='Target coordinates clustering overlap threshold. Default is 0.9.'
    )
    parser.add_argument(
        '-tp',
        '--pair-coverage-threshold',
        default=0.9,
        type=float,
        help='Minimum length coverage between pair of coordinates. Default is 0.9.'
    )
    parser.add_argument(
        '-ta',
        '--fraction-of-aligned-positions',
        default=0.9,
        type=float,
        help='Fraction of aligned positions threshold. Default is 0.9.'
    )
    parser.add_argument(
        '-ti',
        '--peptide-identity-threshold',
        default=0.4,
        type=float,
        help='Local search peptide identity threshold. Default is 0.4.'
    )
    parser.add_argument(
        '-le',
        '--min-exon-length-structural-search',
        default=75,
        type=int,
        help='Minimum exon length for structural search. Default is 75.'
    )
    parser.add_argument(
        '-plddt_th',
        '--plddt-threshold',
        default=70,
        type=float,
        help='Structural search pLDDT (per-residue measure of local confidence) threshold. Default is 0.7.'
    )
    parser.add_argument(
        '-RMSD_th',
        '--RMSD-threshold',
        default=1,
        type=float,
        help='Structural search RMSD (Root Mean Square Deviation) threshold. Default is 0.2.'
    )
    parser.add_argument(
        '-tms_th',
        '--TM-score-threshold',
        default=0.5,
        type=float,
        help='Normalized TM-score threshold. Default is 0.5.'
    )
    # Optional Arguments for Flags
    parser.add_argument(
        '-op',
        '--output_prefix',
        type=str,
        help='Species identifier - used for naming output files.')
    parser.add_argument(
        '--global-search',
        action='store_true',
        default=False,
        help='Exonize will perform a global search'
    )
    parser.add_argument(
        '--local-search',
        action='store_true',
        default=False,
        help='Exonize will perform a local search'
    )
    parser.add_argument(
        '--structural-search',
        action='store_true',
        default=False,
        help='Exonize will perform a structural search'
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        default=False,
        help='Enable DEBUG mode, which saves input and output tblastx files.'
    )
    parser.add_argument(
        '--soft-force',
        action='store_true',
        default=False,
        help='If set, the results database will be overwritten if it already exists.'
    )
    parser.add_argument(
        '--hard-force',
        action='store_true',
        default=False,
        help='If set, all internal files will be overwritten if they already exist.'
    )
    parser.add_argument(
        '--csv',
        action='store_true',
        default=False,
        help='If set, exonize will output a .zip file with a reduced set of the results in CSV format.'
    )
    # Optional Arguments for Numerical Values and Thresholds
    parser.add_argument(
        '-p',
        '--sleep-max-seconds',
        default=5,
        type=int,
        help='Max seconds to sleep. Default is 5.'
    )
    parser.add_argument(
        '-cn',
        '--cpus-number',
        default=os.cpu_count(),  # This is pretty greedy, could be changed and put in a config file
        type=int,
        help='Number of CPUs to use. Default is the number of CPUs available.')
    parser.add_argument(
        '-to',
        '--timeout-database',
        default=160,
        type=int,
        help='Database timeout. Default is 160.'
    )
    parser.add_argument(
        '-odp',
        '--output-directory-path',
        type=Path,
        help='Output directory path. Default is current directory.'
    )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s '+__version__
    )
    args = parser.parse_args()
    return args


def main():
    args = argument_parser()
    exonize_ascii_art_logo()
    exonize_obj = Exonize(
        gff_file_path=args.gff_file_path,
        genome_file_path=args.genome_file_path,
        pdb_ids_mapping_db_path=args.pdb_ids_mapping_db_path,
        pdb_structures_path=args.pdb_structures_path,
        gene_annot_feature=args.gene_annot_feature,
        cds_annot_feature=args.cds_annot_feature,
        transcript_annot_feature=args.transcript_annot_feature,
        sequence_base=args.sequence_base,
        frame_base=args.frame_base,
        min_exon_length=args.min_exon_length,
        evalue_threshold=args.evalue_threshold,
        self_hit_threshold=args.self_hit_threshold,
        query_coverage_threshold=args.query_coverage_threshold,
        exon_clustering_overlap_threshold=args.exon_clustering_overlap_threshold,
        targets_clustering_overlap_threshold=args.targets_clustering_overlap_threshold,
        pair_coverage_threshold=args.pair_coverage_threshold,
        fraction_of_aligned_positions=args.fraction_of_aligned_positions,
        peptide_identity_threshold=args.peptide_identity_threshold,
        min_exon_length_structural_search=args.min_exon_length_structural_search,
        plddt_threshold=args.plddt_threshold,
        RMSD_threshold=args.RMSD_threshold,
        TM_score_threshold=args.TM_score_threshold,
        output_prefix=args.output_prefix,
        csv=args.csv,
        enable_debug=args.debug,
        soft_force=args.soft_force,
        hard_force=args.hard_force,
        global_search=args.global_search,
        structural_search=args.structural_search,
        local_search=args.local_search,
        sleep_max_seconds=args.sleep_max_seconds,
        cpus_number=args.cpus_number,
        timeout_database=args.timeout_database,
        output_directory_path=args.output_directory_path
    )
    exonize_obj.run_exonize_pipeline()


if __name__ == '__main__':
    main()
