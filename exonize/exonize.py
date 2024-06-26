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
          "             Liam M. Longo, Earth-Life Science Institute (ELSI), Tokyo Institute of Technology\n"
          "             Christopher Wheat, Department of Zoology, Stockholm University\n"
          "             Lars Arvestad, Department of Mathematics, Stockholm University\n"

          "maintainers: Marina Herrera Sarrias, Department of Mathematics, Stockholm University,\n"
          "             Lars Arvestad, Department of Mathematics, Stockholm University\n"
          "    Contact: arvestad@math.su.se\n"
          "     GitHub: https://github.com/msarrias/exonize\n"
          "\n",
          file=sys.stderr)


def argument_parser():
    parser = argparse.ArgumentParser(
        description='Exonize: A tool for discovering exon duplications.'
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
    # Optional Arguments for Flags
    parser.add_argument(
        '--output_prefix',
        type=str,
        help='Species identifier - used for naming output files.')
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
        '--multigraphs',
        action='store_true',
        default=False,
        help='Generate expansion graphs figures. Default is False.'
    )
    parser.add_argument(
        '--csv',
        # metavar='csv-output-prefix',
        action='store_true',
        default=False,
        help='File prefix for storing CSV file with identified exon duplicates and their classification.'
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
        '-el',
        '--min-exon-length',
        default=30,
        type=int,
        help='Minimum exon length. Default is 30.'
    )
    parser.add_argument(
        '-et',
        '--evalue-threshold',
        default=1e-2,
        type=float,
        help='E-value threshold. Default is 1e-2.'
    )
    parser.add_argument(
        '-ht',
        '--self-hit-threshold',
        default=0.5,
        type=float,
        help='Self-hit threshold. Default is 0.5.'
    )
    parser.add_argument(
        '-ot',
        '--cds-overlapping-threshold',
        default=0.9,
        type=float,
        help='CDS overlapping threshold. Default is 0.9.'
    )
    parser.add_argument(
        '-qt',
        '--query-overlapping-threshold',
        default=0.9,
        type=float,
        help='tblastx query overlapping threshold. Default is 0.9.'
    )
    parser.add_argument(
        '--cpus_number',
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
        output_prefix=args.output_prefix,
        draw_event_multigraphs=args.multigraphs,
        csv=args.csv,
        enable_debug=args.debug,
        soft_force=args.soft_force,
        hard_force=args.hard_force,
        evalue_threshold=args.evalue_threshold,
        sleep_max_seconds=args.sleep_max_seconds,
        min_exon_length=args.min_exon_length,
        cds_overlapping_threshold=args.cds_overlapping_threshold,
        query_overlapping_threshold=args.query_overlapping_threshold,
        self_hit_threshold=args.self_hit_threshold,
        cpus_number=args.cpus_number,
        timeout_database=args.timeout_database,
        output_directory_path=args.output_directory_path
    )
    exonize_obj.run_exonize_pipeline()


if __name__ == '__main__':
    main()
