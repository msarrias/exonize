import argparse
import cProfile

from exonize.exonize_handler import Exonize
from exonize.profiling import get_run_performance_profile, PROFILE_PATH


def exonize_ascii_art_logo() -> None:
    exonize_ansi_regular = """
    
    ███████╗██╗  ██╗ ██████╗ ███╗   ██╗██╗███████╗███████╗
    ██╔════╝╚██╗██╔╝██╔═══██╗████╗  ██║██║╚══███╔╝██╔════╝
    █████╗   ╚███╔╝ ██║   ██║██╔██╗ ██║██║  ███╔╝ █████╗
    ██╔══╝   ██╔██╗ ██║   ██║██║╚██╗██║██║ ███╔╝  ██╔══╝
    ███████╗██╔╝ ██╗╚██████╔╝██║ ╚████║██║███████╗███████╗
    ╚══════╝╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚═╝╚══════╝╚══════╝
        """
    print(exonize_ansi_regular)
    print("Exonize v1.0\n"
          " Developed by: Marina Herrera Sarrias, Mathematics Department, Stockholm University\n"
          "Supervised by: Lars Arvestad, Mathematics Department, Stockholm University\n"
          "               & Liam Longo, Earth-Life Science Institute (ELSI), Tokyo Institute of Technology\n"
          "      Contact: marina.sarrias@math.su.se\n"
          "       GitHub: https://github.com/msarrias/exonize\n"
          "\n")


def argument_parser():
    parser = argparse.ArgumentParser(description='Exonize: A tool for discovering exon duplications.')
    # Required Arguments
    parser.add_argument('gff_file_path', type=str, help='Path to GFF file.')
    parser.add_argument('genome_path', type=str, help='Path to genome file.')
    parser.add_argument('specie_identifier', type=str, help='Species identifier.')
    # Optional Arguments for Flags
    parser.add_argument('--debug', dest='debug_mode', action='store_true', default=False,
                        help='Enable DEBUG mode, which saves input and output tblastx files.')
    parser.add_argument('--soft-force', dest='soft_force', action='store_true', default=False,
                        help='If set, the results database will be overwritten if it already exists.')
    parser.add_argument('--hard-force', dest='hard_force', action='store_true', default=False,
                        help='If set, all internal files will be overwritten if they already exist.')
    parser.add_argument('--hard-masking', dest='hard_masking', action='store_true', default=False,
                        help='Enable hard masking.')
    # Optional Arguments for Numerical Values and Thresholds
    parser.add_argument('-p', '--sleep-max-seconds', default=5, type=int,
                        help='Max seconds to sleep. Default is 5.')
    parser.add_argument('-el', '--min-exon-length', type=int, default=30,
                        help='Minimum exon length. Default is 30.')
    parser.add_argument('-et', '--evalue-threshold', default=1e-2, type=float,
                        help='E-value threshold. Default is 1e-2.')
    parser.add_argument('-ht', '--self_hit_threshold', default=0.5, type=float,
                        help='Self-hit threshold. Default is 0.5.')
    parser.add_argument('-ot', '--cds_overlapping_threshold', default=0.9, type=float,
                        help='CDS overlapping threshold. Default is 0.9.')
    parser.add_argument('-mt', '--masking_perc_threshold', default=0.8, type=float,
                        help='Masking percentage threshold. Default is 0.8.')
    # Optional Arguments for Parallel and Batch Processing
    parser.add_argument('-bn', '--batch-number', default=100, type=int,
                        help='Number of batches. Default is 100.')
    parser.add_argument('-t', '--threads', default=7, type=int,
                        help='Number of threads to use. Default is 7.')
    # Optional Argument for Timeout
    parser.add_argument('-to', '--timeout-db', default=160, type=int,
                        help='Database timeout. Default is 160.')
    # Optional Argument for saving the parsed genome as a pickle file
    parser.add_argument('--genome_pickled_file_path', default='.parsed_genome.pkl', type=str,
                        help='Parsed genome pickled file path. Default is .parsed_genome.pkl.')
    args = parser.parse_args()
    return args


def main():
    exonize_ascii_art_logo()
    args = argument_parser()
    exonize_obj = Exonize(gff_file_path=args.gff_file_path,
                          genome_path=args.genome_path,
                          specie_identifier=args.specie_identifier,
                          enable_debug=args.debug_mode,
                          soft_force=args.soft_force,
                          hard_force=args.hard_force,
                          hard_masking=args.hard_masking,
                          evalue_threshold=args.evalue_threshold,
                          sleep_max_seconds=args.sleep_max_seconds,
                          min_exon_length=args.min_exon_length,
                          cds_overlapping_threshold=args.cds_overlapping_threshold,
                          masking_perc_threshold=args.masking_perc_threshold,
                          self_hit_threshold=args.self_hit_threshold,
                          batch_number=args.batch_number,
                          threads=args.threads,
                          timeout_db=args.timeout_db,
                          genome_pickled_file_path=args.genome_pickled_file_path,
                          )
    exonize_obj.run_exonize_pipeline()


if __name__ == '__main__':
    main()
