from .exonize_handler import *
import argparse


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
    parser = argparse.ArgumentParser(description='Exonize Description')
    parser.add_argument('gff_file_path', type=str, help='Path to GFF file')
    parser.add_argument('genome_path', type=str, help='Path to genome file')
    parser.add_argument('specie_identifier', type=str, help='Species identifier')
    parser.add_argument('-o', '--results_db_name', default='', type=str, help='Results database name')
    parser.add_argument('-db', '--enable_debug', default=False, help='DEBUG MODE - Saves input and output tblastx files')
    parser.add_argument('-sf', '--soft_force', default=False, help='soft force - if True, the results database will be overwritten if it already exists')
    parser.add_argument('-hf', '--hard_force', default=False, help='soft force - if True, all internal files will be overwritten if they already exist')
    parser.add_argument('-hm', '--hard_masking', default=False, action='store_true', help='Hard masking flag')
    parser.add_argument('-p', '--sleep_max_seconds', default=5, type=int, help='Max sleep seconds')
    parser.add_argument('-el', '--min_exon_length', type=int, default=30, help='Minimum exon length')
    parser.add_argument('-et', '--evalue_threshold', default=1e-2, type=float, help='E-value threshold')
    parser.add_argument('-ht', '--self_hit_threshold', default=0.5, type=float, help='Self-hit threshold')
    parser.add_argument('-ot', '--cds_overlapping_threshold', default=0.9, type=float, help='CDS overlapping threshold')
    parser.add_argument('-mt', '--masking_perc_threshold', default=0.8, type=float, help='Masking percentage threshold')
    parser.add_argument('-bn', '--batch_number', default=100, type=int, help='Batch number')
    parser.add_argument('-t', '--threads', default=7, type=int, help='Number of threads')
    parser.add_argument('-to', '--timeout_db', default=160, type=int, help='Database timeout')
    args = parser.parse_args()
    return args


def main():
    exonize_ascii_art_logo()
    args = argument_parser()
    exonize_obj = Exonize(gff_file_path=args.gff_file_path,
                          genome_path=args.genome_path,
                          specie_identifier=args.specie_identifier,
                          results_db_name=args.results_db_name,
                          enable_debug=args.enable_debug,
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
                          timeout_db=args.timeout_db)
    exonize_obj.run_exonize_pipeline()


if __name__ == '__main__':
    main()
