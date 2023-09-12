from exonize_handler import Exonize
import argparse


def check_requirements():
    req_file_name = "requirements.txt"
    if os.path.exists(req_file_name):
        installed_packages = subprocess.getoutput("pip list").split("\n")[2:]  # The first two lines are headers
        installed = [package.split()[0] for package in installed_packages]
        with open(req_file_name, "r") as f:
            required_packages = f.read().splitlines()
        missing_packages = [package.split("==")[0] for package in required_packages
                            if package.split("==")[0] not in installed]
        if missing_packages:
            print("The following packages are missing:")
            for package in missing_packages:
                print(package)
            print("Please install them before running Exonize.")
            sys.exit(1)


def argument_parser():
    parser = argparse.ArgumentParser(description='Exonize Description')
    parser.add_argument('gff_file_path', type=str, help='Path to GFF file')
    parser.add_argument('genome_path', type=str, help='Path to genome file')
    parser.add_argument('specie_identifier', type=str, help='Species identifier')
    parser.add_argument('--results_db_name', default='', type=str, help='Results database name')
    parser.add_argument('--spec_attribute', default='ID', type=str, help='Feature reference for writting intron annotations with gffutils')
    parser.add_argument('--save_input_files', action='store_true', help='Saves input and output tblastx files')
    parser.add_argument('--verbose', default=False, action='store_true', help='Verbose mode')
    parser.add_argument('--hard_masking', default=False, action='store_true', help='Hard masking flag')
    parser.add_argument('--evalue_threshold', default=1e-2, type=float, help='E-value threshold')
    parser.add_argument('--sleep_max_seconds', default=5, type=int, help='Max sleep seconds')
    parser.add_argument('--min_exon_length', type=int, default=30, help='Minimum exon length')
    parser.add_argument('--self_hit_threshold', default=0.5, type=float, help='Self-hit threshold')
    parser.add_argument('--cds_overlapping_threshold', default=0.9, type=float, help='CDS overlapping threshold')
    parser.add_argument('--masking_perc_threshold', default=0.8, type=float, help='Masking percentage threshold')
    parser.add_argument('--batch_number', default=100, type=int, help='Batch number')
    parser.add_argument('--threads', default=7, type=int, help='Number of threads')
    parser.add_argument('--timeout_db', default=160, type=int, help='Database timeout')
    args = parser.parse_args()
    return args


def main():
    check_requirements()
    args = argument_parser()
    Exonize(gff_file_path=args.gff_file_path,
            genome_path=args.genome_path,
            specie_identifier=args.specie_identifier,
            results_db_name=args.results_db_name,
            spec_attribute=args.spec_attribute,
            save_input_files=args.save_input_files,
            verbose=args.verbose,
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


if __name__ == '__main__':
    main()
