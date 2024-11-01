```
███████╗██╗  ██╗ ██████╗ ███╗   ██╗██╗███████╗███████╗
██╔════╝╚██╗██╔╝██╔═══██╗████╗  ██║██║╚══███╔╝██╔════╝
█████╗   ╚███╔╝ ██║   ██║██╔██╗ ██║██║  ███╔╝ █████╗
██╔══╝   ██╔██╗ ██║   ██║██║╚██╗██║██║ ███╔╝  ██╔══╝
███████╗██╔╝ ██╗╚██████╔╝██║ ╚████║██║███████╗███████╗
╚══════╝╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚═╝╚══════╝╚══════╝
Marina Herrera Sarrias, Department of Mathematics, Stockholm University
Liam M. Longo, Earth-Life Science Institute (ELSI), Institute of Science Tokyo
Christopher Wheat, Department of Zoology, Stockholm University
Lars Arvestad, Department of Mathematics, Stockholm University
```

## Welcome!
`exonize` is an open-source command-line tool and Python package for identifying and classifying coding exon duplications in annotated genomes. `exonize` identifies full exon duplications using local and global alignment methods and implements a graph-based framework to handle clusters of exons formed by repetitive duplication events. In addition, `exonize` categorizes the interdependence between duplicated exons (or groups of exons) across transcripts. For data parsing and downstream analysis, the `exonize_analysis` module is available for Python notebooks.

<div align="center">
	<img src="https://github.com/msarrias/exonize/blob/main/figures/pipeline.png" style="width:70%;">
</div>
</div>

## Citation
If you use `exonize` in a publication, please cite:
```
```

## Dependencies
For running `exonize`, a local installation of [`BLAST+`](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) and [`MUSCLE`](https://www.drive5.com/muscle/) are required.

## Getting started
You are best off installing `exonize` from [PyPI.org](https://pypi.org/project/Exonize/1.0/) using
```
pip install exonize
```
If installing from the repo:
```
$ git clone git@github.com:msarrias/exonize.git
$ cd exonize
$ pip install .
```
You should now be able to run `exonize -h`.

For use in Python notebooks:
```python
import exonize_analysis
```
## Positional arguments
Exonize requires two positional arguments:

- **annotation path**: Path to the genome annotations file. This file should be in GFF3 or GFF format.
- **genome path**: Path to the genome sequence file. The file should be in FASTA format, with a `.zip` version also accepted.
  
## Optional Arguments
- **[-gfeat GENE_ANNOT_FEATURE]**: Specifies the gene feature in the genome annotations. Default is gene.
- **[-cdsfeat CDS_ANNOT_FEATURE]**: Specifies the coding sequence feature in the genome annotations. Default is CDS.
- **[-transfeat TRANSCRIPT_ANNOT_FEATURE]**: Specifies the transcript feature in the genome annotations. Default is transcript.
- **[-sb SEQUENCE_BASE]**: Annotation coordinates base (0 or 1). Default is 1.
- **[-fb FRAME_BASE]**: Frame base (0 or 1). Default is 0.
- **[-el MIN_EXON_LENGTH]**: Minimum length required for the search. Default is 30.
- **[-et EVALUE_THRESHOLD]**: E-value threshold for search sensitivity in the local search. Default is 1e-3.
- **[-ht SELF_HIT_THRESHOLD]**: Self-hit threshold. Default is 0.5.
- **[-qt QUERY_COVERAGE_THRESHOLD]**: Minimum percentage of the aligned query for a match to be considered. Default is 0.9.
- **[-ect EXON_CLUSTERING_OVERLAP_THRESHOLD]**: Overlapping threshold for constructing the set of representative exons. Default is 0.9.
- **[-tct TARGETS_CLUSTERING_OVERLAP_THRESHOLD]**: Threshold for clustering target coordinates in the reconciliation. Default is0.9.
- **[-fap FRACTION_OF_ALIGNED_POSITIONS]**: Threshold for fraction of aligned positions in local search. Default is 0.9.
- **[-pit PEPTIDE_IDENTITY_THRESHOLD]**: Peptide identity threshold for local search. Default is 0.4.
- **[-op OUTPUT_PREFIX]**: Search identifier. Default is the stem of the annotations file.
- **[-cn CPUS_NUMBER]**: Number of CPUs to use. Default is the available CPU count.
- **[-odp OUTPUT_DIRECTORY_PATH]**: Path for output directory. Default is the current directory.
- **[--global-search]**: Enables exclusively the global search mode; it cannot be combined with --local-search.
- **[--local-search]**: Enables only the local search mode; it cannot be combined with --global-search.
- **[--debug]**: Enables debug mode, saving input and output files for the local search.
- **[--soft-force]**: Overwrites the results database if it already exists.
- **[--hard-force]**: Overwrites all internal files if they already exist.
- **[--csv]**: Outputs a .zip file with a reduced set of results in CSV format.
- **[-v, --version]**: show program's version number and exit.
- **[-h, --help]**: show this help message and exit

**Note:** If neither flag `--global-search` or `--local-search` is specified, both searches will run by default.

## Usage

```
exonize <gff_file_path> <genome_file_path>
```
Optional arguments should be added in the command following the positional arguments.


## Output Description:
- **`<output_prefix>_results.db`**: This is the main output database created by `exonize`.

	| Table                             | Description                                                                                                                                                                                   |
	|-----------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
	| **Genes**                         | Lists all protein-coding genes analyzed for exon duplications.                                                                                                                               |
	| **Global_matches**                | Contains all non-reciprocal matches identified by aligning pairs of representative exons using `MUSCLE`. All matches recorded here meet the peptide identity and aligned position fraction criteria.  |
	| **Local_matches**                 | Contains all matches found by querying each representative exon within its gene using `tblastx`. These are the raw, unfiltered matches.                                                                        |
	| **Local_matches_non_reciprocal**  | Contains all the combined non-reciprocal matches found during the global and (or) local search.                                                                                                  |
	| **Expansions**                    | Contains all the group of duplicates (expansions) found within the genes. Each record represents a node in the expansion graph.                                                                   |
	| **Expansions_full**               | Includes expansions formed by full-mode events only.                                                                                                                                           |
	| **Expansions_transcript_interdependence** | Classifies transcript interdependence for each expansion in the `Expansions_full` table.                                                                                       |
	| **Expansions_full_tandem**        | Each record represents a pair of consecutive full events within an expansion, indicating whether they are in tandem or not.                                                                   |


	### Database structure:

	A detailed description of each table's columns can be found in [].

	 <div align="center">
		<img src="https://github.com/msarrias/exonize/blob/main/figures/database.png" style="width:100%;">
	</div>
	</div>
	
 
- **`csvs.tar.gz`** (optional): compressed directory in .tar.gz format, containing individual files in CSV format for each table, except for the `Local_matches` and `Genes` tables in the results database.
- **`<output_prefix>.log`**: Contains the search settings.


# Tutorial

## Example dataset: _Homo sapiens_ Y chromosome
The following steps demonstrate how to use `exonize` on a test dataset for the human Ensembl annotations.

1. Download the test data.
	* **If you have installed the package from the repo** move into the `test_human_chrom_Y` directory and download the test data:
	```
	cd test_human_chrom_Y
	source fetch_data.sh
	```
 	* **If you installed the package via `pip`** download the [`fetch_data.sh`](https://github.com/msarrias/exonize/blob/main/test_human_chrom_Y/fetch_data.sh) script and execute it.

2. Run `exonize` on the human Y chromosome dataset. 
```
exonize Homo_sapiens.GRCh38.111.chromosome.Y.gff3 \
        Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz \
        --gene-annot-feature gene \
        --cds-annot-feature CDS \
        --transcript-annot-feature mRNA \
        --output_prefix Homo_sapiens_chrom_Y \
        --csv
```
This command performs both global and local searches.


