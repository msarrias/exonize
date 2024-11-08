
Usage
============

Exonize requires two positional arguments that can be followed by optional arguments. 

```
exonize <gff_file_path> <genome_file_path> [OPTIONS]
```

Required Arguments
---------------------

- `<gff_file_path>`: Path to the genome annotation file. This file should be in GFF3 or GFF format.
- `<genome_file_path>`: Path to the genome sequence file. The file should be in FASTA format. A `.zip` version is also accepted.
 
Optional Arguments
---------------------
- `[-gfeat]`: Specifies the gene feature in the annotations file. **Default**: `gene`.
- `[-cdsfeat]`: Specifies the coding sequence feature in the annotations file. **Default**: `CDS`.
- `[-transfeat]`: Specifies the transcript feature in the genome annotations. **Default**: `transcript`.
- `[-sb]`: Annotation coordinates base of input annotations, either `0` or `1`. **Default**: `1`.
- `[-fb]`: Frame base of input annotations, either `0` or `1`. **Default**: `0`.
- `[-l]`: Minimum exon length in bases required for the search. **Default**: `30`.
- `[-e]`: E-value threshold (local BLAST search param). **Default**: `1e-3`.
- `[-ts]`: Self-hit overlap threshold. **Default**: `0.5`.
- `[-te]`: Minimum query coverage cut-off. **Default**: `0.9`.
- `[-ce]`: Overlap threshold for constructing the set of representative exons. **Default**: `0.9`.
- `[-ct]`: Overlap threshold for target coordinate reconciliation (local BLAST search param). **Default**: `0.9`.
- `[-tp]`: Minimum length coverage between the query and target (global MUSCL search param). **Default**: `0.9`.
- `[-ta]`: Threshold for the fraction of aligned positions (global MUSCL search param). **Default**: `0.9`.
- `[-ti]`: Peptide identity threshold  (global MUSCL search param). **Default**: `0.4`.
- `[-op]`: Search identifier. **Default**: The stem of the annotations file.
- `[-cn]`: Number of CPUs to use. **Default**: The available CPU count.
- `[-odp]`: Path to the output directory. **Default**: The current directory.
- `[--global-search]`: Enables only the global search mode.
- `[--local-search]`: Enables only the local search mode.
- `[--debug]`: Enables debug mode, saving input and output files for the local search.
- `[--soft-force]`: Overwrites the results database if it already exists.
- `[--hard-force]`: Overwrites all internal files if they already exist.
- `[--csv]`: Outputs a `.zip` file with a reduced set of results in CSV format.
- `[-v, --version]`: Shows the program's version number and exits.
- `[-h, --help]`: Shows this help message and exits.

---
> NOTE:
> If neither flag `--global-search` or `--local-search` is specified, both searches will run by default.
---

Example: Human Y chromosome
---------------------

The following steps demonstrate how to run `exonize` on a test dataset.

**Download the test data**


* If you have installed the package from the repo, move into the `test_human_chrom_Y` directory and download the test data:

```
cd test_human_chrom_Y
source fetch_data.sh
```

  * If you installed the package via `pip`, you can download the script here: [`fetch_data.sh`](https://github.com/msarrias/exonize/blob/main/test_human_chrom_Y/fetch_data.sh).

**Run exonize:**

```
exonize Homo_sapiens.GRCh38.111.chromosome.Y.gff3 \
       Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz \
       --gene-annot-feature gene \
       --cds-annot-feature CDS \
       --transcript-annot-feature mRNA \
       --output_prefix Homo_sapiens_chrom_Y \
       --csv
```

This command will enable the global and local searches.
