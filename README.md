```
███████╗██╗  ██╗ ██████╗ ███╗   ██╗██╗███████╗███████╗
██╔════╝╚██╗██╔╝██╔═══██╗████╗  ██║██║╚══███╔╝██╔════╝
█████╗   ╚███╔╝ ██║   ██║██╔██╗ ██║██║  ███╔╝ █████╗
██╔══╝   ██╔██╗ ██║   ██║██║╚██╗██║██║ ███╔╝  ██╔══╝
███████╗██╔╝ ██╗╚██████╔╝██║ ╚████║██║███████╗███████╗
╚══════╝╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚═╝╚══════╝╚══════╝
```

Exonize is an open-source command-line tool to identify and classify coding exon duplications in annotated genomes. The tool uses both; local and global alignment methods, supported by a graph-based framework to handle clusters of exons formed through repetitive duplication events. Exonize categorizes the interdependence among duplicated exons or exon clusters across transcripts. For data parsing and additional analysis, the `exonize_analysis` module is available for Python notebooks.

## Dependencies
For running Exonize, a local installation of [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) and [MUSCLE](https://www.drive5.com/muscle/) are required.

## Getting started
You are best off installing `exonize` from PyPI.org using
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

## Required arguments
## Optional arguments
## Usage
## Analyzing an example dataset
## Citation
