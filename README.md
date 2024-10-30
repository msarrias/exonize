```
███████╗██╗  ██╗ ██████╗ ███╗   ██╗██╗███████╗███████╗
██╔════╝╚██╗██╔╝██╔═══██╗████╗  ██║██║╚══███╔╝██╔════╝
█████╗   ╚███╔╝ ██║   ██║██╔██╗ ██║██║  ███╔╝ █████╗
██╔══╝   ██╔██╗ ██║   ██║██║╚██╗██║██║ ███╔╝  ██╔══╝
███████╗██╔╝ ██╗╚██████╔╝██║ ╚████║██║███████╗███████╗
╚══════╝╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚═╝╚══════╝╚══════╝
```

Exonize is an open-source command-line tool for identifying and classifying coding exon duplications in annotated genomes. Exonize identifies full exon duplications using local and global alignment methods and implements a graph-based framework to handle clusters of exons formed through repetitive duplication events. In addition, exonize categorizes the interdependence between duplicated exons (or groups of exons) across transcripts. For data parsing and downstream analysis, the `exonize_analysis` module is available for Python notebooks.

<div align="center">
	<img src="https://github.com/msarrias/exonize/blob/msarrias-readme/figures/pipeline.png" style="width:50%;">
</div>
</div>

## Dependencies
For running `exonize`, a local installation of [`BLAST+`](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) and [`MUSCLE`](https://www.drive5.com/muscle/) are required.

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
