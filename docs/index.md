```
███████╗██╗  ██╗ ██████╗ ███╗   ██╗██╗███████╗███████╗
██╔════╝╚██╗██╔╝██╔═══██╗████╗  ██║██║╚══███╔╝██╔════╝
█████╗   ╚███╔╝ ██║   ██║██╔██╗ ██║██║  ███╔╝ █████╗
██╔══╝   ██╔██╗ ██║   ██║██║╚██╗██║██║ ███╔╝  ██╔══╝
███████╗██╔╝ ██╗╚██████╔╝██║ ╚████║██║███████╗███████╗
╚══════╝╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚═╝╚══════╝╚══════╝
Marina Herrera Sarrias, Department of Mathematics, Stockholm University
Christopher Wheat, Department of Zoology, Stockholm University
Liam M. Longo, Earth-Life Science Institute, Institute of Science Tokyo
Lars Arvestad, Department of Mathematics, Stockholm University
```

Welcome!
============

`exonize` is an open-source command-line tool and [Python package](https://pypi.org/project/exonize/) for identifying and classifying coding exon duplications in annotated genomes. `exonize` identifies full exon duplications using local and global alignment methods and implements a graph-based framework to handle clusters of exons formed by repetitive duplication events. `exonize` introduces a classification schema to categorize the interdependence between duplicated exons (or groups of exons) across transcripts. For data parsing and downstream analysis, we introduce `exonize_analysis` a Python module within the `exonize` package.

Installation
---------------------
You are best off installing `exonize` from [PyPI.org](https://pypi.org/project/exonize/) using

```
pip install exonize
```

If installing from the [GitHub](https://github.com/msarrias/exonize) repo

```
git clone git@github.com:msarrias/exonize.git
cd exonize
pip install .
```

You should now be able to run `exonize -h`.

`exonize` is tested  with Python 3.9, 3.10, 3.11, 3.12

Requirements
---------------------

`exonize` requires a local installation of:

* [`BLAST+`](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) \[[download link](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)\]  
* [`MUSCLE (v.5.3)`](https://www.drive5.com/muscle/) \[[download link](https://github.com/rcedgar/muscle/releases)\]  
* [`SQLite`](https://www.sqlite.org/index.html)[[download link](https://www.sqlite.org/download.html)]

