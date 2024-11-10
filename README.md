```
███████╗██╗  ██╗ ██████╗ ███╗   ██╗██╗███████╗███████╗
██╔════╝╚██╗██╔╝██╔═══██╗████╗  ██║██║╚══███╔╝██╔════╝
█████╗   ╚███╔╝ ██║   ██║██╔██╗ ██║██║  ███╔╝ █████╗
██╔══╝   ██╔██╗ ██║   ██║██║╚██╗██║██║ ███╔╝  ██╔══╝
███████╗██╔╝ ██╗╚██████╔╝██║ ╚████║██║███████╗███████╗
╚══════╝╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚═╝╚══════╝╚══════╝
Marina Herrera Sarrias, Department of Mathematics, Stockholm University
Christopher Wheat, Department of Zoology, Stockholm University
Liam M. Longo, Earth-Life Science Institute (ELSI), Institute of Science Tokyo
Lars Arvestad, Department of Mathematics, Stockholm University
```

[![PyPI version](https://badge.fury.io/py/exonize.svg)](https://badge.fury.io/py/exonize)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=msarrias_exonize&metric=alert_status)](https://sonarcloud.io/dashboard?id=msarrias_exonize)

## Welcome!
`exonize` is an open-source command-line tool and [Python package](https://pypi.org/project/exonize/) for identifying and classifying coding exon duplications in annotated genomes. `exonize` identifies full exon duplications using local and global alignment methods and implements a graph-based framework to handle clusters of exons formed by repetitive duplication events. In addition, `exonize` categorizes the interdependence between duplicated exons (or groups of exons) across transcripts. For data parsing and downstream analysis, the `exonize_analysis` module is available for Python notebooks.

## Documentation

Please see the [Documentation](https://msarrias.github.io/exonize/) for a full user guide and an introductory tutorial to the `exonize_analysis` module.

## Installation

You are best off installing `exonize` from [PyPI.org](https://pypi.org/project/Exonize/1.0/) using

```bash
pip install exonize
```

If installing from the [GitHub](https://github.com/msarrias/exonize) repo

```bash
git clone git@github.com:msarrias/exonize.git
cd exonize
pip install .
```

You should now be able to run `exonize -h`.

`exonize` is tested with Python 3.9, 3.10, 3.11, 3.12

## Requirements

`exonize` requires a local installation of:


* [`BLAST+`](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) \[[download link](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)\]: exonize uses the `tblastx` program for conducting the local search.  
* [`MUSCLE`](https://www.drive5.com/muscle/) \[[download link](https://github.com/rcedgar/muscle/releases)\]: used for conducting the global search and correcting the identity of reconciled matches.  
* [`SQLite`](https://www.sqlite.org/index.html)[[download link](https://www.sqlite.org/download.html)] : for storing the search results.



## Support

If you need help with `exonize`, contact us! To report a bug or request a new feature, open an **[Issue](https://github.com/msarrias/exonize/issues)** on the `exonize` repo.

## Citation

If you use `exonize` in a publication, please cite:
```
TBA
```

