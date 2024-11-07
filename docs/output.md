
Output Description
============
- **`<output_prefix>_results.db`**: This is the main output database created by `exonize`.

    | **Table**                                    | **Description**                                                                                                                                                                                           |
    |------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
    | **Genes**                                | Lists all protein-coding genes analyzed <br/>for exon duplications.                                                                                                                                   |
    | **Global_matches_non_reciprocal**        | Contains all non-reciprocal matches identified by aligning<br/> pairs of representative exons using `MUSCLE`.<br/> All matches meet the peptide identity and aligned <br/>position fraction criteria. |
    | **Local_matches**                        | Contains all matches found by querying representative <br/>exons against the genes within which they are situated.<br/> These are raw, unfiltered matches found using `tblastx`.                      |
    | **Local_matches_non_reciprocal**         | Contains non-reciprocal matches found in the local search. <br/>This is a subset of the `Local_matches` table. Matches here do not overlap with<br/> `Global_matches_non_reciprocal` matches.         |
    | **Expansions**                           | Contains all groups of duplicates (expansions) found<br/> within genes. Each record represents a node in the<br/> expansion graph, created by combining local and/or <br/>global matches.                  |
    | **Expansions_full**                      | Includes expansions formed by full-mode events only.                                                                                                                                                  |
    | **Expansions_transcript_interdependence**| Classifies transcript interdependence for each expansion <br/>in the `Expansions_full` table.                                                                                                         |
    | **Expansions_full_tandem**               | Represents pairs of consecutive full events within<br/> an expansion, indicating whether they are in tandem <br/>or not.                                                                              |

- **`csvs.tar.gz`** (optional): compressed directory in .tar.gz format, containing individual files in CSV format for each table in the SQLite database, excluding the `Local_matches` and `Genes` tables.
- **`<output_prefix>.log`**: Search settings.

Database structure
---------------------
![Database](https://github.com/msarrias/exonize/raw/main/figures/database.png)


