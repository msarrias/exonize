
Output Description
============

### Main file: `<output_prefix>_results.db`
This is the main output database created by `exonize`.
    
##### Main tables

If you are only interested in the final results, these are the tables you want to look at:
    
| **Table**                                | **Description**                                                                                                                                                                             |
|------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Expansions**                           | Contains all groups of duplicates (expansions) found<br/> within genes. Each record represents a node in the<br/> expansion graph, created by combining local and (or) <br/>global matches. |
| **Expansions_full**                      | Includes expansions formed by full-mode events only.                                                                                                                                        |
| **Expansions_full_tandem**               | Represents pairs of consecutive full events within<br/> an expansion, indicating whether they are in tandem <br/>or not.                                                                    |
| **Expansions_transcript_interdependence**| Classifies transcript interdependence for each expansion <br/>in the `Expansions_full` table.                                                                                               |


##### Other tables
   
| **Table**                                    | **Description**                                                                                                                                                                                                                                                                   |
    |------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
    | **Genes**                                | Lists all protein-coding genes analyzed for exon duplications.                                                                                                                                                                                                                    |
    | **Local_matches**                        | Contains all matches found by querying representative <br/>exons against the genes within which they are situated.<br/> These are raw, unfiltered matches found using BLAST.                                                                                                      |
    | **Local_matches_non_reciprocal**         | This table is a filtered version of `Local_matches`, containing only <br/>non-reciprocal and reconciled matches identified in the <br/>BLAST local search. These matches meet the search criteria,<br/> including query coverage, in-frame alignment, <br/>and significant non-overlapping. |
    | **Global_matches_non_reciprocal**        | Contains all matches identified by aligning pairs of representative<br/> exons using `MUSCLE`. All matches meet the peptide identity and <br/>fraction of aligned position criteria.                                                                                              |

> **Note:**
> Combining the `Local_matches_non_reciprocal` and `Global_matches_non_reciprocal` tables provides a complete list of all matches identified by `exonize`. When full matches are detected by both methods, they are retained only in the `Local_matches_non_reciprocal` table and omitted from the `Global_matches_non_reciprocal` table.


### Additional files

- **`csvs.tar.gz`** (optional): compressed directory in .tar.gz format, containing individual files in CSV format for each table in the SQLite database, excluding the `Local_matches` and `Genes` tables.
- **`<output_prefix>.log`**: Search settings.

Database structure
---------------------

<a href="https://github.com/msarrias/exonize/blob/main/figures/database.png" target="_blank">
    <img src="https://github.com/msarrias/exonize/raw/main/figures/database.png" alt="Database Structure">
</a>

