
Output Description
============

### Main file: `<output_prefix>_results.db`
This is the main output database created by `exonize`.
    
##### Main tables

If you are only interested in the final results, these are the tables you want to look at:
    
| **Table**                                | **Description**                                                                                                          |
|------------------------------------------|--------------------------------------------------------------------------------------------------------------------------|
| **Expansions**                           | Each record represents a node in an expansion graph,<br/> created by combining local and (or) global matches.            |
| **Expansions_full**                      | Reduced version of `Expansions_full` containing only <br/>expansion graphs constructed from full-length matches.              |
| **Expansions_full_tandem**               | Arranges the nodes of expansions in `Expansions_full`,<br/> in consecutive pairs, indicating whether they are in tandem. |
| **Expansions_transcript_interdependence**| Transcript interdependence classifications for expansions in <br/>`Expansions_full`.                                     |


##### Other tables

If you want to analyze the matches directly, these tables contain the raw data:
   
| **Table**                                    | **Description**                                                                                                                                                                                                                                        |
    |------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
    | **Genes**                                | Lists all protein-coding genes analyzed for exon duplications.                                                                                                                                                                                         |
    | **Local_matches**                        | Contains all matches found by querying representative exons<br/> against the genes within which they are situated. These are raw,<br/> unfiltered matches found using BLAST.                                                                           |
    | **Local_matches_non_reciprocal**         | Filtered version of `Local_matches`, containing only non-reciprocal <br/>and reconciled matches identified in the  local search. These matches <br/>meet the query coverage, in-frame alignment, and significant<br/> non-overlapping search criteria. |
    | **Global_matches_non_reciprocal**        | Contains all matches identified by aligning pairs of representative<br/> exons using `MUSCLE`. All matches meet the peptide identity and <br/>fraction of aligned position criteria.                                                                   |

> **Note:**
> Combining the `Local_matches_non_reciprocal` and `Global_matches_non_reciprocal` tables provides a complete list of all matches identified by `exonize`. When full-length matches are detected by both methods, they are retained only in the `Global_matches_non_reciprocal` table and omitted from the `Local_matches_non_reciprocal` table.


### Additional files

- **`csvs.tar.gz`** (optional): compressed directory in .tar.gz format, containing individual files in CSV format for each table in the SQLite database, excluding the `Local_matches` and `Genes` tables.
- **`<output_prefix>.log`**: Search settings.

Database structure
---------------------

<a href="https://github.com/msarrias/exonize/blob/main/figures/database.png" target="_blank">
    <img src="https://github.com/msarrias/exonize/raw/main/figures/database.png" alt="Database Structure">
</a>

