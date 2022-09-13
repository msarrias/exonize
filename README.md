# Exon duplication analysis

### 1. Can tools for sequence comparison detect exon duplications?
- Use different tools to see different bias.
- Does Augustus allow for exon duplications?
### 2. How to test?
- Take well-annotated region and simulate fake exon duplications to compare tool outcomes.
  * Tools to evaluate:
    * [exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate-manual) - [test cases](#exonerate)
    * [miniprot](https://github.com/lh3/miniprot) - [test cases](#miniprot)
    * [spaln2](https://github.com/ogotoh/spaln)
- Problem: may not exist when predictions are RNA-supported.
  * what does Darwin ToL do?
### 3. How do species differ in % duplicated and divergence of duplicates?
### 4. How can we distinguish b/t TEs and recently duplicated genes? Repeat masking software may have issues.

## References:
- [Gao et al.](https://www.pnas.org/doi/10.1073/pnas.0911093106) 2009.
- [Zhang et al.](https://www.pnas.org/doi/10.1073/pnas.0603042103) 2006.
- [Merkin et al.](https://www.sciencedirect.com/science/article/pii/S2211124715002351) 2015.

## Comparison of the tools' outcomes on simulated genes with coding exons duplications.

### Exonerate
![alt text](https://github.com/msarrias/exon-duplication-analysis/blob/main/Figures/exonerate_alignms_with_is.jpg?raw=true)

### miniprot
![alt text](https://github.com/msarrias/exon-duplication-analysis/blob/main/Figures/miniprot_alignms_with_is.jpg?raw=true)
