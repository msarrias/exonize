from exon_analysis import *
from global_var import *
from exon_duplication_search import *


dirct = FASTA_DIRECT + 'pieris_napi/data/Pieris_napi_brakerProt_rename_agat.gff'
genome_dirct = FASTA_DIRECT + 'pieris_napi/Pieris_napi-GCA_905231885.1-softmasked.fa'
db_direct = DB_DIRECT + 'pieris_napi.db'
gene_hierarchy_path = PKL_DIRECT + 'Pieris_napi_gene_hierarchy_dict.pkl'

Pn_genome_analysis = ExonAnalysis(db_direct,
                                  gene_hierarchy_path,
                                  gene_hierarchy=True)

Pn_genome_analysis.get_annotations_dict()
Pn_genome_analysis.read_genome(genome_dirct)
Pn_genome_analysis.check_overlaps()
Pn_genome_analysis.get_gene_hierarchy_dict_with_coding_exons()
genes_list = list(Pn_genome_analysis.gene_hierarchy_dict.keys())
exon_dup_analysis = ExonDupSearch(ExonAnalysis_obj=Pn_genome_analysis, model='coding2genome')
_ = exon_dup_analysis.search_for_exon_duplicates(genes_list,100,19)
