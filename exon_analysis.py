from Bio import SeqIO # parsing FASTA files
from tqdm import tqdm # progress bar
import re # regular expressions
import matplotlib.pyplot as plt
from genome_analysis import *


class exon_analysis(genome_analysis):
    def __init__(self, db_path, gene_hierarchy_path, gene_hierarchy=True):
        genome_analysis.__init__(self, db_path, gene_hierarchy_path)
    
    def read_genome(self, genome_file_dir):
        genome = SeqIO.parse(open(genome_file_dir), 'fasta')
        self.genome = {
            fasta.id : {'+' : str(fasta.seq),
                        '-': str(fasta.seq.reverse_complement())
                       } 
            for fasta in genome
        }
        
        
    def get_within_gene_real_exon_dupl(self, exons_min_len,  filename=''):
        """
        this might take a long time to run
        """
        self.genes_coding_exons = {}
        for i, (ID, ID_dict) in zip(
            tqdm(range(self.n_genes)),
                 self.gene_hierarchy_dict_with_coding_exons.items()):
            trans_temp = {}
            for transcript_id, transcript_dict in ID_dict.items():
                temp = {
                    intev_dict['id']:
                    [
                        i.span() for i in re.finditer(
                        self.genome[intev_dict['chrom']]['+'][intv.lower:intv.upper],
                        self.genome[intev_dict['chrom']]['+'][self.db[ID].start - 20: 
                                                              self.db[ID].end + 20])
                    ] 
                    for intv, intev_dict in transcript_dict.items() 
                    if intev_dict['type'] == 'coding_exon'
                    and (intv.upper - intv.lower) > exons_min_len
                }
                # we want to get duplicates
                temp = {i:val for i, val in temp.items() if len(val)>1}
                if temp:
                    trans_temp[transcript_id] = temp
            if trans_temp:
                self.genes_coding_exons[ID] = trans_temp
        if filename != '':
            self.dump_pkl_file(filename, self.genes_coding_exons)
    
    
    def plot_histograms_on_real_exon_dup(self):
        exon_dup_per_transctipt = []
        len_of_dup_exon = []
        n_times_dup_exon = []
        for gene, gene_trascpt_dict in self.genes_coding_exons.items():
            for transcpt, trascpt_dict in gene_trascpt_dict.items():
                for exon_id, exons_list in trascpt_dict.items():
                    len_of_dup_exon.append(exons_list[0][1] - exons_list[0][0])
                    n_times_dup_exon.append(len(exons_list))
                exon_dup_per_transctipt.append(len(trascpt_dict))
        data = [n_times_dup_exon, len_of_dup_exon,exon_dup_per_transctipt]
        xlabels = ['Number of distinct exon duplications per transcript',
                   'Length of duplicated exons', 
                  'Single exon duplication frequency']
        plt.subplots(nrows=1, ncols=3, sharex=False, figsize=(13, 5))
        for i in range(3):
            plt.tight_layout()
            plt.subplot(1,3,(i+1)) 
            _, bins, _ = plt.hist(x=data[i], bins='auto',
                                  color='#0504aa',
                                  alpha=0.7,
                                  rwidth=0.85)
            plt.xlabel(xlabels[i], fontsize=14)
            plt.ylabel('count', fontsize=18)
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
        plt.show()

    