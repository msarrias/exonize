from Bio import SeqIO # parsing FASTA files
from tqdm import tqdm # progress bar
import re # regular expressions
import matplotlib.pyplot as plt
from genome_analysis import *


class exon_analysis(genome_analysis):
    def __init__(self, db_path,
                 gene_hierarchy_path,
                 gene_hierarchy=True):
        genome_analysis.__init__(self, 
                                 db_path,
                                 gene_hierarchy_path, 
                                 gene_hierarchy)
    
    
    def read_genome(self, 
                    genome_file_dir,
                    return_=False):
        genome = {
            fasta.id : {'+' : str(fasta.seq),
                        '-': str(fasta.seq.reverse_complement())
                       } 
            for fasta in SeqIO.parse(open(genome_file_dir), 
                             'fasta')
        }
        if not return_:
            self.genome = genome
        else:
            return genome
        
        
    def get_within_gene_real_exon_dupl(self,
                                       exons_min_len, 
                                       filename=''):
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
                        self.genome[intev_dict['chrom']]['+'][self.db[ID].start: 
                                                              self.db[ID].end])
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
        self.exon_dup_per_transctipt = []
        self.len_of_dup_exon = []
        self.n_times_dup_exon = []
        for gene, gene_trascpt_dict in self.genes_coding_exons.items():
            for transcpt, trascpt_dict in gene_trascpt_dict.items():
                for exon_id, exons_list in trascpt_dict.items():
                    self.len_of_dup_exon.append(
                        exons_list[0][1] - exons_list[0][0]
                    )
                    self.n_times_dup_exon.append(
                        len(exons_list)
                    )
                self.exon_dup_per_transctipt.append(
                    len(trascpt_dict)
                )
        data = [self.n_times_dup_exon,
                self.len_of_dup_exon,
                self.exon_dup_per_transctipt]
        xlabels = ['Number of distinct exon duplications per transcript',
                   'Length of duplicated exons', 
                  'Single exon duplication frequency']
        plt.subplots(nrows=1, ncols=3, sharex=False, figsize=(13, 5))
        for i in range(len(xlabels)):
            plt.tight_layout()
            plt.subplot(1,len(xlabels),(i+1)) 
            _, bins, _ = plt.hist(x=data[i], bins='auto',
                                  color='#0504aa',
                                  alpha=0.7,
                                  rwidth=0.85)
            plt.xlabel(xlabels[i], fontsize=14)
            plt.ylabel('count', fontsize=18)
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
        plt.show()
        
        
    def plot_exons_hist_cdf(self):
        temp = []
        for gene, gene_dict in self.gene_hierarchy_dict.items():
            for transcript, transcpt_list in gene_dict.items():
                temp.append(len([
                    i for i in transcpt_list
                    if 'exon' in i['type']
                ])
                           )

        plt.subplots(nrows=1,
                     ncols=2, 
                     sharex=False,
                     figsize=(13, 5))
        plt.subplot(1, 2, 1) 
        _, bins, _ = plt.hist(x=temp, bins='auto',
                              color='#0504aa',
                              alpha=0.7,
                              rwidth=0.85)
        plt.xlabel('exons per mRNA', fontsize=18)
        plt.ylabel('count', fontsize=18)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.subplot(1, 2, 2)
        plt.plot(bins,
                 self.cdf(len(temp), 
                          bins,
                          temp), 
                 '-x',color='#0504aa', 
                 markersize=2)
        plt.ylabel (r'Cumulative distribution',
                    fontsize=18)
        plt.xlabel ('exons per mRNA',
                    fontsize=18)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.show()

    def cdf(self, N,y,v):
        p = np.zeros((np.shape(y)))
        for i in range(len(y)):
            tr = y[i] - v
            p[i] =  len(tr[tr >= 0]) / N
        return p

    