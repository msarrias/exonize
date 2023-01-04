from data_base_op import * # object for handling DBs
import portion as P # for working with intervals
import pickle # saving and loading dictionaries
import copy # deep copy mutable objects
import pandas as pd # stats df
import numpy as np # rounding
# for dumping FASTA files
from Bio import SeqIO, SearchIO
from Bio.SeqRecord import SeqRecord 
from Bio.Seq import Seq
        
    
class GenomeAnalysis(DataBaseOp):
    def __init__(self,
                 db_path,
                 annot_path,
                 gh_path,
                 verbose):
        DataBaseOp.__init__(self, 
                            db_path=db_path,
                            in_file_path=annot_path,
                            verbose = verbose)
        self.chrm_gene_dict = None
        self.n_genes = None
        self.gene_hierarchy_dict_with_coding_exons = None
        self.genes_with_incorrect_intron_exon_overlaps = None
        self.intersections_counter = None
        self.gene_interval_dict = None
        self.basic_stat = None
        self.annot_dict = None
        self.UTR_features = ['five_prime_UTR',
                             'three_prime_UTR']
        self.feat_of_interest = ['CDS',
                                 'exon',
                                 'intron'] + self.UTR_features
        self.gene_hierarchy_path = gh_path
        
        
    def generate_genome_info_for_analysis(self):
        self.create_parse_or_update_database()
        self.db_features = list(self.db.featuretypes())
        if os.path.exists(self.gene_hierarchy_path):
            self.gene_hierarchy_dict = self.read_pkl_file(
                self.gene_hierarchy_path
            )
        else:
            self.create_gene_hierarchy_dict()
        self.chrom_dict = {gene_id: self.db[gene_id].chrom 
                           for gene_id in self.gene_hierarchy_dict.keys()
                          } 
        self.gene_loc = {gene_id: P.open(self.db[gene_id].start, self.db[gene_id].end)
                         for gene_id in self.gene_hierarchy_dict.keys()
                        } 
        self.strand_dict = {gene_id: self.db[gene_id].strand 
                           for gene_id in self.gene_hierarchy_dict.keys()
                           } 
    
    
    @staticmethod
    def read_pkl_file(filepath):
        with open(filepath, 'rb') as handle:
            read_file = pickle.load(handle)
        return read_file
    
    
    @staticmethod
    def dump_pkl_file(out_filepath, obj):
        with open(out_filepath, 'wb') as handle:
            pickle.dump(obj, handle)
    
    
    @staticmethod
    def dump_fasta_file(out_filepath, seq_dict):
        with open(out_filepath, "w") as handle:
            for annot_id, annot_seq in seq_dict.items():
                record = SeqRecord(Seq(annot_seq),
                                   id = str(annot_id),
                                          description= '')
                SeqIO.write(record, handle, "fasta")

                
    @staticmethod
    def flatten(l):
        return [item for sublist in l for item in sublist]

    
    def create_gene_hierarchy_dict(self):
        self.gene_hierarchy_dict = {}
        for gene in self.db.features_of_type('gene'):
            features = {}
            mRNA_transcripts = [mRNA_t for mRNA_t in self.db.children(gene.id,
                                                                      featuretype='mRNA',
                                                                      order_by='start')
                               ]
            if mRNA_transcripts:
                for mRNA_annot in mRNA_transcripts:
                    temp_i = []
                    for child in self.db.children(mRNA_annot.id,
                                                  featuretype=self.feat_of_interest,
                                                  order_by='start'):
                        temp_i += [{'chrom':child.chrom,
                                    'coord': P.open(child.start, child.end),
                                    'id':child.id, 
                                    'strand': child.strand,
                                    'type': child.featuretype
                                   }]
                    # sort first by the start and then by the end
                    temp_j = sorted(temp_i, key = lambda 
                                    item: (item['coord'].lower,
                                           item['coord'].upper))
                    features[mRNA_annot.id] = temp_j
                self.gene_hierarchy_dict[gene.id] = features
                
        self.dump_pkl_file(self.gene_hierarchy_path, 
                           self.gene_hierarchy_dict)
            
  
    def transcript_interval_dict(self, gene_dict):
        transcript_dict = {}
        copy_gene_dict = copy.deepcopy(gene_dict)
        for transcript_id, trans_dict in copy_gene_dict.items():
            interval_dict = {} # create an interval dict
            for idx, feat_annot in enumerate(trans_dict):
                feature_inv = feat_annot['coord']
                temp_dict = {'chrom':feat_annot['chrom'],
                             'strand':feat_annot['strand'],
                             'id': feat_annot['id'],
                             'type': feat_annot['type']}
                if feature_inv not in interval_dict:
                    interval_dict[feature_inv] = copy.deepcopy(temp_dict)
                # CASE 1: duplicated intervals 
                elif feature_inv in interval_dict:
                    # CASE 1.1: if there is an utr, followed by exon,
                    # we keep the utr:
                    if interval_dict[feature_inv]['type'] in self.UTR_features:
                        continue
                    # CASE 1.2: if there is an exon followed by an utr,
                    # we keep the utr:
                    if (
                        interval_dict[feature_inv]['type'] == 'exon' 
                        and feat_annot['type'] in self.UTR_features
                    ):
                        interval_dict[feature_inv] = copy.deepcopy(temp_dict)
                    # CASE 1.3: if there is a CDS, followed by exon we 
                    # register a coding exon:
                    if (
                        interval_dict[feature_inv]['type'] == 'CDS'
                        and feat_annot['type'] == 'exon'
                    ):
                        interval_dict[feature_inv] = copy.deepcopy(temp_dict)
                        interval_dict[feature_inv]['type'] = 'coding_exon'
                    # CASE 1.4: if there is an exon, followed by a CDS we
                    # register a coding exon:
                    if (
                        interval_dict[feature_inv]['type'] == 'exon'
                        and feat_annot['type'] == 'CDS'
                    ):
                        interval_dict[feature_inv]['type'] = 'coding_exon'              
            transcript_dict[transcript_id] = interval_dict
        return transcript_dict
            
        
    def check_overlaps(self):
        inters_type_counter = dict.fromkeys(self.feat_of_interest, 0)
        self.gene_interval_dict = {}
        self.intersections_counter = {key: copy.deepcopy(inters_type_counter)
                                      for key in self.feat_of_interest}
         # do not allow exon and introns annotations 
        no_overlaps = ['intron', 'exon']
        # for collecting those annotations
        self.genes_with_incorrect_intron_exon_overlaps = [] 
        for gene_id, gene_dict in self.gene_hierarchy_dict.items():
            self.gene_interval_dict[gene_id] = self.transcript_interval_dict(gene_dict)
            for transcript_id, transcript_dict in self.gene_interval_dict[gene_id].items():
                overlapping_dict = {}
                intervals_list = list(transcript_dict.keys())
                for idx, (feat_interv, feature_annot) in enumerate(transcript_dict.items()):
                    overlapping_dict[feat_interv] = []
                    if idx != (len(transcript_dict) - 1):
                        for interval_i in intervals_list[idx+1:]:
                            if feat_interv.overlaps(interval_i):
                                overlapping_dict[feat_interv].append(interval_i)
                for interv, interv_overlap in overlapping_dict.items():
                    if interv_overlap:
                        for interval_j in interv_overlap:
                            interval_type = transcript_dict[interv]['type']
                            next_interval_type = transcript_dict[interval_j]['type']
                            if interval_type != next_interval_type:
                                if (interval_type in no_overlaps 
                                    and next_interval_type in no_overlaps):
                                    self.genes_with_incorrect_intron_exon_overlaps.append(gene_id)
                                else:
                                    self.intersections_counter[interval_type][next_interval_type] += 1
        # As instructed by Chris, we neglect those genes with 
        # exon/intron overlapping_dict annotations
        for gene in set(self.genes_with_incorrect_intron_exon_overlaps):
            del self.gene_hierarchy_dict[gene]
            
        if len(self.genes_with_incorrect_intron_exon_overlaps)>0:
#             self.dump_pkl_file(self.gene_hierarchy_path,
#                                self.gene_hierarchy_dict)
            print(f'{len(set(self.genes_with_incorrect_intron_exon_overlaps))} genes \
                  have been removed')
            for gene_i in set(self.genes_with_incorrect_intron_exon_overlaps):
                print(gene_i)
        self.gene_hierarchy_dict = self.read_pkl_file(self.gene_hierarchy_path)
            
            
    @staticmethod
    def get_overlapping_dict(interval_dict):
            copy_interval_dict = copy.deepcopy(interval_dict)
            overlapping_dict = {}
            intervals_list = list(copy_interval_dict.keys())
            for idx, (feat_interv, feature_annot) in enumerate(
                copy_interval_dict.items()
            ):
                overlapping_dict[feat_interv] = []
                if idx != (len(copy_interval_dict) - 1):
                    if feat_interv.overlaps(intervals_list[idx+1]):
                        overlapping_dict[feat_interv] = intervals_list[idx+1]
            return overlapping_dict

        
    @staticmethod
    def sort_interval_dict(interval_dict):
        keys_ = list(interval_dict.keys())
        keys_ = sorted(keys_, key = lambda item: (item.lower, item.upper))
        return {key: copy.deepcopy(interval_dict[key]) for key in keys_}

    
    @staticmethod
    def check_for_short_overlaps(x_intev, y_intev):
        i, j = x_intev.lower, x_intev.upper
        m, n = y_intev.lower, y_intev.upper
        if (m - i) == 1: m = i
        if (i - m) == 1: i = m
        if (n - j) == 1: j = n
        if (j - n) == 1: n = j
        if (j - m) == 1: m = m + 2
        x_intev, y_intev = P.open(i, j), P.open(m, n)
        o_intv = x_intev.intersection(y_intev)
        if o_intv:
            k, l = o_intv.lower, o_intv.upper
            if (k - i) == 1: k = i
            if (n - l) == 1: l = n
            o_intv = P.open(k, l)
        return x_intev, y_intev, o_intv

    
    @staticmethod
    def check_global_range(list_intervals):
        upper =list_intervals[0].upper
        for i in list_intervals[1:]:
            lower = i.lower
            if lower == upper + 1:
                upper = i.upper
            else:
                raise ValueError('set of intervals should be consecutive')
        return True

    
    @staticmethod
    def change_type(dic, feature):
        copy_dict = copy.deepcopy(dic)
        copy_dict['type'] = feature
        return copy_dict


    def get_coords_with_coding_exons(self, interval_dict):
        """
        get_coords_with_coding_exons is a recursive function that ...
        Let e = (i, j) and x = (m, n) be intervals
        Based on the assumption that i <= m and j <= n ==> e_i ∪ x_j = (i, n)
        and i\neq m and j \neq n . The intersections can be of the following type:
        Exon with intron/UTR annotations: 
        - overlaps: o = e_i  ∩ x_j = (k, l)
            - if o ⊆ e_i 
                - if k == i ==> (i, n) = (k, l) ∪ (l + 1, n)
                - if l == j  ==> (i, n) = (i, k-1) ∪ (k, l)
                - if k > i and l < j ==> (i, n) = (i, k-1) ∪ (k, l) ∪ (l + 1, n)
            - else:
                - (i, n) = (i, k-1) ∪ (k, l) ∪ (l+1 n)
        """
        ce_types = ['exon', 'coding_exon']
        introns_utr = ['intron'] + self.UTR_features
        overlapping_dict = self.get_overlapping_dict(interval_dict)
        # if there are no overlaps, it's done
        if not any(list(overlapping_dict.values())):
            final_res = self.sort_interval_dict(interval_dict)
            if self.check_global_range(list(final_res.keys())):
                return final_res
        # otherwise we start solving intersections so that the final transcript
        # annotation do not include overlaping annotations
        new_gene_interv_dict = {}
        intervals_list = list(interval_dict.keys())
        for int_idx, (x_intev, feature_descrip) in enumerate(interval_dict.items()):
            # if there is no overlap do nothing
            y_intev = copy.deepcopy(overlapping_dict[x_intev])
            if not y_intev:
                new_gene_interv_dict[x_intev] = feature_descrip
            else:
                x_dict = copy.deepcopy(feature_descrip)
                x_type = x_dict['type']
                y_dict = copy.deepcopy(interval_dict[y_intev])
                y_type = y_dict['type']
                x_intev, y_intev, o_intv = self.check_for_short_overlaps(x_intev, y_intev)
                if o_intv:
                    i, j = x_intev.lower, x_intev.upper
                    m, n = y_intev.lower, y_intev.upper
                    k, l = o_intv.lower, o_intv.upper
                    if x_intev.contains(y_intev):
                        if x_type not in introns_utr:
                            # overlaps between coding regions        
                            if x_type in ce_types and y_type == 'CDS':
                                new_gene_interv_dict[o_intv] = self.change_type(x_dict, 'coding_exon')
                            elif x_type == 'CDS' and y_type in ce_types:
                                new_gene_interv_dict[o_intv] = self.change_type(y_dict, 'coding_exon')
                            # overlaps between coding and non coding regions
                            elif y_type in introns_utr:
                                new_gene_interv_dict[o_intv] = y_dict
                            elif x_type in ce_types + ['CDS'] and y_type in introns_utr:
                                new_gene_interv_dict[o_intv] = y_dict
                            else: print('case 1.1')
                            if i == k:
                                new_gene_interv_dict[P.open(l+1, n)] = x_dict
                            elif j == n:
                                new_gene_interv_dict[P.open(i, k-1)] = x_dict
                            elif k > i and l < j:
                                new_gene_interv_dict[P.open(i, k-1)] = x_dict
                                new_gene_interv_dict[P.open(l+1, n)] = x_dict
                            # overlap between non coding and coding region
                            else: print('case 1.3')
                        elif x_type in introns_utr:
                                new_gene_interv_dict[x_intev] = x_dict
                        else: print('case 1.2')
                    elif y_intev.contains(x_intev):
                        if y_type not in introns_utr:
                            if y_type == 'CDS' and x_type in ce_types:
                                new_gene_interv_dict[o_intv] = self.change_type(x_dict, 'coding_exon')
                            elif y_type in ce_types and x_type == 'CDS':
                                new_gene_interv_dict[o_intv] = self.change_type(y_dict, 'coding_exon')
                            elif x_type in introns_utr:
                                new_gene_interv_dict[o_intv] = x_dict
                            else: print('case 2.1')
                            if n == l:
                                new_gene_interv_dict[P.open(m, k-1)] = y_dict
                            elif k == m:
                                new_gene_interv_dict[P.open(l+1, n)] = y_dict
                            elif k > m and l < n:
                                new_gene_interv_dict[P.open(m, k-1)] = y_dict
                                new_gene_interv_dict[P.open(l+1, n)] = y_dict
                            else: print((k, l) , (m, n), 'casw 2.3')
                        elif y_type in introns_utr:
                            new_gene_interv_dict[P.open(m, n)] = y_dict
                        else: print('case 2.4')
                    elif not x_intev.contains(y_intev) | y_intev.contains(x_intev):
                        if x_type in ce_types and y_type == 'CDS':
                            new_gene_interv_dict[P.open(i, k-1)] = x_dict
                            new_gene_interv_dict[o_intv] = self.change_type(x_dict, 'coding_exon')
                            new_gene_interv_dict[P.open(l+1, n)] = y_dict
                        elif x_type == 'CDS' and y_type in ce_types:
                            new_gene_interv_dict[P.open(i, k-1)] = x_dict
                            new_gene_interv_dict[o_intv] = self.change_type(y_dict, 'coding_exon')
                            new_gene_interv_dict[P.open(l+1, n)] = y_dict
                        elif y_type in introns_utr:
                            new_gene_interv_dict[P.open(i, k-1)] = x_dict
                            new_gene_interv_dict[P.open(k, n)] = y_dict
                        elif x_type in introns_utr:
                            new_gene_interv_dict[P.open(i, l)] = x_dict
                            new_gene_interv_dict[P.open(l + 1, n)] = y_dict
                        else: print('case 3.1')
                    else: print(x_intev, y_intev, x_intev.contains(y_intev),  y_intev.contains(x_intev), 'case 4.1')
                else:
                    new_gene_interv_dict[x_intev] = x_dict
                    new_gene_interv_dict[y_intev] = y_dict
                new_gene_interv_dict = self.sort_interval_dict(new_gene_interv_dict)
                for intrv in intervals_list[int_idx+2:]:
                    if intrv not in new_gene_interv_dict:
                        new_gene_interv_dict[intrv] = copy.deepcopy(interval_dict[intrv])
                    elif (new_gene_interv_dict[intrv]['type'] in ce_types 
                          and interval_dict[intrv]['type'] in introns_utr
                         ):
                        new_gene_interv_dict[intrv] = copy.deepcopy(interval_dict[intrv])
                return self.get_coords_with_coding_exons(self.sort_interval_dict(new_gene_interv_dict))           
                
                
    def get_gene_hierarchy_dict_with_coding_exons(self):
        self.gene_hierarchy_dict_with_coding_exons = {}
        copy_gene_hierarchy_dict = copy.deepcopy(self.gene_hierarchy_dict)
        for gene_id, gene_hierarchy_dict_ in copy_gene_hierarchy_dict.items():
            temp_transcript_dict = {}
            for transcript_id, transcript_dict in self.gene_interval_dict[gene_id].items():
                temp_transcript_dict[transcript_id] = self.get_coords_with_coding_exons(
                    transcript_dict
                )
            self.gene_hierarchy_dict_with_coding_exons[gene_id] = temp_transcript_dict
        self.n_genes = len(self.gene_hierarchy_dict_with_coding_exons)
        
        
    def get_chrom_gene_dict(self):
        self.chrm_gene_dict = {i:{} for i in self.genome.keys()}
        for gene_id, gene_dict in self.gene_hierarchy_dict_with_coding_exons.items():
            chm = self.db[gene_id].chrom
            self.chrm_gene_dict[chm][gene_id] = copy.deepcopy(gene_dict)
            
            