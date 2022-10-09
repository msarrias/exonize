from data_base_op import * # object for handling DBs
import portion as P # for working with intervals
import pickle # saving and loading dictionries
import copy # deep copy mutable objects
import pandas as pd # stats df
import numpy as np # rounding
# for dumping FASTA files
from Bio import SeqIO 
from Bio.SeqRecord import SeqRecord 
from Bio.Seq import Seq
        
    
class genome_analysis(data_base_op):
    def __init__(self,
                 db_path,
                 gene_hierarchy_path,
                 gene_hierarchy=True):
        data_base_op.__init__(self, db_path)
        self.load_db()
        self.feature_types = list(self.db.featuretypes())
        self.UTR_features = ['five_prime_UTR',
                             'three_prime_UTR']
        self.feat_of_intst = ['CDS',
                              'exon',
                              'intron'] + self.UTR_features
        self.gene_hierarchy = gene_hierarchy
        if self.gene_hierarchy:
            self.gene_hierarchy_dict = self.read_pkl_file(
                gene_hierarchy_path
            )
        else:
            self.gene_hierarchy_path = gene_hierarchy_path
            self.create_gene_hierarchy_dict()
    
    
    def read_pkl_file(self, filepath):
        with open(filepath, 'rb') as handle:
            read_file = pickle.load(handle)
        return read_file
    
    
    def dump_pkl_file(self, out_filepath, obj):
        with open(out_filepath, 'wb') as handle:
            pickle.dump(obj, handle)
    
    
    def dump_fasta_file(self, out_filepath, seq_dict):
        with open(out_filepath, "w") as handle:
            for annot_id, annot_seq in seq_dict.keys():
                record = SeqRecord(Seq(annot_seq),
                                   id = str(annot_id),
                                          description= '')
                SeqIO.write(record, handle, "fasta")

            
    def get_annotations_dict(self):
        self.annot_dict = {annot_type: 
                           {idx:gene for idx, gene in 
                            enumerate(self.db.features_of_type(annot_type))
                           }
                           for annot_type in self.feature_types}
        
        
    def generate_basic_statistics(self):
        self.basic_stat = basic_stat = {"Number": [],
                                        "Size total (kb)":[],
                                        "Size mean (bp)": []}
        for annotation_type, ann_dic in self.annot_dict.items(): 
            lengths = [item.end - item.start 
                       for item in ann_dic.values()]
            self.basic_stat["Number"].append(len(ann_dic))
            self.basic_stat["Size total (kb)"].append(
                np.round(sum(lengths)/1000, 2)
            )
            self.basic_stat["Size mean (bp)"].append(
                np.round(np.mean(lengths), 2)
            )
        for value in self.basic_stat.values():
            value.append(sum(value))
            
        df = pd.DataFrame(self.basic_stat,
                          index=self.feature_types+['Total'])
        return df
    
    
    def create_gene_hierarchy_dict(self):
        self.gene_hierarchy_dict = {}
        for gene in self.db.features_of_type('gene'):
            features = {}
            for mRNA_annot in self.db.children(gene.id,
                                               featuretype='mRNA',
                                               order_by='start'):
                temp_i = []
                for child in self.db.children(mRNA_annot.id,
                                              featuretype=self.feat_of_intst,
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
            
        if not self.gene_hierarchy:
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
                    # CASE 1.2: if there is an exon folowed by an utr,
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
                    # CASE 1.4: if there is a exon, followed by a CDS we
                    # register a coding exon:
                    if (
                        interval_dict[feature_inv]['type'] == 'exon'
                        and feat_annot['type'] == 'CDS'
                    ):
                        interval_dict[feature_inv]['type'] = 'coding_exon'              
            transcript_dict[transcript_id] = interval_dict
        return transcript_dict
            
        
    def check_overlaps(self):
        inters_type_counter = dict.fromkeys(self.feat_of_intst, 0)
        self.gene_interval_dict = {}
        self.intersections_counter = {key: copy.deepcopy(inters_type_counter)
                                      for key in self.feat_of_intst}
         # do not allow exon and introns annotations 
        no_overlaps = ['intron', 'exon']
        # for collecting those annot
        self.genes_with_incorrect_intron_exon_overlaps = [] 
        for gene_id, gene_dict in self.gene_hierarchy_dict.items():
            self.gene_interval_dict[gene_id] = self.transcript_interval_dict(gene_dict)
            for transcript_id, transcript_dict in self.gene_interval_dict[gene_id].items():
                overlaping_dict = {}
                intervals_list = list(transcript_dict.keys())
                for idx, (feat_interv, feature_annot) in enumerate(transcript_dict.items()):
                    overlaping_dict[feat_interv] = []
                    if idx != (len(transcript_dict) - 1):
                        for interval_i in intervals_list[idx+1:]:
                            if feat_interv.overlaps(interval_i):
                                overlaping_dict[feat_interv].append(interval_i)
                for interv, interv_overlap in overlaping_dict.items():
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
        # exon/intron overlaping annotations
        for gene in set(self.genes_with_incorrect_intron_exon_overlaps):
            del self.gene_hierarchy_dict[gene]
            
            
    def get_overlaping_dict(self, interval_dict):
            copy_interval_dict = copy.deepcopy(interval_dict)
            overlaping_dict = {}
            intervals_list = list(copy_interval_dict.keys())
            for idx, (feat_interv, feature_annot) in enumerate(copy_interval_dict.items()):
                overlaping_dict[feat_interv] = []
                if idx != (len(copy_interval_dict) - 1):
                    if feat_interv.overlaps(intervals_list[idx+1]):
                        overlaping_dict[feat_interv] = intervals_list[idx+1]
            return overlaping_dict

        
    def sort_interval_dict(self, interval_dict):
        keys_ = list(interval_dict.keys())
        keys_ = sorted(keys_, key = lambda item: (item.lower, item.upper))
        return {key: copy.deepcopy(interval_dict[key]) for key in keys_}


    def get_coords_with_coding_exons(self, interval_dict):
        """HORRIBLE RECURSIVE FUNCTION - SORRY"""
        overlaping_dict = self.get_overlaping_dict(interval_dict)
        if not any(list(overlaping_dict.values())):
            return self.sort_interval_dict(interval_dict)
        new_gene_interv_dict = {}
        intervals_list = list(interval_dict.keys())
        for int_idx, (interval, feature_descrip) in enumerate(interval_dict.items()):
            overlap_interval = copy.deepcopy(overlaping_dict[interval])
            if not overlap_interval:
                new_gene_interv_dict[interval] = feature_descrip
            else:
                # CASE 1: feature is an exon
                overlap_descript = copy.deepcopy(
                    interval_dict[overlap_interval]
                )
                if feature_descrip['type'] == 'exon':
                    # CASE 1.1: overlap is with a CDS
                    if overlap_descript['type'] == 'CDS':
                        if interval.intersection(overlap_interval):
                            # if the exon is comprised in the CDS region 
                            # - can this happen?
                            if interval.contains(overlap_interval):
                                # the intersection will be a coding exon
                                intersection = interval.intersection(
                                    overlap_interval
                                )
                                new_gene_interv_dict[intersection] = copy.deepcopy(
                                    feature_descrip
                                )
                                new_gene_interv_dict[intersection]['type'] = 'coding_exon'
                                # exon / CDS will be a non coding exon
                                new_nc_exon_coord_list = [
                                    i for i in interval - overlap_interval
                                ]
                                if new_nc_exon_coord_list:
                                    for new_nc_coord in new_nc_exon_coord_list:
                                        if(
                                            (new_nc_coord.upper - new_nc_coord.lower) > 1
                                        ):
                                            new_gene_interv_dict[new_nc_coord] = copy.deepcopy(
                                                feature_descrip
                                            )
                            else:
                                new_gene_interv_dict[overlap_interval] = copy.deepcopy(
                                    feature_descrip
                                )
                                new_gene_interv_dict[overlap_interval]['type'] = 'coding_exon'
                            # pad the rest of the dictionary and repeat 
                            for intrv in intervals_list[int_idx+2:]:
                                if not intrv in new_gene_interv_dict:
                                    new_gene_interv_dict[intrv] = copy.deepcopy(
                                        interval_dict[intrv]
                                    )
                            return self.get_coords_with_coding_exons(
                                self.sort_interval_dict(new_gene_interv_dict)
                            )  
                    # CASE 1.2: overlap is with an UTR
                    if overlap_descript['type'] in self.UTR_features:
                        if interval.intersection(overlap_interval):
                            if interval.contains(overlap_interval):
                                utr_coord = interval.intersection(overlap_interval)
                                new_gene_interv_dict[utr_coord] = copy.deepcopy(
                                    overlap_descript
                                )
                                new_exon_coord_list = [
                                    i for i in interval - overlap_interval
                                ]
                                if new_exon_coord_list:
                                    for new_exon_coord in new_exon_coord_list:
                                        if (
                                            (new_exon_coord.upper - new_exon_coord.lower) > 1
                                        ):
                                            new_gene_interv_dict[exon_coord] = copy.deepcopy(
                                                feature_descrip
                                            )
                            else:
                                new_gene_interv_dict[overlap_interval] = copy.deepcopy(
                                    overlap_descript
                                )
                            for intrv in intervals_list[int_idx+2:]:
                                if not intrv in new_gene_interv_dict:
                                    new_gene_interv_dict[intrv] = copy.deepcopy(
                                        interval_dict[intrv]
                                    )
                            return self.get_coords_with_coding_exons(
                                self.sort_interval_dict(new_gene_interv_dict)
                            )
                # CASE 2: feature is a CDS
                if feature_descrip['type'] == 'CDS':
                    # CASE 2.1: overlap is with an exon:
                    if overlap_descript['type'] == 'exon':
                        if interval.intersection(overlap_interval):
                            # if the CDS region is comprised in the exon region 
                            if overlap_interval.contains(interval):
                                # the intersection will be a coding exon
                                intersection = interval.intersection(
                                    overlap_interval
                                )
                                new_gene_interv_dict[intersection] = copy.deepcopy(
                                    overlap_descript
                                )
                                new_gene_interv_dict[intersection]['type'] = 'coding_exon'
                                # CDS / exon will be a non coding exon
                                new_exon_coord_list = [
                                    i for i in overlap_interval - interval
                                ]
                                if new_exon_coord_list:
                                    for new_exon_coord in new_exon_coord_list:
                                        if (
                                            (new_exon_coord.upper - new_exon_coord.lower) > 1
                                        ):
                                            new_gene_interv_dict[new_exon_coord] = copy.deepcopy(
                                                overlap_descript
                                            )
                            else:
                                new_gene_interv_dict[overlap_interval] = copy.deepcopy(
                                    feature_descrip
                                )
                                new_gene_interv_dict[overlap_interval]['type'] = 'coding_exon'
                        for intrv in intervals_list[int_idx+2:]:
                            if not intrv in new_gene_interv_dict:
                                new_gene_interv_dict[intrv] = copy.deepcopy(
                                    interval_dict[intrv]
                                )
                        return self.get_coords_with_coding_exons(
                            self.sort_interval_dict(new_gene_interv_dict)
                        )
                # CASE 3: feature is an UTR
                if feature_descrip['type'] in self.UTR_features:
                    if interval.intersection(overlap_interval):
                        if overlap_descript['type'] == 'exon':
                            if overlap_interval.contains(interval):
                                new_gene_interv_dict[interval] = copy.deepcopy(
                                    feature_descrip
                                )
                                new_exon_coord_list = [
                                    i for i in overlap_interval - interval
                                ]
                                if new_exon_coord_list:
                                    for new_exon in new_exon_coord_list:
                                        if (
                                            (new_exon.upper - new_exon.lower) > 1
                                        ):
                                            new_gene_interv_dict[new_exon] = copy.deepcopy(
                                                overlap_descript
                                            )
                                for intrv in intervals_list[int_idx+2:]:
                                    if not intrv in new_gene_interv_dict:
                                        new_gene_interv_dict[intrv] = copy.deepcopy(
                                            interval_dict[intrv]
                                        )
                                return self.get_coords_with_coding_exons(
                                    self.sort_interval_dict(new_gene_interv_dict)
                                ) 
                            
                            
    def get_gene_hierarchy_dict_with_coding_exons(self):
        self.gene_hierarchy_dict_with_coding_exons = {}
        copy_gene_hierarchy_dict = copy.deepcopy(self.gene_hierarchy_dict)
        for gene_id, gene_hierarchy_dict_ in copy_gene_hierarchy_dict.items():
            temp_transcript_dict = {}
            for transcript_id, transcript_hierarchy_dict in self.gene_interval_dict[gene_id].items():
                temp_transcript_dict[transcript_id] = self.get_coords_with_coding_exons(
                    transcript_hierarchy_dict
                )
            self.gene_hierarchy_dict_with_coding_exons[gene_id] = temp_transcript_dict
        self.n_genes = len(self.gene_hierarchy_dict_with_coding_exons)
            
            