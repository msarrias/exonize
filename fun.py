import sys, os
import subprocess
import gffutils 
import pandas as pd
import numpy as np
import pybedtools
import pickle
import matplotlib.pyplot as plt
from BCBio import GFF
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Blast import NCBIXML, NCBIWWW
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
pd.options.display.float_format = "{:.2f}".format


def create_database(in_file_path, outfile_path):
    gffutils.create_db(in_file_path, dbfn=outfile_path,
                            force = True,keep_order=True,
                            merge_strategy='merge',
                            sort_attribute_values=True)
    
    
def get_annotations_dict(db):
    annotation_types = list(db.featuretypes())
    annot_dict = {}
    for annot_type in annotation_types:
        annot_dict[annot_type] = {idx:gene for idx,
                                  gene in enumerate(db.features_of_type(annot_type))
                                 }
    return annot_dict


def get_annotations_df(annot_db, annotation_list = ['mRNA', 'transcript']):
    df = pd.DataFrame(columns=['seq id',
                               'source',
                               'feat.',
                               'start',
                               'end',
                               'score',
                               'strand',
                               'frame',
                               'ID'])
    i = 0
    for annot in annotation_list:
        if annot in annot_db:
            for annot_id, annot_dic in annot_db[annot].items():
                #note that there might be different id attributes
                id_ = [id_key for id_key in list(annot_dic.attributes.keys()) if 'id' in id_key.lower()]
                temp = get_line_annotations(annot_dic)[:-1]
                if not id_:
                    temp.append('Na')
                else:
                    temp.append(annot_dic.attributes[id_[0]][0])
                df.loc[i] = temp
                i+=1
    return df


def get_line_annotations(l):
    return [l.seqid,
            l.source,
            l.featuretype,
            l.start,
            l.end,
            l.score,
            l.strand,
            l.frame,
            l.attributes
           ]


def create_parent_child_dic(annot_db_dict, parent_type = ['mRNA'], child_type = ['exon']):
    parent_dic = {}
    parent_child_dic = {}
    for annot_i in parent_type:
        if annot_i in annot_db_dict:
            for annot_id, annot_dic in annot_db_dict[annot_i].items():
                id_ = [id_key for id_key in list(annot_dic.attributes.keys()) if 'id' in id_key.lower()][0]
                annot_ID = annot_dic.attributes[id_][0]
                if annot_ID not in parent_dic.keys():
                    parent_dic[annot_ID] = {'annot': [1, annot_dic.end - annot_dic.start + 1]}
                else:
                    parent_dic[annot_ID][0] += 1
                    parent_dic[annot_ID][1] += annot_dic.end - annot_dic.start + 1
                parent_child_dic[annot_ID] = {'annot':get_line_annotations(annot_dic)}
                for annot_ii in child_type:
                    temp = {}
                    parent_dic[annot_ID][annot_ii] = [0, 0]
                    parent_child_dic[annot_ID][annot_ii] = {}
                    for line_key, line in annot_db_dict[annot_ii].items():
                        id_ = [id_key for id_key in list(line.attributes.keys()) if 'id' in id_key.lower()][0]
                        if line['Parent'][0] in annot_ID:
                            temp[line.attributes[id_][0]] = get_line_annotations(line)
                            parent_dic[line['Parent'][0]][annot_ii][-2] += 1
                            parent_dic[line['Parent'][0]][annot_ii][-1] += line.end - line.start + 1
                    parent_child_dic[annot_dic.attributes['ID'][0]][annot_ii] = temp
    
    return parent_dic, parent_child_dic


def get_exon_overlaps(parent_child_dic, parent_dic):
    parents_exon_coverage = {}
    for parent_idx, (parent_id, parent_dict) in enumerate(parent_child_dic.items()):
        out_file_a = '../temp/temp_file_a.bed'
        seq_id = parent_dict['annot'][0]
        with open(out_file_a, "w") as out_handle:
            out_handle.write(seq_id 
                             + "\t"  + str(parent_dict['annot'][3]) 
                             + "\t"  + str(parent_dict['annot'][4])
                             + "\t"  + str(parent_dict['annot'][6])
                             + "\t"  + str(parent_dict['annot'][-1]['ID'][0]) + '\n')
        parents_exon_coverage[parent_id] = {'parent': {'source': parent_dict['annot'][1],
                                                       'len': parent_dict['annot'][4] - parent_dict['annot'][3],
                                                       'coord':(parent_dict['annot'][3], parent_dict['annot'][4])},
                                            'exons':{}, 'introns':[]}
        out_file_b = '../temp/temp_file_b.bed'
        with open(out_file_b, "w") as out_handle_b:  
            for line_id, line_attrib in parent_dict['exon'].items():
                id_ = line_attrib[-1]['exon_id'][0]
                parents_exon_coverage[parent_id]['exons'][id_] = {'coord':(line_attrib[3], line_attrib[4])}
                out_handle_b.write(line_attrib[0] + "\t"  + str(line_attrib[3]) 
                                   + "\t"  + str(line_attrib[4])
                                   + "\t"  + str(line_attrib[6]) 
                                   + "\t"  + str(id_) + '\n') 
        a = pybedtools.BedTool(out_file_a)
        b = pybedtools.BedTool(out_file_b)
        c = a.subtract(b)
        os.remove(out_file_a)
        os.remove(out_file_b)
        introns_dict = dict(enumerate(list(c)))
        
        total_length = 0
        for line_id, line_interval in introns_dict.items():
            total_length += (line_interval.end - line_interval.start) - 1
            parents_exon_coverage[parent_id]['introns'].append((line_interval.start, line_interval.end))
        parent_dic[parent_id]['intron'] = [len(parents_exon_coverage[parent_id]['introns']), total_length]
        
    return parent_dic, parents_exon_coverage    


def basic_stat(annotation_dict):
    basic_stat = {"Number": [], "Size total (kb)":[], "Size mean (bp)": []}
    for annotation_type, ann_dic in annotation_dict.items():  
        lengths = [item.end - item.start for item in ann_dic.values()]
        basic_stat["Number"].append(len(ann_dic))
        basic_stat["Size total (kb)"].append(np.round(sum(lengths)/1000,2))
        basic_stat["Size mean (bp)"].append(np.round(np.mean(lengths),2))
        
    for value in basic_stat.values():
        value.append(sum(value))
      
    return basic_stat


def get_summary_dict(parent_dic):
    """
    python implementation of: https://darencard.net/blog/2018-01-10-gene-structure-summary/
    1. transcript ID
    2. transcript sequence length
    3. number of exons
    4. total exon sequence length
    5. number of introns
    6. total intron sequence length
    7. number of CDS chunks
    8. total CDS sequence length
    9. number of 5' UTR sequences
    10. total 5' UTR sequence length
    11. number of 3' UTR sequences
    12. total 3' UTR sequence length
    """
    summary_dict = {'transcript_id' : [], 
                    'transcript_len': [],
                    'n_exons': [],
                    'total_exon_len': [],
                    'n_introns': [],
                    'total_intron_len': [],
                    'n_cds': [],
                    'cds_len': []
                   }
    for key, value in parent_dic.items():
        summary_dict['transcript_id'].append(key)
        summary_dict['transcript_len'].append(value['annot'][1])
        summary_dict['n_exons'].append(value['exon'][0])
        summary_dict['total_exon_len'].append(value['exon'][1])
        summary_dict['n_introns'].append(value['intron'])
        summary_dict['total_intron_len'].append(value['exon'][2])
        summary_dict['n_cds'].append(value['CDS'][0])
        summary_dict['cds_len'].append(value['CDS'][1])

    for key, value in summary_dict.items():
        if key == 'transcript_id':
            value.append('Total')
        else:
            value.append(sum(value))
    return summary_dict


def write_gff_file(parent_child_dic, outfile= 'dtemp_file.gff3'):
    for parent_idx, (parent_id, parent_dict) in enumerate(parent_child_dic.items()):
        seq_id = parent_dict['annot'][0]
        rec = SeqRecord(' ', seq_id)
        qualifiers = dict(parent_dict['annot'][-1])
        if '' in qualifiers:
            del qualifiers['']
        qualifiers["source"] = parent_dict['annot'][1]
        qualifiers["score"] = parent_dict['annot'][5]
        if parent_dict['annot'][6] == '+':
            strand_ = 1
        else:
            strand_ = -1
        top_feature = SeqFeature(FeatureLocation(parent_dict['annot'][3],
                                                 parent_dict['annot'][4]),
                                 type=parent_dict['annot'][2],
                                 strand=strand_, qualifiers=qualifiers)
    
        temp_child_list = []
        for exon_id, exon_attrib in parent_dict['exons'].items():
            exon_qualif = dict(exon_attrib[-1])
            if '' in exon_qualif:
                del exon_qualif['']
            if exon_attrib[6] == '+':
                strand_ = 1
            else:
                strand_ = -1
            temp_child_list.append(SeqFeature(FeatureLocation(exon_attrib[3],
                                                              exon_attrib[4]),
                                              type=exon_attrib[2],
                                              strand=strand_,
                                              qualifiers=exon_qualif))
        top_feature.sub_features = temp_child_list
        rec.features = [top_feature]
        with open(out_file, "w") as out_handle:
            GFF.write([rec], out_handle)
            

def dump_parent_child_dict(parent_dic, parent_child_dic, filename):
    merged_dict = {'parent_dic':parent_dic, 'parent_child_dic':parent_child_dic}
    with open(filename, 'wb') as handle:
        pickle.dump(merged_dict, handle)

        
def sort_gff_file(path, infile):
    os.system('sort -k1,1 -k4,4n -k5,5n ' + path + infile + ' > ' + path + 'sorted_' + infile )
    
    
def covert_gtf_gff(gft_in_filename, gff_out_filename):
    subprocess.check_output(['gffread',
                             gft_in_filename,
                             "-o",
                             gff_out_filename,
                             '--t-adopt'
                            ])
    
    
def exons_hist(exons_list, bins, fig_size=(8,7), color='#0504aa'):
    plt.figure(figsize=fig_size)
    _, _, _ = plt.hist(x=exons_list, bins=bins, color=color, alpha=0.7, rwidth=0.85)
    plt.xlabel('exons', fontsize=18)
    plt.ylabel('count', fontsize=18)
    plt.xticks(range(min(exons_list), max(bins)))
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.show()
    
    
def plot_hist_cdf(x):
    plt.subplots(nrows=1, ncols=2, sharex=False, figsize=(13, 5))
    plt.subplot(1, 2, 1) 
    _, bins, _ = plt.hist(x=x, bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    plt.xlabel('exons', fontsize=18)
    plt.ylabel('count', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.subplot(1, 2, 2)
    plt.plot(bins, cdf(len(x), bins, x), '-x',color='#0504aa', markersize=2)
    plt.ylabel (r'Cumulative distribution',fontsize=18)
    plt.xlabel ('exons',fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.show()
    
    
def cdf(N,y,v):
    p = np.zeros((np.shape(y)))
    for i in range(len(y)):
        tr = y[i] - v
        p[i] =  len(tr[tr >= 0]) / N
    return p