import sys, os
import subprocess
import gffutils 
import random
import pandas as pd
import numpy as np
import pybedtools
import pickle, copy, re
import portion as P
from collections import Counter
import matplotlib.pyplot as plt
from BCBio import GFF
from Bio.Seq import Seq
from Bio import SeqIO, SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML, NCBIWWW
from Bio.SeqFeature import SeqFeature, FeatureLocation
from dna_features_viewer import GraphicFeature, GraphicRecord
pd.options.display.float_format = "{:.2f}".format


def create_database(in_file_path, outfile_path):
    gffutils.create_db(in_file_path, dbfn=outfile_path,
                            force = True,keep_order=True,
                            merge_strategy='merge',
                            sort_attribute_values=True)
    

def generate_gene_hierarchy_dict(db):
    gene_hierarchy_dict = {}
    for gene in db.features_of_type('gene'):
        features = {}
        for mRNA_annot in db.children(gene.id, featuretype='mRNA', order_by='start'):
            temp_i = []
            for child in db.children(mRNA_annot.id, featuretype = ['CDS', 'exon',
                                                                   'intron', 
                                                                   'five_prime_UTR',
                                                                   'three_prime_UTR'],
                                     order_by='start'):
                if P.open(child.start, child.end):
                    temp_i += [{'coord': P.open(child.start, child.end),
                              'id':child.id, 
                              'strand': child.strand,
                              'type': child.featuretype}]
            # sort first by the start and then by the end
            temp_j = sorted(temp_i, key = lambda item: (item['coord'].lower, item['coord'].upper))
            features[mRNA_annot.id] = temp_j
        gene_hierarchy_dict[gene.id] = features
    return gene_hierarchy_dict


def get_gene_hierarchy_dict_with_coding_exons(gene_hierarchy_dict):
    gene_hierarchy_dict_with_coding_exons = {}
    copy_gene_hierarchy_dict = copy.deepcopy(gene_hierarchy_dict)
    gene_transcript_dict = {}
    for gene_id, gene_hierarchy_dict_ in copy_gene_hierarchy_dict.items():
        temp_transcript_dict = {}
        gene_transcript_dict[gene_id] = transcript_interval_dict(gene_hierarchy_dict_)
        for transcript_id, transcript_hierarchy_dict in gene_transcript_dict[gene_id].items():
            temp_transcript_dict[transcript_id] = get_coords_with_coding_exons(transcript_hierarchy_dict)
        gene_hierarchy_dict_with_coding_exons[gene_id] = temp_transcript_dict
    return gene_transcript_dict, gene_hierarchy_dict_with_coding_exons


def simulate_genes_with_coding_exons_duplications(duplicate_coding_exon_n,
                                                  insert_intron,
                                                  genes_sample_dict,
                                                  gene_hierarchy_dict,
                                                  genome_seq, db):
    simulated_genes = {}
    copy_genes_sample_dict = copy.deepcopy(genes_sample_dict)
    for gene_id, hierarchy_dict in copy_genes_sample_dict.items():
        transcript_id = next(iter(gene_hierarchy_dict[gene_id].keys()))
        temp_dict = {}
        exon_count = 0
        sim_gene_coord = []
        sim_gene_comp = []
        genome_seq_copy = copy.deepcopy(genome_seq)
        intervals_list = list(hierarchy_dict.keys())
        for idx, (interval, annotation) in enumerate(hierarchy_dict.items()):
            sim_gene_coord.append(interval)
            sim_gene_comp.append(annotation['id'])
            if annotation['type'] == 'coding_exon':
                exon_count += 1
            if exon_count == duplicate_coding_exon_n:
                if insert_intron == 'before':
                    insert_in_interval = intervals_list[idx-1]
                    insert_in = copy.deepcopy(hierarchy_dict[insert_in_interval])
                if insert_intron == 'after':
                    insert_in_interval = intervals_list[idx+1]
                    insert_in = copy.deepcopy(hierarchy_dict[insert_in_interval])
                if insert_in['type'] == 'intron':
                    start, end = insert_in_interval.lower, insert_in_interval.upper
                    center = (end - start) // 2
                    # we want to leave an evolutionary marker so we insert the 
                    # exon duplication leaving 10 bp in the start and end of the intron.
                    if center > 10:
                        temp_list = copy.deepcopy(sim_gene_coord)
                        temp_coords = [P.open(start, start + center), interval, P.open(start + center, end)]
                        if start > temp_list[-1].lower:
                            sim_gene_coord.extend(temp_coords)
                        else:
                            for coord_idx, coord in enumerate(temp_list):
                                if start < coord.lower:
                                    sim_gene_coord = (sim_gene_coord[:(coord_idx -1)] 
                                                      + temp_coords 
                                                      + sim_gene_coord[coord_idx:])
                                    break
                        if insert_intron == 'before':
                            sim_gene_coord.extend(intervals_list[idx+1:])
                            sim_gene_comp.extend([insert_in['id']])
                            sim_gene_comp.extend([hierarchy_dict[intv]['id'] for intv in intervals_list[idx:]])
                        if insert_intron == 'after':
                            sim_gene_coord.extend(intervals_list[idx+2:])
                            sim_gene_comp.extend([insert_in['id'], annotation['id'], insert_in['id']])
                            sim_gene_comp.extend([hierarchy_dict[intv]['id'] for intv in intervals_list[idx+2:]])
                    else:
                        print('the intron is too short to leave an evolutionary marker')
                else:
                    print('two consecutive exons')
                break  
                
        temp_sim_gene_coord = []
        for idx, i in enumerate(sim_gene_coord):
            if i.lower not in [j.upper  for j in sim_gene_coord] and idx != 0:
                temp_sim_gene_coord.append(P.open((i.lower - 1), i.upper))
            else:
                temp_sim_gene_coord.append(P.open(i.lower, i.upper))
                
        if insert_intron == 'before':
            mut_loct = 1
        if insert_intron == 'after':
            mut_loct = 0
        ce_dup_loc = (duplicate_coding_exon_n - mut_loct)
        insertion_in_intron = {}
        ce_exon_counter = 0
        for feat_idx, feat in enumerate(sim_gene_comp):
            if 'exon' in feat:
                if ce_exon_counter == ce_dup_loc:
                    insertion_in_intron['before'] = temp_sim_gene_coord[feat_idx-1]
                    insertion_in_intron['after'] = temp_sim_gene_coord[feat_idx+1]
                    break
                else:
                    ce_exon_counter += 1
        seq = ''            
        structure_breakdown = {}
        for sim_coord, sim_id in zip(temp_sim_gene_coord, sim_gene_comp):
            if 'intron' in sim_id:
                temp = copy.deepcopy(genome_seq_copy[(sim_coord.lower):(sim_coord.upper)])
                if sim_coord == insertion_in_intron['before']:
                    if temp[-2:] != 'AG':
                        temp = temp[:-2] + 'AG'
                if sim_coord == insertion_in_intron['after']:
                    if temp[:2] != 'GT':
                        temp = 'GT' + temp[2:]
            else:
                temp = copy.deepcopy(genome_seq_copy[(sim_coord.lower):(sim_coord.upper)])
            if sim_id not in structure_breakdown:
                structure_breakdown[sim_id] = copy.deepcopy(temp)
            else:
                # since we're inserting an exon in the middle of the intron
                structure_breakdown[sim_id + '_1'] = copy.deepcopy(temp)
            seq += temp
        simulated_genes[gene_id] = {'coord': P.closed(db[transcript_id].start, db[transcript_id].end),
                                   'features_order': sim_gene_comp,
                                   'features_intervals': temp_sim_gene_coord,
                                    'structure_breakdown':structure_breakdown,
                                   'seq': seq}
    return simulated_genes


def generate_new_genome_seq(sorted_gene_sample, gene_hierarchy_dict, genome_seq, db):
    ch_seq_with_dipl = ''
    end = 0
    len_dup = 0
    continue_ = True
    for idx, (sim_gene, sim_trans) in enumerate(sorted_gene_sample.items()):
        if continue_:
            transcript = next(iter(gene_hierarchy_dict[sim_gene].keys()))
            start, e = sim_trans['coord'].lower, sim_trans['coord'].upper
            temp = (copy.deepcopy(genome_seq[end:start]) + sim_trans['seq'])
            if  idx < (len(sorted_gene_sample)-1):
                ch_seq_with_dipl += temp
            else:
                ch_seq_with_dipl += (temp + genome_seq[e:])
                continue_ = False
            end = sim_trans['coord'].upper
        len_dup += len(sim_trans['seq']) - (db[transcript].end - db[transcript].start)  
    return ch_seq_with_dipl, len_dup


def get_simulated_new_intervals(genome_seq, new_genome_seq, sorted_gene_sample, coding_exons_interv):
    sim_genes_interv_in_chr = {gene_id : {'seq': [i.span() for i in re.finditer(gene_dict['seq'],
                                                                            new_genome_seq)],
                                     'coding_exons_dups': {exon_id: [
                                         i.span() 
                                         for i in re.finditer(
                                             genome_seq[intev.lower:intev.upper],
                                             new_genome_seq
                                         )
                                     ] for exon_id, intev in coding_exons_interv[gene_id].items()
                                                          },
                                      'coding_exons': {exon_id: [
                                          i.span() 
                                         for i in re.finditer(
                                             genome_seq[intev.lower:intev.upper],
                                             genome_seq
                                         )
                                     ] for exon_id, intev in coding_exons_interv[gene_id].items()
                                                          }
                                     }
                           for gene_id, gene_dict in sorted_gene_sample.items()}
    return sim_genes_interv_in_chr


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

def transcript_interval_dict(gene_dict):
    UTR_features = ['five_prime_UTR','three_prime_UTR']
    # we look at all genes' transcripts
    transcript_dict = {}
    copy_gene_dict = copy.deepcopy(gene_dict)
    for transcript_id, trans_dict in copy_gene_dict.items():
        interval_dict = {}
        # create an interval dict
        for idx, feature_annot in enumerate(trans_dict):
            feature_inv = feature_annot['coord']
            if feature_inv not in interval_dict:
                interval_dict[feature_inv] = {'id': feature_annot['id'], 'type': feature_annot['type']}
            # CASE 1: duplicated intervals 
            elif feature_inv in interval_dict:
                 # CASE 1.1: if there is an utr, followed by exon we keep the utr:
                if interval_dict[feature_inv]['type'] in UTR_features:
                    continue
                # CASE 1.2: if there is an exon folowed by an utr, we keep the utr:
                if interval_dict[feature_inv]['type'] == 'exon' and feature_annot['type'] in UTR_features:
                    interval_dict[feature_inv] = {'id': feature_annot['id'], 'type': feature_annot['type']}
                # CASE 1.3: if there is a CDS, followed by exon we register a coding exon:
                if interval_dict[feature_inv]['type'] == 'CDS' and feature_annot['type'] == 'exon':
                    interval_dict[feature_inv] = {'id': feature_annot['id'], 'type': 'coding_exon'}
                # CASE 1.4: if there is a exon, followed by a CDS we register a coding exon:
                if interval_dict[feature_inv]['type'] == 'exon' and feature_annot['type'] == 'CDS':
                    interval_dict[feature_inv]['type'] = 'coding_exon'
        transcript_dict[transcript_id] = interval_dict
    return transcript_dict
    
    
def get_overlaping_dict(interval_dict):
    copy_interval_dict = copy.deepcopy(interval_dict)
    overlaping_dict = {}
    intervals_list = list(copy_interval_dict.keys())
    for idx, (feat_interv, feature_annot) in enumerate(copy_interval_dict.items()):
        overlaping_dict[feat_interv] = []
        if idx != (len(copy_interval_dict) - 1):
            if feat_interv.overlaps(intervals_list[idx+1]):
                overlaping_dict[feat_interv] = intervals_list[idx+1] 
    return overlaping_dict


def get_transcript_overlaping_dict(transcript_dict):
    transcript_overlaping_dict = {}
    copy_transcript_dict = copy.deepcopy(transcript_dict)
    for idx, (transcript_id, transcript_info_dict) in enumerate(copy_transcript_dict.items()):
        transcript_overlaping_dict[transcript_id] = get_overlaping_dict(transcript_info_dict)
    return transcript_overlaping_dict


def sort_interval_dict(interval_dict):
    keys_ = list(interval_dict.keys())
    keys_ = sorted(keys_, key = lambda item: (item.lower, item.upper))
    return {key: copy.deepcopy(interval_dict[key]) for key in keys_}


def get_coords_with_coding_exons(interval_dict):
    UTR_features = ['five_prime_UTR','three_prime_UTR']
    overlaping_dict = get_overlaping_dict(interval_dict)
    # if there are no intersections nothing to do
    
#     for item in interval_dict.items():
#         print(item)
#     print('  ')
#     for item in overlaping_dict.items():
#         print(item)
#     print('  ')
#     print(' ---------------------- ')

    if not any(list(overlaping_dict.values())):
#         print('no more intersections')
        return sort_interval_dict(interval_dict)
    new_gene_interv_dict = {}
    intervals_list = list(interval_dict.keys())
    for int_idx, (interval, feature_descrip) in enumerate(interval_dict.items()):
        overlap_interval = copy.deepcopy(overlaping_dict[interval])
        if not overlap_interval:
            new_gene_interv_dict[interval] = feature_descrip
        else:
            # CASE 1: feature is an exon
            overlap_descript = copy.deepcopy(interval_dict[overlap_interval])
            if feature_descrip['type'] == 'exon':
                # CASE 1.1: overlap is with a CDS
                if overlap_descript['type'] == 'CDS':
                    if interval.intersection(overlap_interval):
                        # if the exon is comprised in the CDS region - can this happen?
                        if interval.contains(overlap_interval):
                            # the intersection will be a coding exon
                            intersection = interval.intersection(overlap_interval)
                            new_gene_interv_dict[intersection] = copy.deepcopy(feature_descrip)
                            new_gene_interv_dict[intersection]['type'] = 'coding_exon'
                            # exon / CDS will be a non coding exon
                            new_nc_exon_coord_list = [i for i in interval - overlap_interval]
                            if new_nc_exon_coord_list:
                                for new_nc_coord in new_nc_exon_coord_list:
                                    if (new_nc_coord.upper - new_nc_coord.lower) > 1:
                                        new_gene_interv_dict[new_nc_coord] = copy.deepcopy(feature_descrip)
                        else:
                            new_gene_interv_dict[overlap_interval] = copy.deepcopy(feature_descrip)
                            new_gene_interv_dict[overlap_interval]['type'] = 'coding_exon'
                        # pad the rest of the dictionary and repeat   
                        for intrv in intervals_list[int_idx+2:]:
                            if not intrv in new_gene_interv_dict:
                                new_gene_interv_dict[intrv] = copy.deepcopy(interval_dict[intrv])
#                         for key,value in new_gene_interv_dict.items():
#                             print(key, value)
                        return get_coords_with_coding_exons(sort_interval_dict(new_gene_interv_dict))      
                # CASE 1.2: overlap is with an UTR
                if overlap_descript['type'] in UTR_features:
                    if interval.intersection(overlap_interval):
                        if interval.contains(overlap_interval):
                            utr_coord = interval.intersection(overlap_interval)
                            new_gene_interv_dict[utr_coord] = copy.deepcopy(overlap_descript)
                            new_exon_coord_list = [i for i in interval - overlap_interval]
                            if new_exon_coord_list:
                                for new_exon_coord in new_exon_coord_list:
                                    if (new_exon_coord.upper - new_exon_coord.lower) > 1:
                                        new_gene_interv_dict[exon_coord] = copy.deepcopy(feature_descrip)
                        else:
                            new_gene_interv_dict[overlap_interval] = copy.deepcopy(overlap_descript)
                        for intrv in intervals_list[int_idx+2:]:
                            if not intrv in new_gene_interv_dict:
                                new_gene_interv_dict[intrv] = copy.deepcopy(interval_dict[intrv])
#                         for key,value in new_gene_interv_dict.items():
#                             print(key, value)
                        return get_coords_with_coding_exons(sort_interval_dict(new_gene_interv_dict))    
            # CASE 2: feature is a CDS
            if feature_descrip['type'] == 'CDS':
                # CASE 2.1: overlap is with an exon:
                if overlap_descript['type'] == 'exon':
                    if interval.intersection(overlap_interval):
                        # if the CDS region is comprised in the exon region 
                        if overlap_interval.contains(interval):
                            # the intersection will be a coding exon
                            intersection = interval.intersection(overlap_interval)
                            new_gene_interv_dict[intersection] = copy.deepcopy(overlap_descript)
                            new_gene_interv_dict[intersection]['type'] = 'coding_exon'
                            # CDS / exon will be a non coding exon
                            new_exon_coord_list = [i for i in overlap_interval - interval]
                            if new_exon_coord_list:
                                for new_exon_coord in new_exon_coord_list:
                                    if (new_exon_coord.upper - new_exon_coord.lower) > 1:
                                        new_gene_interv_dict[new_exon_coord] = copy.deepcopy(overlap_descript)
                        else:
                            new_gene_interv_dict[overlap_interval] = copy.deepcopy(feature_descrip)
                            new_gene_interv_dict[overlap_interval]['type'] = 'coding_exon'
                    for intrv in intervals_list[int_idx+2:]:
                        if not intrv in new_gene_interv_dict:
                            new_gene_interv_dict[intrv] = copy.deepcopy(interval_dict[intrv])
                    return get_coords_with_coding_exons(sort_interval_dict(new_gene_interv_dict))
            # CASE 3: feature is an UTR
            if feature_descrip['type'] in UTR_features:
                if interval.intersection(overlap_interval):
                    if overlap_descript['type'] == 'exon':
                        if overlap_interval.contains(interval):
                            new_gene_interv_dict[interval] = copy.deepcopy(feature_descrip)
                            new_exon_coord_list = [i for i in overlap_interval - interval]
                            if new_exon_coord_list:
                                for new_exon in new_exon_coord_list:
                                    if (new_exon.upper - new_exon.lower) > 1:
                                        new_gene_interv_dict[new_exon] = copy.deepcopy(overlap_descript)
                            for intrv in intervals_list[int_idx+2:]:
                                if not intrv in new_gene_interv_dict:
                                    new_gene_interv_dict[intrv] = copy.deepcopy(interval_dict[intrv])
                            return get_coords_with_coding_exons(sort_interval_dict(new_gene_interv_dict)) 


def get_utr_coords(gs, ge, gene_breakdown, genome_seq, mutation=False, dup_genome_seq=''):
    
    utr_coords = [j for i, j in zip(gene_breakdown['features_order'],
                                    gene_breakdown['features_intervals']) if 'UTR' in i]
    features = []
    if mutation:
        for utr_coord in utr_coords:
            s, e = utr_coord.lower, utr_coord.upper
            utr_coords_in_ch_wd = [i.span() for i in re.finditer(genome_seq[s:e],dup_genome_seq)
                                   if P.open(i.span()[0], i.span()[1]).intersection(P.open(gs, ge))]
        if utr_coords_in_ch_wd:
            for coord in utr_coords_in_ch_wd:
                s, e = coord
                features.append(GraphicFeature(start=s, end=e, color='gainsboro'))
    else:
        for utr_coord in utr_coords:
            s, e = utr_coord.lower, utr_coord.upper
            features.append(GraphicFeature(start=s, end=e, color='gainsboro'))
    return features 
        
            
def get_record(idx, gene_id, exons_dict, gs, ge, gene_breakdown, overlaping_hit_dict,
               mutations=False, genome_seq = '', mut_loc='', exon_dup=0, genome_seq_with_mut=''):
    features = []
    if mutations:
        if mut_loc == 'before':
            mut = 0
        if mut_loc == 'after':
            mut = 1
        for exon_id, exon_coords in exons_dict[gene_id]['coding_exons_dups'].items():
            for coord in exon_coords:
                s, e = coord
                color="#ffcccc"
                labl = 'exon ' + str(idx)
                if idx == (exon_dup + mut):
                    color = 'green'
                    labl = 'exon ' + str(idx - mut) + ' dup'
                if idx > (exon_dup + mut):
                    idx_ = idx - 1
                    labl = 'exon ' + str(idx_)
                idx += 1
                features.append(GraphicFeature(start=s, end=e, color=color,label=labl))
    else:
        for exon_id, exon_coords in exons_dict[gene_id]['coding_exons'].items():
            for coord in exon_coords:
                s, e = coord
                color="#ffcccc"
                labl = 'exon ' + str(idx)
                idx += 1
                features.append(GraphicFeature(start=s, end=e, color=color,label=labl))
            
    features += get_utr_coords(gs, ge, gene_breakdown, genome_seq=genome_seq, 
                               mutation=mutations, dup_genome_seq=genome_seq_with_mut)
    
    for idx_j, query_feat_coord in enumerate(overlaping_hit_dict['hit'], 1):
        qs, qe = query_feat_coord.lower, query_feat_coord.upper
        labl = 'exon ' + str(idx_j)
        features.append(GraphicFeature(start=qs, end=qe, color="blue",label=labl))
    record = GraphicRecord(first_index = gs, sequence_length = ge-gs, features=features)
    return record


def parse_exonerate_output(all_qresult):
    hits_feature_intervals = {}
    for query in all_qresult:
        query_dict = {}
        query_id = query.id.rsplit('.')[0]
        if query_id not in hits_feature_intervals:
            hits_feature_intervals[query_id] = []
        hits_dict = {}
        for hit_idx, hit in enumerate(query):
            hsp_dict = {}
            for hsp_idx, hsp in enumerate(hit): # HSP stands for High-scoring Segment Pair
                query_ranges = [frag.query_range for frag in hsp.fragments]
                hit_ranges = [frag.hit_range for frag in hsp.fragments]
                hsp_dict[hsp_idx] = { 'hit': [P.open(i[0], i[1]) for i in hit_ranges],
                                     'query': [P.open(i[0], i[1]) for i in query_ranges]}
            hits_dict[hit_idx] = hsp_dict
        hits_feature_intervals[query_id].append(hits_dict)
    return hits_feature_intervals


def get_overlaping_queries(protein_gene_dict, genes_location, hits_feature_intervals):
    overlaping_queries = {}
    for prot_id, hits_list in hits_feature_intervals.items():
        if prot_id in protein_gene_dict:
            for hit in hits_list:
                for hsp_idx, hsp in hit.items():
                    for key, value in hsp.items():
                        if all([genes_location[protein_gene_dict[prot_id]].contains(i) for i in value['hit']]):
                            if prot_id not in overlaping_queries:
                                overlaping_queries[prot_id] = copy.deepcopy(value)
                            else:
                                print('duplication')
    return overlaping_queries


def get_seq(strand, seq, start, end):
    start, end = sorted([start, end])
    temp_seq = (seq['+'][start-1:end-1])
    if strand == '+':
        return temp_seq
    if strand == '-':
        return str(Seq(temp_seq).reverse_complement())

    
def get_exon_overlaps(parent_child_dic, parent_dic):
    parents_exon_coverage = {}
    for parent_idx, (parent_id, parent_dict) in enumerate(parent_child_dic.items()):
        seq_id = parent_dict['annot'][0]
        parents_exon_coverage[parent_id] = {'parent': {'source': parent_dict['annot'][1], 
                                                       'strand':parent_dict['annot'][6],
                                                       'len': parent_dict['annot'][4] - parent_dict['annot'][3],
                                                       'coord':(parent_dict['annot'][3], parent_dict['annot'][4])},
                                            'structure':[]}
        
        out_file_a = '../temp/temp_file_a.bed'
        with open(out_file_a, "w") as out_handle:
            out_handle.write(seq_id 
                             + "\t"  + str(parent_dict['annot'][3]) 
                             + "\t"  + str(parent_dict['annot'][4])
                             + "\t"  + str(parent_dict['annot'][6])
                             + "\t"  + str(parent_dict['annot'][-1]['ID'][0]) + '\n')
        
        out_file_b = '../temp/temp_file_b.bed'
        with open(out_file_b, "w") as out_handle_b:  
            for exon_idx, (line_id, line_attrib) in enumerate(parent_dict['exon'].items()):
                id_ = line_attrib[-1]['exon_id'][0]
                out_handle_b.write(line_attrib[0] + "\t"  + str(line_attrib[3]) 
                                   + "\t"  + str(line_attrib[4])
                                   + "\t"  + str(line_attrib[6]) 
                                   + "\t"  + str(id_) + '\n')   
                parents_exon_coverage[parent_id]['structure'].append({id_ : {'coord':(line_attrib[3],
                                                                                          line_attrib[4]),
                                                                                 'type':'exon'}})
        a = pybedtools.BedTool(out_file_a)
        b = pybedtools.BedTool(out_file_b)
        c = a.subtract(b)
        os.remove(out_file_a)
        os.remove(out_file_b)
        introns_dict = dict(enumerate(list(c)))
        total_length = 0
        for intron_idx, line_interval in enumerate(introns_dict.values()):
            total_length += (line_interval.end - line_interval.start) - 1
            for idx, comp in enumerate(parents_exon_coverage[parent_id]['structure']):
                if line_interval.start < list(comp.values())[0]['coord'][0]:
                    parents_exon_coverage[parent_id]['structure'].insert(idx,
                                                                         {'intron' + str(intron_idx): {
                                                                             'coord':(line_interval.start,line_interval.end),
                                                                             'type':'intron'}}
                                                                        ) 
                                                                                
                    break
                
            parent_dic[parent_id]['intron'] = [len(introns_dict), total_length]
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
        if 'annot' in value:
            summary_dict['transcript_len'].append(value['annot'][1])
        if 'exon' in value:
            summary_dict['n_exons'].append(value['exon'][0])
            summary_dict['total_exon_len'].append(value['exon'][1])
        if 'intron' in value:
            summary_dict['n_introns'].append(value['intron'][0])
            summary_dict['total_intron_len'].append(value['exon'][1])
        if 'CDS' in value:
            summary_dict['n_cds'].append(value['CDS'][0])
            summary_dict['cds_len'].append(value['CDS'][1])
    
    for key,value in summary_dict.copy().items():
        if not value:
            del summary_dict[key]
            
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
    plt.xlabel('exons per mRNA', fontsize=18)
    plt.ylabel('count', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.subplot(1, 2, 2)
    plt.plot(bins, cdf(len(x), bins, x), '-x',color='#0504aa', markersize=2)
    plt.ylabel (r'Cumulative distribution',fontsize=18)
    plt.xlabel ('exons per mRNA',fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.show()
 
    
def get_overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def flatten(l):
    return [item for sublist in l for item in sublist]


def cdf(N,y,v):
    p = np.zeros((np.shape(y)))
    for i in range(len(y)):
        tr = y[i] - v
        p[i] =  len(tr[tr >= 0]) / N
    return p