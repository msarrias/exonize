from Bio import SearchIO
import subprocess
import portion as P
import copy, re
import os # to remove temp file
from dna_features_viewer import GraphicFeature, GraphicRecord # for plotting 
import matplotlib.pyplot as plt
from global_var import *

class CompareTools:
    def __init__(self):
        self.exonerate_models = ['ungapped',
                                 'ungapped:trans',
                                 'affine:global',
                                 'affine:bestfit',
                                 'affine:local',
                                 'affine:overlap',
                                 'est2genome',
                                 'ner',
                                 'protein2dna',
                                 'protein2dna:bestfit',
                                 'protein2genome',
                                 'protein2genome:bestfit',
                                 'coding2coding',
                                 'coding2genome',
                                 'cdna2genome',
                                 'genome2genome'
                                ]

    @staticmethod
    def create_database(outfile_path):
        gffutils.create_db(in_file_path,
                           dbfn=outfile_path,
                           force = True,
                           keep_order=True,
                           merge_strategy='merge',
                           sort_attribute_values=True)


    @staticmethod
    def flatten_dict(dict_):
        return {i:j for chm, chm_dict in dict_.items()
                for i,j in chm_dict.items()
               }


    def run_exonerate(self,
                      model,
                      query_filename,
                      target_filename,
                      output_file):
        if not model in self.exonerate_models:
            raise Exception('Not a known Exonerate model')
        exonerate_command = ['exonerate',
                             '--showalignment', 'no',
                             '--model', model,
                             '--query', query_filename,
                             '--target', target_filename]
        with open(output_file, 'w') as out:
            subprocess.run(exonerate_command, stdout=out)
            
            
    @staticmethod
    def run_spaln(model,
                  output_format,
                  output_file,
                  genome_direct,
                  prot_queries_direct
                  ):
        spaln_command = [ "spaln",
                         "-Q", str(model), 
                         "-O", str(output_format),
                         "-o", output_file,
                         genome_direct, 
                         prot_queries_direct
                        ]
        subprocess.call(spaln_command)
        
    
    @staticmethod
    def run_miniprot(genome_direct,
                     prot_queries_direct,
                     output_file):
        miniprot_command = [MINIPROT_DIR,
                            "-t8",
                            genome_direct,
                            prot_queries_direct]
        with open(output_file, 'w') as out:
            subprocess.run(miniprot_command, stdout=out)
        
    
    @staticmethod
    def get_align_intervals_from_cigar_string(hit_dict):
        seq = hit_dict['cg']
        i = 0
        cigar_dict = {}
        query_start = hit_dict['query_coord'].lower
        target_start= hit_dict['target_coord'].lower
        while seq:
            find = re.search(r"[0-9]*", seq)
            n = int(find.group())
            op = seq[find.end()]
            # Alignment match. Consuming n*3 nucleotides 
            # and n amino acids 
            if op == 'M': 
                query_intv = P.open(query_start, query_start + n)
                query_start = query_start + n
                target_intv = P.open(target_start, target_start + (n*3))
                target_start = target_start + (n*3)
            # Insertion, Phase-0 intron. Consuming n amino acids                           
            elif op == 'I': 
                query_intv = P.open(query_start, query_start + n)
                query_start = query_start + n
            # Deletion. Consuming n*3 nucleotides
            elif op == 'D': 
                target_intv = P.open(target_start, target_start + (n*3))
                target_start = target_start + (n*3)
            # Frameshift deletion. Consuming n nucleotides
            elif op in ['F', 'N']: 
                target_intv = P.open(target_start, target_start + n)
                target_start = target_start + n
            # Frameshift match, Phase-1 intron, Phase-2 intron. 
            # Consuming n nucleotides and 1 amino acid
            elif op in ['U', 'V', 'G']: 
                target_intv = P.open(target_start, target_start + n)
                target_start = target_start + n
                query_intv = P.open(query_start, query_start + 1)
                query_start = query_start + 1
            cigar_dict[i] = {'pos': n,
                             'op': op, 
                             'query_intv': query_intv,
                             'target_intv': target_intv}
            i +=1
            seq = seq[find.end()+1:]
        return cigar_dict


    def get_cigar_string_overlapping_queries(self,
                                            hits_dit):
        overlapping_queries = {}
        for query_id, query_hit_list in hits_dit.items():
            exon_itv = []
            for hit in query_hit_list:
                if hit['overlaps']:
                    cigar_dict = self.get_align_intervals_from_cigar_string(hit)
                    exon_itv.append(cigar_dict)
            if exon_itv:
                overlapping_queries[query_id] = exon_itv
        return overlapping_queries
    
    
    @staticmethod
    def get_miniprot_matching_intervals(overlapping_queries,
                                        genes_order):
        hits_feature_intervals = {}
        for key, hits_list in overlapping_queries.items():
            temp = {}
            for hit_idx, hit in enumerate(hits_list):
                temp[hit_idx] = {'hit': [hit_dict['target_intv'] 
                                         for hit_id, hit_dict in hit.items()
                                         if hit_dict['op'] == 'M'],
                                'query': [hit_dict['query_intv'] 
                                          for hit_id, hit_dict in hit.items() 
                                          if hit_dict['op'] == 'M']}
            hits_feature_intervals[key] = temp
        return {i : hits_feature_intervals[i][0] for i in genes_order
                if i in hits_feature_intervals}


    def parse_spaln_output(self,
                           spaln_out_fname,
                           genes_loc,
                           genes_order):
        self.create_database(spaln_out_fname,
                             'temp.db')
        db = gffutils.FeatureDB('temp.db',
                                keep_order=True)
        gene_hierarchy_dict = {}
        for gene in db.features_of_type('gene'):
            features = []
            for mRNA_annot in db.children(gene.id,
                                          featuretype='mRNA',
                                          order_by='start'):
                query = []
                hit = []
                for child in db.children(mRNA_annot.id,
                                         featuretype = ['cds'],
                                         order_by='start'):
                    if (child.end - child.start) > 0:
                        attributes = child.attributes['Target'][0].rsplit(' ')
                        query_id = attributes[0]#.rsplit('.')[0]
                        query_start, query_end = map(int, attributes[1:3])
                        query.append(P.open(query_start,
                                            query_end))
                        hit.append(P.open(child.start,
                                          child.end))
                features.append({'query':query,
                                 'hit':hit,
                                 'overlaps': genes_loc[query_id].contains(
                                     P.open(hit[0].lower, hit[-1].upper))})
            if not prot_id in gene_hierarchy_dict:
                gene_hierarchy_dict[query_id] = features
            else:
                print('hit within same region')
        os.remove('temp.db')   
        return {i : gene_hierarchy_dict[i][0] for i in genes_order
                if i in gene_hierarchy_dict.keys() and 
                gene_hierarchy_dict[i][0]['overlaps'] == True} 
                
        
    @staticmethod
    def parse_exonerate_output_only_higher_score_hits(ex_out_fname,
                                                      ex_format,
                                                      genes_order):
        hits_feature_intervals = dict()
        for query in SearchIO.parse(ex_out_fname, 
                                    ex_format):
            query_id = query.id.rsplit('.')[0]
            if query_id not in hits_feature_intervals:
                hits_feature_intervals[query_id] = []
            hits_dict = dict()
            for hit_idx, hit in enumerate(query):
                hsp_dict = dict()
                # HSP stands for High-scoring Segment Pair
                for hsp_idx, hsp in enumerate(hit): 
                    query_ranges = [frag.query_range
                                    for frag in hsp.fragments]
                    hit_ranges = [frag.hit_range
                                  for frag in hsp.fragments]
                    hsp_dict[hsp.score] = {'score':hsp.score,
                                           'hit': [P.open(i[0], i[1])
                                                   for i in hit_ranges],
                                           'query': [P.open(i[0], i[1])
                                                     for i in query_ranges]}
                #let's choose the hsp with higher score
                max_score = max([(key,value['score'])
                                 for key, value in hsp_dict.items()],
                                key=lambda item:item[1])[0]
                hits_dict[hit_idx] = hsp_dict[max_score]
            hits_feature_intervals[query_id].append(hits_dict)
        new_hits_dict = dict()
        for ID, ID_dict in hits_feature_intervals.items():
            temp_score_list = [(hit_idx, hsp_key, hsp_value['score']) 
                               for hit_idx, hit in enumerate(ID_dict) 
                               for hsp_key, hsp_value in hit.items()]
            get_max_score = max(temp_score_list,
                                key=lambda item:item[2])
            hit_n, hsp_n, _ = get_max_score
            new_hits_dict[ID] = copy.deepcopy(ID_dict[hit_n][hsp_n])    

        return {i:new_hits_dict[i] for i in genes_order 
                if i in new_hits_dict}
                    
        
    @staticmethod
    def check_overlaps(target_genes_loc, hits_loc_dict):
        overlaps_dict = dict()
        if target_genes_loc.keys() == hits_loc_dict.keys():
            for query_id, query_dict in hits_loc_dict.items():
                gene_loc = P.closed(target_genes_loc[query_id].lower - 20, 
                                    target_genes_loc[query_id].upper + 20)
                overlaps_dict[query_id] = all([gene_loc.contains(i) 
                                               for i in query_dict['hit']])
        if all((overlaps_dict.values())):
            print('all hits are within the regions of interest')
        else:
            print([(key, value) for key, value in overlaps_dict.items() if value == False])
            
        
    
#     def parse_exonerate_output(self,
#                                ex_out_fname,
#                                ex_format,
#                                target_genes_loc,
#                                genes_order):
#         hits_feature_intervals = dict()
#         overlapping_queries = dict()
#         # our queries are the prot seqs
#         for query in SearchIO.parse(ex_out_fname, 
#                                     ex_format):
#             query_dict = dict()
#             query_id = query.id.rsplit('.')[0]
#             if query_id not in hits_feature_intervals:
#                 hits_feature_intervals[query_id] = []
#             hits_dict = dict()
#             # each query can have many hits
#             for hit_idx, hit in enumerate(query):
#                 hsp_dict = dict()
#                 # and each hit can have many HSPs
#                 # HSP stands for High-scoring Segment Pair
#                 for hsp_idx, hsp in enumerate(hit): 
#                     query_ranges = [frag.query_range 
#                                     for frag in hsp.fragments]
#                     hit_ranges = [frag.hit_range
#                                   for frag in hsp.fragments]
#                     hsp_dict[hsp_idx] = { 'hit':
#                                          [P.open(i[0], i[1]) 
#                                           for i in hit_ranges],
#                                          'query': 
#                                          [P.open(i[0], i[1]) 
#                                           for i in query_ranges]}
#                     gene_loc = P.closed(target_genes_loc[query_id].lower - 20, 
#                                         target_genes_loc[query_id].upper + 20)
#                     if all([gene_loc.contains(i) for i in hsp_dict[hsp_idx]['hit']]):
#                         if query_id not in overlapping_queries:
#                             overlapping_queries[query_id] = copy.deepcopy(
#                                 hsp_dict[hsp_idx]
#                             )
#                         else:
#                             print('extra overlapping hit:',
#                                   query_id,
#                                   hsp_dict[hsp_idx]['hit'])
#                 hits_dict[hit_idx] = hsp_dict
#             hits_feature_intervals[query_id].append(hits_dict)
#         return {i:overlapping_queries[i] for i in genes_order
#                 if i in overlapping_queries}
    
    
    @staticmethod
    def get_utr_coords(gs, ge,
                       gene_breakdown,
                       genome_seq,
                       mutation=False,
                       dup_genome_seq=''):
        utr_coords = [j for i, j in zip(gene_breakdown['features_order'],
                                        gene_breakdown['features_intervals']) 
                      if 'UTR' in i
                     ]
        features = []
        if mutation:
            for utr_coord in utr_coords:
                s, e = utr_coord.lower, utr_coord.upper
                utr_coords_in_ch_wd = [i.span() 
                                       for i in re.finditer(genome_seq[s:e],
                                                            dup_genome_seq)
                                       if P.open(i.span()[0],
                                                 i.span()[1]).intersection(P.open(gs, ge))
                                      ]
                if utr_coords_in_ch_wd:
                    for coord in utr_coords_in_ch_wd:
                        s, e = coord
                        features.append(GraphicFeature(start=s,
                                                       end=e,
                                                       color='gainsboro'))
        else:
            for utr_coord in utr_coords:
                s, e = utr_coord.lower, utr_coord.upper
                features.append(GraphicFeature(start=s,
                                               end=e,
                                               color='gainsboro'))
        return features 

    
    def get_record(self,
                   idx,
                   gene_id,
                   exons_dict,
                   gs, ge,
                   gene_breakdown,
                   overlapping_hit_dict,
                   mutations=False,
                   genome_seq = '',
                   mut_loc='',
                   exon_dup=0,
                   genome_seq_with_mut=''):
        features = []
        if mutations:
            if mut_loc == 'before':
                mut = 0
            elif mut_loc == 'after':
                mut = 1
            else:
                raise Exception('the mutation location should be either \
                before or after')
            for exon_id, exon_coords in exons_dict[gene_id]['coding_exons_dups'].items():
                for coord in exon_coords:
                    s, e = coord
                    color = "#ffcccc"
                    label = 'exon ' + str(idx)
                    if idx == (exon_dup + mut):
                        color = 'green'
                        label = 'exon ' + str(idx - mut) + ' dup'
                    if idx > (exon_dup + mut):
                        idx_ = idx - 1
                        label = 'exon ' + str(idx_)
                    idx += 1
                    features.append(GraphicFeature(start=s,
                                                   end=e,
                                                   color=color,
                                                   label=label))
        else:
            for exon_id, exon_coords in exons_dict[gene_id]['coding_exons'].items():
                for coord in exon_coords:
                    s, e = coord
                    color = "#ffcccc"
                    label = 'exon ' + str(idx)
                    idx += 1
                    features.append(GraphicFeature(start=s,
                                                   end=e,
                                                   color=color,
                                                   label=label))
        features += self.get_utr_coords(gs,
                                        ge,
                                        gene_breakdown,
                                        genome_seq=genome_seq, 
                                        mutation=mutations,
                                        dup_genome_seq=genome_seq_with_mut)

        for idx_j, query_feat_coord in enumerate(overlapping_hit_dict['hit'], 1):
            qs, qe = query_feat_coord.lower, query_feat_coord.upper
            label = 'exon ' + str(idx_j)
            features.append(GraphicFeature(start=qs,
                                           end=qe,
                                           color="blue",
                                           label=label))
        record = GraphicRecord(first_index=gs, 
                               sequence_length=ge-gs,
                               features=features)
        return record


    def visualize_seq_alignments(self,
                                 overlapping_queries_list,
                                 genes_locations_list,
                                 sim_genes_breakdown_list,
                                 mut_loc_list,
                                 genomes_list,
                                 exon_dup_n,
                                 id_gene_dict,
                                 exons_intervals):
        fig, ax = plt.subplots(len(max(overlapping_queries_list,
                                       key=len)),
                               3,
                               figsize=(13, 32))
        plt.tight_layout()
        genome_seq = genomes_list[0]
        queries_ref = overlapping_queries_list[0].keys()
        for col_idx, (overlapping_queries,
                      genes_location,
                      genes_breakdown,
                      mut_loc,
                      genome,
                      exon_intv) in enumerate(zip(overlapping_queries_list,
                                                  genes_locations_list,
                                                  sim_genes_breakdown_list,
                                                  mut_loc_list,
                                                  genomes_list,
                                                  exons_intervals)):
            overlapping_queries = {i: overlapping_queries[i] for i in queries_ref
                                  if i in overlapping_queries
                                 }
            mutations = True
            genome_seq_w_mut = genome
            if mut_loc == '':
                mutations = False
                genome_seq_w_mut = ''
            for idx_, (ID, overlapping_hit_dict) in enumerate(overlapping_queries.items()):
                coords_gene = genes_location[ID]
                breakdown_dict = self.flatten_dict(genes_breakdown)[ID]
                chm = breakdown_dict['chrom']
                gs, ge = coords_gene.lower, coords_gene.upper
                idx = 1
                record = self.get_record(idx,
                                         ID,
                                         exon_intv[chm],
                                         gs,
                                         ge, 
                                         breakdown_dict,
                                         overlapping_hit_dict,
                                         mutations=mutations,
                                         genome_seq = genome_seq,
                                         mut_loc=mut_loc,
                                         exon_dup=exon_dup_n,
                                         genome_seq_with_mut=genome_seq_w_mut)
                ax[idx_, col_idx].set_title(ID + ':' + id_gene_dict[ID], loc='center')
                record.plot(ax=ax[idx_, col_idx],
                            figure_width=5,
                            with_ruler=False)
