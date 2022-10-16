from Bio import SearchIO
import subprocess
import portion as P
import copy, re
import os # to remove temp file
from dna_features_viewer import GraphicFeature, GraphicRecord # for plotting 
import matplotlib.pyplot as plt

class compare_tools():
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
        
        
    def create_database(self, in_file_path, outfile_path):
        gffutils.create_db(in_file_path,
                           dbfn=outfile_path,
                           force = True,
                           keep_order=True,
                           merge_strategy='merge',
                           sort_attribute_values=True)
        
        
    def flatten_dict(self,
                     dict_):
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
            
            
    def run_splan(self,
                  model,
                  output_format,
                  output_file,
                  genome_direct,
                  prot_queries_direct
                  ):
        splan_command = [ "spaln",
                         "-Q", str(model), 
                         "-O", str(output_format),
                         "-o", output_file,
                         genome_direct, 
                         prot_queries_direct
                        ]
        subprocess.call(splan_command)
        
    
    def run_miniprot(self,genome_direct,
                     prot_queries_direct,
                     output_file):
        miniprot_command = ["/domus/h1/msarrias/bin/miniprot/miniprot",
                            "-t8",
                            genome_direct,
                            prot_queries_direct]
        with open(output_file, 'w') as out:
            subprocess.run(miniprot_command, stdout=out)
        
    
    def get_align_intervals_from_CIGAR_string(self, 
                                              hit_dict):
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
            # Delection. Consuming n*3 nucleotides 
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


    def get_CIGAR_string_overlaping_queries(self,
                                            hits_dit):
        overlaping_queries = {}
        for prot_id, query_hit_list in hits_dit.items():
            exon_itv = []
            for hit in query_hit_list:
                if hit['overlaps'] == True:
                    cigar_dict = self.get_align_intervals_from_CIGAR_string(hit)
                    exon_itv.append(cigar_dict)
            if exon_itv:
                overlaping_queries[prot_id] = exon_itv
        return overlaping_queries
    
    
    def get_miniprot_matching_intervals(self,
                                        overlaping_queries, 
                                        genes_order):
        hits_feature_intervals = {}
        for key, hits_list in overlaping_queries.items():
            temp = {}
            for hit_idx, hit in enumerate(hits_list):
                temp[hit_idx] = {'hit': [hit_dict['target_intv'] 
                                         for hit_id, hit_dict in hit.items()
                                         if hit_dict['op'] == 'M'],
                                'query': [hit_dict['query_intv'] 
                                          for hit_id, hit_dict in hit.items() 
                                          if hit_dict['op'] == 'M']}
            hits_feature_intervals[key] = temp
        return {key : hits_feature_intervals[i][0] for i in genes_order 
                if i in hits_feature_intervals} 


    def parse_splan_output(self, 
                           splan_out_fname,
                           genes_loc,
                           pg_dict,
                           genes_order):
        self.create_database(splan_out_fname,
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
                        prot_id = attributes[0].rsplit('.')[0]
                        query_start, query_end = map(int, attributes[1:3])
                        query.append(P.open(query_start,
                                            query_end))
                        hit.append(P.open(child.start,
                                          child.end))
                features.append({'query':query,
                                 'hit':hit,
                                 'overlaps': genes_loc[prot_id].contains(
                                     P.open(hit[0].lower, hit[-1].upper))})
            if not prot_id in gene_hierarchy_dict:
                gene_hierarchy_dict[prot_id] = features
            else:
                print('hit within same region')
        os.remove('temp.db')   
        return {i : gene_hierarchy_dict[i][0] for i in genes_order
                if i in gene_hierarchy_dict.keys() and 
                gene_hierarchy_dict[i][0]['overlaps'] == True} 
                
        
    def parse_exonerate_output(self,
                               ex_out_fname,
                               ex_format,
                               target_genes_loc,
                               genes_order):
        hits_feature_intervals = dict()
        overlaping_queries = dict()
        # our queries are the prot seqs
        for query in SearchIO.parse(ex_out_fname, 
                                    ex_format):
            query_dict = dict()
            query_id = query.id.rsplit('.')[0]
            if query_id not in hits_feature_intervals:
                hits_feature_intervals[query_id] = []
            hits_dict = dict()
            # each query can have many hits
            for hit_idx, hit in enumerate(query):
                hsp_dict = dict()
                # and each hit can have many HSPs
                # HSP stands for High-scoring Segment Pair
                for hsp_idx, hsp in enumerate(hit): 
                    query_ranges = [frag.query_range 
                                    for frag in hsp.fragments]
                    hit_ranges = [frag.hit_range
                                  for frag in hsp.fragments]
                    hsp_dict[hsp_idx] = { 'hit':
                                         [P.open(i[0], i[1]) 
                                          for i in hit_ranges],
                                         'query': 
                                         [P.open(i[0], i[1]) 
                                          for i in query_ranges]}
                    gene_loc = P.closed(target_genes_loc[query_id].lower - 20, 
                                        target_genes_loc[query_id].upper + 20)
                    if all([gene_loc.contains(i) for i in hsp_dict[hsp_idx]['hit']]):
                        if query_id not in overlaping_queries:
                            overlaping_queries[query_id] = copy.deepcopy(
                                hsp_dict[hsp_idx]
                            )
                        else:
                            print('extra overlaping hit:',
                                  query_id,
                                  hsp_dict[hsp_idx]['hit'])
                hits_dict[hit_idx] = hsp_dict
            hits_feature_intervals[query_id].append(hits_dict)
        return {i:overlaping_queries[i] for i in genes_order 
                if i in overlaping_queries}
    
    
    def get_utr_coords(self,
                       gs,
                       ge,
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
                   overlaping_hit_dict,
                   mutations=False,
                   genome_seq = '',
                   mut_loc='',
                   exon_dup=0,
                   genome_seq_with_mut=''):
        features = []
        if mutations:
            if mut_loc == 'before':
                mut = 0
            if mut_loc == 'after':
                mut = 1
            for exon_id, exon_coords in exons_dict[gene_id]['coding_exons_dups'].items():
                for coord in exon_coords:
                    s, e = coord
                    color = "#ffcccc"
                    labl = 'exon ' + str(idx)
                    if idx == (exon_dup + mut):
                        color = 'green'
                        labl = 'exon ' + str(idx - mut) + ' dup'
                    if idx > (exon_dup + mut):
                        idx_ = idx - 1
                        labl = 'exon ' + str(idx_)
                    idx += 1
                    features.append(GraphicFeature(start=s,
                                                   end=e,
                                                   color=color,
                                                   label=labl))
        else:
            for exon_id, exon_coords in exons_dict[gene_id]['coding_exons'].items():
                for coord in exon_coords:
                    s, e = coord
                    color = "#ffcccc"
                    labl = 'exon ' + str(idx)
                    idx += 1
                    features.append(GraphicFeature(start=s,
                                                   end=e,
                                                   color=color,
                                                   label=labl))
        features += self.get_utr_coords(gs,
                                        ge,
                                        gene_breakdown,
                                        genome_seq=genome_seq, 
                                        mutation=mutations,
                                        dup_genome_seq=genome_seq_with_mut)

        for idx_j, query_feat_coord in enumerate(overlaping_hit_dict['hit'], 1):
            qs, qe = query_feat_coord.lower, query_feat_coord.upper
            labl = 'exon ' + str(idx_j)
            features.append(GraphicFeature(start=qs,
                                           end=qe,
                                           color="blue",
                                           label=labl))
        record = GraphicRecord(first_index=gs, 
                               sequence_length=ge-gs,
                               features=features)
        return record


    def visualize_seq_alignments(self,
                                 overlaping_queries_list,
                                 genes_locations_list,
                                 sim_genes_bkdown_list,
                                 mut_loc_list,
                                 genomes_list,
                                 exon_dup_n,
                                 id_gene_dict,
                                 exons_intervals,
                                 fig_filename):
        fig, ax = plt.subplots(len(max(overlaping_queries_list,
                                       key=len)),
                               3,
                               figsize=(13, 32))
        plt.tight_layout()
        genome_seq = genomes_list[0]
        queries_ref = overlaping_queries_list[0].keys()
        for col_idx, (overlaping_queries,
                      genes_location,
                      genes_bkdown,
                      mut_loc,
                      genome,
                      exon_intv) in enumerate(zip(overlaping_queries_list,
                                                  genes_locations_list,
                                                  sim_genes_bkdown_list,
                                                  mut_loc_list,
                                                  genomes_list,
                                                  exons_intervals)):
            overlaping_queries = {i: overlaping_queries[i] for i in queries_ref
                                  if i in overlaping_queries
                                 }
            muts = True
            genome_seq_w_mut = genome
            if mut_loc == '':
                muts = False
                genome_seq_w_mut = ''
            for idx_, (ID, overlaping_hit_dict) in enumerate(overlaping_queries.items()):
                coords_gene = genes_location[ID]
                bkdown_dict = self.flatten_dict(genes_bkdown)[ID]
                chm = bkdown_dict['chrom']
                gs, ge = coords_gene.lower, coords_gene.upper
                features = []
                idx = 1
                record = self.get_record(idx,
                                         ID,
                                         exon_intv[chm],
                                         gs,
                                         ge, 
                                         bkdown_dict,
                                         overlaping_hit_dict,
                                         mutations=muts,
                                         genome_seq = genome_seq,
                                         mut_loc=mut_loc,
                                         exon_dup=exon_dup_n,
                                         genome_seq_with_mut=genome_seq_w_mut)
                ax[idx_, col_idx].set_title(ID + ':' + id_gene_dict[ID], loc='center')
                record.plot(ax=ax[idx_, col_idx],
                            figure_width=5,
                            with_ruler=False)
