# for multithreading the search of exon duplicates
from multiprocessing.pool import ThreadPool
import time, random
from tqdm import tqdm # progress bar
# for writing exonerate input/output as temporal files
import tempfile 
import pickle, copy
import subprocess # for running Exonerate
from Bio import SeqIO, SearchIO # for parsing exonerate output
import portion as P # for working with intervals


class ExonDupSearch:
    def __init__(self,
                 ExonAnalysis_obj, 
                 model = 'coding2genome',
                 percentage = 70,
                 sleep_max_seconds = 5,
                 min_exon_length = 50):
        self.exon_analysis_obj = ExonAnalysis_obj
        self.model = model
        self.percent = percentage
        self.secs = sleep_max_seconds
        self.min_len = min_exon_length

        
    @staticmethod
    def batch(iterable, n=1):
        l = len(iterable)
        for ndx in range(0, l, n):
            yield iterable[ndx:min(ndx + n, l)]
        
        
    @staticmethod
    def parse_exonerate_output_only_higher_score_hits(ex_out_fname):
        hits_feature_intervals = dict()
        for query in SearchIO.parse(ex_out_fname, 
                                    'exonerate-vulgar'):
            query_dict = dict()
            query_id = query.id #.rsplit('.')[0]
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
                    hsp_dict[hsp_idx] = {'query_id':query_id,
                                           'hit_id':hit.id,
                                           'score':hsp.score,
                                           'hit': [P.open(i[0], i[1])
                                                   for i in hit_ranges],
                                           'query': [P.open(i[0], i[1])
                                                     for i in query_ranges]}
                    #let's choose the hsp with higher score 
                    if len(set([value['hit_id'] for key, value in hsp_dict.items()])) < 2:
                        max_score = max([(key, value['score'])
                                         for key, value in hsp_dict.items()
                                        ],
                                        key=lambda item:item[1]
                                       )[0]
                        hits_dict[hit_idx] = hsp_dict[max_score]
                    else:
                        hits_dict[hit_idx] = hsp_dict

            hits_feature_intervals[query_id].append(hits_dict)
        if len(hits_feature_intervals) == 1:
            return hits_feature_intervals[query_id]   
        return hits_feature_intervals
    
    
    def get_query_and_hits_seqs(self, 
                                exon_id,
                                gene_id,
                                coding_exons_introns_dict,
                                exon_coord):
        chrm = self.exon_analysis_obj.chrom_dict[gene_id]
        strand = self.exon_analysis_obj.strand_dict[gene_id]
        query_seq = self.exon_analysis_obj.genome[chrm][strand][exon_coord.lower:exon_coord.upper]
        
        hits_seqs = {feat_dict['id']: 
                     self.exon_analysis_obj.genome[chrm][strand][coord.lower:coord.upper]
                     for coord, feat_dict in coding_exons_introns_dict.items() 
                     if feat_dict['id'] != exon_id and (coord.upper - coord.lower) > 1
                    }
        return query_seq ,hits_seqs
        
        
    @staticmethod    
    def sort_interval_dict(interval_dict):
        """
        Sort intervals dictionary
        the intervals objects are created with the portion library
        https://github.com/AlexandreDecan/portion
        """
        keys_ = list(interval_dict.keys())
        keys_ = sorted(keys_, key = lambda item: (item.lower, item.upper))
        return {key: copy.deepcopy(interval_dict[key]) for key in keys_}

    
    @staticmethod
    def get_overlapping_dict(interval_dict):
        """
        gets the next overlap to the right
        for each interval in the interval
        dictionary keys `interval_dict`
        """
        copy_interval_dict = copy.deepcopy(interval_dict)
        overlapping_dict = {intv:[] for intv in copy_interval_dict}
        intervals_list = list(copy_interval_dict.keys())
        for idx, (feat_interv, feature_annot) in enumerate(copy_interval_dict.items()):
            if (idx != (len(copy_interval_dict) - 1) 
                and feat_interv.overlaps(intervals_list[idx+1])
               ):
                overlapping_dict[feat_interv] = intervals_list[idx+1]
        return overlapping_dict

    
    @staticmethod
    def complete_dict(int_idx, 
                      new_gene_interv_dict, 
                      interval_dict,
                      step = 2):
        """
        `complete_dict` populates the `new_gene_interv_dict` dictionary
        with items in `interval_dict`. `complete_dict` is a help function 
        to the `get_coding_exons_across_transcripts` function, where after
        resolving a  not allowed intersection on `interval_dict`  
        the `complete_dict` function updates the new interval dictionary
        that will be passed on the next recursive step of 
        `get_coding_exons_across_transcripts`
        """
        intervals_list = list(interval_dict.keys())
        new_gene_interv_dict_copy = copy.deepcopy(new_gene_interv_dict)
        for intrv in intervals_list[int_idx+step:]:
            if not intrv in new_gene_interv_dict_copy:
                new_gene_interv_dict_copy[intrv] = copy.deepcopy(interval_dict[intrv])
        return new_gene_interv_dict_copy

    
    @staticmethod
    def get_coding_exons_annotations(coords_dict):
        return {i:y for key, v in coords_dict.items()
                for i, y in v.items() if 'coding_exon' in y['type'] 
               }

    
    @staticmethod
    def get_coverage_perc(a,b):
        return (a.upper - a.lower)/(b.upper - b.lower)

    
    @staticmethod
    def get_small_large_interv(a,b):
        len_a = a.upper - a.lower
        len_b = b.upper - b.lower
        if len_a < len_b:
            return a, b
        else:
            return b, a

        
    def get_coding_exons_across_transcripts(self, interval_dict):
        """
        `get_coding_exons_across_transcripts` is a recursive function that 
        takes as an input `interval_dict`: a nested dictionary of a 
        gene, its transcripts and a breakdown of each transcript's architecture. 
        input e.g.,
        {geneID:{transcript1:{exon1:{}, intron1{},...}, transcript2:{...}...}}
        and following certain criteria, it returns a sorted dictionary
        `new_gene_interv_dict` of the most inclusive set of coding exons:
        return e.g.,
        {exonID: feature_dictionary, ...}
        Cases of intersections are the following:
        - Case 1: exon_i is contained in exon_j 
        - Case 2: exon_i intersects exon_j but exon_i is not contained in exon_j 
            - Case 2.1: exon_i intersects exon_j in at least .7 of the length of each exon.
            - Case 2.2: opposite to Case2.1
        `get_coding_exons_across_transcripts` resolves what can me categorazied in 
        allowed and not allowed intersections:
        Not allowed: Case 1, Case 2.1, Allowed: Case 2.2
        Each overlap is handled at a recursion step, so there are many 
        recursion steps as overlaps.
        """
        overlapping_dict = self.get_overlapping_dict(interval_dict)
        #if sorted interval dictionary has no overlaps - nothing to do
        if not any(list(overlapping_dict.values())):
            return self.sort_interval_dict(interval_dict)
        new_gene_interv_dict = {}
        intervals_list = list(interval_dict.keys())
        for int_idx, (interval, feature_descrip) in enumerate(interval_dict.items()):
            overlap_interval = copy.deepcopy(overlapping_dict[interval])
            # if interval has no overlaps
            if not overlap_interval:
                new_gene_interv_dict[interval] = feature_descrip
            # resolve interval overlaps:
            else:
                small_itv, large_itv = self.get_small_large_interv(interval, overlap_interval)
                # Case 1: exon_i is contained in exon_j but exon_i is not contained in exon_j 
                if large_itv.contains(small_itv):
                    new_gene_interv_dict[large_itv] = interval_dict[large_itv]
                    new_gene_interv_dict = self.complete_dict(int_idx, new_gene_interv_dict, interval_dict)
                    return self.get_coding_exons_across_transcripts(self.sort_interval_dict(new_gene_interv_dict))
                # Case 2: exon_i intersects exon_j but exon_i is not contained in exon_j 
                else:
                    intersection = large_itv & small_itv
                    perc_long = self.get_coverage_perc(intersection, large_itv)
                    # Case 2.1: If the intersection covers more than 70% of the longest exon
                    if round(perc_long) > 0.7:
                        new_gene_interv_dict[large_itv] = interval_dict[large_itv]
                        new_gene_interv_dict = self.complete_dict(int_idx,new_gene_interv_dict,interval_dict)
                        return self.get_coding_exons_across_transcripts(self.sort_interval_dict(new_gene_interv_dict))
                    # Case 2.2: If the intersection covers less than 70% of the longest exon
                    # we allow the overlap
                    else:
                        new_gene_interv_dict[interval] = copy.deepcopy(interval_dict[interval])
                        exclude_overlap_dict = {}
                        for intrv in intervals_list[int_idx+1:]:
                            if not intrv in new_gene_interv_dict:
                                exclude_overlap_dict[intrv] = copy.deepcopy(interval_dict[intrv])
                        new_block = self.get_coding_exons_across_transcripts(self.sort_interval_dict(exclude_overlap_dict))
                        new_gene_interv_dict.update(new_block)
                        return new_gene_interv_dict

                    
    def trans_exons_and_introns_coords(self, intv_dict, gene_id):
        new_dict, i = {}, 0
        copy_intv_dict = copy.deepcopy(intv_dict)
        keys_list = list(intv_dict.keys())
        first_itv, last_itv = keys_list[0], keys_list[-1]
        gene_coord = self.exon_analysis_obj.gene_loc[gene_id]
        s, e = first_itv.lower, first_itv.upper
        gc_s, gc_e = gene_coord.lower, gene_coord.upper
        if P.open(gc_s, s):
            copy_intv_dict[P.open(gc_s, s)] = {'id': f'intron_{i}'}
            i+=1
        for coord in keys_list[1:]:
            if P.open(e, coord.lower):
                copy_intv_dict[P.open(e, coord.lower)] = {'id': f'intron_{i}'}
                i+=1
                e = coord.upper
        if P.open(last_itv.upper, gc_e):
            copy_intv_dict[P.open(last_itv.upper, gc_e)] = {'id': f'intron_{i}'}
        return self.sort_interval_dict(copy_intv_dict)

    
    def get_duplicates(self, gene_id):
        time.sleep(random.randrange(0, self.secs))
        temp_ce_dict = dict()
        gene_dict = copy.deepcopy(self.exon_analysis_obj.gene_hierarchy_dict_with_coding_exons[gene_id])
        if gene_dict:
            sorted_exons = self.sort_interval_dict(self.get_coding_exons_annotations(gene_dict))
            coding_exons_dict = self.get_coding_exons_across_transcripts(sorted_exons)
            
            if coding_exons_dict:
                coding_exons_with_introns = self.trans_exons_and_introns_coords(coding_exons_dict, gene_id)
                # if exon x is a duplicate of exon y, we skip looking for 
                # x duplicates, since they will be comprised in the y hits
                pass_exons = []
                for exon_coord, coding_exon_dict in coding_exons_dict.items():
                    exon_id = coding_exon_dict['id']
                    if exon_id not in pass_exons:
                        if (exon_coord.upper - exon_coord.lower) > self.min_len:
                            query_seq, hits_seqs = self.get_query_and_hits_seqs(exon_id,
                                                                                gene_id,
                                                                                coding_exons_with_introns,
                                                                                exon_coord)
                            if hits_seqs:
                                with tempfile.TemporaryDirectory() as tmpdirname:
                                    self.exon_analysis_obj.dump_fasta_file(f'{tmpdirname}/query.fa',
                                                                           {exon_id:query_seq})
                                    self.exon_analysis_obj.dump_fasta_file(f'{tmpdirname}/target.fa',
                                                                           hits_seqs)
                                    output_file = f'{tmpdirname}/output.txt'
                                    exonerate_command = ['exonerate',
                                                         '--showalignment', 'no',
                                                         '--model', self.model,
                                                         '--percent', str(self.percent),
                                                         '--query', 'query.fa',
                                                         '--target', 'target.fa']
                                    with open(output_file, 'w') as out:
                                        subprocess.run(exonerate_command,
                                                       cwd=tmpdirname,
                                                       stdout=out)
                                    temp = self.parse_exonerate_output_only_higher_score_hits(output_file)
                                    for hit in temp:
                                        for hit_id, hit_dict in hit.items():
                                            hit_dict['query_coord'] = {exon_coord:query_seq}
                                            hit_dict['hit_coord'] = {coord : hits_seqs[coord_dict['id']]
                                                                     for coord, coord_dict in coding_exons_with_introns.items()
                                                                     if coord_dict['id'] != hit_dict['query'] 
                                                                     and coord_dict['id'] == hit_dict['hit_id']}

                                            hit_dict['hit_seqs'] = {coord: 
                                                                    hits_seqs[hit_dict['hit_id']][coord.lower:coord.upper]
                                                                    for coord in hit_dict['hit']
                                                                   }
                                            hit_dict['query_seqs'] = {coord: 
                                                                      query_seq[coord.lower:coord.upper]
                                                                      for coord in hit_dict['query']
                                                                     }
                                            del hit_dict['query']
                                            del hit_dict['hit']
                                            if 'exon' in hit_dict['hit_id']:
                                                pass_exons.append(hit_dict['hit_id'])
                                    if temp:
                                        temp_ce_dict[exon_id] = temp
        return temp_ce_dict

    
    def search_for_exon_duplicates(self, 
                                   args_list,
                                   batch_n,
                                   threads=10):
        res = []
        tic = time.time()
        batches_list = [i for i in self.batch(args_list, batch_n)]
        with tqdm(total=len(args_list)) as progress_bar:
            for arg_batch in batches_list:              
                t = ThreadPool(processes=threads)
                res += list(t.map(self.get_duplicates, arg_batch))
                t.close()
                t.join()
                progress_bar.update(len(arg_batch))
        tac = time.time()
        res = [i for i in res if i]
        res = {next(iter(exon_hit_dict)).rsplit('.')[0]:exon_hit_dict for exon_hit_dict in res}
        print(f'process took: {(tac - tic)/60} minutes')
        self.exon_analysis_obj.dump_pkl_file(f'exon_dups_analysis{str(time.time())}.pkl', res)
        return res
        ## alternative paralellizing - this should be done 
        ## outside the class
#         pool = Pool(processes=10)
#         mapped_values = list(tqdm(pool.imap_unordered(self.get_duplicates, args_list), total=len(args_list)))
#         pool.close()
#         return mapped_values