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
from datetime import date

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
                                transcpt_dict,
                                exon_coord):
        chrm = self.exon_analysis_obj.chrom_dict[gene_id]
        strand = self.exon_analysis_obj.strand_dict[gene_id]
        query_seq = self.exon_analysis_obj.genome[chrm][strand][exon_coord.lower:exon_coord.upper]
        hits_seqs = {feat_dict['id']: 
                     self.exon_analysis_obj.genome[chrm][strand][coord.lower:coord.upper]
                     for coord, feat_dict in transcpt_dict.items() 
                     if feat_dict['id'] != exon_id
                    }
        return query_seq ,hits_seqs
        
        
    def get_duplicates(self, gene_id):
        time.sleep(random.randrange(0, self.secs))
        temp_ce_dict = dict()
        gene_dict = copy.deepcopy(self.exon_analysis_obj.gene_hierarchy_dict_with_coding_exons[gene_id])
        if gene_dict:
            #might change later
            transcript_id = next(iter(gene_dict))
            transcpt_dict = copy.deepcopy(gene_dict[transcript_id])
            coding_exons_ids = [(coord, feat_dict['id'])
                                for coord, feat_dict in transcpt_dict.items() 
                                if feat_dict['type'] == 'coding_exon'
                               ]
            if coding_exons_ids:
                # if exon x is a duplicate of exon y, we skip looking for 
                # x duplicates, since they will be comprised in the y hits
                pass_exons = []
                for coding_exon in coding_exons_ids:
                    exon_coord, exon_id = coding_exon
                    if exon_id not in pass_exons:
                        if (exon_coord.upper - exon_coord.lower) > self.min_len:
                            query_seq, hits_seqs = self.get_query_and_hits_seqs(exon_id,
                                                                                gene_id,
                                                                                transcpt_dict,
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
        print(f'process took: {(tac - tic)/60} minutes')
        self.exon_analysis_obj.dump_pkl_file(f'exon_dups_analysis{str(date.today())}.pkl', {'res':res})
        return res
        ## alternative paralellizing - this should be done 
        ## outside the class
#         pool = Pool(processes=10)
#         mapped_values = list(tqdm(pool.imap_unordered(self.get_duplicates, args_list), total=len(args_list)))
#         pool.close()
#         return mapped_values