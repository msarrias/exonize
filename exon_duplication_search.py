from exon_analysis import *
# for multithreading the search of exon duplicates
from multiprocessing.pool import ThreadPool
import time
import random
# for writing exonerate input/output as temporal files
import tempfile 
from tqdm import tqdm # progress bar
import sqlite3

class ExonDupSearch(ExonAnalysis):
    def __init__(self,
                 specie,
                 db_path, 
                 gff_file_path,
                 genome_file_path,
                 gene_hierarchy_path,
                 model = 'coding2genome',
                 percentage = 70,
                 sleep_max_seconds = 5,
                 min_exon_length = 50,
                 cutoff = 0.7,
                 verbose = True):
        ExonAnalysis.__init__(self,
                              db_path,
                              gff_file_path,
                              gene_hierarchy_path,
                              verbose
                             )
        self.specie = specie
        self.timeout_db = 160.0
        self.stop_codons = ["TAG","TGA","TAA"]
        self.results_df = f'{self.specie}_results.db'
        self.verbose = verbose
        self.model = model
        self.percent = percentage
        self.secs = sleep_max_seconds
        self.min_len = min_exon_length
        self.cutoff = cutoff
        self.genome_file_path = genome_file_path
        self.most_inclusive_transcript_dict = dict()
        self.exons_across_trascpts_dict = dict()
        
        
    def generate_data_for_analysis(self):
        self.generate_genome_info_for_analysis()
        try:
            if self.verbose: print(f"Reading genome:", end = " ")
            self.read_genome(self.genome_file_path)
            if self.verbose: print("Done!")
        except ValueError:
            print("Wrong genome path")
        if self.verbose: print("Checking overlaps:", end = " ")
        self.check_overlaps()
        if self.verbose:
            print("Done!")
            print(f"Computing genome-wide coding exons coordinates:",end = " ")
        self.get_gene_hierarchy_dict_with_coding_exons()
        if self.verbose: 
            print("Done!")
            print(f"Generate most inclusive exon set across transcripts:",end = " ")
        self.generate_most_inclusive_transcript_dict()
        if self.verbose:
            print("Done!")   
            
            
    def get_summary_statistics(self):
        # raw annotations - gff file features count
        self.raw_annot_dict = {annot_type: len(list(self.db.features_of_type(annot_type)))
                               for annot_type in self.db_features}
        self.filtered_annot_dict = copy.deepcopy(self.raw_annot_dict)
        # filtered coding annotations count  
        self.filtered_annot_dict['gene'] = len(self.gene_hierarchy_dict_with_coding_exons)
        self.filtered_annot_dict['mRNA'] = sum([len(gene_dict)
                                                for gene, gene_dict
                                                in self.gene_hierarchy_dict.items()
                                               ])
        self.filtered_annot_dict['intron'] = sum([1 for i, j in self.most_inclusive_transcript_dict.items() 
                                                   for k, l in j['transcpt'].items()
                                                   if l['type'] == 'intron' 
                                                 ])
        self.filtered_annot_dict['exon'] = sum([len(j['ce_dict']) 
                                                 for i, j in self.most_inclusive_transcript_dict.items()
                                               ])
        #count coding regions only once
        CDS_count = 0
        for gene_id, gene_dict in self.gene_hierarchy_dict.items():
            temp_dict = {transcpt:{i['coord']:{'type':i['type']} for i in list_annot}
                         for transcpt, list_annot in gene_dict.items()}
            if temp_dict:
                sorted_cds = self.sort_interval_dict(self.get_annotations(temp_dict, 'CDS'))
                temp = self.get_annotations_across_transcripts(sorted_cds)
                CDS_count+= len(temp)
        self.filtered_annot_dict['CDS'] = CDS_count
        #create df
        df = pd.DataFrame({"Raw": list(self.raw_annot_dict.values()),
                           'Filtered': list(self.filtered_annot_dict.values())},
                          index=self.db_features) 
        df["Raw"] = df["Raw"].map("{:,}".format)
        df["Filtered"] = df["Filtered"].map("{:,}".format)
        return df 
    
    
    @staticmethod
    def batch(iterable, n=1):
        l = len(iterable)
        for ndx in range(0, l, n):
            yield iterable[ndx:min(ndx + n, l)]
        
        
    def get_fragment_annotation(self, frag):
        """
        Collect attributes of a single HSPFragment
        """
        return {'query_coord': P.open(frag.query_start, frag.query_end),
                'hit_coord':P.open(frag.hit_start, frag.hit_end),
                'query_annotation': frag.aln_annotation['query_annotation'],
                'hit_annotation':frag.aln_annotation['hit_annotation'],
                'query_seq':str(frag.query.seq),
                'hit_seq':str(frag.hit.seq)}


    def get_HSP_annotations(self, HSP):
        """ 
        get_HSP_annotations returns a dictionary with the
        fragments annotation for a High-scoring Segment Pair
        :param HSP: object of the class HSP
        """
        fragments_dict = {}
        for frag_idx, frag in enumerate(HSP.fragments):
            fragments_dict[frag_idx] = self.get_fragment_annotation(frag)
        return fragments_dict

    
    def get_query_and_hits_seqs(self, 
                                annot_id,
                                gene_id,
                                coding_exons_introns_dict,
                                exon_coord):
        """
        get_query_and_hits_seqs is a function given a gene and exon id, and
        a features dictionary with coordinates, returns the dna sequences.
        :param annot_id: as in the GFF file.
        :param gene_id: as in the GFF file.
        :param coding_exons_introns_dict: target annotations.
        param exon_coord: coordinates of the hit annotation.
        :return: the query sequence and a dictionary with the target sequences
        """
        chrm = self.chrom_dict[gene_id]
        strand = self.strand_dict[gene_id]
        query_seq = self.genome[chrm][strand][exon_coord.lower:exon_coord.upper]
        
        hits_seqs = {feat_dict['id']: 
                     self.genome[chrm][strand][coord.lower:coord.upper]
                     for coord, feat_dict in coding_exons_introns_dict.items() 
                     if (feat_dict['id'] != annot_id 
                         and (coord.upper - coord.lower) >= round(0.9*(exon_coord.upper-exon_coord.lower)))}
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
        to the `get_annotations_across_transcripts` function, where after
        resolving a  not allowed intersection on `interval_dict`  
        the `complete_dict` function updates the new interval dictionary
        that will be passed on the next recursive step of 
        `get_annotations_across_transcripts`
        """
        intervals_list = list(interval_dict.keys())
        new_gene_interv_dict_copy = copy.deepcopy(new_gene_interv_dict)
        for intrv in intervals_list[int_idx+step:]:
            if not intrv in new_gene_interv_dict_copy:
                new_gene_interv_dict_copy[intrv] = copy.deepcopy(interval_dict[intrv])
        return new_gene_interv_dict_copy

    
    @staticmethod
    def get_annotations(coords_dict, type_annot = "coding_exon"):
        """
        get_coding_exons_annotations returns a dictionary
        with annotations of the type "coding_exon".
        :param coords_dict: gene annotations dict.
        """
        return {i:y for key, v in coords_dict.items()
                for i, y in v.items() if type_annot in y['type']}

    
    @staticmethod
    def intervals_proportion(a,b):
        """Computes the proportion between two intervals"""
        return (a.upper - a.lower)/(b.upper - b.lower)

    
    @staticmethod
    def get_small_large_interv(a,b):
        """
        Given two intervals, the function 
        get_small_large_interv returns the smaller 
        and the larger interval in length.
        """
        len_a = a.upper - a.lower
        len_b = b.upper - b.lower
        if len_a < len_b:
            return a, b
        else:
            return b, a

        
    def get_annotations_across_transcripts(self, interval_dict):
        """
        `get_annotations_across_transcripts` is a recursive function that 
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
        `get_annotations_transcripts` resolves what can be categorazied in 
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
                    return self.get_annotations_across_transcripts(self.sort_interval_dict(new_gene_interv_dict))
                # Case 2: exon_i intersects exon_j but exon_i is not contained in exon_j 
                else:
                    intersection = large_itv & small_itv
                    perc_long = self.intervals_proportion(intersection, large_itv)
                    # Case 2.1: If the proportion is more than 70% of the longest exon
                    if round(perc_long) > self.cutoff:
                        new_gene_interv_dict[large_itv] = interval_dict[large_itv]
                        new_gene_interv_dict = self.complete_dict(int_idx,new_gene_interv_dict,interval_dict)
                        return self.get_annotations_across_transcripts(self.sort_interval_dict(new_gene_interv_dict))
                    # Case 2.2: If the proportion is less than 70% of the longest exon
                    # we allow the overlap
                    else:
                        new_gene_interv_dict[interval] = copy.deepcopy(interval_dict[interval])
                        exclude_overlap_dict = {}
                        for intrv in intervals_list[int_idx+1:]:
                            if not intrv in new_gene_interv_dict:
                                exclude_overlap_dict[intrv] = copy.deepcopy(interval_dict[intrv])
                        new_block = self.get_annotations_across_transcripts(self.sort_interval_dict(exclude_overlap_dict))
                        new_gene_interv_dict.update(new_block)
                        return new_gene_interv_dict

                    
    @staticmethod                
    def trans_exons_and_introns_coords(intv_dict, gene_loc):
        """
        trans_exons_and_introns_coords is a function that for a given gene 
        annotations dictionary, it creates intron annotations on the 
        non-annotated regions.
        :param intv_dict: gene annotations breakdown.
        :param gene_id: gene identifier.
        :return: new dictionary with the new and sorted annotations.
        """
        new_dict, i = {}, 0
        keys_list = list(intv_dict.keys())
        first_itv, last_itv = keys_list[0], keys_list[-1]
        s, e = first_itv.lower, first_itv.upper
        gc_s, gc_e = gene_loc.lower, gene_loc.upper
        if P.open(gc_s, s-1):
            new_dict[P.open(gc_s, s-1)] = {'id': f'intron_{i}', 'type': 'intron'}
            i +=1
        new_dict[first_itv] = intv_dict[first_itv]
        prev_intv = first_itv
        for intv_idx, intv in enumerate(keys_list[1:]):
            if  P.open(prev_intv.upper+1, intv.lower-1):
                new_dict[P.open(prev_intv.upper+1, intv.lower-1)] = {'id': f'intron_{i}', 'type': 'intron'}
                i +=1
            new_dict[intv] = intv_dict[intv]
            prev_intv = intv
        if P.open(last_itv.upper+1, gc_e):
            new_dict[P.open(last_itv.upper+1, gc_e)] = {'id': f'intron_{i}', 'type': 'intron'}
        return new_dict

        
    def generate_most_inclusive_transcript_dict(self):
        for gene_id, gene_dict in self.gene_hierarchy_dict_with_coding_exons.items():
            if gene_dict:
                sorted_exons = self.sort_interval_dict(self.get_annotations(gene_dict))
                coding_exons_dict = self.get_annotations_across_transcripts(sorted_exons)
                if coding_exons_dict:
                    temp = self.trans_exons_and_introns_coords(coding_exons_dict, self.gene_loc[gene_id])
                    self.most_inclusive_transcript_dict[gene_id] = {'trans_exons':sorted_exons,
                                                                    'transcpt': temp, 
                                                                    'ce_dict': coding_exons_dict }
                    
                    
    def parse_exonerate_output_only_higher_score_hits(self, ex_out_fname):
        """
        `parse_exonerate_output_only_higher_score_hits` parses single query
        exonerate-text outputs for coding2genome model.
        Each query is composed by a hit(s), each hit by HSP(s), and each HSP
        by fragments (since ungappend). Some criteria:
        - In case we have multiple HSps within a single hit, we keep the one
        with the highest score. 
        - If we have repeated hits on the same sequence, we keep the hit with 
        the highest HSP score.
        """
        #we take the first element of the list since there's only one query 
        all_results = list(SearchIO.parse(ex_out_fname, 'exonerate-text'))
        if all_results:
            hit_dict = {}
            for hit_idx, hit in enumerate(all_results[0]):
                hsp_dict = {'score': 0, 'frag_dict': {}}
                for idx_hsp, HSP in enumerate(hit):
                    #We keep the HSP with the highest score
                    if (not hsp_dict['frag_dict'] or HSP.score > hsp_dict['score']):
                        try:
                            hsp_dict = {'score':HSP.score, 'frag_dict': self.get_HSP_annotations(HSP)} 
                        except:
                            print('not hsp', hit.id, hsp_dict)
                            pass
                # if there's more than one hit on the same target seq, 
                # we keep the one with the highest score
                if (hit.id not in hit_dict or HSP.score > hit_dict[hit.id]['score']):
                    hit_dict[hit.id] = hsp_dict
            return hit_dict
        return {}
    
    
    def run_exonerate(self, annotation_id, query_seq, hits_seqs):
        temp = {}
        with tempfile.TemporaryDirectory() as tmpdirname:
            self.dump_fasta_file(f'{tmpdirname}/query.fa', {annotation_id:query_seq})
            self.dump_fasta_file(f'{tmpdirname}/target.fa', hits_seqs)
            output_file = f'{tmpdirname}/output.txt'
            exonerate_command = ['exonerate',
                                     '--showalignment', 'yes',
                                     '--showvulgar', 'no',
                                     '--model', self.model,
                                     '--percent', str(self.percent),
                                     '--query', 'query.fa',
                                     '--target', 'target.fa']
            with open(output_file, 'w') as out:
                subprocess.run(exonerate_command, cwd=tmpdirname, stdout=out)
                temp = self.parse_exonerate_output_only_higher_score_hits(output_file)
        return temp
           

    def get_duplicates(self, gene_id):
        """
        get_duplicates is a function that given a gene_id, and a set of constraints it
        performs a search of exon duplications. The constraints on the exon to query are
        the following:
        1. The exon must be larger than a minimum length, default is 50 bps.
        2. If exon_i has a hit on exon_j we exclude the exon_j search.
        :param gene_id: gene in question.
        :return temp_ce_dict: hits dictionary.
        """
        time.sleep(random.randrange(0, self.secs))
        temp_ce_dict = dict()
        # this gives a list of "distinct" coding exons across transcripts, i.e., no duplicates
        coding_exons_dict = copy.deepcopy(self.most_inclusive_transcript_dict[gene_id]['ce_dict'])
        # this creates intron annotations for the regions in between exons.
        coding_exons_with_introns = copy.deepcopy(self.most_inclusive_transcript_dict[gene_id]['transcpt'])
        # if exon x is a duplicate of exon y, we skip looking for 
        # x duplicates, since they will be comprised in the y hits
        introns_list = []
        for coord, annotation_dict in coding_exons_dict.items():
            annotation_id = annotation_dict['id']
            if (coord.upper - coord.lower) > self.min_len:
                query_seq, hits_seqs = self.get_query_and_hits_seqs(annotation_id, gene_id, coding_exons_with_introns, coord)
                if hits_seqs:
                    temp = self.run_exonerate(annotation_id, query_seq, hits_seqs)
                    if temp:
                        temp_ce_dict[annotation_id] = temp
        if temp_ce_dict:
            self.insert_fragments_table(gene_id, temp_ce_dict)
        else:
            gene_arg_tuple = (gene_id,
                              self.chrom_dict[gene_id],
                              self.strand_dict[gene_id],
                              self.gene_loc[gene_id].lower,
                              self.gene_loc[gene_id].upper,
                              0)
            self.insert_gene_ids_table(gene_arg_tuple)
            
            
    def connect_create_results_db(self):
        db = sqlite3.connect(self.results_df, timeout=self.timeout_db) 
        c = db.cursor()
        c.execute('''
                  CREATE TABLE IF NOT EXISTS GeneIDs (
                  gene_id INTEGER PRIMARY KEY AUTOINCREMENT,
                  gene_name VARCHAR(100) NOT NULL,
                  gene_chrom VARCHAR(100) NOT NULL,
                  gene_strand VARCHAR(1) NOT NULL,
                  gene_start INTEGER NOT NULL,
                  gene_end INTEGER NOT NULL,
                  has_duplicated_exon BINARY(1) NOT NULL,
                  UNIQUE(gene_name)
                  )
                  ''')
        c.execute('''
                  CREATE TABLE IF NOT EXISTS Fragments (
                  fragment_id INTEGER PRIMARY KEY AUTOINCREMENT,
                  gene_id  VARCHAR(100) NOT NULL REFERENCES GeneIDs(gene_id),
                  query_id VARCHAR(100) NOT NULL,
                  query_start INTEGER NOT NULL,
                  query_end INTEGER NOT NULL,
                  target_id VARCHAR(100) NOT NULL,
                  target_start INTEGER NOT NULL,
                  target_end INTEGER NOT NULL,
                  category BINARY(1) NOT NULL,
                  is_pseudo_exon BINARY(1) NOT NULL,
                  score INTEGER NOT NULL,
                  frag_query_start INTEGER NOT NULL,
                  frag_query_end INTEGER NOT NULL,
                  frag_target_start INTEGER NOT NULL,
                  frag_target_end INTEGER NOT NULL,
                  query_dna_sequence VARCHAR NOT NULL,
                  target_dna_sequence VARCHAR NOT NULL,
                  query_prot_sequence VARCHAR NOT NULL,
                  target_prot_sequence VARCHAR NOT NULL)
                  ''')        
        db.commit()
        db.close()
    

    def insert_gene_ids_table(self, gene_args_tuple):
        db = sqlite3.connect(self.results_df, timeout=self.timeout_db) 
        cursor = db.cursor()
        insert_gene_table_param = """  
        INSERT INTO GeneIDs 
        (gene_name,
        gene_chrom,
        gene_strand,
        gene_start, 
        gene_end,  
        has_duplicated_exon) 
        VALUES (?, ?, ?, ?, ?, ?)
        """
        cursor.execute(insert_gene_table_param, gene_args_tuple)
        db.commit()
        db.close()
   
    
    @staticmethod
    def is_novel(target_id):
        dup_type = 0
        if 'intron' in target_id:
            dup_type = 1
        return dup_type
        
           
    def is_pseudo_exon(self, fragment):
        is_pseudo_exon = 0
        if sum([True for stop_codon in self.stop_codons
                if stop_codon in fragment]) > 1:
            is_pseudo_exon = 1
        return is_pseudo_exon
        
    def get_fragment_annotation_coord(self, gene_id, annotation_id):
        return [i for i, j in self.most_inclusive_transcript_dict[gene_id]['transcpt'].items()
                if j['id'] == annotation_id][0]
        
        
    def insert_fragments_table(self, gene_id, temp_ce_dict):
        tuple_list = []
        for query_id, gene_query_dict in temp_ce_dict.items():
            annot_coord = self.get_fragment_annotation_coord(gene_id, query_id)
            for target_id, target_dict in gene_query_dict.items():
                for fragment_id, fragment_dict in target_dict['frag_dict'].items():
                    fragment_coord = self.get_fragment_annotation_coord(gene_id, target_id)
                    tuple_list.append(
                        (gene_id,                                        #gene_id
                         query_id,                                       #query_id
                         annot_coord.lower,                              #query_start
                         annot_coord.upper,                              #query_end
                         target_id,                                      #target_id
                         fragment_coord.lower,                           #target_start
                         fragment_coord.upper,                           #target_end       
                         self.is_novel(target_id),                       #category
                         self.is_pseudo_exon(fragment_dict['hit_seq']),  #is_pseudo_exon
                         target_dict['score'],                           #score
                         fragment_dict['query_coord'].lower,             #query_start
                         fragment_dict['query_coord'].upper,             #query_end
                         fragment_dict['hit_coord'].lower,               #hit_start
                         fragment_dict['hit_coord'].upper,               #hit_end
                         fragment_dict['query_seq'],                     #query_dna_sequence
                         fragment_dict['hit_seq'],                       #hit_dna_sequence
                         fragment_dict['query_annotation'],              #query_prot_sequence
                         fragment_dict['hit_annotation']                 #hit_prot_sequence
                        )
                    )
        insert_fragments_table_param = """
        INSERT INTO Fragments 
        (gene_id,
        query_id, 
        query_start,
        query_end,
        target_id,
        target_start,
        target_end,
        category,
        is_pseudo_exon,
        score,
        frag_query_start,
        frag_query_end,
        frag_target_start,
        frag_target_end,
        query_dna_sequence,
        target_dna_sequence,
        query_prot_sequence,
        target_prot_sequence)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """
        insert_gene_table_param = """  
        INSERT INTO GeneIDs 
        (gene_name,
        gene_chrom,
        gene_strand,
        gene_start, 
        gene_end,  
        has_duplicated_exon) 
        VALUES (?, ?, ?, ?, ?, ?)
        """
        gene_arg_tuple = (gene_id,
                          self.chrom_dict[gene_id],
                          self.strand_dict[gene_id],
                          self.gene_loc[gene_id].lower,
                          self.gene_loc[gene_id].upper,
                          1)
        db = sqlite3.connect(self.results_df, timeout=self.timeout_db) 
        cursor = db.cursor()
        cursor.execute(insert_gene_table_param, gene_arg_tuple)
        cursor.executemany(insert_fragments_table_param, tuple_list)
        db.commit()
        db.close()   

            
    def query_gene_ids_in_res_db(self):
        db = sqlite3.connect(self.results_df, timeout=self.timeout_db) 
        cursor = db.cursor()
        cursor.execute("SELECT gene_name FROM GeneIDs")
        rows = cursor.fetchall()
        db.close()
        return [i[0] for i in rows]

        
    def search_for_exon_duplicates(self, 
                                   args_list,
                                   batch_n,
                                   threads=6):
        """
        search_for_exon_duplicates is a function that performs a 
        genome-wide search of exons duplications, multithreading 
        is performed in barches of genes.
        :param args_list: list of gene ids.
        :param batch_n: batches size
        :param threads: number of threads WARNING: Mind that a large 
        selection of threads number can yield to a SQLite3.OperationalError
        since many processes trying to access the db will potentially exhaust
        the db timeout. Hence, choose a large enough value that will not 
        result in locking the db.
        """
        self.connect_create_results_db()
        processed_gene_ids = self.query_gene_ids_in_res_db()
        if processed_gene_ids:
            args_list = [i for i in args_list if i not in processed_gene_ids] 
        if args_list:
            batches_list = [i for i in self.batch(args_list, batch_n)]
            tic = time.time()
            with tqdm(total=len(args_list)) as progress_bar:
                for arg_batch in batches_list:              
                    t = ThreadPool(processes=threads)
                    t.map(self.get_duplicates, arg_batch)
                    t.close()
                    t.join()
                    progress_bar.update(len(arg_batch))
            tac = time.time()
            print(f'process took: {(tac - tic)/60} minutes')

        ## alternative paralellizing - this should be done 
        ## outside the class
#         pool = Pool(processes=10)
#         mapped_values = list(tqdm(pool.imap_unordered(self.get_duplicates, args_list), total=len(args_list)))
#         pool.close()
#         return mapped_values