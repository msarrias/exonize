from global_var import *
from exon_analysis import *
import random # for sampling and making sampling reproducible

class genome_constructor(exon_analysis):
    def __init__(self, 
                 db_path,
                 gene_hierarchy_path,
                 genome_path,
                 seed,
                 verbose=True):
        exon_analysis.__init__(self,
                               db_path,
                               gene_hierarchy_path,
                               gene_hierarchy = True
                              )
        self.verbose = verbose
        self.seed = seed
        self.read_genome(genome_path)
        self.chm_order = list(self.genome.keys())
        self.check_overlaps()
        self.get_gene_hierarchy_dict_with_coding_exons()
        self.intron_signal = 20
        
        
    def collect_genes(self,
                     n_coding_exons,
                     max_exon_len,
                     min_intron_len):
        """
        to make the sampling easier, we will set some criteria:
        1. genes in the positive strand
        2. genes with a single transcript
        :par n_coding_exons: number of coding exons
        :max_exon_len: CHECK ON THIS - how does it depend on the 
        neighbouring intron length
        """
        random.seed(self.seed)
        genes_sample = dict()
        cp_gene_hierarchy_dict = copy.deepcopy(
            self.gene_hierarchy_dict_with_coding_exons
        )
        for gene_id, gene_hierarchy_dict in cp_gene_hierarchy_dict.items():
            if self.db[gene_id].strand == '+' and len(gene_hierarchy_dict) == 1:
                transcript_temp = dict()
                trapt_hier_dict = copy.deepcopy(
                    gene_hierarchy_dict[next(iter(gene_hierarchy_dict.keys()))]
                )
                coding_exons = [(value['type'],
                                 key.upper - key.lower <= max_exon_len
                                ) 
                                for key, value in trapt_hier_dict.items() 
                                if value['type'] == 'coding_exon'
                               ]
                if (len(coding_exons) == n_coding_exons
                    and all([i[1] for i in coding_exons])
                   ):
                    introns_dist = [key.upper - key.lower >= min_intron_len 
                                    for key, value in trapt_hier_dict.items() 
                                    if 'intron' in value['type']]
                    if all(introns_dist):
                        genes_sample[gene_id] = copy.deepcopy(trapt_hier_dict)
        return genes_sample

    
    def get_sample(self, genes_set, sample_size):
        count = 0
        sample_dict = dict()
        random.seed(self.seed)
        while count < sample_size:
            overlaps = True
            # It might happen that the coords of distinct
            # gene annotations overlap. We want a sample
            # of non-overlaping genes
            while overlaps:
                key = random.sample(list(genes_set),
                                    1)[0]
                genes_ids = list(sample_dict.keys())
                if not genes_ids:
                    sample_dict[key] = copy.deepcopy(genes_set[key])
                    count += 1
                    continue

                if genes_ids:
                    if any([P.closed(self.db[key].start,
                                   self.db[key].end).overlaps(k) 
                            for k in [
                                P.closed(self.db[i].start, self.db[i].end)
                                for i in genes_ids
                            ]
                           ]):
                        continue
                    else:
                        sample_dict[key] = copy.deepcopy(genes_set[key])
                        overlaps = False
                        count += 1
        return sample_dict
        
        
    def sort_gene_sample(self, genes_dict):
        sorted_gene_sample = {chrm :
                              {gene[0]: 
                               copy.deepcopy(genes_dict[gene[0]]) 
                               for gene in sorted(
                                   [
                                       (key, value['coord']) 
                                       for key, value in genes_dict.items()
                                       if value['chrom'] == chrm
                                   ],
                                   key = lambda item: (item[1].lower, item[1].upper)
                               )
                              }
                              for chrm in self.chm_order
                             }
        return sorted_gene_sample
            
            
    def simulate_genes_with_exons_dup(self,
                                      duplicate_coding_exon_n,
                                      insert_location,
                                      genes_sample_dict):

        simulated_genes = dict()
        copy_genes_sample_dict = copy.deepcopy(genes_sample_dict)
        for gene_id, hierarchy_dict in copy_genes_sample_dict.items():
            transcript_id = next(iter(
                self.gene_hierarchy_dict[gene_id].keys()
            ))
            temp_dict = dict()
            exon_count = 0
            sim_gene_comp, sim_gene_coord = list(), list()
            intervals_list = list(hierarchy_dict.keys())
            for idx, (interval, annotation) in enumerate(hierarchy_dict.items()):
                sim_gene_coord.append(interval)
                sim_gene_comp.append(annotation['id'])
                if annotation['type'] == 'coding_exon':
                    exon_count += 1
                if exon_count == duplicate_coding_exon_n:
                    if insert_location == 'before':
                        ins = -1
                    elif insert_location == 'after':
                        ins = 1
                    insert_in_interval = intervals_list[idx+ins]
                    insert_in = copy.deepcopy(
                        hierarchy_dict[insert_in_interval]
                    )
                    if insert_in['type'] == 'intron':
                        start, end = (insert_in_interval.lower,
                                      insert_in_interval.upper)
                        center = (end - start) // 2
                        # we want to leave an evolutionary marker so we insert the 
                        # exon duplication leaving 10 bp in the start and end of 
                        # the intron.
                        temp_list = copy.deepcopy(sim_gene_coord)
                        temp_coords = [
                            P.open(start, start + center),
                            P.open((interval.lower - 10),
                                   (interval.upper + 10)),
                            P.open(start + center, end)
                                      ]
                        if start > temp_list[-1].lower:
                            sim_gene_coord.extend(temp_coords)
                        else:
                            for coord_idx, coord in enumerate(temp_list):
                                if start < coord.lower:
                                    sim_gene_coord = (
                                        sim_gene_coord[:(coord_idx -1)] 
                                        + temp_coords 
                                        + sim_gene_coord[coord_idx:]
                                    )
                                    break
                        if insert_location == 'before':
                            sim_gene_coord.extend(intervals_list[idx+1:])
                            sim_gene_comp.extend([insert_in['id']])
                            sim_gene_comp.extend([hierarchy_dict[intv]['id']
                                                  for intv in intervals_list[idx:]])
                        if insert_location == 'after':
                            sim_gene_coord.extend(intervals_list[idx+2:])
                            sim_gene_comp.extend([insert_in['id'],
                                                  annotation['id'],
                                                  insert_in['id']])
                            sim_gene_comp.extend([hierarchy_dict[intv]['id']
                                                  for intv in intervals_list[idx+2:]])
                    else:
                        raise Exception('two consecutive exons')
                    break
            temp_sim_gene_coord = list()
            for idx, i in enumerate(sim_gene_coord):
                if (i.lower not in [j.upper  for j in sim_gene_coord]
                    and idx != 0
                   ):
                    temp_sim_gene_coord.append(P.open((i.lower - 1),
                                                      i.upper))
                else:
                    temp_sim_gene_coord.append(P.open(i.lower,
                                                      i.upper))
            if insert_location == 'before':
                mut_loct = 1
            elif insert_location == 'after':
                mut_loct = 0
            ce_dup_loc = (duplicate_coding_exon_n - mut_loct)
            insertion_in_intron = dict()
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
            structure_breakdown = dict()
            for sim_coord, sim_id in zip(temp_sim_gene_coord,
                                         sim_gene_comp):
                chrm = annotation['chrom']
                strand = annotation['strand']
                if 'intron' in sim_id:
                    temp = copy.deepcopy(
                        self.genome[chrm][strand][sim_coord.lower:sim_coord.upper]
                    )
                    if sim_coord == insertion_in_intron['before']:
                        if temp[-2:] != 'AG':
                            temp = temp[:-2] + 'AG'
                    if sim_coord == insertion_in_intron['after']:
                        if temp[:2] != 'GT':
                            temp = 'GT' + temp[2:]
                else:
                    temp = copy.deepcopy(
                        self.genome[chrm][strand][sim_coord.lower:sim_coord.upper])
                if sim_id not in structure_breakdown:
                    structure_breakdown[sim_id] = copy.deepcopy(temp)
                else:
                    # since we're inserting an exon in the middle of the intron
                    structure_breakdown[sim_id + '_1'] = copy.deepcopy(temp)
                seq += temp
            simulated_genes[gene_id] = {'chrom': chrm,
                                        'coord': P.closed(self.db[transcript_id].start,
                                                          self.db[transcript_id].end),
                                        'features_order': sim_gene_comp,
                                        'features_intervals': temp_sim_gene_coord,
                                        'structure_breakdown':structure_breakdown,
                                        'seq': seq}

        return self.sort_gene_sample(simulated_genes)

    
    def simulate_genome(self,
                        sorted_sim_genes,
                        insert_location,
                        out_filename):
        
        padding = len(max(self.chm_order, key=len))
        if insert_location == 'before':
            exon_idx = 2
        if insert_location == 'after':
            exon_idx = 1
        new_genome_dict = dict()
        sanity_check = list()
        for chrm in self.chm_order:
            if not sorted_sim_genes[chrm]:
                new_genome_dict[chrm] = copy.deepcopy(self.genome[chrm]['+'])
            else:
                end = 0
                len_dup = 0
                continue_ = True
                ch_seq_with_dipl = ''
                temp = copy.deepcopy(sorted_sim_genes[chrm])
                for idx, (sim_gene, sim_trans) in enumerate(temp.items()):
                    if continue_:
                        transcript = next(iter(self.gene_hierarchy_dict[sim_gene].keys()))
                        start, e = (sim_trans['coord'].lower,
                                    sim_trans['coord'].upper)
                        temp_seq = (copy.deepcopy(self.genome[chrm]['+'][end:start])
                                    + sim_trans['seq'])
                        temp_sb = copy.deepcopy(sim_trans['structure_breakdown'])
                        len_dup += (len(temp_sb[
                            [i for i in temp_sb.keys() if 'exon' in i][exon_idx]
                        ]) 
                                    + self.intron_signal)
                        if  idx < (len(sorted_sim_genes[chrm])-1):
                            ch_seq_with_dipl += temp_seq
                        else:
                            ch_seq_with_dipl += (temp_seq 
                                                 + copy.deepcopy(self.genome[chrm]['+'][e:])
                                                )
                            continue_ = False
                        end = sim_trans['coord'].upper
                #sanity check
                sanity_check.append([len(ch_seq_with_dipl) - len_dup  == len(self.genome[chrm]['+'])])
                new_genome_dict[chrm] = ch_seq_with_dipl
        if all(sanity_check):
            print('sanity check: passed')
        else:
            print('sanity check: not passed - check lengths')
                
        if out_filename != '':
            with open(out_filename, "w") as handle:
                for chrom in self.chm_order:
                    if chrom not in new_genome_dict:
                        temp_seq = Seq(self.genome[chrom]['+'])
                    else:
                        temp_seq = Seq(new_genome_dict[chrom])
                        
                    record = SeqRecord(temp_seq,
                                       id=chrom,
                                       description='')
                    SeqIO.write(record,
                                handle,
                                "fasta")
                    
        return new_genome_dict
    
    
    def get_coding_exons_interv(self, sorted_sim_genes):
        ce_dict = {chrom_id :
                    {gene_id : {
                        i:j for i, j in zip(transcr_dict['features_order'],
                                            transcr_dict['features_intervals']) 
                        if 'exon' in i
                    } 
                     for gene_id, transcr_dict in chrom_dict.items() 
                    }
                    for chrom_id, chrom_dict in sorted_sim_genes.items() 
                   }
        return ce_dict
    
    
    def get_new_exon_intervals(self,
                               new_genome_dict,
                               coding_exons_interv,
                               sorted_gene_sample):
        sim_genes_interv_in_chr = {chrm:
                                   {
                                       gene_id : 
                                       {
                                           'seq': [
                                               P.closed(i.span()[0], i.span()[1]) for i in 
                                               re.finditer(transcr_dict['seq'],
                                                           new_genome_dict[chrm])
                                           ],
                                          'coding_exons_dups': {
                                              exon_id:[
                                                  i.span() for i in re.finditer(
                                                      self.genome[chrm]['+'][intev.lower:intev.upper],
                                                      new_genome_dict[chrm])
                                              ]
                                              for exon_id, intev in 
                                              coding_exons_interv[chrm][gene_id].items()
                                          },
                                           'coding_exons': {
                                               exon_id: [i.span() for i in re.finditer(
                                                   self.genome[chrm]['+'][intev.lower:intev.upper],
                                                   self.genome[chrm]['+']
                                               )
                                                        ] 
                                               for exon_id, intev in
                                               coding_exons_interv[chrm][gene_id].items()
                                                          }
                                         }
                                       for gene_id, transcr_dict in chrom_dict.items()
                                   }
                                   for chrm, chrom_dict in sorted_gene_sample.items()
                                  }
        return sim_genes_interv_in_chr
    
    
    def intron_signal_sanity_check(self,
                                   sample_size,
                                   insert_location,
                                   new_genome_dict, 
                                   sim_genes_interv_in_chr):
        padding = max(self.flatten([[len(gene_id) for gene_id,
                                gene_dict in chrm_dict.items()] 
                               for chrm, chrm_dict 
                               in sim_genes_interv_in_chr.items()]))
        # sanity check - introns signal - 
        # All introns start with an - AG and end with a GT
        if insert_location == 'before':
            loc = 0
        if insert_location == 'after': 
            loc = 1
        sanity_check = list()
        for chrm, chrm_dict in sim_genes_interv_in_chr.items():
            for gene_id, gene_dict in chrm_dict.items():
                for id_, annot in gene_dict['coding_exons_dups'].items():
                    if len(annot) == 2:
                        s, e = annot[loc]
                        sanity_check.append([gene_id, 
                                             new_genome_dict[chrm][(s-2):(s)],
                                             new_genome_dict[chrm][e:(e+2)]])
        if len(sanity_check) > sample_size:
            print('Warning: real exon duplications within a gene')
        # if all intron signals are the same; "AG", "GT"
        if all([(i[1].upper(), i[2].upper()) == ("AG", "GT") for i in sanity_check]):
            print('sanity check: passed')
        else:
            for i in [i for i in sanity_check 
                      if (i[1].upper(), i[2].upper()) != ("AG", "GT")
                     ]:
                print(i[0] + ' ' * (padding - len(i[0])), 
                      i[1].upper(), i[2].upper())

                
    def collect_protein_seqs(self,
                         file_direct,
                         split_ref):
        prot_seq = dict()
        fasta_sequences = SeqIO.parse(open(file_direct),
                                      'fasta')
        for fasta in fasta_sequences:
            gene = [i.rsplit(split_ref)[1] 
                    for i in fasta.description.rsplit(' ') 
                    if 'gene' in i]
            if gene:
                if gene[0] not in prot_seq:
                    prot_seq[gene[0]] = dict()
                prot_seq[gene[0]][fasta.id] = str(fasta.seq)
        return prot_seq

    
    def get_alignment_queries(self,
                              prot_seq_dict,
                              genes_ids_list,
                              out_file_direct):
        proteins_simulated_genes = dict()
        protein_gene_dict = dict()
        with open(out_file_direct, "w") as handle:
             for ID in genes_ids_list:
                gene_id = self.db[ID].attributes['gene_id'][0]
                gene_id_prot_dict = [i for i in list(prot_seq_dict.keys())
                                     if gene_id in i]
                if gene_id_prot_dict:
                    id_ = gene_id_prot_dict[0]
                    temp = copy.deepcopy(prot_seq_dict[id_])
                    for prot_id, seq in temp.items():  
                        record = SeqRecord(Seq(seq),
                                           id = str(prot_id),
                                          description = id_)
                        SeqIO.write(record,
                                    handle,
                                    "fasta")
                        protein_gene_dict[prot_id.rsplit('.')[0]] = id_
        return protein_gene_dict

                               