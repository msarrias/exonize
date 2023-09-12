from utils import *
import sqlite3


class ExonizeEvents(object):
    def __init__(self, exonize_res_db_path):
        self.classified_events = {}
        self.results_db = exonize_res_db_path
        self.skip_duplicated_pairs = list()
        self.all_events = dict()
        self.coverage_threshold = 0.9
        self.pair_id_idx = -2
        self.evalue_idx = -3
        self.fragment_id_idx = 0
        self.event_type_idx = -1
        self.optional_group_category_idx = 15
        self.query_start_idx = 5
        self.query_end_idx = 6
        self.target_start_idx = 7
        self.target_end_idx = 8
        self.dna_perc_identity_idx = 9
        self.prot_perc_identity_idx = 10
        self.coding_event_dict = {'1001': 'obligate_events',
                                  '0111': 'MXEs_events',
                                  '0101': 'deactivated_pairs',
                                  '1111': 'flexible_pairs',
                                  '1101': 'flexible_pairs',
                                  '1011': 'flexible_pairs'}
        self.event_coding_dict = {'obligate_events': ['1001'],
                                  'MXEs_events': ['0111'],
                                  'deactivated_pairs': ['0101'],
                                  'flexible_pairs': ['1111', '1101', '1011']}
        self.deactiv_features = ['DEACTIVATED',
                                 'OUT_OF_MRNA',
                                 'INS_UTR',
                                 'TRUNC',
                                 'NEITHER']
# b =

    def fragments_count_sanity_check(self) -> None:
        pairs_ids_len = 0
        split_pairs_len = 0
        all_full_events = self.get_all_full_events()
        for mut_class, mut_dict in self.all_events.items():
            pairs_ids_len += len([*mut_dict['pairs_fragment_ids'], *mut_dict['orphan_fragment_ids']])
            split_pairs_len += len(mut_dict['splits_fragment_ids'])
        if not (split_pairs_len + pairs_ids_len + len(self.skip_duplicated_pairs) == len(all_full_events)):
            print('WARNING: some events are doubled counted or missing from the analysis')

    def check_for_duplications(self) -> None:
        db = sqlite3.connect(self.results_db)
        cursor = db.cursor()
        cursor.execute("""
        SELECT
            pair_id 
        FROM (SELECT pair_id, COUNT(*) as count 
        FROM Full_length_events_cumulative_counts
        GROUP BY pair_id) WHERE count < 2""")
        self.skip_duplicated_pairs = [i[0] for i in cursor.fetchall()]
        db.close()
        if self.skip_duplicated_pairs:
            print('WARNING: overlapping events in the Full_length_events_cumulative_counts table,'
                  ' skipping said events')

    def organize_event_dict(self, events_list) -> dict:
        pair_ids_dict = self.get_count_pair_ids()
        records = dict(pairs_fragment_ids=list(), pairs_unique_record=list(),
                       splits_fragment_ids=list(), split_records=list(),
                       orphan_fragment_ids=list(), orphan_records=list())
        processed_ids = list()
        for record in events_list:
            if (record[self.pair_id_idx] is not None
                    and record[self.pair_id_idx] not in processed_ids
                    and record[self.pair_id_idx] not in self.skip_duplicated_pairs):
                records_list = [r for r in events_list if r[self.pair_id_idx] == record[self.pair_id_idx]]
                if len(records_list) == pair_ids_dict[record[self.pair_id_idx]]:
                    records['pairs_fragment_ids'].extend([i[self.fragment_id_idx] for i in records_list])
                    # we just want one event per pair, so we take the one with the lowest evalue
                    exclude_reciprocals = [i for i in records_list if 'TRUNC' not in i[self.event_type_idx]]
                    if not exclude_reciprocals:
                        exclude_reciprocals = records_list
                    records['pairs_unique_record'].append(min(exclude_reciprocals, key=lambda x: x[self.evalue_idx]))
                    processed_ids.append(record[self.pair_id_idx])
                else:
                    records['split_records'].extend(records_list)
                    records['splits_fragment_ids'].extend([i[self.fragment_id_idx] for i in records_list])
                    processed_ids.append(record[self.pair_id_idx])
            elif record[self.pair_id_idx] is None:
                records['orphan_records'].append(record)
                records['orphan_fragment_ids'].append(record[self.fragment_id_idx])
        return records

    def search_for_split_events(self):
        new_dict = dict()
        pair_ids_dict = self.get_count_pair_ids()
        all_split_records = [record for event_dict in self.all_events.values() for record in event_dict['split_records']]
        if all_split_records:
            new_dict = {
                i[self.pair_id_idx]: [j for j in all_split_records if i[self.pair_id_idx] == j[self.pair_id_idx]]
                for i in all_split_records}
            # ### CHECK THERE ARE NOT UNRECONCILED ANNOTATIONS
            for key, value in new_dict.items():
                if len(value) != pair_ids_dict[key]:
                    print('WARNING: Unreconciled annotations in the split events')
        return new_dict

    def reassign_split_pairs(self):
        split_events_dict = self.search_for_split_events()
        skip_frag_id, skip_pair_id = list(), list()
        for mut_type in self.all_events.keys():
            new_split_records = list(self.all_events[mut_type]['split_records'])
            new_split_fragments = list(self.all_events[mut_type]['splits_fragment_ids'])
            new_split_records_ids = [i[self.fragment_id_idx] for i in new_split_records]
            for split_event in self.all_events[mut_type]['split_records']:
                if split_event[self.fragment_id_idx] not in skip_frag_id:
                    pairs = [i for i in split_events_dict[split_event[self.pair_id_idx]]
                             if i[self.fragment_id_idx] != split_event[self.fragment_id_idx]
                             and split_event[self.pair_id_idx] not in skip_pair_id]
                    if len(pairs) > 1:
                        id_min_evalue = None
                        skip_id, annot_grouping, list_events = list(), list(), list()
                        events_list = list(split_events_dict[split_event[self.pair_id_idx]])
                        for annot in events_list:
                            if annot[self.fragment_id_idx] not in skip_id:
                                temp = [i for i in events_list if
                                        i != annot and P.open(i[3], i[4]).overlaps(P.open(annot[3], annot[4]))]
                                if temp:
                                    skip_id.extend([annot[0], *[i[0] for i in temp]])
                                annot_grouping.extend([[annot, *temp]])
                        if len(annot_grouping) == 1:
                            # this means that there is not a reciprocal hit
                            id_min_evalue = min([(i[self.fragment_id_idx], i[self.evalue_idx])
                                                 for i in annot_grouping[0]], key=lambda x: x[1])[0]
                            list_events = annot_grouping[0]
                        elif len(annot_grouping) == 2:  # pair of hit and reciprocal hit - or the other way around
                            hit, recip_hit = annot_grouping
                            if (all('INS' in i[self.event_type_idx] for i in hit)
                                    and any('TRUNC' in i[self.event_type_idx] for i in recip_hit) or
                                    all('INS' in i[self.event_type_idx] for i in recip_hit)
                                    and any('TRUNC' in i[self.event_type_idx] for i in hit)):
                                id_min_evalue = min([(i[self.fragment_id_idx], i[self.evalue_idx])
                                                     for i in [*hit, *recip_hit]], key=lambda x: x[1])[0]
                                list_events = [*hit, *recip_hit]
                        if split_event[self.fragment_id_idx] == id_min_evalue:
                            self.all_events[mut_type]['pairs_unique_record'].append(split_event)
                            for event in list_events:
                                self.all_events[mut_type]['pairs_fragment_ids'].append(event[self.fragment_id_idx])
                                if event[self.fragment_id_idx] in new_split_records_ids:
                                    new_split_records.remove(event)
                                    new_split_fragments.remove(event[self.fragment_id_idx])
                                else:
                                    skip_frag_id.append(event[self.fragment_id_idx])
                                    skip_pair_id.append(event[self.pair_id_idx])
                    if len(pairs) == 1:
                        pair = pairs[0]
                        if pair[self.fragment_id_idx] not in skip_frag_id:
                            ins_event = False
                            group_categ = get_unmatched_events(split_event[self.event_type_idx],
                                                               pair[self.event_type_idx])
                            if group_categ:
                                ins_event = (any(element in self.deactiv_features for element in group_categ)
                                             and any('INS' in element for element in group_categ)
                                             and split_event[self.evalue_idx] < pair[self.evalue_idx])
                            # check if the split event and the pair overlap
                            pair_t = P.open(pair[self.target_start_idx], pair[self.target_end_idx])
                            pair_q = P.open(pair[self.query_start_idx], pair[self.query_end_idx])
                            split_t = P.open(split_event[self.target_start_idx], split_event[self.target_end_idx])
                            split_q = P.open(split_event[self.query_start_idx], split_event[self.query_end_idx])
                            overlap_se = ((get_average_overlapping_percentage(pair_t, split_t) > self.coverage_threshold
                                           and get_average_overlapping_percentage(pair_q, split_q) > self.coverage_threshold)
                                          and split_event[self.evalue_idx] < pair[self.evalue_idx])
                            if ins_event or overlap_se or 'TRUNC' in pair[self.event_type_idx]:
                                self.all_events[mut_type]['pairs_fragment_ids'].extend([i[self.fragment_id_idx]
                                                                                        for i in [split_event, pair]])
                                self.all_events[mut_type]['pairs_unique_record'].append(split_event)
                                # handle the reciprocal event

                                if split_event[self.fragment_id_idx] in new_split_records_ids:
                                    new_split_records.remove(split_event)
                                    new_split_fragments.remove(split_event[self.fragment_id_idx])
                                else:
                                    skip_frag_id.append(split_event[self.fragment_id_idx])
                                    skip_pair_id.append(split_event[self.pair_id_idx])
                                # handle the event that is not reciprocal
                                if pair in new_split_records:
                                    new_split_records.remove(pair)
                                    new_split_fragments.remove(pair[self.fragment_id_idx])
                                else:
                                    skip_frag_id.append(pair[self.fragment_id_idx])
                                    skip_pair_id.append(pair[self.pair_id_idx])
            self.all_events[mut_type]['split_records'] = list(new_split_records)
            self.all_events[mut_type]['splits_fragment_ids'] = list(new_split_fragments)
        # If there are still events pending to re-assign
        if skip_frag_id:
            new_remove_frag_id = list(skip_frag_id)
            for mut_type in self.all_events.keys():
                split_records = list(self.all_events[mut_type]['split_records'])
                split_fragments = list(self.all_events[mut_type]['splits_fragment_ids'])
                records_split_ids = [i[self.fragment_id_idx] for i in self.all_events[mut_type]['split_records']]
                remove_ids = list(set(new_remove_frag_id) & set(records_split_ids))
                if remove_ids:
                    self.all_events[mut_type]['split_records'] = [i for i in split_records if i[self.fragment_id_idx]
                                                                  not in remove_ids]
                    self.all_events[mut_type]['splits_fragment_ids'] = [i for i in split_fragments if i not in remove_ids]
                    new_remove_frag_id = list(set(new_remove_frag_id) - set(remove_ids))
        self.fragments_count_sanity_check()

    def reorganize_schema_according_to_type(self) -> None:
        for key, value in self.all_events.items():
            all_events = [*value['orphan_records'], *value['pairs_unique_record']]
            if key != 'optional_pairs':
                self.classified_events[key] = {i: [j for j in all_events if j[self.event_type_idx] == i]
                                               for i in set([i[self.event_type_idx] for i in all_events])}
            else:
                self.classified_events[key] = {self.coding_event_dict[i]: {
                    k: [j for j in all_events if j[self.event_type_idx] == k
                        and j[self.optional_group_category_idx] == i]
                    for k in set([i[self.event_type_idx] for i in [j for j in all_events
                                                                          if j[self.optional_group_category_idx] == i]])}
                    for i in set([i[self.optional_group_category_idx] for i in all_events])}

    def get_classification_schema(self) -> None:
        self.all_events = dict(flexible_pairs=self.fetch_flexible_events(),
                               optional_pairs=self.fetch_optional_events(),
                               deactivated_pairs=self.fetch_deactivated_event(),
                               MXEs_events=self.fetch_exclusive_events(),
                               obligate_events=self.fetch_obligate_pairs())
        self.fragments_count_sanity_check()
        self.reassign_split_pairs()
        self.reorganize_schema_according_to_type()

    # ###### QUERY TABLES ########
    def get_all_full_events(self):
        db = sqlite3.connect(self.results_db)
        cursor = db.cursor()
        cursor.execute(""" 
        SELECT 
            DISTINCT(fragment_id)
        FROM Full_length_events_cumulative_counts;""")
        all_events = [i[0] for i in cursor.fetchall()]
        db.close()
        return all_events

    def fetch_obligate_pairs(self) -> dict:
        db = sqlite3.connect(self.results_db)
        cursor = db.cursor()
        cursor.execute("""
        SELECT 
            f.fragment_id,
            f.gene_id,
            f.mrna_count,
            f.CDS_start - f.gene_start as CDS_start,
            f.CDS_end - f.gene_start as CDS_end,
            f.query_start,
            f.query_end,
            f.target_start,
            f.target_end,
            fr.dna_perc_identity,
            fr.prot_perc_identity,
            f.cum_both AS cb,
            f.cum_query AS cq,
            f.cum_target AS ct,
            f.cum_neither AS cn,
            f.evalue,
            f.pair_id,
            f.event_type
        FROM Full_length_events_cumulative_counts AS f 
        INNER JOIN Fragments AS fr ON fr.fragment_id = f.fragment_id
        WHERE f.mrna_count = f.cum_both
        ORDER BY f.pair_id;
            """)
        obligate_events_dict = self.organize_event_dict(cursor.fetchall())
        return obligate_events_dict

    def fetch_exclusive_events(self) -> dict:
        db = sqlite3.connect(self.results_db)
        cursor = db.cursor()
        cursor.execute("""
        SELECT 
            f.fragment_id,
            f.gene_id,
            f.mrna_count,
            f.CDS_start - f.gene_start as CDS_start,
            f.CDS_end - f.gene_start as CDS_end,
            f.query_start,
            f.query_end,
            f.target_start,
            f.target_end,
            fr.dna_perc_identity,
            fr.prot_perc_identity,
            f.cum_both AS cb,
            f.cum_query AS cq,
            f.cum_target AS ct,
            f.cum_neither AS cn,
            f.evalue,
            f.pair_id,
            f.event_type
        FROM Full_length_events_cumulative_counts AS f 
        INNER JOIN Fragments AS fr ON fr.fragment_id = f.fragment_id
        WHERE f.mrna_count = (f.cum_query + f.cum_target) 
        AND f.cum_target < f.mrna_count 
        AND f.cum_target > 0
        ORDER BY f.pair_id;
            """)
        MXEs_events_dict = self.organize_event_dict(cursor.fetchall())
        db.close()
        return MXEs_events_dict

    def fetch_deactivated_event(self) -> dict:
        db = sqlite3.connect(self.results_db)
        cursor = db.cursor()
        cursor.execute("""
        SELECT 
            f.fragment_id,
            f.gene_id,
            f.mrna_count,
            f.CDS_start - f.gene_start as CDS_start,
            f.CDS_end - f.gene_start as CDS_end,
            f.query_start,
            f.query_end,
            f.target_start,
            f.target_end,
            fr.dna_perc_identity,
            fr.prot_perc_identity,
            f.cum_both AS cb,
            f.cum_query AS cq,
            f.cum_target AS ct,
            f.cum_neither AS cn,
            f.evalue,
            f.pair_id,
            f.event_type
        FROM Full_length_events_cumulative_counts AS f 
        INNER JOIN Fragments AS fr ON fr.fragment_id = f.fragment_id
        WHERE f.cum_query=f.mrna_count
        ORDER BY f.pair_id;
            """)
        deactivated_events_dict = self.organize_event_dict(cursor.fetchall())
        db.close()
        return deactivated_events_dict

    def fetch_optional_events(self) -> dict:
        db = sqlite3.connect(self.results_db)
        cursor = db.cursor()
        cursor.execute("""
        SELECT
            gs.fragment_id,
            gs.gene_id,
            gs.mrna_count,
            gs.CDS_start - gs.gene_start as CDS_start,
            gs.CDS_end - gs.gene_start as CDS_end,
            gs.query_start,
            gs.query_end,
            gs.target_start,
            gs.target_end,
            gs.dna_perc_identity,
            gs.prot_perc_identity,
            gs.cum_both AS cb,
            gs.cum_query AS cq,
            gs.cum_target AS ct,
            gs.cum_neither AS cn,
            gs.cb || gs.cq || gs.ct || gs.cn AS optional_group_category,
            gs.evalue,
            gs.pair_id,
            gs.event_type
        FROM (
        SELECT
            subquery.*,
        CASE WHEN cum_query <> 0 THEN 1 ELSE 0 END AS cq,
        CASE WHEN cum_target <> 0 THEN 1 ELSE 0 END AS ct,
        CASE WHEN cum_both <> 0 THEN 1 ELSE 0 END AS cb,
        CASE WHEN cum_neither <> 0 THEN 1 ELSE 0 END AS cn
        FROM (
        SELECT
            fle.fragment_id,
            fle.gene_id,
            fle.gene_start,
            fle.gene_end,
            fle.mrna_count,
            fle.CDS_start,
            fle.CDS_end,
            fle.query_start,
            fle.query_end,
            fle.target_start,
            fle.target_end,
            fr.dna_perc_identity,
            fr.prot_perc_identity,
            fle.cum_both,
            fle.cum_query,
            fle.cum_target,
            fle.cum_neither,
            fle.evalue,
            fle.pair_id,
            fle.event_type
        FROM Full_length_events_cumulative_counts AS fle
        INNER JOIN Fragments AS fr ON fr.fragment_id = fle.fragment_id
        WHERE (fle.cum_neither > 0 AND (fle.cum_query > 0 OR fle.cum_target > 0 OR fle.cum_both > 0))
        ) AS subquery
        ) AS gs
        ORDER BY gs.pair_id;
        """)
        optional_events_dict = self.organize_event_dict(cursor.fetchall())
        db.close()
        return optional_events_dict

    def fetch_flexible_events(self) -> dict:
        db = sqlite3.connect(self.results_db)
        cursor = db.cursor()
        cursor.execute("""
        SELECT
            gs.fragment_id,
            gs.gene_id,
            gs.mrna_count,
            gs.CDS_start - gs.gene_start as CDS_start,
            gs.CDS_end - gs.gene_start as CDS_end,
            gs.query_start,
            gs.query_end,
            gs.target_start,
            gs.target_end,
            gs.dna_perc_identity,
            gs.prot_perc_identity,
            gs.cum_both AS cb,
            gs.cum_query AS cq,
            gs.cum_target AS ct,
            gs.cum_neither AS cn,
            gs.cb || gs.cq || gs.ct || gs.cn AS optional_group_category,
            gs.evalue,
            gs.pair_id,
            gs.event_type
        FROM (
        SELECT
            subquery.*,
        CASE WHEN cum_query <> 0 THEN 1 ELSE 0 END AS cq,
        CASE WHEN cum_target <> 0 THEN 1 ELSE 0 END AS ct,
        CASE WHEN cum_both <> 0 THEN 1 ELSE 0 END AS cb,
        CASE WHEN cum_neither <> 0 THEN 1 ELSE 0 END AS cn
        FROM (
        SELECT
            fle.fragment_id,
            fle.gene_id,
            fle.gene_start,
            fle.gene_end,
            fle.mrna_count,
            fle.CDS_start,
            fle.CDS_end,
            fle.query_start,
            fle.query_end,
            fle.target_start,
            fle.target_end,
            fr.dna_perc_identity,
            fr.prot_perc_identity,
            fle.cum_both,
            fle.cum_query,
            fle.cum_target,
            fle.cum_neither,
            fle.evalue,
            fle.pair_id,
            fle.event_type
        FROM Full_length_events_cumulative_counts AS fle
        INNER JOIN Fragments AS fr ON fr.fragment_id = fle.fragment_id
        WHERE (fle.cum_both > 0 AND (fle.cum_query > 0 OR fle.cum_target > 0) 
        AND fle.cum_neither = 0)
        ) AS subquery
        ) AS gs
        ORDER BY gs.pair_id;
        """)
        flexible_events_dict = self.organize_event_dict(cursor.fetchall())
        db.close()
        return flexible_events_dict

    def get_count_pair_ids(self) -> dict:
        db = sqlite3.connect(self.results_db)
        cursor = db.cursor()
        cursor.execute("""
        SELECT 
            pair_id, 
            COUNT(*) 
        FROM Full_length_events_cumulative_counts GROUP BY pair_id""")
        pair_ids_dict = {i[0]: i[1] for i in cursor.fetchall()}
        return pair_ids_dict
