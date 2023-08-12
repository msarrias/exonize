from utils import *
import sqlite3


class ExonizeEvents(object):
    def __init__(self, exonize_res_db_path):
        self.results_db = exonize_res_db_path
        self.skip_duplicated_pairs = []

    def check_for_duplications(self):
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

    def get_n_pair_ids(self):
        db = sqlite3.connect(self.results_db)
        cursor = db.cursor()
        cursor.execute("""
        SELECT 
            pair_id, 
            COUNT(*) 
        FROM Full_length_events_cumulative_counts GROUP BY pair_id""")
        res = {i[0]: i[1] for i in cursor.fetchall()}
        return res

    def organize_event_dict(self, events_list):
        pair_id_idx = -2
        pair_ids_dict = self.get_n_pair_ids()
        records = dict(pairs=dict(), orphans=[], split_pairs=[])
        processed_ids = []
        for record in events_list:
            if record[pair_id_idx] is not None:
                if record[pair_id_idx] not in processed_ids:
                    ids_list = [r for r in events_list if r[pair_id_idx] == record[pair_id_idx]]
                    if record[pair_id_idx] not in self.skip_duplicated_pairs:
                        if len(ids_list) == pair_ids_dict[record[pair_id_idx]]:
                            records['pairs'][record[pair_id_idx]] = ids_list
                            processed_ids.append(record[pair_id_idx])
                        else:
                            records['split_pairs'][record[pair_id_idx]] = ids_list
                            processed_ids.append(record[pair_id_idx])
            else:
                records['orphans'].append(record)
        return records

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
            f.cum_both AS cb,
            f.cum_query AS cq,
            f.cum_target AS ct,
            f.cum_neither AS cn,
            f.evalue,
            f.concat_event_type,
            f.pair_id
        FROM Full_length_events_cumulative_counts AS f 
        WHERE f.mrna_count = f.cum_both
        ORDER BY f.pair_id;
            """)
        obligate_events_list = generate_event_list(cursor.fetchall())
        obligate_events_dict = self.organize_event_dict(obligate_events_list)
        return obligate_events_dict

    def fetch_exclusive_events(self):
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
            f.cum_both AS cb,
            f.cum_query AS cq,
            f.cum_target AS ct,
            f.cum_neither AS cn,
            f.evalue,
            f.concat_event_type,
            f.pair_id
        FROM Full_length_events_cumulative_counts AS f 
        WHERE f.mrna_count = (f.cum_query + f.cum_target) 
        AND f.cum_target < f.mrna_count 
        AND f.cum_target > 0
        ORDER BY f.pair_id;
            """)
        MXEs_pairs = cursor.fetchall()
        MXEs_events_dict = generate_event_list(MXEs_pairs)
        MXEs_events_dict = self.organize_event_dict(MXEs_events_dict)
        db.close()
        return MXEs_events_dict

    def fetch_deactivated_event(self):
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
            f.cum_both AS cb,
            f.cum_query AS cq,
            f.cum_target AS ct,
            f.cum_neither AS cn,
            f.evalue,
            f.concat_event_type,
            f.pair_id
        FROM Full_length_events_cumulative_counts AS f 
        WHERE f.cum_query==f.mrna_count
        ORDER BY f.pair_id;
            """)
        deactivated_events_list = generate_event_list(cursor.fetchall())
        deactivated_events_dict = self.organize_event_dict(deactivated_events_list)
        db.close()
        return deactivated_events_dict

    def fetch_optional_events(self):
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
            gs.cum_both AS cb,
            gs.cum_query AS cq,
            gs.cum_target AS ct,
            gs.cum_neither AS cn,
            gs.cq || gs.ct || gs.cb || gs.cn AS optional_group_category,
            gs.evalue,
            gs.concat_event_type,
            gs.pair_id
        FROM (
        SELECT
            subquery.*,
        CASE WHEN cum_query <> 0 THEN 1 ELSE 0 END AS cq,
        CASE WHEN cum_target <> 0 THEN 1 ELSE 0 END AS ct,
        CASE WHEN cum_both <> 0 THEN 1 ELSE 0 END AS cb,
        CASE WHEN cum_neither <> 0 THEN 1 ELSE 0 END AS cn
        FROM (
        SELECT
            fle.*
        FROM Full_length_events_cumulative_counts AS fle
        WHERE (fle.cum_neither > 0 AND (fle.cum_query > 0 OR fle.cum_target > 0 OR fle.cum_both > 0))
        ) AS subquery
        ) AS gs
        ORDER BY gs.pair_id;
        """)
        optional_events_list = generate_event_list(cursor.fetchall())
        optional_events_dict = self.organize_event_dict(optional_events_list)
        db.close()
        return optional_events_dict

    def fetch_flexible_events(self):
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
            gs.cum_both AS cb,
            gs.cum_query AS cq,
            gs.cum_target AS ct,
            gs.cum_neither AS cn,
            gs.cq || gs.ct || gs.cb || gs.cn AS optional_group_category,
            gs.evalue,
            gs.concat_event_type,
            gs.pair_id
        FROM (
        SELECT
            subquery.*,
        CASE WHEN cum_query <> 0 THEN 1 ELSE 0 END AS cq,
        CASE WHEN cum_target <> 0 THEN 1 ELSE 0 END AS ct,
        CASE WHEN cum_both <> 0 THEN 1 ELSE 0 END AS cb,
        CASE WHEN cum_neither <> 0 THEN 1 ELSE 0 END AS cn
        FROM (
        SELECT
            fle.*
        FROM Full_length_events_cumulative_counts AS fle
        WHERE (fle.cum_both > 0 AND (fle.cum_query > 0 OR fle.cum_target > 0) 
        AND fle.cum_neither = 0)
        ) AS subquery
        ) AS gs
        ORDER BY gs.pair_id;
        """)
        flexible_events_list = generate_event_list(cursor.fetchall())
        flexible_events_dict = self.organize_event_dict(flexible_events_list)
        db.close()
        return flexible_events_dict
