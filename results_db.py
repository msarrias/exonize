import sqlite3
import re # regular expressions


class SQlite_res():
    def __init__(self, db_path):
        self.db_path = db_path
        self.colnames = []
        self.create_table_concatenate_fragments()
        self.get_table_col_names('concatenate_fragments')
        self.insert_conservation_percentage_columns()
          
            
    def get_tables_name(self):
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables_names = cursor.fetchall()
        conn.close()
        return tables_names
                              
        
    def get_table_col_names(self,table_name):
        conn = sqlite3.connect(self.db_path)
        info = conn.execute(f"PRAGMA table_info('{table_name}')").fetchall()
        colnames = [col_name[1] for col_name in info]
        conn.close()
        return colnames
      
    
    @staticmethod        
    def hd_dna_seqs(seq_a, seq_b):
        return sum([i != j for i, j in zip(seq_a, seq_b)])

    
    @staticmethod
    def hd_amino_seqs(seq_a, seq_b):
        return sum([seq_a[i:i+3] != seq_b[i:i+3] for i in range(0,len(seq_a), 3)])
    
    @staticmethod
    def ungap_sequence(seq):
        return re.sub(r'[^\w]', '', seq)

    
    def create_table_concatenate_fragments(self):
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        new_table = """
        CREATE TABLE IF NOT EXISTS concatenate_fragments AS
        SELECT * 
        FROM (
        SELECT 
        gene_id,
        query_id,
        target_id, 
        query_start,
        query_end,
        target_start,
        target_end,
        query_end - query_start as query_length,
        target_end - target_start as target_length,
        HSP_n,
        MAX(score) AS score,
        MAX(category) AS category,
        MAX(strand) AS strand,
        GROUP_CONCAT(query_frame, ';') AS query_frame,
        GROUP_CONCAT(hit_frame, ';') AS hit_frame,
        SUM(n_stop_codons) AS n_stop_codons,
        GROUP_CONCAT(frag_query_start, ';') AS frag_query_start,
        GROUP_CONCAT(frag_query_end, ';') AS frag_query_end,
        GROUP_CONCAT(frag_target_start, ';') AS frag_hit_start,
        GROUP_CONCAT(frag_target_end, ';') AS frag_hit_end,
        GROUP_CONCAT(query_aln_dna_sequence, '') AS query_aln_dna_sequence, 
        GROUP_CONCAT(target_aln_dna_sequence, '') AS target_aln_dna_sequence,
        GROUP_CONCAT(query_aln_prot_sequence, '') AS query_aln_prot_sequence,
        GROUP_CONCAT(target_aln_prot_sequence, '') AS target_aln_prot_sequence
        FROM (
        SELECT * 
        FROM Fragments
        ORDER BY 
        gene_id ASC,
        query_id ASC,
        target_id ASC,
        HSP_n ASC,
        query_start ASC,
        query_end ASC)
        GROUP BY 
        gene_id, 
        query_id, 
        target_id,
        HSP_n);
        """
        #create table
        cursor.execute(new_table)
        #create index
        cursor.execute("""CREATE INDEX IF NOT EXISTS concat_frag_idx ON concatenate_fragments (gene_id, query_id, target_id);""")
        cursor.execute("""CREATE INDEX IF NOT EXISTS concat_frag_gene_idx ON concatenate_fragments (gene_id);""")
        conn.close()
        
        
    def insert_conservation_percentage_columns(self):
        concatenate_fragments_colnames = self.get_table_col_names('concatenate_fragments')
        new_colnames = ['dna_perc_conservation','prot_perc_conservation',
                        'query_dna_sequence','target_dna_sequence',
                        'query_prot_sequence','target_prot_sequence', 'query_aligned_fraction']
        
        if sum([True for i in new_colnames if i in concatenate_fragments_colnames]) == 0:
            conn = sqlite3.connect(self.db_path)
            cur = conn.cursor()
            cur.execute("ALTER TABLE concatenate_fragments ADD COLUMN dna_perc_conservation REAL;")
            cur.execute("ALTER TABLE concatenate_fragments ADD COLUMN prot_perc_conservation REAL;")
            cur.execute("ALTER TABLE concatenate_fragments ADD COLUMN query_dna_sequence VARCHAR;")
            cur.execute("ALTER TABLE concatenate_fragments ADD COLUMN target_dna_sequence VARCHAR;")
            cur.execute("ALTER TABLE concatenate_fragments ADD COLUMN query_prot_sequence VARCHAR;")
            cur.execute("ALTER TABLE concatenate_fragments ADD COLUMN target_prot_sequence VARCHAR;")
            cur.execute("ALTER TABLE concatenate_fragments ADD COLUMN query_aligned_fraction REAL;")
            #populate column
            cur.execute("""
            SELECT 
            gene_id,
            query_id,
            target_id,
            query_length,
            query_aln_dna_sequence,
            target_aln_dna_sequence,
            query_aln_prot_sequence,
            target_aln_prot_sequence,
            HSP_n,
            score
            FROM
            concatenate_fragments
            """)
            for row in cur:
                cursor_insert = conn.cursor()
                (gene_id, 
                 query_id, 
                 target_id,
                 query_length,
                 query_aln_dna_sequence,
                 target_aln_dna_sequence,
                 query_aln_prot_sequence,
                 target_aln_prot_sequence,
                 HSP_n,
                 score) = row
                dna_consev = 1 - (self.hd_dna_seqs(query_aln_dna_sequence, target_aln_dna_sequence) / len(query_aln_dna_sequence))
                prot_consev = 1 - (self.hd_amino_seqs(query_aln_prot_sequence, target_aln_prot_sequence) / len(query_aln_prot_sequence))
                query_dna_sequence = self.ungap_sequence(query_aln_dna_sequence)
                target_dna_sequence = self.ungap_sequence(target_aln_dna_sequence)
                query_prot_sequence = self.ungap_sequence(query_aln_prot_sequence)
                target_prot_sequence = self.ungap_sequence(target_aln_prot_sequence)
                query_aligned_fraction = len(query_dna_sequence)/query_length
                update_string = """
                UPDATE 
                concatenate_fragments
                SET 
                dna_perc_conservation = ?,
                prot_perc_conservation = ?,
                query_dna_sequence = ?,
                target_dna_sequence = ?,
                query_prot_sequence = ?,
                target_prot_sequence = ?,
                query_aligned_fraction = ?
                WHERE 
                (gene_id = ? AND query_id = ? AND target_id = ? AND HSP_n = ? AND score = ?)
                """
                cursor_insert.execute(update_string, (dna_consev, prot_consev,
                                                      query_dna_sequence, target_dna_sequence,
                                                      query_prot_sequence, target_prot_sequence,
                                                      query_aligned_fraction,
                                                      gene_id, query_id, target_id,
                                                      HSP_n, score))
            conn.commit()
            conn.close()
