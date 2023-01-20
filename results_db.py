import sqlite3


class SQlite_res():
    def __init__(self, db_path):
        self.db_path = db_path
        self.colnames = []
        self.create_table_concatenate_fragments()
        self.get_concatenate_fragments_col_names()
        self.insert_conservation_percentage_columns()
          
            
    def get_tables_name(self):
        con = sqlite3.connect(self.db_path)
        cursor = con.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        return cursor.fetchall()
                              
        
    def get_concatenate_fragments_col_names(self):
        conn = sqlite3.connect(self.db_path)
        info = conn.execute("PRAGMA table_info('concatenate_fragments')").fetchall()
        self.concatenate_fragments_colnames = [col_name[1] for col_name in info]
        conn.close()
      
    
    @staticmethod        
    def hd_dna_seqs(seq_a, seq_b):
        return sum([i != j for i, j in zip(seq_a, seq_b)])

    
    @staticmethod
    def hd_amino_seqs(seq_a, seq_b):
        return sum([seq_a[i:i+3] != seq_b[i:i+3] for i in range(0,len(seq_a), 3)])

    
    def create_table_concatenate_fragments(self):
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        new_table = """
        CREATE TABLE IF NOT EXISTS concatenate_fragments AS
        SELECT *,
        CAST(query_ungapped_dna_seq_length as REAL)/query_length AS query_aligned_fraction
        FROM(
        SELECT *, 
        LENGTH(query_ungapped_dna_sequence) AS query_ungapped_dna_seq_length
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
        MAX(category) AS category,
        MAX(is_pseudo_exon) AS is_pseudo_exon,
        GROUP_CONCAT(frag_query_start, ';') AS frag_query_start,
        GROUP_CONCAT(frag_query_end, ';') AS frag_query_end,
        GROUP_CONCAT(frag_target_start, ';') AS frag_hit_start,
        GROUP_CONCAT(frag_target_end, ';') AS frag_hit_end,
        GROUP_CONCAT(query_dna_sequence, '') AS query_dna_sequence, 
        replace(replace(replace(replace(replace(query_dna_sequence,'<',''),
        '>',''), '*',''), '-',''), '#', '') 
        AS query_ungapped_dna_sequence,
        GROUP_CONCAT(target_dna_sequence, '') AS target_dna_sequence,
        GROUP_CONCAT(query_prot_sequence, '') AS query_prot_sequence,
        GROUP_CONCAT(target_prot_sequence, '') AS target_prot_sequence
        FROM (
        SELECT * 
        FROM Fragments
        ORDER BY 
        gene_id ASC,
        query_id ASC,
        target_id ASC,
        query_start ASC,
        query_end ASC)
        GROUP BY 
        gene_id, 
        query_id, 
        target_id));
        """
        #create table
        cursor.execute(new_table)
        #create index
        cursor.execute("""CREATE INDEX concat_frag_idx ON concatenate_fragments (gene_id, query_id, target_id);""")
        cursor.execute("""CREATE INDEX concat_frag_gene_idx ON concatenate_fragments (gene_id);""")
        conn.close()
        
        
    def insert_conservation_percentage_columns(self):
        conn = sqlite3.connect(self.db_path)
        cur = conn.cursor()
        if 'dna_perc_conservation' and 'prot_perc_conservation' not in self.concatenate_fragments_colnames:
            cur.execute("ALTER TABLE concatenate_fragments ADD COLUMN dna_perc_conservation REAL;")
            cur.execute("ALTER TABLE concatenate_fragments ADD COLUMN prot_perc_conservation REAL;")
            #populate column
            cur.execute("""
            SELECT 
            gene_id,
            query_id,
            target_id,
            query_dna_sequence,
            target_dna_sequence,
            query_prot_sequence,
            target_prot_sequence
            FROM
            concatenate_fragments
            """)
            for row in cur:
                cursor_insert = conn.cursor()
                (gene_id, 
                query_id, 
                target_id,
                query_dna_sequence,
                target_dna_sequence,
                query_prot_sequence,
                target_prot_sequence) = row
                dna_consev = 1 - (self.hd_dna_seqs(query_dna_sequence, target_dna_sequence) / len(query_dna_sequence) )
                prot_consev = 1 - (self.hd_amino_seqs(query_prot_sequence, target_prot_sequence) / len(query_prot_sequence))
                update_string = """
                UPDATE 
                concatenate_fragments
                SET 
                dna_perc_conservation = ?,
                prot_perc_conservation = ?
                WHERE 
                (gene_id = ? AND query_id = ? AND target_id = ?)
                """
                cursor_insert.execute(update_string, (dna_consev, prot_consev, gene_id, query_id, target_id))
            conn.commit()
            conn.close()
