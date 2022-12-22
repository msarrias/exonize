import gffutils # for creating/loading DBs
import subprocess # for calling gffread
import os

class DataBaseOp(object):
    def __init__(self,
                 db_path,
                 in_file_path,
                 verbose
                ):
        self.verbose = verbose
        self.db = None
        self.db_path = db_path
        
        
    def create_parse_or_update_database(self):
        if not os.path.exists(self.db_path):
            self.in_file_path = in_file_path
            if 'gtf' in self.in_file_path:
                self.old_filename = self.in_file_path
                self.in_file_path = f"{self.old_filename.rsplit('.gtf')[0]}.gff"
                self.convert_gtf_to_gff()
                print('the GTF file has been converted into a GFF3 file')
                print(f'with filename: {self.in_file_path}')
            self.create_database()
        if not self.db:
            if self.verbose: print("Reading annotations database:", end = " ")
            self.load_db()
            if self.verbose: print("Done!")
        self.db_features = list(self.db.featuretypes())
        self.create_intron_annotations()
            
    
    def convert_gtf_to_gff(self):
        gffread_command = [ "gffread",
                           self.old_filename, 
                           "-O",
                           "-o", 
                           self.in_file_path
                          ]
        subprocess.call(gffread_command)
        
        
    def create_database(self):
        try:
            if self.verbose: print("Creating annotations database", end=" ")
            self.db = gffutils.create_db(self.in_file_path,
                                         dbfn=self.db_path,
                                         force=True, 
                                         keep_order=True,
                                         merge_strategy='create_unique',
                                         sort_attribute_values=True,
                                         disable_infer_genes=True,
                                         disable_infer_transcripts=True)
            if self.verbose: print("Done!")
        except ValueError:
            print("Wrong infile path")

            
    def load_db(self):
        try:
            self.db = gffutils.FeatureDB(self.db_path,
                                         keep_order=True)
        except ValueError:
            print("Wrong db file path")
        
        
    def create_intron_annotations(self):
        """
        run only once otherwise duplicated 
        annotations will be created
        """
        if 'intron' not in self.db_features:
            if self.verbose: print(f"Writing intron annotations in database:", end=" ")
            introns = list(self.db.create_introns())
            self.db.update(introns, 
                           id_spec={'intron': [lambda f : ','.join(f['ID'])]},
                           make_backup=False
                          )
            if self.verbose: print("Done!")
                
        
        