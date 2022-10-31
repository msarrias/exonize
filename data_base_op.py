import gffutils # for creating/loading DBs
import subprocess # for calling gffread

class DataBaseOp(object):
    def __init__(self,
                 db_path,
                 in_file_path='',
                 create_db=False,
                 create_introns=False):
        self.db = None
        self.db_path = db_path
        self.create_introns = create_introns
        if create_db:
            self.in_file_path = in_file_path
            if 'gtf' in self.in_file_path:
                self.old_filename = self.in_file_path
                self.in_file_path = (self.old_filename.rsplit('.gtf')[0]
                                     + '.gff')
                self.convert_gtf_to_gff()
                print('the GTF file has been converted into a GFF3 file')
                print(f'with filename: {self.in_file_path}')
            self.create_database()
            if self.create_introns:
                self.load_db()
                self.get_introns()
               
            
    def convert_gtf_to_gff(self):
        gffread_command = [ "gffread",
                           self.old_filename, 
                           "-O", "-o", self.in_file_path
                          ]
        subprocess.call(gffread_command)
        
        
    def create_database(self):
        try:
            gffutils.create_db(self.in_file_path,
                               dbfn=self.db_path,
                               force = True, keep_order=True,
                               merge_strategy='merge',
                               sort_attribute_values=True,
                               )
        except ValueError:
            print("Wrong infile path")

            
    def load_db(self):
        try:
            self.db = gffutils.FeatureDB(self.db_path,
                                     keep_order=True)
        except ValueError:
            print("Wrong db file path")
        
        
    def get_introns(self):
        """run only once otherwise duplicated 
        annotations will be created"""
        introns = list(self.db.create_introns())
        self.db.update(introns, keep_order=True)
        
        