import gffutils # for creating/loading DBs


class data_base_op(object):
    def __init__(self,
                 db_path,
                 in_file_path='',
                 create_db=False,
                 create_introns=False):
        self.db_path = db_path
        self.create_introns = create_introns
        if create_db:
            self.in_file_path = in_file_path
            self.create_database()
            if self.create_introns == True:
                self.load_db()
                self.get_introns()
                
                
    def create_database(self):
        try:
            gffutils.create_db(self.in_file_path,
                               dbfn=self.db_path,
                               force = True, keep_order=True,
                               merge_strategy='merge',
                               sort_attribute_values=True)
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
        
        