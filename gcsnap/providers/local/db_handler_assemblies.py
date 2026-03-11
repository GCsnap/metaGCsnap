import os
import sqlite3

class AssembliesDBHandler:
    def __init__(self, db_path: str ,db_name: str = 'assemblies.db'):
        self.db = os.path.join(db_path, db_name)
        self.db_name = db_name
        
    def create_tables(self) -> None:
        self.create_assembly_table()
        self.create_mapping_table()
        self.disable_indices()
    
    def create_mapping_table(self) -> None:
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute('DROP TABLE IF EXISTS mappings')
        cursor.execute('''
            CREATE TABLE mappings (
                ncbi_code TEXT PRIMARY KEY,
                assembly_accession TEXT
            )
        ''')
        conn.commit()
        conn.close()
        
    def create_assembly_table(self) -> None:
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute('DROP TABLE IF EXISTS assemblies')
        cursor.execute('''
            CREATE TABLE assemblies (
                assembly_accession TEXT PRIMARY KEY,
                url TEXT,
                taxid INTEGER,                 
                species TEXT
            )
        ''')
        conn.commit()
        conn.close()    
        
    def disable_indices(self) -> None:
        conn = sqlite3.connect(self.db)
        # Disabling indices and other performance-related settings
        conn.execute('PRAGMA synchronous = OFF')
        conn.execute('PRAGMA journal_mode = OFF')
        conn.execute('PRAGMA temp_store = MEMORY')
        conn.execute('PRAGMA cache_size = -1048576')  # 1 GB cache size (1,048,576 KB) 
        conn.commit()
        conn.close()
        
    def enable_indices(self) -> None:
        conn = sqlite3.connect(self.db)
        # Re-enabling indices and other settings
        conn.execute('PRAGMA synchronous = NORMAL')
        conn.execute('PRAGMA journal_mode = WAL')
        conn.execute('PRAGMA temp_store = DEFAULT')        
        conn.execute('PRAGMA cache_size = -2000')  # Reset to default cache size  
        conn.commit()
        conn.close()  
        
    def reindex(self) -> None:
        self.reindex_mappings()
                
    def reindex_mappings(self) -> None:
        self.enable_indices()        
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute('REINDEX mappings')
        conn.commit()
        conn.close()    
        
    def reindex_assemblies(self) -> None:
        self.enable_indices()
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute('REINDEX assemblies')
        conn.commit()
        conn.close()         

    def batch_insert_mappings(self, mappings: list[tuple[str,str,str]]) -> None:
        self.disable_indices()
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.executemany('INSERT OR REPLACE INTO mappings (ncbi_code, assembly_accession) VALUES (?, ?)', mappings)
        conn.commit()
        conn.close()
            
    def batch_insert_assemblies(self, assemblies: list[tuple[str,str,int,str]]) -> None:
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.executemany('INSERT OR REPLACE INTO assemblies (assembly_accession, url , taxid, species) VALUES (?, ?, ?, ?)', assemblies)
        conn.commit()
        conn.close()                       
        
    def parse_and_insert_assemblies(self, file_path: str) -> None:     
        with open(file_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()
            
        assembly_list = []    
        
        for line in lines:
            if not line.startswith('#'):
                data = line.strip().split('\t')
                # the first column has the Assemlby id, column 19 the url to it
                assembly_accession = data[0]
                taxid = data[1]
                species = data[2]
                url = data[3]
                
                assembly_list.append((assembly_accession, url, taxid, species))        

                if len(assembly_list) >= 100000:  # Adjust batch size as needed
                    self.batch_insert_assemblies(assembly_list)
                    assembly_list = []

        # write the last part not complete 1000 assemblies
        if assembly_list:
            self.batch_insert_assemblies(assembly_list)     
            
        # index assemblies table and set PRAGMA back to fast wrtiting
        self.reindex_assemblies()
        self.disable_indices()
        
    def insert_assemblies_from_summary(self, file_path: str) -> None:
        self.parse_and_insert_assemblies(file_path)
        
    def insert_mappings(self, mappings: list[tuple[str,str,str]]) -> None:
        self.batch_insert_mappings(mappings)

    def select(self, codes: list[str], return_fields: list[str] = None, table: str = 'mappings', request_size: int = 5000) -> list[tuple]:
        # split the codes into chunks of size request_size, this is only of effect if the workload for each process is large
        # for instance with 1M sequences and just 10 MPI ranks, each rank has to process 100k sequences which is not efficient for DB requests
        # For smaller workloads, the loop reduces to just one iteration
        chunks = [codes[i:i + request_size] for i in range(0, len(codes), request_size)]
        
        # get the results for each chunk
        results = []
        for codes_chunk in chunks:
            results += self.select_chunk(codes_chunk, return_fields, table)
        
        return results
                
    def select_chunk(self, codes: list[str], return_fields: list[str] = None, table: str = 'mappings') -> list[tuple]:
        # combine query
        
        # which fields to return
        if return_fields is None:
            select_fields = '*'
        else:
            select_fields = ', '.join(return_fields)

        # determin correct colum
        if table == 'mappings':
            column = 'ncbi_code'
        elif table == 'assemblies':
            column = 'assembly_accession'
        else:
            raise ValueError(f'Unknown table {table}')

        # which records based on ncbi_codes
        # the query contains ?,?,?,?, with the cursor.execute(query, codes) the values are inserted
        records = f"({','.join(['?']*len(codes))})" # (?,?,?,?) for each entry in ncbi_code       
        query = 'SELECT {} FROM {} WHERE {} IN {}'.format(select_fields,table,column,records)
        
        # execute query
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute(query, codes)
        result = cursor.fetchall()  # fetchall() gets all, fetchone() just the next in the result list
        conn.close()
        
        return result
    
    def select_all(self, table: str = 'mappings') -> list[tuple]:
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute(f'SELECT * FROM {table}')
        result = cursor.fetchall()  # fetchall() gets all, fetchone() just the next in the result list
        conn.close()
        
        return result
    
    def select_number_of_entries(self, table: str = 'mappings') -> int:
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute('SELECT COUNT(*) FROM {}'.format(table))
        result = cursor.fetchone()  # fetchall() gets all, fetchone() just the next in the result list
        conn.close()   
        
        # extract the count from the tuple
        return result[0]   
