import os
import gzip
import sqlite3

class SequenceDBHandler:
    def __init__(self, db_path: str , db_name: str = 'sequences.db'):
        self.db = os.path.join(db_path, db_name)
        self.db_name = db_name
        
    def create_table(self) -> None:
        self.create_sequence_table()
        self.disable_indices()
    
    def create_sequence_table(self) -> None:
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute('DROP TABLE IF EXISTS sequences')
        cursor.execute('''
            CREATE TABLE sequences (
                ncbi_code TEXT PRIMARY KEY,
                sequence TEXT
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
        self.reindex_sequences()        
        
    def reindex_sequences(self) -> None:
        self.enable_indices()
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute('REINDEX sequences')
        conn.commit()
        conn.close()         

    def batch_insert_sequences(self, sequences: list[tuple[str,str]]) -> None:
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()          
        cursor.executemany('INSERT OR REPLACE INTO sequences (ncbi_code, sequence) VALUES (?, ?)', sequences)
        conn.commit()
        conn.close()         

    def batch_update_sequences(self, sequences: list[tuple[str,str]]) -> None:
        # we need (new_sequence, ncbi_code)   
        updates = [(sequence[1], sequence[0]) for sequence in sequences]        
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()    
        # Disable journaling and set performance-related settings in the same connection
        # Perform the bulk update using executemany
        sql = 'UPDATE sequences SET sequence = ? WHERE ncbi_code = ?'
        cursor.executemany(sql, updates)        
        conn.commit()
        conn.close()         

    def read_gzip_file(self, file_path: str) -> str:
        with gzip.open(file_path, 'rt', encoding='utf-8') as file:
            content = file.read()
        return content.splitlines()

    def read_txt_file(self, file_path: str) -> str:
        with open(file_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()
        return lines
            
    def parse_sequences(self, file_paths: list[str]) -> tuple[list[tuple[str,str]],list[tuple[str,str]]]:
        sequence_list = []  
        mapping_list = []
        for file_path in file_paths:
            # accession taken from file name: GCF_000247695.1_HetGla_female_1.0_protein.faa.gz
            # all but the last are the accession
            assembly_accession = '_'.join(os.path.basename(file_path).split('_')[:-1])

            # read the file in once
            if file_path.endswith('.gz'):
                lines = self.read_gzip_file(file_path)
            else:
                lines = self.read_txt_file(file_path)
            
            # parse the lines, the challange, the sequence can be several lines long
            # make one string of it with a splitabe character
            content = '$%'.join(lines)
            #ORIGINAL WORNG: content = ''.join(lines) --> can't be splitted anymore and the sequence is empty
            # split that string to extract each sequence id
            # EFB12766.1 hypothetical protein PANDA_022614, partial [Ailuropoda melanoleuca]WSDGHLIYYDDQTRQSVEDKVHMPVDCINIRTGHECRGT
            # the first is an empty result
            entries = content.split('>')[1:]
            
            for entry in entries:
                # split the info from the acutal sequence str
                entry_split = entry.split('$%')
                sequence = ''.join(entry_split[1:])
                # split the organism name in []
                info_split = entry_split[0].split('[') 
                # orgname is same for the assembly, put it into assembly table
                #orgname = info_split[-1].replace(']','')
                ncbi_code = info_split[0].split(' ')[0]
                
                sequence_list.append((ncbi_code, sequence))
                mapping_list.append((ncbi_code, assembly_accession))   
        
        # return what is needed from the .faa file
        return (sequence_list, mapping_list)      
    
    def insert_sequences(self, sequence_list: list) -> None:
        # insert all as batch
        self.batch_insert_sequences(sequence_list)  

    def parse_sequences_from_faa_files(self, file_paths: list[str]) -> list[tuple[str,str]]:
        # new use to do this in parallel, as the database writing is the bottleneck
        # SQLite does not support parallel write, so just batches to reduce write calls          
        return self.parse_sequences(file_paths) 
                    
    def select(self, ncbi_codes: list[str], return_fields: list[str] = None) -> list[tuple]:
        # combine query
        
        # which fields to return
        if return_fields is None:
            select_fields = '*'
        else:
            select_fields = ', '.join(return_fields)

        # which records based on ncbi_codes
        records = f"({','.join(['?']*len(ncbi_codes))})" # (?,?,?,?) for each entry in ncbi_code       
        query = 'SELECT {} FROM sequences WHERE ncbi_code IN {}'.format(select_fields, records)
        
        # execute query
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute(query, ncbi_codes)
        result = cursor.fetchall()  # fetchall() gets all, fetchone() just the next in the result list
        conn.close()  
        
        return result
    
    def select_all(self) -> list[tuple]:
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute('SELECT * FROM sequences')
        result = cursor.fetchall()  # fetchall() gets all, fetchone() just the next in the result list
        conn.close()   
        
        return result
    
    def select_as_dict(self, ncbi_codes: list[str], return_fields: list[str] = None) -> dict:
        # get the result as a dictionary
        result = self.select(ncbi_codes, return_fields)
        return {record[0]: record[1] for record in result}

    def select_number_of_entries(self) -> int:
        conn = sqlite3.connect(self.db)
        cursor = conn.cursor()
        cursor.execute('SELECT COUNT(*) FROM sequences')
        result = cursor.fetchone()  # fetchall() gets all, fetchone() just the next in the result list
        conn.close()   
        
        # extract the count from the tuple
        return result[0]        
    
    def get_db_size(self) -> int:
        return os.path.getsize(self.db)