import os
import sys
import glob
import time
import datetime

# import database handlers
from db_handler_assemblies import AssembliesDBHandler
from db_handler_sequences import SequenceDBHandler

from gcsnap.utils import processpool_wrapper
import argparse


def split_into_batches(lst: list, batch_size: int = 1000):
    for i in range(0, len(lst), batch_size):
        # yield produces a generator
        yield lst[i:i + batch_size]
        
        
def split_into_parts(lst: list, n: int) -> list[list]:
    # k, how many times n fits in len(list), m the reminder
    q, r = divmod(len(lst), n)
    return [lst[i * q + min(i, r):(i + 1) * q + min(i + 1, r)] for i in range(n)]
    

def execute_handler(args: tuple):
    handler, batch = args
    # each handler returns a list of tuples
    sequences, mapping = handler.parse_sequences_from_faa_files(batch)
    
    return sequences, mapping


def execute_handlers(handler: SequenceDBHandler, batch: list, n_processes: int) -> tuple[list[tuple[str,str]], list[tuple[str,str]]]:
    # split the batch in equal sized parts
    batches = split_into_parts(batch, n_processes) 
    parallel_args = [(handler, subbatch) for subbatch in batches]
    
    # retrurns a list of tuples containing lists of tuples
    result = processpool_wrapper(n_processes, parallel_args, execute_handler)
    
    # flatten the list of tuples containing lists of tuples
    mappings = [item for sublist in result for item in sublist[1]]
    sequences = [item for sublist in result for item in sublist[0]]
    
    return (sequences, mappings)


def reindex(handler) -> None:
    handler.reindex()


def print_assemblies(assemblies: list[str]) -> None:
    for assembly in assemblies:
        assembly_accession = '_'.join(os.path.basename(assembly).split('_')[:2])
        print('Assembly {} done'.format(assembly_accession))
  

def create_dbs(path: str, n_processes: int, sources: list[str]) -> None:
    
    # where to create the databases
    db_dir = os.path.join(path, 'db')        
    
    if not os.path.exists(db_dir):
        os.makedirs(db_dir)

    # open assembly database handler and create tables
    asse_db_handler = AssembliesDBHandler(db_dir, 'assemblies.db')
    asse_db_handler.create_tables()
    seq_db_handler = SequenceDBHandler(db_dir, 'sequences.db')    
    seq_db_handler.create_table()
    
    print('Assembly database created: {}'.format(os.path.join(db_dir, 'assemblies.db')))
    print('Sequence database created: {}'.format(os.path.join(db_dir, 'sequences.db')))
    # number of databases to write to in parallel
    batch_size = n_processes * 500
    # keep track of sequences and assemblies
    n_assemblies = 0 
    n_sequences = 0

    for loop_var in sources:

        ML = []
        db_type = loop_var
                
        # 1. fill assemblies table from assembly_summary_{db_type}.txt
        # (assembly_accession, url, taxid, species)
        summary_file = os.path.join(path, db_type, 'assembly_summary_{}.txt'.format(db_type))
        asse_db_handler.insert_assemblies_from_summary(summary_file)
        print('Summary {} done'.format(os.path.basename(summary_file)))

        # 2. Extract information from .faa files       
        # a) (ncbi_code, assembly_accession and seq_db) go into the mapping table of assemblies.db
        # b) (ncbi_code, sequence) go into sepearte databases (as this will get huge)
            # We use new databses for next file when reaching the limit
        # Directory containing the .faa files
        faa_data_dir = os.path.join(path, db_type, 'data')
        # list of all files to parse
        # file_paths = [os.path.join(faa_data_dir, file_name) for file_name 
        #               in os.listdir(faa_data_dir)]        
        file_paths = glob.glob(os.path.join(faa_data_dir,'*_protein.faa.gz'))
        # file_paths = glob.glob(os.path.join(faa_data_dir,'*_protein.faa'))
                
        # we loop over all those files in batches, each database takes 
        # indexing is switched off to speed up
        for batch in split_into_batches(file_paths, batch_size):
            
            # start the parsing for each handler, in parallel 
            sequence_list, mapping_list = execute_handlers(seq_db_handler, batch, n_processes)
            # keep track of done things
            n_sequences += len(mapping_list)
            n_assemblies += len(batch)
                    
            # add sequences
            seq_db_handler.insert_sequences(sequence_list)
            
            # add mappings to assembly db
            asse_db_handler.insert_mappings(mapping_list)           
                        
            # Format numbers with thousand separators
            formatted_assemblies = "{:,}".format(n_assemblies)
            formatted_sequences = "{:,}".format(n_sequences)
            
            #print_assemblies(batch)
            print('{} assemblies and {} sequences done so far'.format(formatted_assemblies, formatted_sequences))
            
    parallel_args = [asse_db_handler, seq_db_handler]
    processpool_wrapper(2, parallel_args, reindex)
    
    print('Assembly indexings done')        
    print('Sequence indexings done')
    
    # get number of unique sequences
    n_unique = seq_db_handler.select_number_of_entries()

    return n_assemblies, n_sequences, n_unique                   

# Example usage

if __name__ == "__main__":    
   
    parser = argparse.ArgumentParser(description="Create assemblies and sequences databases.")
    parser.add_argument( "-p", "--processes", type=int, required=True, help="Number of parallel processes to use." )
    parser.add_argument( "--path", type=str, required=True, help="Path to the database directory." )
    parser.add_argument( "--sources", nargs='+', type=str, required=True, help="Source databases (one or more)." )
    args = parser.parse_args()

    n_processes = args.processes
    path = args.path
            
    st = time.time()
    n_assemblies, n_sequences, n_unique = create_dbs(path, n_processes, args.sources)
    
    elapsed_time = time.time() - st
    formatted_time = str(datetime.timedelta(seconds=round(elapsed_time)))

    # Format numbers with thousand separators
    formatted_assemblies = "{:,}".format(n_assemblies)
    formatted_sequences = "{:,}".format(n_sequences)
    formatted_unique = "{:,}".format(n_unique)
    
    print('{} assemblies with {} sequences ({} unique sequences found) done in {}'.
          format(formatted_assemblies, formatted_sequences, formatted_unique, formatted_time))
  
    # some selecting
    # ncbi_codes = ['AAK02085.1','AAK02086.1','AAC70070.1','AAC70069.1']
    # res = db_handler.select(ncbi_codes)
    
    # # select all assemblies
    # res = db_handler.select_all('assemblies')
    
    # # select all sequences (DEFAULT)
    # res = db_handler.select_all()
    
    # # Retrieve a sample sequence
    # sequence = db_handler.get_sequence('AAV61620.1')
    # print(sequence)