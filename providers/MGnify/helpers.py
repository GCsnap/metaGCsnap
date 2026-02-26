from itertools import islice

##########################################################################
############################################################### parameters

# MMSEQS PARAMETERS
class MMSeqsParams:
    #threads = 32
    split = 4
    #min_seq_id = 0.15
    #cov_mode = 1
    #coverage = 0.85
    #e_value = 1e-5
    #sensitivity = 5

    tsv_columns = [
    "query",           # 0 - query sequence name
    "target",          # 1 - target sequence name
    "pident",          # 2 - percent identity
    "alnlen",          # 3 - alignment length
    "mismatch",        # 4 - number of mismatches
    "gapopen",         # 5 - number of gap openings
    "qstart",          # 6 - start of alignment in query
    "qend",            # 7 - end of alignment in query
    "tstart",          # 8 - start of alignment in target
    "tend",            # 9 - end of alignment in target
    "evalue",          # 10 - E-value of the alignment
    ]

mmseqs = MMSeqsParams()

##########################################################################
####################################################### files organization

# DEFAULT TMP DIRECTORY
tmp= '/tmp'

# DEFAULT PATHS INSIDE THE MGNIFY DATABASE

# database for MGnify clusters
mgyp_database = 'mmseqsDBs/mgyc/DB' # final
# database for MGnify clusters

taxonomy_database = 'GTDB/GTDB' # final
#mgyp_database='mmseqsDBs/mgy_clusters_10M/DB'# test

# DEFAULT PATHS INSIDE OUTPUT DIRECTORY

# hits for sequence search on mgnify clusters
mgyp_search_out = 'MGnify_mgyp_hits.tsv'
mgyp_search_ids = 'MGnify_mgyp_hits.txt'
mgyp_search_fasta = 'MGnify_mgyp_hits.fasta'
shiflanking_search_out = 'CG_mgyp_hits.tsv'
# results of taxonomic assignment of metagenomics contigs
taxonomic_search_out = 'taxonomic_assignments.tsv'

query_db = 'mmseq/query/DB' 

# mgyp metadata file
mgyp_metadata = 'mgyp_metadata.csv'

##########################################################################
################################################################ constants

# MGNIFY PIPELINE VERSIONS INCLUDED IN THE WORKFLOW
PIPELINES = ['4.1','5.0']

# tuple dict separator
TSEP='_,_'

##########################################################################
################################################################# scraping

# MGNIFY API URL
MGnify_API = "https://www.ebi.ac.uk/metagenomics/api/v1"

# pipeline, v5 targets
pipelineV5_targets = [ 'Processed contigs','Predicted CDS (aa)','Complete GO annotation','InterPro matches']
pipelineV4_targets = [ 'Processed contigs','Predicted CDS without annotation','Predicted CDS with annotation','Complete GO annotation','InterPro matches']

contigs_description = [ 'Processed contigs']
cds_description = [ 'Predicted CDS (aa)','Predicted CDS without annotation','Predicted CDS with annotation']
fannotation_description = ['InterPro matches'] # , 'Complete GO annotation' ]

##########################################################################
########################################################### tiny functions

def chunked_iterable(iterable, size):
    it = iter(iterable)
    while True:
        batch = list(islice(it, size))
        if not batch:
            break
        yield batch

##########################################################################
################################################################# taxonomy

def tax_ranks(binning: str = 'binning') -> tuple:
    
    b = (binning == 'binning')

    if b:
        return ['root','metagenomic bin']
    else:
        return ['root','domain','phylum','class','order','family','genus','species','strain']

def tax_ranks_dict(binning: str = 'binning') -> dict:
    
    b = (binning == 'binning')

    if b:
        return {'r': 'root', 'm': 'metagenomic bin'}
    else:
        return {'r': 'root',
                'd': 'domain',
                'p': 'phylum',
                'c': 'class',
                'o': 'order',
                'f': 'family',
                'g': 'genus',
                's': 'species'
                }