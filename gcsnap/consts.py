
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