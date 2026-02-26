"""
cs.py - A module for entity search using pre-indexed Parquet files and Dask.

The module provides functions to:
  1. Load delimiters from a CSV file.
  2. Group query entity IDs by file and block (storing block boundaries).
  3. Process each file/block by reading the corresponding Parquet file with Dask,
     filtering for the block, and using binary search (on a Pandas DataFrame)
     to locate each query.
"""

import os
import time
import bisect
import pandas as pd
import dask.dataframe as dd

def load_delimiters(delimiters_path,entity):
    """
    Load delimiters from a CSV file.
    
    The CSV must have the following columns:
      - filename
      - block_number (integer)
      - first_entity_id
      - last_entity_id

    Returns:
      tuple: (delimiters, first_ids)
        - delimiters: a list of dictionaries (one per block), sorted by first_entity_id.
        - first_ids: a list of first_entity_id strings for binary search.
    """
    df = pd.read_csv(delimiters_path)
    df['block_number'] = df['block_number'].astype(int)
    df.sort_values(f'first_{entity}_id', inplace=True)
    delimiters = df.to_dict(orient='records')
    first_ids = [rec[f'first_{entity}_id'] for rec in delimiters]
    return delimiters, first_ids

def find_block_for_query(query_id, delimiters, first_ids, entity):
    """
    Find the delimiter record corresponding to a given query_id.

    Args:
      query_id (str): The entity id to search for.
      delimiters (list): List of delimiter records (dicts).
      first_ids (list): List of first_entity_id values from delimiters.

    Returns:
      dict or None: The delimiter record (with filename, block_number,
                    first_entity_id, last_entity_id) if query_id is within range,
                    otherwise None.
    """
    index = bisect.bisect_right(first_ids, query_id) - 1
    if index >= 0:
        rec = delimiters[index]
        if rec[f'first_{entity}_id'] <= query_id <= rec[f'last_{entity}_id']:
            return rec
    return None

def group_queries_by_block(queries, delimiters, first_ids, entity):
    """
    Group queries by file and block, storing for each block its boundary values.

    Args:
      queries (iterable): List/array of entity id strings.
      delimiters (list): List of delimiter records.
      first_ids (list): List of first_entity_id values.

    Returns:
      tuple: (queries_by_file, results)
        queries_by_file is a dict structured as:
          { filename: { block_number: {'queries': [q1, q2, ...],
                                        'low_bound': first_entity_id,
                                        'high_bound': last_entity_id } } }
        results is a dict mapping each query to None (to be updated later).
    """
    queries_by_file = {}
    results = {}
    for query in queries:
        block_info = find_block_for_query(query, delimiters, first_ids, entity)
        if block_info is None:
            results[query] = None
            continue
        filename = block_info['filename']
        block_number = block_info['block_number']
        low_bound = block_info[f'first_{entity}_id']
        high_bound = block_info[f'last_{entity}_id']
        if filename not in queries_by_file:
            queries_by_file[filename] = {}
        if block_number not in queries_by_file[filename]:
            queries_by_file[filename][block_number] = {
                'queries': [],
                'low_bound': low_bound,
                'high_bound': high_bound
            }
        queries_by_file[filename][block_number]['queries'].append(query)
        results[query] = None
    return queries_by_file, results

def binary_search_block(pdf, query_id, entity):
    """
    Perform a binary search for query_id in a sorted Pandas DataFrame.

    Args:
      pdf (pandas.DataFrame): DataFrame with columns 'entity_id' and 'data'. It must be sorted by 'entity_id'.
      query_id (str): The entity id to search for.

    Returns:
      The associated data if found; otherwise, None.
    """
    entity_ids = pdf[f'{entity}_id'].tolist()
    data_list = pdf['data'].tolist()
    lo, hi = 0, len(entity_ids) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        if entity_ids[mid] == query_id:
            return data_list[mid]
        elif entity_ids[mid] < query_id:
            lo = mid + 1
        else:
            hi = mid - 1
    return None

def process_file_blocks(parquet_base_path, filename, blocks, results, entity):
    """
    For each block in the given file, read the corresponding Parquet file using Dask,
    filter for rows within the block's entity_id boundaries, convert to a Pandas DataFrame,
    and perform a binary search for each query in that block.

    Args:
      parquet_base_path (str): Base folder where the Parquet files are stored.
      filename (str): The original file name (e.g., "mgy_entity_map_1.tsv.gz").
      blocks (dict): A dictionary mapping block_number to a dict with keys 'queries', 'low_bound', and 'high_bound'.
      results (dict): A dictionary mapping query id to its result (updated in-place).
    """
    # Convert the filename to a corresponding Parquet file name.
    parquet_file = os.path.join(parquet_base_path, filename.replace('.tsv.gz', '.parquet'))
    
    # Read the entire Parquet file with Dask.
    ddf = dd.read_parquet(parquet_file)
    
    for block_number, block_data in blocks.items():
        low_bound = block_data['low_bound']
        high_bound = block_data['high_bound']
        query_list = block_data['queries']
        ql = ', '.join(query_list)
        
        # Filter Dask DataFrame to rows within this block's boundaries.
        filtered_ddf = ddf[(ddf[f'{entity}_id'] >= low_bound) & (ddf[f'{entity}_id'] <= high_bound)]
        
        # Compute the filtered Dask DataFrame into Pandas.
        pdf = filtered_ddf.compute()
        pdf.sort_values(f'{entity}_id', inplace=True)
        
        # For each query in this block, perform binary search and update the results.
        for query in query_list:
            result = binary_search_block(pdf, query, entity)
            results[query] = result
