import numpy as np
import pandas as pd
from collections import OrderedDict, defaultdict
from gcsnap.providers.MGnify.helpers import tax_ranks
import sys
from typing import List, Tuple

###################################################################
# Expand contig distance to cds distance matrix
###################################################################


def expand_contig_matrix_to_cds( cds_contig_pairs: List[Tuple[str, str]],  contig_distance_matrix: pd.DataFrame) -> pd.DataFrame:
    """
    Expands a contig-level distance matrix to a CDS-level distance matrix.

    This function takes a list of (contig_id, cds_id) tuples
    and a contig-contig distance matrix. It then "broadcasts" the
    contig-level distances to all their associated CDS IDs, resulting
    in a new, larger matrix indexed by CDS ID.

    Args:
        cds_contig_pairs: 
            A list of tuples, where each tuple is (contig_id, cds_id).
            Example: [('contig_A', 'cds_A1'), ('contig_A', 'cds_A2'), ...]
        contig_distance_matrix: 
            A pandas DataFrame with contig_ids as both index and columns, 
            containing the distances.

    Returns:
        A new pandas DataFrame with all cds_ids as both index and columns,
        populated with the distances of their parent contigs.
        
    Raises:
        KeyError: If one or more contigs from the `cds_contig_pairs`
                  are not found in the `contig_distance_matrix` index/columns.
    """
    
    # 1. Create the reverse mapping (CDS ID -> Contig ID)

    cds_to_contig = {}
    for contig_id, cds_id in cds_contig_pairs:
        # This check warns you if your metadata has a CDS ID
        # associated with more than one contig.
        if cds_id in cds_to_contig and cds_to_contig[cds_id] != contig_id:
            print(f"Warning: CDS ID {cds_id} is in multiple contigs. Using last found: {contig_id}")
        cds_to_contig[cds_id] = contig_id

    # 2. Get the complete, sorted list of all relevant CDS IDs
    all_cds_ids = sorted(cds_to_contig.keys())
    num_cds_ids = len(all_cds_ids)

    # 3. Handle edge case
    if num_cds_ids == 0:
        print("No CDS IDs found. Returning an empty DataFrame.")
        return pd.DataFrame()

    # 4. Create a helper Series for mapping
    cds_to_contig_series = pd.Series(cds_to_contig)

    # 5. Get the corresponding contig for each CDS ID, in the sorted order
    contig_labels_for_indexing = cds_to_contig_series.loc[all_cds_ids].values

    # 6. Robustness Check: Ensure all contigs are in the matrix
    contigs_in_map = set(contig_labels_for_indexing)
    contigs_in_matrix = set(contig_distance_matrix.index)
    
    if not contigs_in_map.issubset(contigs_in_matrix):
        missing_contigs = contigs_in_map - contigs_in_matrix
        raise KeyError(
            f"{len(missing_contigs)} contigs are in the metadata but "
            f"missing from the distance matrix. Missing: {missing_contigs}"
        )

    # 7. Build the new matrix using efficient pandas .loc indexing
    temp_matrix = contig_distance_matrix.loc[contig_labels_for_indexing, contig_labels_for_indexing]

    # 8. Create the final DataFrame with the correct CDS ID labels
    adjusted_cds_distance_matrix = pd.DataFrame(
        temp_matrix.values,
        index=all_cds_ids,
        columns=all_cds_ids
    )
    
    return adjusted_cds_distance_matrix

###################################################################
# Function to parse the taxonomy tree and build the distance matrix
###################################################################

def get_target_taxonomy_paths(taxonomy_data,genome_classification):
    
    target_paths_map = {}
    ranks = list(tax_ranks(genome_classification))

    if not ranks:
        return target_paths_map

    root_rank = ranks[0]

    root_node = taxonomy_data.get(root_rank)

    if root_node is None:
        raise KeyError(f"Taxonomy data must contain a '{root_rank}' section.")

    _traverse_and_collect(root_node, [root_rank], 1, target_paths_map, ranks)
    return target_paths_map


def _traverse_and_collect(current_node, current_path, depth, target_paths_map, ranks):
    if 'target_members' in current_node:
        for target_id in current_node['target_members']:
            path_dict = {rank: '' for rank in ranks}
            for idx, name in enumerate(current_path):
                path_dict[ranks[idx]] = name
            target_paths_map.setdefault(target_id, path_dict)

    if depth >= len(ranks):
        return

    for child_name, child_node in current_node.items():
        if child_name in {'target_members', 'ncbi_codes'} or not isinstance(child_node, dict):
            continue
        _traverse_and_collect(child_node, current_path + [child_name], depth + 1, target_paths_map, ranks)
        
'''def get_target_taxonomy_paths(taxonomy_data):
    """
    Traverses a nested taxonomy dictionary and maps every target ID (leaf)
    to its full taxonomic path.

    Args:
        taxonomy_data (dict): The complete taxonomy dictionary loaded from JSON.

    Returns:
        dict: A dictionary where keys are target_ids and values are
              dictionaries of their taxonomic path.
              e.g., { 'target_1': {'superkingdom': 'Bacteria', 'phylum': 'B...', ...},
                      'target_2': ... }
    """
    
    # Define the ranks corresponding to the depth, as you specified.
    # Level 1 (e.g., 'Bacteria') is 'superkingdom'
    # Level 2 (e.g., 'Bacteroidota') is 'phylum'
    # ...and so on.
    
    
    # This will be our final dictionary: {target_id: path_dict}
    target_paths_map = {}

    # Get the top-level 'root' node
    root_node = taxonomy_data.get('root')
    if not root_node:
        print("Error: 'root' key not found in taxonomy data.", file=sys.stderr)
        return {}

    # Start the recursive search.
    # We iterate through the children of 'root' (e.g., 'Bacteria', 'Archaea'),
    # which are at 'superkingdom' level (depth 0 of our tax_ranks list).
    for node_name, sub_node in root_node.items():
        # Check if the item is a sub-tree (a dict) and not a data key
        if node_name not in ['target_members', 'ncbi_codes'] and isinstance(sub_node, dict):
            # Start the recursion, passing the first part of the path
            _traverse_and_collect(sub_node, [node_name], target_paths_map, tax_ranks)

    return target_paths_map
'''
'''def _traverse_and_collect(current_node, current_path_names, target_paths_map, ranks_list):
    """
    A recursive helper function to find all target_members below a node.

    Args:
        current_node (dict): The current node we are inspecting.
        current_path_names (list): The list of taxa names to get here
                                   (e.g., ['Bacteria', 'Bacteroidota']).
        target_paths_map (dict): The master dictionary to add results to.
        ranks_list (list): The list of rank names ['superkingdom', 'phylum', ...].
    """
    
    # 1. Check for targets at the current level
    if 'target_members' in current_node:
        for target_id in current_node['target_members']:
            
            # This is a leaf. Let's build its path dictionary.
            # Initialize with all ranks as empty strings
            path_dict = {rank: '' for rank in ranks_list}
            
            # Fill in the path names we have collected so far
            for i, name in enumerate(current_path_names):
                if i < len(ranks_list): # Safety check
                    rank_name = ranks_list[i]
                    path_dict[rank_name] = name
            
            # Save the complete path dict for this target ID
            # If a target appears in multiple places, this will
            # keep the path for the first one found.
            if target_id not in target_paths_map:
                target_paths_map[target_id] = path_dict

    # 2. Recurse into children nodes, if we haven't reached max depth
    current_depth_level = len(current_path_names)
    if current_depth_level >= len(ranks_list):
        return  # We are at the 'species' level, stop recursing

    for key, sub_node in current_node.items():
        # Check if this item is a child taxon (a dict, not a data key)
        if key not in ['target_members', 'ncbi_codes'] and isinstance(sub_node, dict):
            # Recurse deeper, adding the current key to the path
            _traverse_and_collect(sub_node, current_path_names + [key], target_paths_map, ranks_list)
'''

###################################################################
# Function to parse the taxonomy tree and build the distance matrix
###################################################################

def parse_taxonomy_tree(node, current_path, all_targets_list, target_to_path_map):
    """
    Recursively traverses the taxonomy JSON tree.
    
    This function builds two objects:
    1. all_targets_list: An ordered list of all unique target IDs.
    2. target_to_path_map: A dictionary mapping each target ID to its full
                           taxonomic path (a list of strings).
    """
    # Check if the current node has targets assigned to it
    if 'target_members' in node:
        for target_id in node['target_members']:
            # Only add a target the first time we see it
            # This handles cases where an ID might be in multiple lists
            if target_id not in target_to_path_map:
                all_targets_list.append(target_id)
                target_to_path_map[target_id] = current_path
    
    # Recurse into child nodes
    for key, sub_node in node.items():
        # Ignore the leaf-node keys and only look at nested dicts (taxa)
        if key not in ['target_members', 'ncbi_codes']:
            if isinstance(sub_node, dict):
                parse_taxonomy_tree(sub_node, current_path + [key], all_targets_list, target_to_path_map)

def get_taxonomic_distance(path_a, path_b):
    """
    Calculates the path distance between two taxonomic paths.
    
    Example:
    path_a = ['root', 'Bacteria', 'Phylum_A', 'Class_A']
    path_b = ['root', 'Bacteria', 'Phylum_B', 'Class_B']
    
    LCA is 'Bacteria' at depth 2.
    dist = (len(a) - lca_depth) + (len(b) - lca_depth)
    dist = (4 - 2) + (4 - 2) = 2 + 2 = 4
    """
    lca_depth = 0
    # Find the depth of the Lowest Common Ancestor (LCA)
    for i in range(min(len(path_a), len(path_b))):
        if path_a[i] == path_b[i]:
            lca_depth = i + 1
        else:
            break # Stop as soon as the paths diverge
            
    depth_a = len(path_a)
    depth_b = len(path_b)
    
    # Distance is the sum of steps from each node up to the LCA
    distance = (depth_a - lca_depth) + (depth_b - lca_depth)
    return distance

def get_dist(taxonomy_data):
    """
    Main function to load, process, and save the distance matrix.
    """
        
    all_targets_list = []
    target_to_path_map = OrderedDict() # Use OrderedDict to maintain insertion order
    
    print("Parsing taxonomic tree and mapping targets to paths...")
    # Start the recursive parse. We use 'root' as the starting path.
    # The JSON is already inside the "root" key, so we pass that node.
    parse_taxonomy_tree(taxonomy_data['root'], ['root'], all_targets_list, target_to_path_map)
    
    num_targets = len(all_targets_list)
    if num_targets == 0:
        print("Error: No targets found in the 'target_members' keys. Exiting.")
        return
        
    print(f"Found {num_targets} unique targets.")
    print("Building all-vs-all distance matrix...")
    
    # Initialize an empty matrix with zeros
    distance_matrix = np.zeros((num_targets, num_targets), dtype=int)
    
    # Pre-fetch all paths for faster lookup
    paths = [target_to_path_map[target] for target in all_targets_list]
    
    # Iterate through the upper triangle of the matrix
    for i in range(num_targets):
        
        path_a = paths[i]
        for j in range(i + 1, num_targets):
            path_b = paths[j]
            
            # Calculate and store the symmetric distance
            dist = get_taxonomic_distance(path_a, path_b)
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist
            
    print("Matrix calculation complete.")
    
    # Save the matrix to a user-friendly CSV file
    
    df = pd.DataFrame(distance_matrix, index=all_targets_list, columns=all_targets_list)
    df.index=df.index.rename('targets')

    return df