import os
from typing import Union

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole 

class Target():
    """ 
    Methods and attributes to parse the targets from the input files.

    Attributes:
        arguments (dict): The dictionary with the arguments.
        targets (Union[str, list]): The targets to parse.
        console (RichConsole): The RichConsole object to print messages.
        targets_dict (dict): The dictionary with the parsed targets.
    """

    def __init__(self, configuration: Configuration):
        """
        Initialize the Target object.

        Args:
            configuration (Configuration): The Configuration object containing the arguments.
        """
        self.arguments = configuration.arguments
        self.targets = configuration.targets
        self.console = RichConsole('base')

        # empty dictionary to store the targets
        # shape: {out_label: {provider: [ids]}}
        self.targets_dict = {}
        self.all_ids = []

    # ── per-provider API ──────────────────────────────────────────────────────

    def add_provider_targets(self, out_label: str, provider: str, ids: list) -> None:
        """
        Register a list of target IDs for a specific provider.

        Args:
            out_label (str): The run label (top-level key).
            provider (str): Provider name, e.g. ``'ncbi'``, ``'mgnify'``, ``'local'``.
            ids (list): List of target ID strings.
        """
        if out_label not in self.targets_dict:
            self.targets_dict[out_label] = {}
        self.targets_dict[out_label][provider] = ids
        self.all_ids.extend(ids)

    def get_provider_targets(self, provider: str) -> list:
        """
        Return the target IDs registered for *provider* under the current
        ``out_label``.

        Args:
            provider (str): Provider name, e.g. ``'ncbi'``.

        Returns:
            list: Target IDs, or an empty list when the provider was not
                  registered.
        """
        out_label = self.arguments['out_label']['value']
        return self.targets_dict.get(out_label, {}).get(provider, [])

    def run(self) -> None:
        """
        Run the parsing of the targets.
        """        
        # ignore --out-label if multiple files are given
        if len(self.targets) > 1 and os.path.isfile(self.targets[1]):
            self.arguments['out_label']['value'] = 'default'
            
        with self.console.status('Parsing targets'):
            self.parse_targets(self.targets)

    def get_targets_dict(self) -> dict:
        """
        Getter for the targets_dict attribute.

        Returns:
            dict: The dictionary with the parsed targets.
        """        
        return self.targets_dict

    def parse_targets(self, targets: Union[str, list]) -> dict:
        """
        Parse the targets from the input files or the list of ids.

        Args:
            targets (Union[str, list]): The targets to parse. Either a string or a list of strings.

        Returns:
            dict: The dictionary with the parsed targets.
        """      
        for target in targets:
            # target is a file
            if os.path.isfile(target):
                if target.endswith('.clans'):
                    self.get_clusters_from_clans(target)
                else:
                    self.get_targets_from_file(target)
                
                # restet outlabel, as it is only used for the first file

            elif target.endswith('.fasta') or target.endswith('.txt') or target.endswith('.clans'):
                # the error case when file path was missspelled, hence file not found
                self.console.print_error('Target file {} does not exist, maybe the path is misspelled?'.format(target))
                self.console.stop_execution()
            else:
                # target is a list of ids
                # the name of this list (the label) by default is 'default' 
                label = self.arguments['out_label']['value']
                if label not in self.targets_dict:
                    self.targets_dict[label] = []
                
                self.targets_dict[label].append(target)

    def get_targets_from_file(self, target_file: str) -> None:
        """
        Get the targets from a file.

        Args:
            target_file (str): The path to the file with the targets.
        """        
        # set label
        if self.arguments['out_label']['value'] == 'default':
            curr_label = os.path.basename(target_file).split('.')[0]
        else:
            curr_label = self.arguments['out_label']['value']

        if curr_label not in self.targets_dict:
            self.targets_dict[curr_label] = []

        is_fasta = False
        with open(target_file, 'r') as f:
            content = f.readlines()

        for line in content:
            if line.startswith('>'):
                is_fasta = True
                curr_target = line.split(' ')[0].split(':')[0].split('|')[0].split('_#')[0].replace('>','').strip()
            elif not is_fasta:
                curr_target = line.strip().split()[0]

            if curr_target not in self.targets_dict[curr_label]:
                self.targets_dict[curr_label].append(curr_target)        

    def get_clusters_from_clans(self, clans_file: str) -> dict:
        """
        Get the clusters from a clans file.

        Args:
            clans_file (str): The path to the clans file.

        Returns:
            dict: The dictionary with the clusters.
        """        
        ncbis_ordered = self.get_ncbicodes_order_in_clans(clans_file)
        cluster_codes = self.arguments['clans_patterns']['value']

        # split cluster codes into list

        
        if cluster_codes != None:
            found_seqgroup = False
            with open(clans_file, 'r') as f:
                content = f.readlines()

            for line in content:
                for cluster_code in cluster_codes:
                    if '<seqgroup' in line:
                        found_seqgroup = True
                        found_allowed_cluster = False

                    elif found_seqgroup:
                        if 'name=' in line and cluster_code in line:

                            current_cluster = line.split('=')[-1].strip()
                            found_allowed_cluster = True
                        
                        elif 'numbers=' in line and found_allowed_cluster:
                            numbers = line.split('=')[-1].split(';')[:-1]
                            numbers = [int(i) for i in numbers]
                            numbers = [ncbis_ordered[i] for i in numbers]
                            self.targets_dict[current_cluster] = numbers

                            found_allowed_cluster = False

        else:
            label = os.path.basename(clans_file).split('.')[0]
            self.targets_dict[label] = ncbis_ordered
    
    def get_ncbicodes_order_in_clans(self, clans_file: str) -> list:   
        """
        Get the NCBI codes in the order they appear in the clans file.

        Args:
            clans_file (str): The path to the clans file.

        Returns:
            list: The list with the NCBI codes in the order they appear in the clans file.
        """        
        ncbids_ordered = []

        found_seq = False
        with open(clans_file, 'r') as clans:
            for line in clans:
                if '<seq>' in line:
                    found_seq = True
                elif '</seq>' in line:
                    found_seq = False
                elif found_seq and line.startswith('>'):
                    line = line[1:]
                    ncbi_code = line.split(' ')[0].split(':')[0].split('|')[0].split('_#')[0].replace('>','').strip()
                    ncbids_ordered.append(ncbi_code)
                    
        return ncbids_ordered
 