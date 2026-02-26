from gcsnap.sequence_mapping import SequenceMapping
from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole
from gcsnap.genomic_context import GenomicContext
from gcsnap.apis import SwissProtAPI, AlphaFoldAPI, EbiAPI

from gcsnap.utils import processpool_wrapper
from gcsnap.utils import split_dict_chunks, split_list_chunks

class FamiliesFunctionsStructures:
    """ 
    Methods and attributes to get the functional annotations and structures of the families.
    
    Attributes:
        config (Configuration): The Configuration object containing the arguments.
        cores (int): The number of CPU cores to use.
        get_pdb (bool): The flag to get the PDB structures.
        get_annotations (bool): The flag to get the functional annotations.
        families (dict): The dictionary with the families and their members.
        console (RichConsole): The RichConsole object to print messages.
        annotations_and_structures (dict): The dictionary with the functional annotations and structures
    """

    def __init__(self, config: Configuration, gc: GenomicContext):
        """
        Initialize the FamiliesFunctionsStructures object.

        Args:
            config (Configuration): The Configuration object containing the arguments.
            gc (GenomicContext): The GenomicContext object containing all genomic context information.
        """        
        self.config = config
        self.cores = config.arguments['n_cpu']['value']
        self.get_pdb = config.arguments['get_pdb']['value']
        self.get_annotations = config.arguments['get_functional_annotations']['value']
        self.annotations_and_structures = {}

        # set parameters
        self.families = gc.get_families()

        self.console = RichConsole()

    def get_annotations_and_structures(self) -> dict:
        """
        Getter for the annotations_and_structures attribute.

        Returns:
            dict: The dictionary with the functional annotations and structures.
        """        
        return self.annotations_and_structures
    
    def run(self) -> None:
        """
        Run the retrieval of the functional annotations and structures:
            - Map every family member to UniProtKB-AC.
            - Request all annotation information from EbiAPI in parallel.
            - Extract the functional annotations and structures for each family in parallel.
        Usese parallel processing with the processpool_wrapper from utils.py.
        Only done when conserved families were found.
        """        
        if (list(self.families.keys()) == [0] or list(self.families.keys()) == [1] or
                    list(self.families.keys()) == [-1, 0]):
            self.console.print_warning('No conserved family found and no functional annotation possible')
        elif self.get_pdb or self.get_annotations:
            self.console.print_step('Functional annotation needs mapping first')
            # find all family members
            # always exclude non-conserved families and pseudogenes as uniprot mapping
            all_members = [member_code for family in self.families.keys() 
                           for member_code in self.families[family]['members']
                           if 0 < family]
            
            # map all members to UniProtKB-AC
            mapping = SequenceMapping(self.config, all_members, 
                                            'UniProtKB-AC', 'family members')
            mapping.run()
            mapping_dict = mapping.get_target_to_result_dict()
            mapping.log_failed()

            # request all annotation information from EbiAPI for all members with uniprot code
            uniprot_codes = [mapping_dict.get(member,'nan') for member in all_members]
            
            # split them into chunks of at most 90 (limit is at most 100)
            n_chunks = len(uniprot_codes) // 90 + 1
            parallel_args = split_list_chunks(uniprot_codes, n_chunks)
            with self.console.status('Get functional annotation from EBI'):
                result_list = processpool_wrapper(self.cores, parallel_args, self.run_get_functional_annotations)
                # combine results
                all_uniprot_data = {k: v for dict_ in result_list for k, v in dict_.items()}

            # # Version 1: Parallel processing with chunks
            # # split the dictionary into chunks
            # dict_list = split_dict_chunks(self.families, self.cores)
            # # create parallel arguments
            # parallel_args = [(dict_, mapping_dict, all_uniprot_data, self.get_pdb, self.get_annotations) 
            #                  for dict_ in dict_list]
            
            # Version 2: Parallel processing each family
            # As this relieas on APIs, we try to minimies load imbalance not
            # know apriori how long each family will take (as depends on the number of members)
            parallel_args = [({k: v}, mapping_dict, all_uniprot_data, self.get_pdb, self.get_annotations)
                             for k,v in self.families.items()]

            with self.console.status('Get functional annotations and structures'):
                result_list = processpool_wrapper(self.cores, parallel_args, self.run_each)
                # combine results
                self.annotations_and_structures = {k: v for dict_ in result_list for k, v in dict_.items()}
        else:
            self.console.print_skipped_step('No annotations or structures retrieval requested')

    def run_each(self, args: tuple[dict,dict,dict,bool,bool]) -> dict:
        """
        Run the retrieval of the functional annotations and structures for each family
        used in parallel processing.

        Args:
            args (tuple[dict,dict,dict,bool,bool]): The arguments for the retrieval.
                First element is a dictionary with the family and its members.
                Second element is a dictionary with the mapping of the members to UniProtKB-AC.
                Third element is a dictionary with the UniProtKB-AC codes and their annotations.
                Fourth element is a flag to get the PDB structures.
                Fifth element is a flag to get the functional annotations.

        Returns:
            dict: The dictionary with the family and its members with the functional annotations and structures.
        """        
        families, mapping_dict, all_uniprot_data, get_pdb, get_annotations = args

        # all families that are not pseudogenes
        # a list of tuples (family, member)
        family_keys = [(key, sorted(families[key]['members'])) 
                        for key in sorted(list(families.keys()))
                        if 'function' not in families[key]]
        
        family_function = {}
        family_structure = ''
        family_uniprot_code = ''
        family_model_state = 'Model does not exist'

        for family, members in family_keys:

            # pseudo genes and non-conserved families
            if family == 0 or family == -1:
                family_function = {}
                family_structure = 'n.a.'
                family_uniprot_code = 'n.a.'
                family_model_state = 'n.a.'

            # conserved families
            else:
                family_function = {}
                family_structure = ''
                family_uniprot_code = ''
                family_model_state = 'Model does not exist'

                # it is done until a result is found
                for member in members:
                    # if not found, it will be nan
                    uniprot_code = mapping_dict.get(member,'nan')

                    if (family_uniprot_code == '' or family_structure == '' or 
                        (family_function == '' or family_function['Function_description'] == '')):
                    
                        if get_pdb and family_structure == '':
                            # check either in SwissModel or AlphaFold
                            curr_pdb = SwissProtAPI.find_uniprot_in_swissmodel(uniprot_code)
                            if 'nan' in curr_pdb:
                                curr_pdb = AlphaFoldAPI.find_uniprot_in_alphafold_database(uniprot_code)

                            if 'nan' not in curr_pdb:
                                family_structure = curr_pdb
                                family_uniprot_code = uniprot_code
                                family_model_state = 'Model exists'
                            elif uniprot_code != 'nan' and curr_pdb == 'nan*':
                                family_uniprot_code = uniprot_code

                        if (get_annotations and (family_function == {} 
                                                 or family_function['Function_description'] == '')):        
                            # get functional annotations
                            uniprot_accession = uniprot_code.split('_')[0]
                            uniprot_data = all_uniprot_data.get(uniprot_accession, {})
                            curr_uniprot_annotations = EbiAPI.parse_uniprot_data(uniprot_data, 
                                                                    previous_annotations = family_function)
                            if curr_uniprot_annotations != 'nan':
                                family_function = curr_uniprot_annotations

                            if uniprot_code != 'nan' and family_structure == '' and family_uniprot_code == '':
                                family_uniprot_code = uniprot_code

                if family_uniprot_code == '':
                    family_model_state = 'Not possible to map'

            families[family]['function'] = 'n.a.' if family_function == {} else family_function
            families[family]['structure'] = family_structure
            families[family]['uniprot_code'] = family_uniprot_code
            families[family]['model_state'] = family_model_state

        return families
    
    def run_get_functional_annotations(self, uniprot_codes: list) -> dict:
        """
        Run the retrieval of the functional annotations for a list of UniProtKB-AC codes
        used in parallel processing.

        Args:
            uniprot_codes (list): The list of UniProtKB-AC codes.

        Returns:
            dict: The dictionary with the UniProtKB-AC codes and their annotations.
        """        
        return EbiAPI.get_uniprot_annotations_batch(uniprot_codes)