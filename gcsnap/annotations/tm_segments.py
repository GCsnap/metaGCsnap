import os
import subprocess
import shutil

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole
from gcsnap.genomic_context import GenomicContext
from gcsnap.sequence_mapping import SequenceMapping
from gcsnap.apis import EbiAPI

from gcsnap.utils import processpool_wrapper
from gcsnap.utils import split_list_chunks

import logging
logger = logging.getLogger('iteration')

class TMsegments:
    """ 
    Methods and attributes to annotate transmembrane segments and signal peptides of flanking genes.

    Attributes:
        config (Configuration): The Configuration object containing the arguments.
        cores (int): The number of CPU cores to use.
        annotate_TM (bool): The boolean to decide whether to annotate TM segments.
        annotate_mode (str): The mode to annotate TM segments.
        annotate_file (str): The path to the annotation file.
        gc (GenomicContext): The GenomicContext object containing all genomic context information.
        syntenies (dict): The dictionary with the syntenies of the target genes.
        out_label (str): The label of the output.
        out_dir (str): The path to store the output.
        fasta_file (str): The path to the fasta file.
        ncbi_code_order (list): The order of the ncbi codes in the fasta file.
        annotation_out_file (str): The path to the annotation output file.
        console (RichConsole): The RichConsole object to print messages.
        annotations (dict): The dictionary with the annotations of the flanking genes.
    """

    def __init__(self, config: Configuration, gc: GenomicContext, out_label: str):
        """
        Initialize the TMsegments object.

        Args:
            config (Configuration): The Configuration object containing the arguments.
            gc (GenomicContext): The GenomicContext object containing all genomic context information.
            out_label (str): The label of the output.
        """        
        self.config = config
        self.cores = config.arguments['n_cpu']['value']
        self.annotate_TM = config.arguments['annotate_TM']['value']
        self.annotate_mode = config.arguments['annotation_TM_mode']['value']
        self.annotate_file = config.arguments['annotation_TM_file']['value']

        # set parameters
        self.gc = gc
        self.syntenies = gc.get_syntenies()
        self.out_label = out_label  
        self.out_dir = os.path.join(os.getcwd(), f'{out_label}_TM_annotations')   
        if not os.path.isdir(self.out_dir) and self.annotate_TM:
            os.mkdir(self.out_dir)    
        self.annotations = {}                          

        self.console = RichConsole()

    def get_annotations(self) -> dict:
        """
        Getter for the annotations attribute.

        Returns:
            dict: The dictionary with the annotations of the flanking genes.
        """        
        return self.annotations

    def run(self) -> None:
        """
        Run the annotation of transmembrane segments and signal peptides of flanking genes:
            - Prepare data for annotation
            - Annotate TM segments
            - Parse annotation file
            - Update flanking genes
        """        
        # abort if no annotation requested
        if not self.annotate_TM:
            msg = 'TM annotation set to False. Transmembrane segments and signal peptides will not be searched'
            self.console.print_skipped_step(msg)
            return
        
        # create fasta file and get order if sequences in fasta file
        self.fasta_file = self.gc.write_to_fasta('{}_flanking_sequences.fasta'.format(self.out_label), 
                                                self.out_dir, exclude_pseudogenes = True)  
        self.ncbi_code_order = self.gc.get_fasta_order(False)      

        # outfile if not specified in config.yaml
        if self.annotate_file is None:
            # set as the name of fasta file without extension            
            self.annotation_out_file = os.path.join(self.out_dir, '{}_{}.out'.format(
                    os.path.basename(self.fasta_file)[:-6], self.annotate_mode))             
            self.signal_annotation()
        else:
            # set as the name of the input file            
            self.annotation_out_file = os.path.join(self.out_dir, '{}'.format(
                    os.path.basename(self.annotate_file)))              
            # copy the specified file to the output directory
            try:
                shutil.copy(self.annotate_file, self.annotation_out_file)
            # catch when it is not there
            except FileNotFoundError:
                self.console.print_error('Your specifed file {} could not be found.'. \
                                         format(self.annotate_file)) 
                msg =  'Check file is at speciefed path {} ans run GCsnap by specifying --annotation-TM-file.'. \
                        format(self.annotate_file)
                self.console.print_hint(msg)
                self.console.stop_execution()             
        
        # at this point, the file might not exist when no phobius or tmhmm installation was found
        if os.path.isfile(self.annotation_out_file):
            with self.console.status('Parsing annotation file and updating flanking genes'):
                self.parse_annotation_file()
                self.update_flanking_genes()
                # save in annotations
                self.annotations = self.syntenies

        # The code finishes upon errors in phobius or tmhmm
        # self.annotations remains empty and hence mergable to syntenies

    def signal_annotation(self) -> None:
        """
        Annotate signal peptides of flanking genes using either Phobius, TMHMM or UniProt.
        """        
        # do single tm annotation depending on mode
        if self.annotate_mode == 'phobius' or self.annotate_mode == 'tmhmm':
            with self.console.status('Annotating TM segments with {}'.format(self.annotate_mode)):                
                self.tool_annotation()
                # if phobius or tmhmm could not be executed, e warning instead of error
                # but no file is written.
        else:
            self.console.print_step('TM annotation with "uniprot" mode needs mapping first')
            # map all members to UniProtKB-AC
            mapping = SequenceMapping(self.config, self.ncbi_code_order, 
                                        to_type = 'UniProtKB-AC', msg = 'TM segments')
            mapping.run()
            mapping_dict = mapping.get_target_to_result_dict()
            mapping.log_failed()

            # request all annotation information from EbiAPI for all ncbi_code
            uniprot_codes = [mapping_dict.get(ncbi_code,'nan') for ncbi_code in self.ncbi_code_order]
            
            # split them into chunks of at most 90 (limit is at most 100)
            n_chunks = len(uniprot_codes) // 90 + 1
            parallel_args = split_list_chunks(uniprot_codes, n_chunks)
            with self.console.status('Get functional annotation from EBI'):
                result_list = processpool_wrapper(self.cores, parallel_args, self.run_get_functional_annotations)
                # combine results
                all_uniprot_data = {k: v for dict_ in result_list for k, v in dict_.items()}        

            with self.console.status('Annotating TM segments with {}'.format(self.annotate_mode)):
                # create parallel arguments
                parallel_args = [(sub_list, mapping_dict, all_uniprot_data) 
                                for sub_list in split_list_chunks(self.ncbi_code_order, self.cores)]

                result_lists = processpool_wrapper(self.cores, parallel_args, self.uniprot_annotation)
                # combine results
                result_list = [result for sub_list in result_lists for result in sub_list]   

                # create writing list
                write_list = ['{}\t{}'.format(ncbi_code, tm_annotation) 
                                for ncbi_code, tm_annotation in result_list]   
                self.write_to_file('\n'.join(write_list)) 

    def tool_annotation(self) -> None:
        """
        Annotate transmembrane segments and signal peptides of flanking genes using Phobius or TMHMM.

        Raises:
            FileNotFoundError: If Phobius or TMHMM is not installed.
        """        
        if self.annotate_mode == 'phobius':
            tool = 'phobius.pl'
        elif self.annotate_mode == 'tmhmm':
            tool = 'tmhmm'

        if not os.path.isfile(self.annotation_out_file):
            try:
                stdout, stderr = self.run_command(tool)
                if len(stderr) > 0 or len(stdout) == 0:
                    raise FileNotFoundError
                # if no error raised
                self.write_to_file(stdout)
            except FileNotFoundError:
                self.console.print_warning('No {} installation was found.'.format(self.annotate_mode)) 
                msg =  'Run {} online and run GCsnap by specifying --annotation-TM-file.'. \
                        format(self.annotate_mode)
                self.console.print_hint(msg)
                self.extended_instructions()      

    def run_command(self, tool: str) -> tuple:
        """
        Run Phobius or TMHMM command to execute.

        Args:
            tool (str): Either 'phobius.pl' or 'tmhmm'.

        Returns:
            tuple: The stdout and stderr of the Phobius or TMHMM command.
        """        
        # returns stdout,stderr
        command = [tool, 
                    self.fasta_file,
                    '-short']
        
        result = subprocess.run(command, capture_output=True, text=True, shell=True)        
        return result.stdout, result.stderr   
    
    def uniprot_annotation(self, args: tuple[list,dict,dict]) -> list[tuple]:
        """
        Annotate transmembrane segments and signal peptides of flanking genes using UniProt
        used in parallel processing.

        Args:
            args (tuple[list,dict,dict]): The arguments to annotate TM segments.
                First element is a list of ncbi codes.
                Second element is a dictionary with the mapping of ncbi codes to uniprot codes.
                Third element is a dictionary with the retrieved uniprot data.

        Returns:
            list[tuple]: The list of tuples with the ncbi code and the annotation.
        """         
        ncbi_codes, mapping_dict, all_uniprot_data = args

        result_list = []

        for ncbi_code in ncbi_codes:
            uniprot_code = mapping_dict.get(ncbi_code, 'nan')
            # if not found, it will be nan
            if uniprot_code == 'nan':
                result_list.append((ncbi_code, 'nan'))
                continue

            # get uniprot annotations, returns a dictionary
            curr_uniprot_annotations = all_uniprot_data.get(uniprot_code)

            if curr_uniprot_annotations != 'nan':
                tm_annotation = curr_uniprot_annotations['TM_topology']
                if tm_annotation == '':
                    tm_annotation = 'nan'
            else:
                tm_annotation = 'nan'

            result_list.append((ncbi_code, tm_annotation))
        return result_list
    
    def parse_annotation_file(self) -> None:
        """
        Parse the annotation file to extract the transmembrane segments and signal peptides of flanking
        genes and store them in a dictionary.
        """        
        # read full file to list
        with open(self.annotation_out_file, 'r') as file:
            lines = file.readlines()

        # there would be no need to check for any other mode, as the parser
        # with choices took care of that
        if self.annotate_mode == 'phobius':
            self.protein_annotations = self.parse_phobius_output(lines)
        elif self.annotate_mode == 'tmhmm':
            self.protein_annotations = self.parse_tmhmm_output(lines)
        elif self.annotate_mode == 'uniprot':
            self.protein_annotations = self.parse_tm_uniprot_output(lines)
        else:
            self.console.print_warning('Annotation mode not recognized.')
            self.protein_annotations = {}

    def parse_phobius_output(self, lines: list) -> dict:
        """
        Parse the Phobius output file.

        Args:
            lines (list): The Phobius output file lines.

        Returns:
            dict: The dictionary with the annotations of the flanking genes.
        """        
        annotations = {}
        for line in lines:
            if 'PREDICTION' not in line:
                ncbi_code = line.split('|')[0]
                tm_pred = int(line.split()[1])
                sp_pred = line.split()[2]

                if tm_pred > 0:
                    annotations[ncbi_code] = 'TM' # it has a transmembrane helix -> it is transmembrane
                elif sp_pred == 'Y':
                    annotations[ncbi_code] = 'SP' # it has a signal peptide
                else:
                    annotations[ncbi_code] = '' # does not have any special feature
        return annotations

    def parse_tmhmm_output(self, lines: list) -> dict:
        """
        Parse the TMHMM output file.

        Args:
            lines (list): The TMHMM output file lines.

        Returns:
            dict: The dictionary with the annotations of the flanking genes.
        """        
        annotations = {}
        for line in lines:
            if 'PREDICTION' not in line:
                ncbi_code = line.split('|')[0]
                tm_pred = line.split('Topology=')[-1].strip()

                if tm_pred != 'o' and tm_pred != 'i':
                    annotations[ncbi_code] = 'TM' # it has a transmembrane helix -> it is transmembrane
                else:
                    annotations[ncbi_code] = '' # does not have any special feature (signal peptides are not identified)
        return annotations     

    def parse_tm_uniprot_output(self, lines: list) -> dict:
        """
        Parse the UniProt output file.

        Args:
            lines (list): The UniProt output file lines.

        Returns:
            dict: The dictionary with the annotations of the flanking genes.
        """        
        annotations = {}
        for line in lines:
            ncbi_code = line.split()[0]
            tm_pred = line.split()[-1].strip()

            if tm_pred != 'nan':
                annotations[ncbi_code] = 'TM' 
            else:
                annotations[ncbi_code] = '' 
        return annotations       
    
    def update_flanking_genes(self) -> None:
        """
        Update the flanking genes with the transmembrane segments and signal peptides.
        """        
        for curr_target in self.syntenies:
            flanking_genes = self.syntenies[curr_target]['flanking_genes']
            flanking_genes['TM_annotations'] = []

            for ncbi_code in flanking_genes['ncbi_codes']:
                if ncbi_code in self.protein_annotations:
                    flanking_genes['TM_annotations'].append(self.protein_annotations[ncbi_code])
                else:
                    flanking_genes['TM_annotations'].append('')

            self.syntenies[curr_target]['flanking_genes'] = flanking_genes        

    def write_to_file(self, data: str) -> None:
        """
        Write the annotation data to a file.

        Args:
            data (str): The annotation data to write to the file.
        """        
        with open(self.annotation_out_file, 'w') as file:
            file.write(data)

    def extended_instructions(self) -> None:
        """
        Print extended instructions for the user.
        """        
        msg_output_file = 'Output file: copy-paste output text below the line to a .txt file'
        msg_save_folder = 'Provide path to file via --annotation-TM-file either in CLI or config.yaml'         
        if self.annotate_mode == 'phobius':
            msg_server = 'Phobius webserver: http://phobius.sbc.su.se/'
            msg_run_mode = 'Run mode: "short"'
        elif self.annotate_mode == 'tmhmm':
            msg_server = 'TMHMM webserver: https://services.healthtech.dtu.dk/service.php?TMHMM-2.0'
            msg_run_mode = 'Run mode: "One line per protein"'

        for msg in [msg_server, msg_run_mode, msg_output_file, msg_save_folder]:
            self.console.print_hint(msg)

    def run_get_functional_annotations(self, uniprot_codes: list) -> dict:
        """
        Get functional annotations from EBI for a list of UniProt codes using
        the EbiAPI.

        Args:
            uniprot_codes (list): The list of UniProt codes.

        Returns:
            dict: The dictionary with the functional annotations of the flanking genes.
        """        
        return EbiAPI.get_uniprot_annotations_batch(uniprot_codes, with_parsing = True)
        

