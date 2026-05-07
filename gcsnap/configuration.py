import os
import importlib.util
import argparse
import psutil

# pip install rich
from rich.console import Console
# pip install pyyaml
import yaml
# pip install ruamel.yaml
from ruamel.yaml import YAML
from ruamel.yaml.constructor import ConstructorError

from gcsnap.rich_console import RichConsole 

class CustomArgumentParser(argparse.ArgumentParser):
    """
    Methods to handle the argument parser and the arguments from the CLI.

    Attributes:
        console (type): 

    Inheritance:
        argparse.ArgumentParser: The ArgumentParser class from the argparse module.
    """
    
    def __init__(self, usage: str = None, epilog: str = None):
        """
        Initialize the CustomArgumentParser object.

        Args:
            usage (str, optional): The usage message discribing the program's usage. Defaults to None.
            epilog (str, optional): The epilog message discribing the program's epilog. Defaults to None.
        """   
        super().__init__(usage=usage, epilog=epilog, add_help=False)
        self.console = RichConsole('base')

    def add_argument_from_config(self, config: dict) -> None:
        """
        Add arguments to the parser from a configuration dictionary.

        Args:
            config (dict): Dictionary with arguments and their details.
        """        
        for key, details in config.items():
            value_type = eval(details['type'])
            help_text = details['help']
            help_choices = details.get('choices')            
            default_value = self.str_to_none(details['value'])
            
            if value_type == bool:
                # argparse does not support bool type in config, so we need to convert it
                self.add_argument(f'--{key}', type=self.str_to_bool, default=default_value, 
                                  help=help_text)
            else:
                self.add_argument(f'--{key}', type=value_type, default=default_value, 
                                  help=help_text, choices=help_choices)

    def error(self, message: str) -> None:
        """
        Overwrite error method from argparse.ArgumentParser to print error message and hint.

        Args:
            message (str): Error message from argparse.
        """        
        self.console.print_error(f'{message}')
        self.console.print_hint('Use --help to see all supported arguments.')
        self.console.stop_execution()

    
    # @staticmethod: Used when you need a method inside a class but don't need to access or modify the instance (self) or class (cls) state.
    # Avoids unnecessary self or cls (class) parameters.
    # Ideal for utility functions related to the class's purpose but not dependent on instance or class data.    
    @staticmethod
    def str_to_bool(value: str) -> bool:
        """
        Convert a string to a boolean.

        Args:
            value (str): The string to convert to a boolean.

        Raises:
            argparse.ArgumentTypeError: Error message if the argument value is not supported.

        Returns:
            bool: The boolean value corresponding to the string.
        """        
        if isinstance(value, bool):
            return value
        if value.lower() in ('yes', 'true', 'True' , 't', 'y', '1'):
            return True
        elif value.lower() in ('no', 'false', 'False', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    @staticmethod
    def str_to_type(type_str: str) -> type:
        """
        Convert a type-string to its type.

        Args:
            type_str (str): The string representing the type.

        Raises:
            ValueError: Error message if the type is not supported.

        Returns:
            type: The type corresponding to the string.
        """        
        if type_str == 'bool':
            return bool
        elif type_str == 'int':
            return int
        elif type_str == 'float':
            return float
        elif type_str == 'str':
            return str
        else:
            raise ValueError(f'Unsupported type: {type_str}')   

    @staticmethod
    def str_to_none(value: str) -> any:
        """
        Convert a string to None if the string is 'None'.

        Args:
            value (str): The string to convert to None.

        Returns:
            any: The value or None.
        """        
        if value == 'None':
            return None
        return value                     


class Configuration:    
    """
    Methods and attributes to parse the arguments from the CLI and the configuration file.
    While the arguments from the CLI and the config.yaml use hyphens, 
    the internal dictionary uses underscores.

    Attributes:
        path (str): The path to the configuration file.
        console (RichConsole): The RichConsole object to print messages.
        default_write (bool): Boolean statement to write the default configuration to file.
        arguments_hyphen (dict): Dictionary with arguments from the configuration file.
        arguments (dict): Dictionary with arguments from the configuration file with underscores.
    """

    def __init__(self, config_path: str = None):       
        """
        Initialize the Configuration object.
        """                      
        self.console = RichConsole('base')

        # path to the configuration file
        self.set_configuration_path(config_path)

        # in case config.yaml is not found, use default configuration
        # and by default write to file after all parsing
        self.default_write = False

        self.read_configuration_yaml()
        with self.console.status('Parsing CLI arguments and config.yaml'):
            self.create_argument_parser()

    def __getstate__(self) -> dict:
        """
        Used by pickle to convert the object into a dictionary for serialization.
        Parser object is not picklable, so it's removed from the state.

        Returns:
            dict: The state of the object without the console.
        """         
        state = self.__dict__.copy()
        del state['parser']
        return state

    def __setstate__(self, state: dict) -> None:
        """
        Used by pickle to restore the object from a dictionary after deserialization.
        Parser object is recreated.

        Args:
            state (dict): The state of the object.
        """              
        self.__dict__.update(state)
        self.create_argument_parser()

    def hyphen_to_underscore(self, argument: str) -> str:
        """
        Convert a string with hyphens to a string with underscores.

        Args:
            argument (str): String with hyphens.

        Returns:
            str: String with underscores.
        """        
        return argument.replace('-', '_')

    def underscore_to_hyphen(self, argument: str) -> str:
        """
        Convert a string with underscores to a string with hyphens.

        Args:
            argument (str): String with underscores.

        Returns:
            str: String with hyphens.
        """        
        return argument.replace('_', '-') 

    def set_configuration_path(self, config_path: str = None) -> None:
        """
        Set the path to the configuration file.
        If GCsnap is executed from a different directory, the path is updated.
        In that case, the configuration file is not likely to be not found
        and a default configuration is created.
        """        

        if config_path is not None:
            self.path = config_path
            return

        # Locate the GCsnap package directory
        spec = importlib.util.find_spec('gcsnap')
        if spec is None or spec.origin is None:
            self.console.print_error('GCsnap package seems not to be installed')
            self.console.print_hint('Please install GCsnap package and try again')
        package_dir = os.path.dirname(spec.origin)   

        # Check if GCsnap is executed from a different directory
        if os.path.samefile(package_dir, os.getcwd()):
            self.path = os.path.join(package_dir, 'config.yaml')
        else:
            self.path = os.path.join(os.getcwd(), 'config.yaml')
        
    def read_configuration_yaml(self) -> None:  
        """
        Read the configuration file and load the arguments from it.
        If the file is not found, create a default configuration and set
        boolean statement to write it to file at the end.
        Additional check for an empty configuration file.
        """         
        self.is_yaml_valid()

        if not os.path.isfile(self.path):
            self.console.print_warning('Configuration file config.yaml not found')
            self.arguments_hyphen = self.get_default_configuration()
            # set default write to True to write the default configuration to file
            self.default_write = True
            self.console.print_done('Default config.yaml created')
        else:
            with open(self.path, 'r') as file:
                self.arguments_hyphen = yaml.load(file, Loader=yaml.FullLoader)
            self.console.print_done('Configuration file config.yaml loaded')

    def is_yaml_valid(self) -> None:
        """
        Check if a YAML file is valid and not empty. If not, delete the broken file
        and return False.
        """
        yaml = YAML()
        try:
            with open(self.path, 'r') as file:
                data = yaml.load(file)
                if data is None:
                    raise Exception
        except FileNotFoundError:
            pass                
        except (ConstructorError, Exception):
            os.remove(self.path)
            
    def write_configuration_yaml(self) -> None:
        """
        Write the configuration dictionary to the configuration file.
        """        
        out = {self.underscore_to_hyphen(key): value for key, value in self.arguments.items()}
        header_comment = (
            'Configuration file to set arguments for GCsnap.\n'
            'To change argument, change: value: entry.\n'
            'E.g. value: 1 to value: 2\n'
            '---------------------------------------\n')
        
        with open(self.path, 'w') as file:
            file.write('# ' + header_comment.replace('\n', '\n# ') + '\n')
            # 4 is the indentation space
            yaml = YAML()
            yaml.indent(mapping=4, sequence=4, offset=2)  # Set indent and block sequence indent
            yaml.dump(out, file)

    def write_configuration_yaml_log(self, file_name: str, file_path: str = None) -> None:
        """
        Write the configuration dictionary to a log file.

        Args:
            file_name (str): The name of the log file.
            file_path (str, optional): The path to the log file. Defaults to None.
        """        
        if file_path is None:
            file_path = os.getcwd()
        lines_to_write = ['{}:\t{}\n'.format(self.underscore_to_hyphen(key), value['value']) 
                          for key, value in self.arguments.items()]

        with open(os.path.join(file_path, file_name), 'w') as file:
            file.writelines(lines_to_write)
            
    def parse_arguments(self) -> None:
        """
        Parse and handle the arguments from the CLI and update the configuration dictionary.
        If --help is detected, call print method from console attribute and exit.
        """        
        args = self.parser.parse_args()

        # show help message (--help) and exit
        if hasattr(args, 'help') and args.help:
            self.console.print_help(self.parser)
            exit(0)

        # Update self.arguments dictionary with parsed values
        for arg in vars(args):
            # in arguments, we use underscores, in parser hyphens
            config_key = self.hyphen_to_underscore(arg)
            if config_key in self.arguments:
                self.arguments[config_key]['value'] = getattr(args, arg)
            elif arg == 'targets':
                self.targets = getattr(args, 'targets')
        
        # handle arguments that require special treatment
        #self.handle_special_arguments()

        # Write updated configuration to file
        if self.arguments['overwrite_config']['value'] or self.default_write:
            self.write_configuration_yaml()
        
    def create_argument_parser(self) -> None:
        """
        Create the argument parser with usage and epilog messages.
        """        
        usage = (
            'GCsnap --ncbi-targets <file> [--mgnify-targets <file>] '
            '[--local-targets <file>] [Optional arguments]'
        )
        epilog = (
            'Example: GCsnap --ncbi-targets ncbi_ids.txt --mgnify-targets mgnify_ids.txt\n'
            'At least one provider flag is required.\n'
            'Legacy positional targets (--targets) are still accepted for '
            'notebook / scripting use.'
        )

        self.parser = CustomArgumentParser(usage=usage, epilog=epilog)

        # help argument (replaces default --help with rich printing)
        self.parser.add_argument('-h', '--help', action='store_true',
                                 help='Show this help message and exit')

        # legacy targets argument – kept for notebook / Target() compatibility
        self.parser.add_argument('-t', '--targets', nargs='+', type=str,
                                 help='(Legacy) List of input targets for use with the '
                                      'Target class or direct scripting.')

        # Add arguments from YAML
        self.parser.add_argument_from_config(self.arguments_hyphen)

        # Convert hyphen arguments to underscore
        self.arguments = {self.hyphen_to_underscore(key): value for key, value
                          in self.arguments_hyphen.items()}

    def handle_special_arguments(self) -> None:
        """
        Handle special arguments that require additional processing.
        """               
        # clanse file
        targets = self.targets
        clans_file = self.arguments['clans_file']['value']
        if len(targets) == 1 and targets[0].endswith('.clans') and clans_file is None:
            clans_file = targets[0]
        if clans_file is not None:
            clans_file = os.path.abspath(clans_file)

        self.arguments['clans_file']['value'] = clans_file

    def get_default_configuration(self) -> dict:
        """
        Define the default configuration dictionary if no config.yaml is found.
        This can be used to create a new config.yaml file if it was messed up.
        Just delete the config.yaml file and run the GCsnap again.

        Returns:
            dict: The default configuration dictionary.
        """        
        return {
            "out-label": {
                "value": "default",
                "type": "str",
                "help": "Name of output directory. If default, name of the input file."
            },
            "ncbi-targets": {
                "value": None,
                "type": "str",
                "help": "Path to a plain-text file of NCBI / UniProt target IDs (one per line)."
            },
            "mgnify-targets": {
                "value": None,
                "type": "str",
                "help": "Path to a plain-text file of MGnify (MGYP…) target IDs (one per line)."
            },
            "local-targets": {
                "value": None,
                "type": "str",
                "help": "Path to a plain-text file of local-database target IDs (one per line)."
            },
            "gff-path": {
                "value": None,
                "type": "str",
                "help": "Path to the flat folder containing all .gff.gz annotation files (used at pipeline runtime)."
            },
            "db-path": {
                "value": None,
                "type": "str",
                "help": "Path to the GCsnap database folder containing assemblies.db, sequences.db and the summary files."
            },
            "tmp-mmseqs-folder": {
                "value": None,
                "type": "str",
                "help": "The temporary folder to store mmseqs files. May be changed so that intermediary mmseqs files are saved somewhere else then the automatic 'out-label' directory."
            },
            "mmseqs-executable-path": {
                "value": None,
                "type": "str",
                "help": "Path of MMseqs executable (i.e., mmseqs.bat) if not installed in Conda environment."
            },            
            "collect-only": {
                "value": False,
                "type": "bool",
                "help": "Boolean statement to make GCsnap collect genomic contexts only, without comparing them."
            },
            "n-flanking5": {
                "value": 4,
                "type": "int",
                "help": "Number of flanking sequences to take on 5' end."
            },
            "n-flanking3": {
                "value": 4,
                "type": "int",
                "help": "Number of flanking sequences to take on 3' end."
            },
            "exclude-partial": {
                "value": False,
                "type": "bool",
                "help": "Exclude partial operon/genomic_context blocks. If turned off, partial cases will still be ignored to get the most common genomic features."
            },
            "max-evalue": {
                "value": 0.001,
                "type": "float",
                "help": "Max e-value at which two sequences are considered to be homologous. Required to define protein families."
            },
            "default-base": {
                "value": 10,
                "type": "int",
                "help": "Artificial distance value for two sequences that do not match with an E-value better than --max-evalue."
            },
            "min-coverage": {
                "value": 0.7,
                "type": "float",
                "help": "Minimum coverage of target and subject a match needs to be so that two sequences are considered to be homologous. Required to define protein families."
            },
            "num-iterations": {
                "value": 1,
                "type": "int",
                "help": "Number of iterations for all-against-all searches. Required to define protein families."
            },
            "operon-cluster-advanced": {
                "value": False,
                "type": "bool",
                "help": "Boolean statement to use the operon clustering advanced mode using PacMAP."
            },
            "max-family-freq": {
                "value": 20,
                "type": "int",
                "help": "Maximum frequency of a family in the set of genomic contexts found to be considered for advanced operon clustering."
            },
            "min-family-freq": {
                "value": 2,
                "type": "int",
                "help": "Minimum frequency of a family in the set of genomic contexts found to be considered for advanced operon clustering."
            },
            "min-family-freq-accross-contexts": {
                "value": 30,
                "type": "int",
                "help": "Minimum frequency of a family in a conserved genomic context type to be considered as a member."
            },
            "n-max-operons": {
                "value": 30,
                "type": "int",
                "help": "Maximum number of top most populated operon/genomic_context block types."
            },
            "interactive": {
                "value": True,
                "type": "bool",
                "help": "Boolean statement to make the interactive html output."
            },
            "genomic-context-cmap": {
                "value": "Spectral",
                "type": "str",
                "help": "Color map (as of matplotlib) to assign colors to and plot the syntenic blocks."
            },
            "gc-legend-mode": {
                "value": "species",
                "type": "str",
                "help": "Mode of the genomic context legend.",
                "choices": ["species", "ncbi_code"]
            },
            "out-format": {
                "value": "png",
                "type": "str",
                "help": "Output format of the core figures.",
                "choices": ["png", "svg", "pdf"]
            },
            "min-coocc": {
                "value": 0.30,
                "type": "float",
                "help": "Minimum maximum co-occurrence of two genes to be connected in the graphs."
            },
            "sort-mode": {
                "value": "taxonomy",
                "type": "str",
                "help": "Mode to sort the genomic contexts.",
                "choices": ["taxonomy", "as_input", "tree", "operon", "operon cluster"]
            },
            "overwrite-config": {
                "value": False,
                "type": "bool",
                "help": "Overwrite the argument value in config file with CLI value."
            }
        }