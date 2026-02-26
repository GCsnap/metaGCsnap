import os
import argparse
# pip install rich
from rich.console import Console
from rich.text import Text
from rich.progress import Progress, BarColumn, TextColumn, TimeRemainingColumn
from contextlib import contextmanager
from gcsnap.utils import CustomLogger

import logging

class RichConsole():    
    """
    Mehtods to print messages in the console with different colors and styles.
 
    Attributes:
        console: Console object from rich library.
        color_grey (str): Color for grey messages.
        color_blue (str): Color for blue messages.
        color_purple (str): Color for purple messages.
        color_gold (str): Color for gold messages.
        color_red (str): Color for red messages.
        color_green (str): Color for green messages.
        out_label (str): Output label for the task.
        logger (logging.Logger): Current logger object from logging module.
    """    
    # outlable set for each task, shared among all instances of the class
    out_label = ''  # output label

    def __init__(self, logger_name : str = 'iteration') -> None:
        """
        Initialize the RichConsole object with the console object and colors to use.

        Args:
            logger_name (str): The name of the logger object ('base' or 'iteration'). 
                               Defaults to 'iteration'.
        """            
        self.console = Console()
        self.logger = logging.getLogger(logger_name)

        # colors to use
        self.color_grey = 'gray78'
        self.color_blue = 'steel_blue1'
        self.color_purple = 'plum1'
        self.color_gold = 'light_goldenrod2'
        self.color_red = 'indian_red1'
        self.color_green = 'pale_green1'

    def __getstate__(self) -> dict:
        """
        Used by pickle to convert the object into a dictionary for serialization.
        Console object is not picklable, so it's removed from the state.

        Returns:
            dict: The state of the object without the console.
        """        
        state = self.__dict__.copy()
        del state['console']  # Remove the console as it's not picklable
        return state

    def __setstate__(self, state: dict) -> None:
        """
        Used by pickle to restore the object from a dictionary after deserialization.
        Console object is recreated.

        Args:
            state (dict): The state of the object.
        """        
        self.__dict__.update(state)
        self.console = Console()  # Recreate the Console object      

    @classmethod
    def set_out_label(cls, out_label: str) -> None:
        """
        Sets the out_label attribute to the specified value.

        Args:
            out_label (str): The value to set the out_label attribute to.
        """        
        cls.out_label = out_label

    def print_title(self) -> None:
        """
        Prints the title of the program in the console at the start of the program.
        """        
        color = self.color_blue
        self.print_line(color)
        self.console.print(Text('metaGCsnap', style=f'bold {color}'))
        self.console.print(Text('metaGCsnap is a helps you explore and analyze metagenomic data from a protein-centric perspective.', style=color))
        self.console.print(Text('Thanks for using it! 🇨🇭 ❤️🤍❤️ 🇨🇭', style=color))        
        self.print_line(color)  
        self.logger.info('GCsnap started')

    def print_final(self) -> None:
        """
        Prints the final message of the program in the console at the end of the program.
        """        
        color = self.color_blue
        self.print_line(color)
        self.console.print(Text('🏁  GCsnap finished successfully! 🥳🎆💫', style=color))
        self.print_line(color)  
        CustomLogger.log_to_base('GCsnap finished successfully')                         

    def print_line(self, color: str) -> None:
        """
        Prints a line in the console with the specified color.

        Args:
            color (str): The color of the line.
        """        
        self.console.print(Text('---------------------------------------------', style=color))

    def print_error(self, message: str) -> None:
        """
        Prints an error message in the console with the specified message.

        Args:
            message (str): The error message to print.
        """        
        self.console.print(Text('❌  Error {}.'.format(message), style=self.color_red))   
        self.logger.error(f'{message}')

    def print_warning(self, message: str) -> None:
        """
        Prints a warning message in the console with the specified message.

        Args:
            message (str): The warning message to print.
        """        
        self.console.print(Text('  ⚠️  {} Check gcsnap_{}.log'.format(message, self.out_label), style=self.color_gold))   
        self.logger.warning(f'{message}')        

    def print_skipped_step(self, message: str) -> None:
        """
        Prints a message in the console that a step was skipped.

        Args:
            message (str): The message to print.
        """        
        self.console.print(Text(message, style=self.color_grey))
        self.logger.info(f'{message}')

    def print_step(self, message: str) -> None:
        """
        Prints a message in the console for the current step.

        Args:
            message (str): The message to print.
        """        
        self.console.print(Text('  {}'.format(message), style=self.color_grey))
        self.logger.info(f'{message}')        

    def print_working_on(self, message: str) -> None:
        """
        Prints a message in the console that the program is working on a specific task.

        Args:
            message (str): The message to print.
        """        
        self.console.print(Text('🔨 Working on {}'.format(message), style=self.color_grey))   
        self.logger.info(f'Working on {message}')

    def print_info(self, message: str) -> None:
        """
        Prints an information message in the console.

        Args:
            message (str): The information message to print.
        """        
        self.console.print(Text('\t{}'.format(message), style=self.color_grey))  
        self.logger.info(f'{message}')        
         
    def print_hint(self, message: str) -> None:
        """
        Prints a hint message in the console.

        Args:
            message (str): The hint message to print.
        """        
        self.console.print(Text('  💡  {}'.format(message), style=self.color_gold))   
        self.logger.warning(f'{message}')          

    def print_done(self, message: str) -> None:
        """
        Prints a done message in the console.

        Args:
            message (str): The message to print.
        """        
        # print done message, make first character of message lowercase            
        self.console.print(Text('✅  Done {}'.format(message[0].lower() + 
                                                      message[1:]), style=self.color_green))     
        self.logger.info(f'Done {message}')   

    def print_stop(self) -> None:
        """
        Prints stop message in the console.
        """        
        message = 'Analysis stopped. Check partial output and log file in {}'.format(os.getcwd())
        self.console.print(Text('⏹️  {}'.format(message), style=self.color_blue))  
        self.logger.info(f'{message}')       

    def stop_execution(self, msg: str = None) -> None:
        """
        Stop the execution of the program and logs message to both loggers. 
        If no message is provided, it only logs a default message and stops the execution.

        Args:
            msg (str): The message to print/log. Defaults to None.
        """    
        if msg is not None:
            # log to iteration log
            self.logger.info(msg)
            self.logger.info('Analysis stopped.')
            # print to console and log to base
            self.logger = logging.getLogger('base')
            self.print_warning(msg)
            self.print_stop()            
        else:
            msg = 'Aborted due to error'
            CustomLogger.log_to_base(msg)
        exit(1)     

    @contextmanager
    def status(self, message: str):
        """
        Context manager to print a status a dynamic message with a working icon in the console.
        The completion message is printed when the context manager is exited.

        Args:
            message (str): The message to print.

        Yields:
            None: Yields nothing.
        """        
        self.logger.info(f'{message}')  
        with self.console.status(Text('{} ...'.format(message), style=self.color_grey)):
            yield  
        self.print_done(message)
        
    @contextmanager
    def progress(self, message: str, total: int):    
        """
        Context manager to print a progress bar in the console.

        Args:
            message (str): The message to print.
            total (int): Total amount of work to be done. E.g., total size of file to donwload.

        Yields:
            progress (Progress): The progress object to update the progress bar.
            task_id (int): The task id associated with the progress bar.
        """          
        self.logger.info(f'{message}')  
        with Progress(
            TextColumn("{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeRemainingColumn(),
            console=self.console
        ) as progress:
            task_id = progress.add_task('[{}]\t{}'.format(self.color_grey, message), total=total)
            yield progress, task_id        
        self.print_done(message)

    def print_help(self, parser: argparse.PARSER) -> None:
        """
        Prints the help message for all arguments as defined in the configuration file
        with the usage, optional arguments and epilog. Shown when using GCsnap --help.

        Args:
            parser (CustomArgumentParser): The argument parser object.
        """              
        usage = 'GCsnap --targets <targets> [Optional arguments]'
        epilog = 'Example: GCsnap --targets PHOL_ECOLI A0A0U4VKN7_9PSED A0'

        usage = Text(parser.format_usage(), style=self.color_grey)
        epilog = Text('\n' + parser.epilog, style=self.color_grey)

        # Print usage
        self.console.print(usage)
        # print config.yaml message
        self.console.print(Text('If CLI arguments not sepcified, default value from config.yaml', style=self.color_green))
        self.console.print(Text('--overwrite-config to replace values with CLI Arguments', style=self.color_green))
        self.console.print(Text('Or change them manually in config.yaml\n', style=self.color_green))

        # Print --targets separately
        for action in parser._actions:
            if '--targets' in action.option_strings:
                self.print_argument(action)
                break

        # Optional arguments
        self.console.print(Text('\nOptional arguments:', style=self.color_grey))
        # Print arguments
        for action in parser._actions:
            # skip --help and --targets
            if '--targets' in action.option_strings or '--help' in action.option_strings:
                continue
            self.print_argument(action)

        # Print epilog
        self.console.print(epilog)
        
    def print_argument(self, action: argparse.Action) -> None:
        """
        Prints the formated help message for an argument in the console based
        on the argument action object including the option strings, default value and help message
        as defined in the configuration file and read by the CustomArgumentParser.

        Args:
            action (argparse.Action): The argument action object.
        """            
        option_strings = ", ".join(action.option_strings)
        metavar = action.metavar or ""
        argument_help = Text(f"{option_strings} {metavar}", style=self.color_blue)
        if action.help:
            argument_help.append(f"\n    {action.help}", style=self.color_grey)

        if action.choices:
            choices_str = ', '.join(map(str, action.choices))
            argument_help.append(f'\n    Supported values: [{choices_str}]', style=self.color_grey)        

        if action.default != argparse.SUPPRESS:
            default_value = f"\n    [config.yaml value: {action.default}]"
            argument_help.append(default_value, style=self.color_green)

        self.console.print(argument_help)

