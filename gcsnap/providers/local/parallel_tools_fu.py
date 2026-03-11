from typing import Callable

from concurrent.futures import as_completed, ProcessPoolExecutor
from mpi4py.futures import MPIPoolExecutor
from mpi4py import MPI

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole

class ParallelTools:

    """
    Methods and attributes to handle parallel processing. The class is initialized with a Configuration object.

    Attributes:
        config (Configuration): The Configuration object with the parsed arguments.
        n_nodes (int): The number of nodes to use for parallel processing.
        n_cpu (int): The number of CPUs per node to use for parallel processing.
        workers (int): The number of workers to use for parallel processing with mpoi4py.
        memory_per_node (int): The memory per node to use for parallel processing.
        console (RichConsole): The RichConsole object.
    """

    def __init__(self, config: Configuration):
        """
        Initialize the ParallelTools object.

        Args:
            config (Configuration): The Configuration object with the parsed arguments.
        """        
        self.config = config
        self.n_nodes = config.arguments['n_nodes']['value']
        self.n_cpu = config.arguments['n_ranks_per_node']['value']
        self.cluster = None

         # for mpi
        self.workers = (self.n_nodes * self.n_cpu)
            
        self.console = RichConsole()

        # Assign this instance to the class-level variable
        # Now self is available in the static method as on class level
        # usually instance is not on class level
        ParallelTools._instance = self

    @staticmethod
    def parallel_wrapper(parallel_args: list[tuple], func: Callable) -> list:
        """
        A static method that calls the process_wrapper method of the stored instance.
        This is called from the modules using parallel processing.

        Args:
            parallel_args (list[tuple]):A list of tuples, where each tuple contains the arguments for the function.
            func (Callable): The function to apply to the arguments.

        Returns:
            list: A list of results from the function applied to the arguments in the order they are provided.            
        """
        if ParallelTools._instance is None:
            raise RuntimeError("ParallelTools instance has not been initialized.")
        
        return ParallelTools._instance.mpiprocess_wrapper(parallel_args, func)       

    def mpiprocess_wrapper(self, parallel_args: list[tuple], func: Callable) -> list:
        """
        Apply a function to a list of arguments using ProcessPoolExecutor. The arguments are passed as tuples
        and are unpacked within the function. As completed is used to get the results in the order they finish.

        Args:
            parallel_args (list[tuple]): A list of tuples, where each tuple contains the arguments for the function.
            func (Callable): The function to apply to the arguments.

        Returns:
            list: A list of results from the function applied to the arguments in the order they finish.
        """    
        # Same as with ProcessPoolExecutor from cuncurrent.futures
        # https://mpi4py.readthedocs.io/en/stable/mpi4py.futures.html#parallel-tasks
        # with MPIPoolExecutor(max_workers = workers) as executor:
        with MPIPoolExecutor(max_workers = self.workers) as executor:        
            # get numbers of worker
            #print('Number of workers given: {}'.format(executor.num_workers))
            #print('Number of workers asked: {}'.format(workers))
            futures = [executor.submit(func, arg) for arg in parallel_args]
            result_list = [future.result() for future in as_completed(futures)]

        return result_list