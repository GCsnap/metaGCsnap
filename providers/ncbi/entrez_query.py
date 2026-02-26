import requests
import time
import datetime
import random
import numpy as np
import xml.etree.ElementTree as ET
import copy

from gcsnap.configuration import Configuration
from gcsnap.rich_console import RichConsole 

from gcsnap.utils import processpool_wrapper
from gcsnap.utils import WarningToLog

import logging
logger = logging.getLogger('iteration')

class EntrezQuery:
    """
    Methods and attributes to query the NCBI Entrez database for information about proteins or taxonomies.

    Attributes:
        cores (int): The number of CPU cores to use.
        api_key (str): The NCBI API key.
        ncbi_email (str): The email address for the NCBI API key.
        console (RichConsole): The RichConsole object to print messages.
        db (str): The database to query.
        rettype (str): The return type of the query.
        retmode (str): The return mode of the query.
        search_codes (list): The list of search codes.
        logging (bool): If True, log the results of the query.
        msg (str): The message to display during the query.
        function (function): The function to run depending on the rettype.
        chunk_size (int): The size of the chunks to query.
    """

    def __init__(self, config: Configuration, search_codes: list, db: str = 'protein', 
                 rettype: str = 'ipg', retmode: str = 'xml', logging: bool = True):
        """
        Initialize the EntrezQuery object.

        Args:
            config (Configuration): The Configuration object containing the arguments.
            search_codes (list): The list of search codes.
            db (str, optional): The database to query. Defaults to 'protein'.
            rettype (str, optional): The return type of the query. Defaults to 'ipg'.
            retmode (str, optional): The return mode of the query. Defaults to 'xml'.
            logging (bool, optional): If True, log what is not found. Defaults to True.
        """        
        # get necessary configuration arguments        
        self.cores = config.arguments['n_cpu']['value']
        self.api_key = config.arguments['ncbi_api_key']['value']
        self.ncbi_email = config.arguments['ncbi_user_email']['value']
        # Information for the user regarding api keys: 
        # https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.chapter2_table1
        # Pargaraph: API Keys

        self.console = RichConsole()

        # set arguments
        self.db = db
        self.rettype = rettype
        self.retmode = retmode
        self.search_codes = search_codes
        self.logging = logging

        # set correct function depending on rettype
        if db == 'protein' and rettype == 'ipg':
            self.msg = 'Find identical protein groups'
            self.function = self.run_accession
            self.chunk_size = 50
        elif db == 'protein' and rettype == 'fasta':
            self.msg = 'Download sequences'
            self.function = self.run_sequence
            self.chunk_size = 100
        elif db == 'taxonomy':
            self.msg = 'Download taxonomies'
            self.function = self.run_taxonomy
            self.chunk_size = 100            

    def run(self) -> dict:
        """
        Run the query to the NCBI Entrez database.

        Returns:
            dict: The dictionary with the results of the query.
        """        
        with self.console.status(self.msg):
            chunk_list = self.create_chunks(self.chunk_size)
            # returns a result tuple with (found_dict, None) or ({}, (WarningToLog, chunk))
            results_list = processpool_wrapper(self.cores, chunk_list, self.function)
            # combine all dictionaries
            found = {k: v for result_tuple in results_list for k, v 
                     in result_tuple[0].items() if result_tuple[1] is None}
            # get those with failed entrez requests:[(error, chunk)]
            errors = [(result_tuple[1][0], result_tuple[1][1]) 
                      for result_tuple in results_list 
                      if result_tuple[1] is not None]
 
        # log errors if desired
        failed_request_ids = []
        if errors and self.logging:
            failed_request_ids = self.log_request_error(errors)

        # log not found if desired
        not_found = [id for id in self.search_codes if id not in found.keys() 
                     and id not in failed_request_ids]
        if not_found and self.logging:
            self.log_not_found(not_found)

        return found
        
    def run_accession(self, chunk: list) -> tuple[dict,tuple]:
        """
        Run the accession query to the NCBI Entrez database.

        Args:
            chunk (list): The list of search codes to query.

        Returns:
            tuple[dict,tuple]: The result dictionary and the error tuple.
        """        
        try:
            root = self.query_chunk(chunk)
            return (self.extract_accessions(chunk, root) , None)       
        except WarningToLog as e:
            return ({}, (e, chunk))      
    
    def run_sequence(self, chunk: list) -> tuple[dict,tuple]:
        """
        Run the sequence query to the NCBI Entrez database

        Args:
            chunk (list): The list of search codes to query.

        Returns:
            tuple[dict,tuple]: The result dictionary and the error tuple.
        """        
        try:
            root = self.query_chunk(chunk)
            return (self.extract_sequences(root), None)
        except WarningToLog as e:
            return ({}, (e, chunk))
        
    def run_taxonomy(self, chunk: list) -> tuple[dict,tuple]:
        """
        Run the taxonomy query to the NCBI Entrez database

        Args:
            chunk (list): The list of search codes to query.

        Returns:
            tuple[dict,tuple]: The result dictionary and the error tuple.
        """        
        try:
            root = self.query_chunk(chunk)
            return (self.extract_taxonomy(root), None)
        except WarningToLog as e:
            return ({}, (e, chunk))

    def create_chunks(self, chunk_size: int = 50) -> list:
        """
        Create chunks of the search codes to query.

        Args:
            chunk_size (int, optional): The size of the chunks. Defaults to 50.

        Returns:
            list: The list of chunks.
        """        
        return [self.search_codes[i:i + chunk_size] for i in range(0, len(self.search_codes), chunk_size)]

    def get_ncbi_base_url(self) -> str:
        """
        Get the base URL for the NCBI Entrez database.

        Returns:
            str: The base URL for the NCBI Entrez database.
        """        
        return 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db='+self.db

    def query_chunk(self, chunk: list) -> ET.Element:
        """
        Query the NCBI Entrez database with a chunk of search codes.
        The API limit is 3 request per second without an API key and 10 request per second with an API key.
        Information about api keys: https://www.ncbi.nlm.nih.gov/books/NBK25497/

        Args:
            chunk (list): The list of search codes to query.

        Raises:
            WarningToLog: If the query fails.

        Returns:
            ET.Element: The root element of the XML response.
        """        
        # information about eutils email and api key
        # https://www.ncbi.nlm.nih.gov/books/NBK25497/

        ids = '&id=' + ','.join(chunk) # works for single ids as well
        rt = '&rettype=' + self.rettype
        rm = '&retmode=' + self.retmode

        # add email if provided
        if self.api_key is not None:
            email = '&email=' + self.ncbi_email
        else:
            email = ''

        # add api key if provided
        if self.api_key is not None:
            api_key = '&api_key=' + self.api_key
        else:
            api_key = ''

        # set timeout depending on chunk size
        #timeout = len(chunk_list) * 0.2
        timeout = None

        url = self.get_ncbi_base_url() + ids + rt + rm + email + api_key

        # Entrez has an API limit of 3 request per second without an API key
        # and 10 request per second with an API key.
        initial_wait_time = 1  # Initial wait time in seconds
        max_wait_time = 10  # Maximum wait time in seconds
        factor = 2  # Exponential backoff factor
        retries = 3  # Number of retries on failure
        retry_delay = 2  # Delay between retries

        attempt = 0  # Track the number of attempts

        while True:
            response = requests.get(url, timeout=None)                
            if response.status_code == 200 and attempt < retries:  # Successful request
                attempt += 1
                # Check content type for XML
                content_type = response.headers.get('Content-Type', '').lower()
                if 'xml' not in content_type:
                    time.sleep(retry_delay)
                    continue  # Retry      
                # empty reply or other parsing error          
                xml_string = response.text.strip()
                if not xml_string:
                    time.sleep(retry_delay)
                    continue  # Retry
                try:
                    # Try parsing the XML
                    return ET.fromstring(xml_string)
                except ET.ParseError as e:
                    time.sleep(retry_delay)
                    continue  # Retry
            elif response.status_code == 429:  # Too Many Requests
                # don't count this as attempt
                # extract the time with wait suggestion from the response
                retry_after = response.headers.get('Retry-After')
                if retry_after:
                    try:
                        retry_after = int(retry_after)
                    except ValueError:
                        # If Retry-After is not a number, it might be a date
                        retry_after_date = requests.utils.parse_date(retry_after)
                        if retry_after_date:
                            wait_time = (retry_after_date - datetime.datetime.now()).total_seconds()
                    wait_time = int(retry_after)
                else:
                    # Decorrelated Jitter
                    wait_time = min(random.uniform(initial_wait_time, initial_wait_time * factor), max_wait_time)
                time.sleep(wait_time)
            else:
                attempt += 1
                if attempt < retries:
                    time.sleep(retry_delay)
                    continue
                else:
                    raise WarningToLog('Entrez request failed in {} attempts with status code {}'.format(retries, response.status_code))
  
    def extract_accessions(self, chunk_list: list, root: ET.Element) -> dict:
        """
        Extract the accessions from the XML response.

        Args:
            chunk_list (list): The list of search codes.
            root (ET.Element): The root element of the XML response.

        Returns:
            dict: The dictionary with the accessions.
        """        
        reports = root.findall('IPGReport') # one report for each id
        # ProteinList: full identical protein group for each id
        ipg_list = [(reports[i].find('ProteinList'), chunk_list[i]) for i in range(len(reports))
                    if reports[i].find('ProteinList') is not None]

        return_dict = {}
        for ipg, ncbi_code in ipg_list: # iterate though all ipgs
            protein_list = ipg.findall('Protein') # all proteins from the ipg group
            
            # all individual CDS of the full sequnce of the protein
            cds_lists = [(protein.find('CDSList').findall('CDS'), protein.get('accver')) 
                         for protein in protein_list if protein.find('CDSList') is not None]
            # flatten and exclude those without assembly accession
            cds_list = [(item,sublist[1]) for sublist in cds_lists for item in sublist[0]
                        if item.get('assembly') is not None]

            assembly_acc = None
            for cds, accver in cds_list:    
                if ncbi_code == accver: 
                    # the assembly accession for the desired ncbi_code
                    assembly_acc_tmp = cds.get('assembly')
                    if assembly_acc_tmp.startswith('GCF'):
                        # if the assembly accession is a GCF (RefSeq), use it directly
                        return_dict[ncbi_code] = assembly_acc_tmp
                        break
                    # if the assembly accession is a GCA (GenBank), continue search
                    assembly_acc = assembly_acc_tmp

            # if still None, don't add and catch later by logger
            if assembly_acc is not None:
                return_dict[ncbi_code] = assembly_acc
        return return_dict
    
    def extract_sequences(self, root: ET.Element) -> dict:
        """
        Extract the sequences from the XML response.

        Args:
            root (ET.Element): The root element of the XML response.

        Returns:
            dict: The dictionary with the sequences.
        """        
        # create dictionary for all found sequences
        sequences = {report.find('TSeq_accver').text : {'seq': report.find('TSeq_sequence').text,
                                'taxid': report.find('TSeq_taxid').text,
                                'species': report.find('TSeq_orgname').text} 
                    for report in root.findall('TSeq')}        
        return sequences
    
    def extract_taxonomy(self, root: ET.Element) -> dict:
        """
        Extract the taxonomy from the XML response.

        Args:
            root (ET.Element): The root element of the XML response.

        Returns:
            dict: The dictionary with the taxonomy in nested form representing the hierarchy.
        """        
        taxonomy_template = {'superkingdom' :  None,
                    'clade': None,
                    'phylum': None,
                    'class': None,
                    'order': None,
                    'family': None,
                    'genus': None}  

        # copy.deepcopy(taxonomy_template) ensures a new independent copy of the 
        # taxonomy_template for each taxon
        # {**copy.deepcopy(taxonomy_template), **updated_ranks} combines the 
        # copied taxonomy_template with the ranks dictionary using dictionary unpacking (**)      
        return {
            taxon.find('TaxId').text: {
                **copy.deepcopy(taxonomy_template),
                **{rank.find('Rank').text: rank.find('ScientificName').text 
                for rank in taxon.find('LineageEx') if rank.find('Rank').text != 'no rank'}
            }
            for taxon in root.findall('Taxon')
        }
    
    def log_not_found(self, not_found: list) -> None:
        """
        Log the search codes that were not found in the NCBI Entrez database.

        Args:
            not_found (list): The list of search codes that were not found.
        """        
        message = '{} entrez queries returned no results for db "{}" with rettype "{}"'.format(len(not_found), self.db, self.rettype)
        self.console.print_warning(message)
        for id in not_found:
            logger.warning('No entrez entry found for {} with db "{}" and rettype "{}"'.format(id, self.db, self.rettype)) 

    def log_request_error(self, errors: list) -> list:
        """
        Log the errors of an NCBI Entrez requests.

        Args:
            errors (list): The list of errors.

        Returns:
            list: The list of search codes that failed.
        """        
        message = '{} entrez requests failed for db "{}" with rettype "{}"'.format(len(errors), self.db, self.rettype)
        self.console.print_warning(message)
        # all filed ncbi_codes
        all_failed = []
        for error, chunk in errors:
            logger.warning(error)
            for each in chunk:
                logger.warning('Request failed for {}'.format(each))
                all_failed.append(each)
        return all_failed