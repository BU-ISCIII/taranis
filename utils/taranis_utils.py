#!/usr/bin/env python3
import logging
#from logging.config import fileConfig
#from logging.handlers import RotatingFileHandler
import os
import re
import glob
import shutil
import subprocess
from Bio import SeqIO
from Bio import Seq
from Bio.Alphabet import generic_dna

from openpyxl import load_workbook
import pandas as pd

from taranis_configuration import *



def check_prerequisites (pre_requisite_list):
    '''
    Description:
        The function check if  the external software
        has the right version 
    Functions:
        check_program_is_exec_version # located at utils.taranis_utils 
    Variable:
        pre_requisite_list  # tupla list containing software and version
        experiment_name # contains the experiment name from the sample sheet
        library_name  # contains the library name from the sample sheet
    Return:
        True if all checking are successful False if any of the check fails
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function check_prerequisites')
    # check if software is installed and has the minimum version
    for program, version in pre_requisite_list :
        if not check_program_is_exec_version (program , version):
            logger.debug ('End function check_prerequisites with error')
            return False
    logger.info('')
    logger.debug ('End function check_prerequisites')
    return True


def check_program_is_exec_version (program, version):
    '''
    Description:
        The function will check if the program is installed in your
        system and if the version installed matched or it is higher
        with pre-requisites
    Input:
        program    # Is the program name 
        version    # version of the software that was tested
    Return:
        False is return in case that version is below 
        True  if equal version or higher
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function check_program_is_exec_version')
    if shutil.which(program) is not None :
        # check version
        version_str= str(subprocess.check_output([program , '-version']))
        if version_str == "b''" :
            version_str = subprocess.getoutput( str (program + ' -version'))
        if not re.search(version, version_str):
            v_str = re.search(".*:\s+(\d\.\d).*", version_str)
            v_float = float(v_str.groups(1)[0])
            if v_float > float(version) :
                logger.info('Found a higher version in the system')
                logger.debug ('End function check_program_is_exec_version')
                return True
            else:
                string_message = program + ' require version ' + version + 'but get ' + version_str
                logging_errors(string_message, False, True)
                return False
        logger.debug ('End function check_program_is_exec_version')
        return True
    else:
        logger.info('Cannot find %s installed on your system', program)
        logger.debug ('End function check_program_is_exec_version with error')
        return False



def logging_errors(string_text, showing_traceback , print_on_screen ):
    '''
    Description:
        The function will log the error information to file.
    Input:
        print_on_screen # Boolean to print warning on screen
        string_text # information text to include in the log

    Variables:
        subject # text to include in the subject email
    '''
    logger = logging.getLogger(__name__)
    logger.error('-----------------    ERROR   ------------------')
    logger.error(string_text )
    if showing_traceback :
        logger.error('Showing traceback: ',  exc_info=True)
    logger.error('-----------------    END ERROR   --------------')
    if print_on_screen :
        from datetime import datetime
        print('********* ERROR **********')
        print(string_text)
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        print('When processing run . Check log for detail information')
        print('******* END ERROR ********')
    return ''
    
def logging_warnings(string_text, print_on_screen ):
    '''
    Description:
        The function will log the error information to file.
        Optional can send an email to inform about the issue
    Input:
        print_on_screen # Boolean to print warning on screen
        string_text # information text to include in the log
    '''
    logger = logging.getLogger(__name__)
    logger.warning('-----------------    WARNING   ------------------')
    logger.warning(string_text )
    logger.warning('-----------------    END WARNING   --------------')
    if print_on_screen :
        from datetime import datetime
        print('******* WARNING ********')
        print(string_text)
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        print('When processing run . Check log for detail information')
        print('**** END WARNING *******')
    return ''


def parsing_fasta_file_to_dict (fasta_file):
    '''
    Description:
        The function get the fasta file and converts it to dictionary.
    Input:
        fasta_file  # fasta file to convert to dictionary
    Import:
        SeqIO.parse
    Variable:
        fasta_dict  # dictionary containing fasta file
    Return:
        fasta_dict 
    '''
    logger = logging.getLogger(__name__)
    logger.debug('Starting the function validate_sample_sheet' )
    
    fasta_dict = {}
    logger.info('Converting  %s to dictionary', fasta_file)
    for contig in SeqIO.parse(fasta_file, "fasta", generic_dna):
        fasta_dict[contig.id] = str(contig.seq.upper())
    logger.debug('End the function validate_sample_sheet' )
    return fasta_dict

def read_xls_file (in_file, logger):
    '''
    Description:
        This function open the Excel file enter by the user in the xlsfile parameter
        Once the excel is read the column information of the gene and protein is
        stored on the gene_list that it is returned back
    Input:
        logger      # Is the logging object
        in_file     # It is the excel file which contains the information to parse
    Variables:
        wb      # Contains the excel workbook object
        ws      # Contains the working sheet object of the workbook
        gene    # Used in the interaction row to get the gene name
        protein # Used in the interaction row to get the protein name
        gene_prot   # Used to get the tupla (gene/protein) for each or excel row
        genes_prots_list  # Is a list containing tuplas of gene, protein
    Return:
        'Error message' is returned in case excel file does not exists
        genes_prots_list is returned as a successful execution
    '''
    logger = logging.getLogger(__name__)
    logger.debug ('Starting function read_xls_file')
    logger.debug('opening the excel file : %s', in_file)
    try:
        wb = load_workbook(in_file)
        logger.info('Excel file has been read and starts processing it.')
    except Exception as e:
        logger.error('-----------------    ERROR   ------------------')
        logger.error('Unable to open the excel file.  %s ', e )
        logger.error('Showing traceback: ',  exc_info=True)
        logger.error('-----------------    END ERROR   --------------')
        #raise
        return 'Error: Unable to open excel file'
    # Only fetch the first working sheet
    ws = wb[wb.sheetnames[0]]

    genes_prots_list = []
    ## Get the content block from A2 : B(latest row in the excel)
    for row in ws.iter_rows(min_row=2, min_col=1, max_row=ws.max_row, max_col=2) :
        gene_prot = []
        for index in range(len(row)) :
            gene_prot.append(row[index].value)
        genes_prots_list.append(gene_prot)
    logger.info('Exiting the function ---read_xls_file-- ')
    logger.info('Returning back the gene/protein list' )
    logger.debug ('End function read_xls_file')
    return genes_prots_list

def download_fasta_locus (locus_list, output_dir, logger):
    '''
    Description:
        This function will download the protein sequence.
        Then it will be translated to nucleotide and saved
        in the output directory specified by the users.
    Input:
        gene_list
        filename    # Is the name of the file to be checked
        logger      # is the logging object to logging information
    Return:
        Error is return in case that file does not exists
        True  if file exists
    '''
    download_counter = 0
    for loci in locus_list :
        tmp_split = loci.split('/')
        loci_name = tmp_split[-1]
        r = requests.get(loci + '/alleles_fasta')
        if r.status_code != 200 :
            logger.error('Unable to download the fasta file  for allele %s ', loci_name)

        else :
            fasta_alleles = r.text
            fasta_file =  os.path.join(output_dir, str(loci_name + '.fasta'))
            with open (fasta_file , 'w') as fasta_fh :
                fasta_fh.write(fasta_alleles)
            download_counter += 1
    if download_counter == len(locus_list) :
        return True
    else :
        logger.info('All alleles have been successfully downloaded and saved on %s', output_dir)
        return False



def check_if_file_exists (filename, logger):
    '''
    Description:
        This function will check if the file exists
    Input:
        filename    # Is the name of the file to be checked
        logger      # is the logging object to logging information
    Return:
        Error is return in case that file does not exists
        True  if file exists
    '''
    if not os.path.isfile(filename):
        logger.info('File  %s , does not exists', filename)
        return 'Error'
    return True



def create_blastdb (file_name, db_name,db_type, logger ):
    f_name = os.path.basename(file_name).split('.')
    db_dir = os.path.join(db_name,f_name[0])
    output_blast_dir = os.path.join(db_dir, f_name[0])
    if not os.path.exists(db_dir):
        try:
            os.makedirs(db_dir)
            logger.debug(' Created local blast directory for Core Gene %s', f_name[0])
        except:
            logger.info('Cannot create directory for local blast database on Core Gene file %s' , f_name[0])
            print ('Error when creating the directory %s for blastdb. ', db_dir)
            exit(0)

        blast_command = ['makeblastdb' , '-in' , file_name , '-parse_seqids', '-dbtype',  db_type, '-out' , output_blast_dir]
        blast_result = subprocess.run(blast_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if blast_result.stderr:
            logger.error('cannot create blast db for %s ', f_name[0])
            logger.error('makeblastdb returning error code %s', blast_result.stderr)
            return False

    else:
        logger.info('Skeeping the blastdb creation for %s, as it is already exists', f_name[0])
    return True

def check_blast (reference_allele, sample_files, db_name, logger) :
    for s_file in sample_files:
        f_name = os.path.basename(s_file).split('.')
        dir_name = os.path.dirname(s_file)
        blast_dir = os.path.join(dir_name, db_name,f_name[0])
        blast_db = os.path.join(blast_dir,f_name[0])
        if not os.path.exists(blast_dir) :
            logger.error('Blast db folder for sample %s does not exist', f_name)
            return False
        cline = NcbiblastnCommandline(db=blast_db, evalue=0.001, outfmt=5, max_target_seqs=10, max_hsps=10,num_threads=1, query=reference_allele)
        out, err = cline()

        psiblast_xml = StringIO(out)
        blast_records = NCBIXML.parse(psiblast_xml)

        for blast_record in blast_records:
            locationcontigs = []
            for alignment in blast_record.alignments:
                # select the best match
                for match in alignment.hsps:
                    alleleMatchid = int((blast_record.query_id.split("_"))[-1])
    return True





def junk ():
    '''
    Description:
        The function will check if the program is installed in your
        system and if the version installed matched or it is higher
        with pre-requisites
    Input:
        program    # Is the program name 
        version    # version of the software that was tested
    Return:
        False is return in case that version is below 
        True  if equal version or higher
    '''
    
    
    
    
    AA_codon = {
            'C': ['TGT', 'TGC'],
            'A': ['GAT', 'GAC'],
            'S': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
            'G': ['CAA', 'CAG'],
            'M': ['ATG'], #Start
            'A': ['AAC', 'AAT'],
            'P': ['CCT', 'CCG', 'CCA', 'CCC'],
            'L': ['AAG', 'AAA'],
            'Q': ['TAG', 'TGA', 'TAA'], #Stop
            'T': ['ACC', 'ACA', 'ACG', 'ACT'],
            'P': ['TTT', 'TTC'],
            'A': ['GCA', 'GCC', 'GCG', 'GCT'],
            'G': ['GGT', 'GGG', 'GGA', 'GGC'],
            'I': ['ATC', 'ATA', 'ATT'],
            'L': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
            'H': ['CAT', 'CAC'],
            'A': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
            'T': ['TGG'],
            'V': ['GTA', 'GTC', 'GTG', 'GTT'],
            'G': ['GAG', 'GAA'],
            'T': ['TAT', 'TAC'] }
    return True




def is_fasta_file (file_name):
    '''
    Description:
        The function will check if file has the fasta format 
    Input:
        file_name    # file name to check
    Return:
        False in case that does not have a fasta format 
        True  if file is fasta
    '''
    logger = logging.getLogger(__name__)
    logger.debug('Starting the function is_fasta_file' )
    logger.info('Reading file %s', file_name)
    with open (file_name, 'r') as fh:
        fasta = SeqIO.parse(fh, 'fasta')
    
    logger.debug('End the function is_fasta_file' )
    return any(fasta)

def get_fasta_file_list (check_directory):
    '''
    Description:
        The function will get the list of the fasta files in the
        directory
    Input:
        check_directory    # Is the program name 
    Functions:
        valid_files     # located at this file
    Variable:
        list_filtered_files # list of files that has the extension fasta
        valid_files     # list containing all fasta files in directory 
    Return:
        False if director does not exist or not fasta files in directory.
        valid_files list with fasta files in the directory 
    '''
    
    
    logger = logging.getLogger(__name__)
    logger.debug('Starting the function get_fasta_file_list' )
    if not os.path.isdir(check_directory):
        string_message = 'directory ' + check_directory + ' does not exists'
        logging_errors(string_message, False, False)
        logger.debug('End the function get_fasta_file_list with error' )
        return False
    filter_files = os.path.join(check_directory, '*.fasta')
    list_filtered_files =  glob.glob(filter_files)
    list_filtered_files.sort()
    if len (list_filtered_files) == 0 :
        string_message = 'directory ' + check_directory + ' does not have any fasta file'
        logging_errors(string_message, False, True)
        logger.debug('End the function get_fasta_file_list with error' )
        return False
    valid_files = []
    for file_name in list_filtered_files:
        if is_fasta_file( file_name):
            valid_files.append(file_name)
        else:
            logger.info('Ignoring file  %s .Does not have a fasta format', file_name)
    if len(valid_files) == 0:
        logger.info('There are not valid fasta files in the directory %s', check_directory)
        logger.debug('Files in the directory are:  $s', list_filtered_files)
        return False
    else:
        logger.debug('Starting the function get_fasta_file_list' )
        return valid_files

def check_sequence_order(allele_sequence, logger) :
    start_codon_forward= ['ATG','ATA','ATT','GTG', 'TTG']
    start_codon_reverse= ['CAT', 'TAT','AAT','CAC','CAA']
    # check forward direction
    if allele_sequence[0:3] in start_codon_forward :
        return 'forward'
    if allele_sequence[len(allele_sequence) -3: len(allele_sequence)] in start_codon_reverse :
        return 'reverse'
    return "Error"

def hamming_distance (pd_matrix):
    '''
    The function is used to find the hamming distance matrix
    Input:
        pd_matrix    # Contains the panda dataFrame
    Variables:
        unique_values   # contains the array with the unique values in the dataFrame
        U   # Is the boolean matrix of differences
        H   # It is accumulative values of U
    Return:
       H where the number of columns have been subtracted
    '''

    unique_values = pd.unique(pd_matrix[list(pd_matrix.keys())].values.ravel('K'))
    # Create binary matrix ('1' or '0' ) matching the input matrix vs the unique_values[0]
    # astype(int) is used to transform the boolean matrix into integer
    U = pd_matrix.eq(unique_values[0]).astype(int)
    # multiply the matrix with the transpose
    H = U.dot(U.T)

    # Repeat for each unique value
    for unique_val in range(1,len(unique_values)):
        U = pd_matrix.eq(unique_values[unique_val]).astype(int)
        # Add the value of the binary matrix with the previous stored values
        H = H.add(U.dot(U.T))

    return len(pd_matrix.columns) - H


def create_distance_matrix (input_dir, input_file):
    '''
    Description:
        The function will check if the program is installed in your
        system and if the version installed matched or it is higher
        with pre-requisites
    Input:
        program    # Is the program name 
        version    # version of the software that was tested
    Return:
        False is return in case that version is below 
        True  if equal version or higher
    '''
    logger = logging.getLogger(__name__)
    logger.debug('Starting the function create_distance_matrix' )
    
    try:
        result_file = os.path.join(input_dir, input_file)
        pd_matrix = pd.read_csv(result_file, sep='\t', header=0, index_col=0)
    except Exception as e:
        string_message = 'Unable to open the matrix distance file'
        logging_errors (string_message, False, True)
        logger.debug('End the function create_distance_matrix with error' )
        return 'Error'

    distance_matrix = hamming_distance (pd_matrix)
    out_file = os.path.join(input_dir, 'matrix_distance.tsv')
    try:
        distance_matrix.to_csv(out_file, sep = '\t')
    except Exception as e:
        string_message = 'Unable to create the matrix distance file'
        logging_errors (string_message, False, True)
        logger.debug('End the function create_distance_matrix with error' )
        return 'Error'
    logger.debug('End the function create_distance_matrix' )
    return True

def open_log(log_name):
    '''
    Description:
        This function open the log file with the configuration defined
        on the config file (loging_config.ini)
        The path for the logging config is defined on the application
        configuration file.
    Input:
        log_name    # Is the name that will be written inside the logfile
    Variables:
        log_folder  # directory extracted from log_name to create the folder
        
    Return:
        Error is return in case that config file does not exists
        logger # containing the logging object
    '''
    logger = logging.getLogger(__name__) 
    logging.basicConfig(filename=log_name, format='%(asctime)s %(funcName)-12s %(levelname)-8s %(lineno)s %(message)s')
    
    handler = logging.StreamHandler()
    #logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)
    try:
        
        log_folder = os.path.dirname(log_name)
        if not os.path.isdir(log_folder) :
            os.makedirs(log_folder)
        
        logger.info('--------------- LOG FILE -----------------')
        logger.info('Log file has been created for process %s', log_name)
    except:
        print('------------- ERROR --------------')
        print('Unable to create the logging file')
        print('check that ', log_folder , ' path has write permissions')
        print('------------------------------------------')
        raise 
    
    return logger

def write_first_allele_seq(fasta_file, full_path_first_allele):
    '''
    Description:
        The function get the fasta file to save the first allele in
        the first allele temporary directory 
    Input:
        fasta_file      # fasta file to be processed
        full_path_first_allele  # full path for saving first allele
    Import:
        SeqIO
    Variable:
        f_name      # fasta file name without extension
        fasta_file  # full path to store the file
    Return:
        fasta_file
    '''
    # split file_sequence into directory and filename
    f_name = os.path.basename(fasta_file)

    first_record = SeqIO.parse(file_sequence, "fasta").__next__()
    # build the fasta file name to store under first_allele_firectory
    fasta_file = os.path.join(full_path_first_allele, f_name)
    SeqIO.write(first_record, fasta_file, "fasta")

    return fasta_file

