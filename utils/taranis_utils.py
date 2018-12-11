#!/usr/bin/env python3
import logging
from logging.config import fileConfig
#from logging.handlers import RotatingFileHandler
import os
import re
import glob
import shutil
#import subprocess
from Bio import SeqIO
from Bio import Seq
from taranis_configuration import *
def open_log(log_name):
    '''
    Description:
        This function open the log file with the configuration defined
        on the config file (loging_config.ini)
        The path for the logging config is defined on the application
        configuration file.
    Input:
        log_name    # Is the name that will be written inside the logfile
        LOGGIN_CONFIGURATION # is the constant value defined on the configuration
                            file of the application
    Return:
        Error is return in case that config file does not exists
        logger # containing the logging object
    '''
    #working_dir = os.getcwd()
    

    #fileConfig('/srv/taranis/logging_config.ini')
    #log_name=os.path.join(working_dir, log_name)
    #def create_log ():
    #logger = logging.getLogger(__name__)
    #logger.setLevel(logging.DEBUG)
    #create the file handler
    #handler = logging.handlers.RotatingFileHandler('pepe.log', maxBytes=4000000, backupCount=5)
    #handler.setLevel(logging.DEBUG)

    #create a Logging format
    #formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    #handler.setFormatter(formatter)
    #add the handlers to the logger
    #logger.addHandler(handler)
    try:
        logging.config.fileConfig(LOGGING_CONFIGURATION)
        logger = logging.getLogger(log_name)
        logger.info('--------------- LOG FILE -----------------')
        logger.info('Log file has been created for process %s', log_name)
    except:
        print('------------- ERROR --------------')
        print(' Unable to create the logging file')
        print('Check in the logging configuration file')
        print('that the path to store the log file exists')
        print('------------------------------------------')
        return 'Error'
    return logger

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


def check_program_is_exec_version (program, version, logger):
    # The function will check if the program is installed in your system and if the version
    # installed matched with the pre-requisites
    if shutil.which(program) is not None :
        # check version
        version_str= str(subprocess.check_output([program , '-version']))
        if version_str == "b''" :
            version_str = subprocess.getoutput( str (program + ' -version'))
        if not re.search(version, version_str):
            logger.info('%s program does not have the right version ', program)
            print ('Exiting script \n, Version of ' , program, 'does not fulfill the requirements')
            return False
        return True
    else:
        logger.info('Cannot find %s installed on your system', program)
        return False

def is_fasta_file (file_name):
    with open (file_name, 'r') as fh:
        fasta = SeqIO.parse(fh, 'fasta')
        return any(fasta)

def get_fasta_file_list (check_directory,  logger):
    if not os.path.isdir(check_directory):
        logger.info('directory %s does not exists', check_directory)
        return False
    filter_files = os.path.join(check_directory, '*.fasta')
    list_filtered_files =  glob.glob(filter_files)
    list_filtered_files.sort()
    if len (list_filtered_files) == 0 :
        logger.info('directory %s does not have any fasta file ', check_directory)
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
