#!/usr/bin/env python3
from utils.taranis_utils import *



def processing_create_schema(arguments):
    '''
    Description:
    
    
    Input:
        
    Variables:
    
    
    Return:
    
    '''
    xls_file = arguments.xlsfile
    output_dir = arguments.outputdir
    logger = open_log('create_schema')
    if logger != 'Error' :
        gene_list = read_xls_file(xls_file, logger)
        if  'Error' not in gene_list :
            for gene, protein in gene_list :
                #curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore;id=NP_214515.1;retmode=text;rettype=fasta"
                pass
        else:
            print('There was an error when accessing the excel file ')
            print(gene_list )
            return gene_list
            #return 'Error when reading the excel file'
        logger.info('test')
        logger.error('error_test ')
    else :
        print ('Exiting the create schema utility because the log file')
        print ('could not be created')
        return 'Error'
    
    return True
