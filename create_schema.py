#!/usr/bin/env python3
from utils.taranis_utils import *

def prueba () :
    print ('hola')
    return


def processing_create_schema(arguments):
    '''
    Description:
    
    
    Input:
        
    Variables:
    
    
    Return:
    
    '''
    logger = open_log('create_schema')
    if logger != 'Error' :
        logger.info('test')
        logger.error('error_test ')
    else :
        print ('Exiting the create schema utility because the log file')
        print ('could not be created')
        return 'Error'
    
    return True