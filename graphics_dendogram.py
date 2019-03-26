import pandas as pd
import os
import logging
from logging.handlers import RotatingFileHandler

from utils.taranis_utils import *
from taranis_configuration import *


def create_dendogram_graphic (out_file, label_list, clustering_matrix, method, metric):
    '''
    Description:
        This function create the dendogram graphic.
    Input:
        out_file
        label_list
        clustering_matrix
        method
        metric
        pd_matrix   # panda dataframe to calculate the hamming distance 
    Variables:
        H   # Distance matrix
        U   # Binary matrix to used to calculate the distance
        unique_values   # unique values to iterate for the distance algorithm
    Return:
        True if graphic creation file is sucessful, exception if error when creating the file
    '''
    logger = logging.getLogger(__name__)
    logger.debug('Starting the function create_graphic' )
    
    
    logger.debug('End the function create_graphic' )
    return True

def hamming_distance (pd_matrix):
    '''
    Description:
        This function will perform the distance matrix for a given panda dataframe.
    Input:
        pd_matrix   # panda dataframe to calculate the hamming distance 
    Variables:
        H   # Distance matrix
        U   # Binary matrix to used to calculate the distance
        unique_values   # unique values to iterate for the distance algorithm
    Return:
        H   # matrix containg the hamming matrix
    '''
    logger = logging.getLogger(__name__)
    logger.debug('Starting the function hamming_distance' )
    unique_values = pd.unique(pd_matrix[list(pd_matrix.keys())].values.ravel('K'))
    # Create binary matrix ('1' or '0' ) matching the input matrix vs the unique_values[0]
    # astype(int) is used to transform the boolean matrix into integer
    U = pd_matrix.eq(unique_values[0]).astype(int)
    # multiply the matrix with the transpose
    H = U.dot(U.T)
    logger.info('Start distance calculation')
    
    # Repeat for each unique value
    for unique_val in range(1,len(unique_values)):
        U = pd_matrix.eq(unique_values[unique_val]).astype(int)
        # Add the value of the binary matrix with the previous stored values
        H = H.add(U.dot(U.T))
    
    logger.debug('End the function hamming_distance' )
    return len(pd_matrix.columns) - H  #,  H #H.div(len(pd_matrix.columns))

def create_dendogram_from_distance (arguments):
    '''
    Description:
        This is the main function for create_dendogram_from_distance.
    Input:
        arguments   # Input arguments given on command line 
    IMPORT:
        LOGGING_NAME
        LOGGING_FOLDER
    Functions:
        
         # located at utils.taranis_utils
        open_log # located at utils.taranis_utils
          # located at this file
            # located at this file
        
    Variables:
        clustering_matrix   # matrix contains the cluster 
        distance_matrix     # matrix having the distance calculation
        
        Return:
        
    '''
    start_time = datetime.now()
    print('Start the execution at :', start_time )
    # check if input files exists
    if not os.path.isfile(arguments.file_distance):
        print ('File ',arguments.file_distance, ' does not exist')
        return 'ERROR'
    if not os.path.isdir(arguments.outputdir) :
        try :
            os.makedirs(arguments.outputdir)
        except :
            string_message = 'Unable to create ouput directory'
            logging_errors(string_message, True, True)
            return 'ERROR'
    
    # open log file
    taranis_log = os.path.join(arguments.outputdir, LOGGING_FOLDER, LOGGING_NAME)
    logger = open_log (taranis_log)
    
    pd_matrix = pd.read_csv(arguments.file_distance, sep='\t', header=0, index_col=0)
    distance_matrix = hamming_distance (pd_matrix)
    logger.info('Saving the matrix distance ')
    out_file = os.path.join(arguments.outputdir, 'matrix_distance.tsv')
    matrix_for_distance.to_csv(out_file, sep = '\t')

    #logger.info('Converting distance matrix to cluster')
    #clustering_matrix = linkage(distance_matrix)
    
    logger.info('Prepare graphic')
    out_file = os.path.join(arguments.outputdir, 'dendogram')
    label_list =  list(distance_matrix.index)
    
    try:
        create_dendogram_graphic (out_file, label_list, distance_matrix, method, metric)
    except:
        print('error')
    
    end_time = datetime.now()
    print('completed execution at :', end_time )
    return True