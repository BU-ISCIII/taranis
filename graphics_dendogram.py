import pandas as pd
import os
from datetime import datetime
import logging
from logging.handlers import RotatingFileHandler

#import plotly.graph_objs as go
import plotly.io as pio

import plotly.figure_factory as ff


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
    #import plotly.plotly as py
    #import plotly.figure_factory as ff
    #import plotly.graph_objs as go
    #import plotly.io as pio
    #import plotly
    #fig = go.Figure()
    
    dendro = ff.create_dendrogram(clustering_matrix, labels=label_list)
    dendro['layout'].update({'width':800, 'height':500})
    #plotly.offline.iplot(dendro, filename='simple_dendrogram')
    #py.iplot(dendro, filename='simple_dendrogram')

    
    
    try:
        pio.write_image(dendro, file=out_file, format = 'png')
    except Exception as e :
        string_message = 'Unable to create ouput picture file'
        logging_errors(string_message, True, True)
        raise
    import pdb; pdb.set_trace()
    logger.debug('End the function create_graphic' )
    return True


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
    # create temporary folder to store the log
    if os.path.isdir(os.path.join(arguments.outputdir,LOGGING_FOLDER)) :
        print( 'Reusing the log folder for a previous execution \n')
    else:
        try:
            os.makedirs(os.path.join(arguments.outputdir,LOGGING_FOLDER))
        except :
            print ('Unable to create folder log \n')
            return 'ERROR'
    
    # open log file
    taranis_log = os.path.join(arguments.outputdir, LOGGING_FOLDER, LOGGING_NAME)
    logger = open_log (taranis_log)
    
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
    
    
    
    logger.info('Opening the input file %s', arguments.file_distance)
    try:
        pd_matrix = pd.read_csv(arguments.file_distance, sep='\t', header=0, index_col=0)
    except Exception as e:
        string_message = 'Unable to open the matrix distance file'
        logging_errors (string_message, False, True)
        logger.debug('End the function create_distance_matrix with error' )
        return 'ERROR'
    
    
    if arguments.exclude :
        indexes = pd_matrix.index.tolist()
        for be_excluded in arguments.exclude :
            #if be_excluded in indexes :
            try:
                index_to_remove = indexes.index(be_excluded)
                logger.info('Deleting index %s in panda matrix', be_excluded)
                pd_matrix = pd_matrix.drop(pd_matrix.index[index_to_remove])
            except ValueError as e:
                import pdb; pdb.set_trace()
                string_message = 'Index ' + be_excluded + ' does not found in file distance'
                logging_errors (string_message, False, True)
                return 'ERROR'
    import pdb; pdb.set_trace()
    try:
        distance_matrix = create_distance_matrix (pd_matrix, arguments.outputdir)
    except Exception as e:
        return 'ERROR'

    #logger.info('Converting distance matrix to cluster')
    #clustering_matrix = linkage(distance_matrix)
    
    logger.info('Prepare graphic')
    out_file = os.path.join(arguments.outputdir, 'dendogram.png')
    label_list =  list(distance_matrix.index)
    method = 'complete'
    metric = 'euclidean'
    try:
        create_dendogram_graphic (out_file, label_list, distance_matrix, method, metric)
    except:
        print('error')
    
    end_time = datetime.now()
    print('completed execution at :', end_time )
    return True