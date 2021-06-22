#!/usr/bin/env python3

import argparse
import sys
import io
import os
import logging
from logging.handlers import RotatingFileHandler
from datetime import datetime
from io import StringIO
import pandas as pd
import shutil
import csv 
import plotly.graph_objects as go
from utils.taranis_utils import *


# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · *  #  
# Filter and remove samples with missing percentage above certain specified threshold from allele calling comparison table for distance matrix calculation #
# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · *  #

def missing_filter_row(pd_matrix, missing_values_dict, sample_missing_threshold, logger):

    logger.info('Filtering samples with missing values in more than %s loci', sample_missing_threshold)

    num_cols = len(pd_matrix.columns)
    num_rows = len(pd_matrix.index) 

    missing_values_all_samples_dict = {} 
    sample_index_to_ignore = []
    

    ## Count unique tags to get missing values tag percentage for each sample or row
    for index, row in pd_matrix.iterrows(): 

        missing_values_per_sample_dict = dict(missing_values_dict)
        row_elements = list(row) 
        row_elements_unique = list(set(row_elements)) 

        for tag in row_elements_unique:

            for missing_tag in missing_values_per_sample_dict:
                if missing_tag in tag:
                    missing_values_per_sample_dict[missing_tag] += row_elements.count(tag)       

        missing_percent = (sum(missing_values_per_sample_dict.values())/num_cols) * 100

        if sample_missing_threshold < missing_percent : 

            sample_index_to_ignore.append(index)
            missing_values_all_samples_dict[index] = missing_values_per_sample_dict

    ## Remove samples with missing value percentage above specified missing threshold
    for sample_index in sample_index_to_ignore:
        
        pd_matrix = pd_matrix.drop(sample_index, axis = 0)
        
    return pd_matrix, missing_values_all_samples_dict, num_rows, num_cols


# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · *  #  
# Filter and remove loci with missing percentage above certain specified threshold from allele calling comparison table for distance matrix calculation #
# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · *  #

def missing_filter_col(pd_matrix, missing_values_dict, locus_missing_threshold, logger):

    logger.info('Filtering loci with missing values in more than %s samples', locus_missing_threshold)

    num_rows = len(pd_matrix.index)
    num_cols = len(pd_matrix.columns)

    pd_matrix_cols = list(pd_matrix) 
    index = 0

    missing_values_all_loci_dict = {}
    locus_index_to_ignore = []

    ## Count unique tags to get missing values tag percentage for each locus or column
    for col in pd_matrix_cols:
        
        missing_values_per_locus_dict = dict(missing_values_dict)
        col_list = []

        for row in range(0, num_rows):
        
            col_list.append(pd_matrix[col][row])  
        
        col_elements_count_dict = {}
        col_elements_unique = list(set(col_list))

        for tag in col_elements_unique:
            for missing_tag in missing_values_per_locus_dict:
                if missing_tag in tag:
                    missing_values_per_locus_dict[missing_tag] += col_list.count(tag)

        missing_percent = (sum(missing_values_per_locus_dict.values())/num_rows) * 100

        if locus_missing_threshold < missing_percent : 
        
            locus_index_to_ignore.append(pd_matrix_cols[index])
            missing_values_all_loci_dict[list(pd_matrix.columns)[index]] = missing_values_per_locus_dict
        
        index += 1

    ## Remove loci with missing value percentage above specified missing threshold
    for locus_index in locus_index_to_ignore:
        
        pd_matrix = pd_matrix.drop(locus_index, axis = 1)

    return pd_matrix, missing_values_all_loci_dict, num_rows, num_cols


# · * · * · * · * · * · * · *  #  
# Find hamming distance matrix #
# · * · * · * · * · * · * · *  #

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


# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * ·  #  
# Create samples distance matrix from filtered allele calling comparison table and get samples and loci filtering report #
# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * ·  #

def create_distance_matrix (alleles_matrix, outputdir, locus_missing_threshold, sample_missing_threshold, paralog_filter, lnf_filter, plot_filter, logger):

    ## Read the allele calling comparison table file
    try:
        pd_matrix = pd.read_csv(alleles_matrix, sep='\t', header=0, index_col=0)

    except Exception as e:
        logger.info('Unable to open the distance matrix file %s', result_file)

        print('------------- ERROR --------------')
        print('Unable to open the allele calling comparison table file')
        print('Check in the logging configuration file')
        print('------------------------------------------')
        return 'Error'

    ## Get list of tags to be considered missing values if required
    if paralog_filter != "false" or lnf_filter != "false" or plot_filter != "false":

        missing_values_dict = {}
        if lnf_filter == "true":
            missing_values_dict["LNF"] = 0

        if paralog_filter == "true":
            missing_values_dict["NIPH"] = 0 

        if plot_filter == "true":
            missing_values_dict["PLOT"] = 0


        ## Filter samples and loci with missing values percentage obove specified threshold if required
        if sample_missing_threshold < 100 or locus_missing_threshold < 100:

            filter_report_file = os.path.join(outputdir, 'matrix_distance_filter_report.tsv')

            if sample_missing_threshold < 100: ### ESTO MIRAR CÓMO LO VOY A INDICAR AL FINAL #####

                # Filter samples from allele calling comparison table file
                pd_matrix, missing_values_all_samples_dict, total_samples, total_loci = missing_filter_row(pd_matrix, missing_values_dict, sample_missing_threshold, logger)

                # Get samples filtering report 
                logger.info('Saving distance matrix filter information to file after samples filtering..')

                header_samples_filter = ["Sample", "LNF", "NIPH + NIPHEM", "PLOT", "Total missing values %"]

                #with open (filter_report_file, 'a') as filter_fh:

                if len(missing_values_all_samples_dict) > 0:

                    filter_fh.write(str(len(missing_values_all_samples_dict)) + "/" + str(total_samples) + " samples with missing values in more than " + str(sample_missing_threshold) + "% loci removed: " + '\n')
                    filter_fh.write('\n' + '\t'.join(header_samples_filter) + '\n')

                    for sample in sorted(missing_values_all_samples_dict):

                        missing_values = [str(missing_values_all_samples_dict[sample][missing_tag]) + "/" + str(total_loci) for missing_tag in sorted(missing_values_all_samples_dict[sample])]
                        filter_fh.write(sample + '\t' + '\t'.join(missing_values) + '\t' + str(sum(missing_values_all_samples_dict[sample].values())/total_loci*100) + '\n')
                    filter_fh.write('\n' + '\n' + '\n')

            
            if locus_missing_threshold < 100: 
                
                # Filter loci from allele calling comparison table file
                pd_matrix, missing_values_all_loci_dict, total_samples, total_loci = missing_filter_col(pd_matrix, missing_values_dict, locus_missing_threshold, logger)

                # Get samples filtering report
                logger.info('Saving distance matrix filter information to file after loci filtering..')

                header_loci_filter = ["Locus", "LNF", "NIPH + NIPHEM", "PLOT", "Total missing values %"]

                with open (filter_report_file, 'a') as filter_fh:
                    
                    #if len(missing_values_all_loci_dict) > 0:

                    filter_fh.write(str(len(missing_values_all_loci_dict)) + "/" + str(total_loci) + " loci with missing values in more than " + str(locus_missing_threshold) + "% of samples removed: " + '\n')

                    if len(missing_values_all_loci_dict)/total_loci > 0.1:
                    
                        filter_fh.write("WARNING! Less than 90% of cgMLST schema loci kept to create distance matrix." + '\n')             
                    
                    filter_fh.write('\n' + '\t'.join(header_loci_filter) + '\n')
                        
                    for locus in sorted(missing_values_all_loci_dict):

                        missing_values = [str(missing_values_all_loci_dict[locus][missing_tag]) + "/" + str(total_samples) for missing_tag in sorted(missing_values_all_loci_dict[locus])]
                        filter_fh.write(locus + '\t' + '\t'.join(missing_values) + '\t' + str(sum(missing_values_all_loci_dict[locus].values())/total_samples*100) + '\n')

            ## Save filtered allele calling results
            out_filtered_results_file = os.path.join(outputdir, 'filtered_result.tsv')
            pd_matrix.to_csv(out_filtered_results_file, sep = '\t')

        else:
            print("WARNING: Samples or loci filtering not specified. Cannot filter missing values")
            logger.info('WARNING: Samples or loci filtering not specified. Cannot filter missing values..')

    else:
        if sample_missing_threshold < 100 or locus_missing_threshold < 100:
            print("WARNING: Missing values to be filtered not specified. Cannot filter missing values")
            logger.info('WARNING: Missing values to be filtered not specified. Cannot filter missing values..')            

    ## Get and keep hamming distance matrix
    distance_matrix = hamming_distance (pd_matrix)
    out_distance_file = os.path.join(outputdir, 'matrix_distance.tsv')
    
    try:
        distance_matrix.to_csv(out_distance_file, sep = '\t')
    except Exception as e:
        logger.info('Unable to create the distance matrix file %s', out_file)

        print('------------- ERROR --------------')
        print('Unable to create the distance matrix file')
        print('Check in the logging configuration file')
        print('------------------------------------------')
        return 'Error'

    return True



# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · *  #  
# Get samples distance matrix from allele calling comparison table #
# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · *  #

def processing_distance_matrix (arguments) :
    '''
    Description:
        This is the main function for getting distance matrix.
        With the support of additional functions it will obtain the samples distance matrix from allele calling comparison table.
    Input:
        arguments   # Input arguments given on command line 
    Functions:
        
    Variables:  ?????

    Return: ?????
    '''

    start_time = datetime.now()
    print('Start the execution at :', start_time )
    # Open log file
    logger = open_log ('taranis_distance_matrix.log')

    #############################
    ## Create output directory ##
    #############################
    try:
        os.makedirs(arguments.outputdir)
    except:
        logger.info('Deleting the output directory for a previous execution without cleaning up')
        shutil.rmtree(arguments.outputdir)
        try:
            os.makedirs(arguments.outputdir)
            logger.info ('Output folder %s has been created again', arguments.outputdir)
        except:
            logger.info('Unable to create again the output directory %s', arguments.outputdir)
            print('Cannot create output directory on ', arguments.outputdir)
            exit(0)

    ################################
    ## Create the distance matrix ## 
    ################################
    try:
        logger.info('Creating distance matrix')
        print ('Creating distance matrix\n')
        create_distance_matrix(arguments.alleles_matrix, arguments.outputdir, int(arguments.locus_missing_threshold), int(arguments.sample_missing_threshold), str(arguments.paralog_filter).lower(), str(arguments.lnf_filter).lower(), str(arguments.plot_filter).lower(), logger)

    except:
        logger.info('There was an error when creating distance matrix')
        print('There was an error when creating distance matrix\n')
        shutil.rmtree(os.path.join(arguments.outputdir))
        exit(0)

    end_time = datetime.now()
    print('Completed execution at :', end_time )

    return True