#!/usr/bin/env python3
import argparse
import os
import shutil
import sys
import glob
from datetime import datetime
import statistics
import matplotlib.pyplot as plt
import numpy as np
#import logging
#from logging.handlers import RotatingFileHandler
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import Seq

#from Bio.Blast.Applications import NcbiblastnCommandline
from io import StringIO
#from Bio.Blast import NCBIXML
#from BCBio import GFF
from utils.taranis_utils import *

def check_arg(args=None):
    
    parser = argparse.ArgumentParser(prog = 'analyze_schema.py', description="This program will analyze the schema that is in schemadir parameter or it will compare 2 schemas ")
    #group = parser.add_mutually_exclusive_group()
    #group.add_argument ('-a', help = 'Interactive locus download.')
    #group.add_argument ('-b' , help = 'opcion b')
    parser.add_argument('-output_dir', help = 'Directory where the result files will be stored')
    subparser = parser.add_subparsers(help = 'analyze schema has 2 available options: (evaluate/compare) Evaluate 1 schema or compare 2 different schemas', dest = 'chosen_option')
    
    evaluate_parser = subparser.add_parser('evaluate', help = 'Evaluate the schema ')
    evaluate_parser.add_argument('-input_dir', help = 'Directory where are the schema files.')
    
    compare_parser = subparser.add_parser('compare', help = 'Compare 2 schema')
    compare_parser.add_argument('-scheme1', help = 'Directory where are the schema files for the schema 1')
    compare_parser.add_argument('-scheme2', help = 'Directory where are the schema files for the schema 2')
    
    return parser.parse_args()
'''
def open_log(log_name):
    working_dir = os.getcwd()
    log_name=os.path.join(working_dir, log_name)
    #def create_log ():
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    #create the file handler
    handler = logging.handlers.RotatingFileHandler(log_name, maxBytes=200000, backupCount=5)
    handler.setLevel(logging.DEBUG)

    #create a Logging format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    #add the handlers to the logger
    logger.addHandler(handler)

    return logger

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
    return False
'''    
def analyze_schema (schema_files,  logger) :
    not_cds_dict = {}
    schema_sequence_dict ={}
    
    schema_info_dict = {}
    reverse_alleles_dict = {}
    protein_dict = {}
    allele_duplicated = {}
    for schema_file in schema_files :
        schema_fasta_dict ={}
        tmp_gene_name = os.path.basename(schema_file).split('.')
        gene_name = tmp_gene_name[0]
        print('analyzing : ' ,gene_name)
        protein_dict[gene_name] = {}
        schema_info_dict[gene_name] = {}
        schema_sequence_dict[gene_name] = {}
        for contig in SeqIO.parse(schema_file, "fasta", generic_dna):
            #check if allele id contain string characters . If yes get only the nummber
            if '_' in contig.id:
                tmp_id = contig.id.split('_')
                contig_id = int(tmp_id[-1])
            else:
                contig_id = int(contig.id)
            schema_fasta_dict[contig_id] = str(contig.seq.upper())
        
        for allele_id in sorted(schema_fasta_dict) :
            query_direction = check_sequence_order(schema_fasta_dict[allele_id], logger)
            sequence = Seq.Seq(schema_fasta_dict[allele_id])
            sequence_str = str(sequence)
            if query_direction == 'reverse' :
                if not gene_name in reverse_alleles_dict :
                    reverse_alleles_dict[gene_name] = {}
                if not allele_id in reverse_alleles_dict[gene_name] :
                    reverse_alleles_dict[gene_name][allele_id] = schema_fasta_dict[allele_id]
                
                sequence = sequence.reverse_complement()
            
            try:
                protein = str(sequence.translate(cds=True))
                protein_dict[gene_name][allele_id] = protein
                coding_cds = 'Yes'
                error_description = 'No error'
            except Exception as error:
                logger.error('Not CDS for gene %s in the allele %s ', gene_name, contig.id)
                if not gene_name in not_cds_dict :
                    not_cds_dict[gene_name] = {}
                coding_cds = 'No'
                error_description = str(error)
                not_cds_dict[gene_name][allele_id] = [error_description , schema_fasta_dict[allele_id]]
                protein_dict[gene_name][allele_id] = 'NOT CDS'
                #print('allele ', allele_id, 'error : ', error, ' seq', str(sequence))
                #print( 'Next')
            schema_info_dict[gene_name][allele_id] = [str(len(sequence_str)), coding_cds, error_description, query_direction]
            if not sequence_str in schema_sequence_dict[gene_name] :
                schema_sequence_dict [gene_name][sequence_str]= [allele_id]
            else:
                schema_sequence_dict [gene_name][sequence_str].append(str(allele_id))
            for allele_found in schema_sequence_dict [gene_name] :
                if len(schema_sequence_dict [gene_name][allele_found]) > 1 :
                    if not gene_name in allele_duplicated [gene_name] :
                        allele_duplicated[gene_name] = []
                    allele_duplicated [gene_name].append([gene_name][allele_found])
            
    return not_cds_dict , reverse_alleles_dict, protein_dict, schema_info_dict , allele_duplicated

def create_bar_graphic (x_data, y_data, x_label, y_label, title , rotation, file_name) :
    
    index = np.arange(len(x_data))
    plt.bar(index, y_data)
    plt.xlabel(x_label, fontsize=5)
    plt.ylabel(y_label, fontsize=5)
    plt.xticks(index, x_data, fontsize=10, rotation=rotation)
    plt.title(title)
    #plt.show()
    plt.savefig(file_name)
    plt.close()
    return True

def summary_schema( schema_info, output_dir , logger) :
    summary_info = {}
    variability_length = {}
    coding_cds = {}
    error_type = {}
    gene_length = {}
    direction = {}
    
    # join all individual information to one item per gene
    for gene in sorted(schema_info) :
        g_length = []
        coding_cds[gene] = {}
        error_type[gene] = {}
        direction[gene] = {}
        logger.debug('dumping g_length for gene  %s ' ,gene)
        for allele in schema_info[gene] :
            values = schema_info[gene][allele]
            g_length.append(int(values[0]))
            #g_coding.append(values[1])
            if not values[1] in coding_cds[gene] :
                coding_cds[gene][values[1]] = 0
            coding_cds[gene][values[1]] += 1
            
            if not values[2] in error_type[gene] :
               error_type[gene][values[2]] = 0
            error_type[gene][values[2]] += 1
            
            if not values[3] in direction [gene]:
                direction[gene][values[3]] = 0
            direction[gene][values[3]] += 1
        try:  
            mode_length=statistics.mode(g_length)
        except:
            import pdb; pdb.set_trace()
        min_length = min(g_length)
        max_length = max(g_length)
        gene_length[gene] = mode_length
        variability_length[gene]=format(max((mode_length-min_length), (max_length-mode_length))/mode_length, '.2f')
    
    # combine the length information to create the graphic to show the number of the lenght gene and the number of times that gene has the same length in the schema     
    summary_length = {}    
    set_of_length = []
    number_of_set_length = []
    
    for value  in gene_length.values() :
        if not value in summary_length :
            summary_length[value] = 0
        summary_length[value] += 1   
    for index, value in sorted(summary_length.items()) :
        set_of_length.append(index)
        number_of_set_length.append(value)
       
    length_graphic_file = os.path.join(output_dir, 'graphic_length_relation.png')
    rotation = 30
    create_bar_graphic (set_of_length, number_of_set_length, 'length of gene', 'Number of gene with the same length', 'Length of the sequence for each gene defined in the schema ' , rotation,  length_graphic_file) 
    
    variation_lenght = {}
    index_variation = []
    value_varation = []
    for gene, v_length in variability_length.items() :
        if not v_length in variation_lenght :
            variation_lenght[v_length] = 0
        variation_lenght [v_length] += 1
    for index, value in sorted(variation_lenght.items()):
        index_variation.append(index)
        value_varation.append(value)
    
    varation_length_graphic_file = os.path.join(output_dir, 'graphic_varation_length.png')
    rotation = 30
    create_bar_graphic (index_variation, value_varation, 'length variability of gene', 'Numbers of gene variability', 'Variability length of the sequence for each gene defined in the schema ' , rotation,  varation_length_graphic_file) 
    
    
    # combine the number of times that an allele is not protein coding
    summary_coding_cds = {}
    count_conting_cds = {}
    percents = []
    percent_value = []
    for gene in coding_cds :
        if 'Yes' in coding_cds[gene] :
            allele_coding_cds = coding_cds[gene]['Yes']
        else:
            allele_coding_cds = 0
        if 'No' in coding_cds[gene] :
            allele_no_coding_cds = coding_cds[gene]['No']
        else:
            allele_no_coding_cds = 0
        percent_not_coding = format(allele_no_coding_cds/(allele_no_coding_cds + allele_coding_cds), '.2f')
        summary_coding_cds[gene] = percent_not_coding
        if not percent_not_coding in count_conting_cds :
            count_conting_cds[percent_not_coding] = 0
        count_conting_cds[percent_not_coding] += 1
    for index, value in sorted(count_conting_cds.items()) :
        percents.append(index)
        percent_value.append(value)
    # create the plot file for the (cdc/non cds) percent relation 
    percent_graphic_file = os.path.join(output_dir, 'graphic_percent_relation.png')
    rotation = 30
    create_bar_graphic (percents, percent_value, 'Percent of non coding CDS', 'Number no coding CDS', 'Percent of the alleles in the schema that are not coding CDS ' , rotation, percent_graphic_file) 
    
    
    # combine the number of times that the error codo arise when trying to conver to cds
    summary_error_type = {}
    error_name = []
    error_value = []
    for gene, errors in error_type.items() :
        for error_code , value_error in errors.items() :
            if error_code != 'No error' :
                
                if 'start codon' in error_code :
                    error_code = 'not start codon'
                elif 'Extra in frame stop' in error_code :
                    error_code = 'extra stop codon'
                elif 'not a stop codon' in error_code :
                    error_code = 'not stop codon'
                else:
                    pass
                if not error_code in summary_error_type :
                    summary_error_type[error_code] = 0
                summary_error_type[error_code] += value_error
    for error , value in summary_error_type.items():
        error_name.append(error)
        error_value.append(value)
    
    #create the plot file for error types when trying to convert to cds
    error_type_graphic_file = os.path.join(output_dir, 'graphic_error_type_cds.png')
    rotation = 0
    create_bar_graphic (error_name, error_value, 'Error type when converting to CDS', 'Number of errors', 'Type of errors that are generated when trying to convert to CDS ' , rotation , error_type_graphic_file) 
    
    
    
    return variability_length, gene_length, coding_cds , error_type, direction


def evaluate_schema (inputdir, outputdir, logger) :

    header_allele_no_cds = ['Gene name', 'Allele id' , 'error description', 'sequence']
    header_reverse_alleles = ['Gene name', 'allele id' , 'sequence']
    header_proteins = ['Gene name', 'allele id' , 'protein']
    header_alleles_duplicated = ['Gene name', 'Duplicated alleles id' ]
    header_schema_info = ['Gene name', 'Allele id' , 'length', 'Coding(Yes/No)' , 'Error description','direction']
    schema_files = get_fasta_file_list(inputdir, logger)
    allele_no_cds , reverse_alleles, proteins , schema_info , allele_duplicated = analyze_schema (schema_files,  logger)
    
    logger.info('Saving alleles not coding to protein to file..')
    
    for schema in sorted (allele_no_cds) :
        allele_no_cds_file =  os.path.join(outputdir, str(schema + '_allele_no_cds.tsv'))
        with open (allele_no_cds_file , 'w') as allele_no_cds_fh :
            allele_no_cds_fh.write('\t'.join(header_allele_no_cds) + '\n')
            for allele in sorted (allele_no_cds[schema], key=int):
                allele_no_cds_fh.write(schema + '\t' + str(allele) + '\t' + '\t'.join(allele_no_cds[schema][allele]) + '\n')
    
    logger.info('Saving dulicate alleles to file..')
    for gene in sorted (allele_duplicated) :
        allele_duplicated_file =  os.path.join(outputdir, str(gene + '_alleles_duplicated.tsv'))
        with open (allele_duplicated_file , 'w') as allele_duplicated_fh :
            allele_duplicated_fh.write('\t'.join(header_allele_duplicated) + '\n')
            for duplication in (allele_duplicated[gene]):
                allele_duplicated_fh.write(gene + '\t'  + '\t'.join(allele_duplicated[gene][duplication]) + '\n')
    
    logger.info('Saving schema info  to file..')
    schema_info_file =  os.path.join(outputdir,  'schema_information.tsv')
    with open (schema_info_file , 'w') as schema_info_fh :
        schema_info_fh.write('\t'.join(header_schema_info) + '\n')
        for gene in sorted (schema_info) :
            for allele in (schema_info[gene]):
                schema_info_fh.write(gene + '\t' + str(allele) + '\t' + '\t'.join(schema_info[gene][allele]) + '\n')


    logger.info('Saving alleles not coding to protein to file..')
    for schema in sorted (reverse_alleles) :
        reverse_alleles_file =  os.path.join(outputdir, str(schema + '_reverse_alleles.tsv'))
        with open (reverse_alleles_file , 'w') as reverse_alleles_fh :
            reverse_alleles_fh.write('\t'.join(header_reverse_alleles) + '\n')
            for allele in sorted (reverse_alleles[schema], key=int):
                reverse_alleles_fh.write(schema + '\t' + str(allele) + '\t' + reverse_alleles[schema][allele] + '\n')
                
    logger.info('Saving proteins to file..')
    os.makedirs(os.path.join(outputdir, 'proteins'))
    for schema in sorted (proteins) :
        proteins_file =  os.path.join(outputdir, 'proteins', str(schema + '_proteins.tsv'))
        with open (proteins_file , 'w') as proteins_fh :
            proteins_fh.write('\t'.join(header_proteins) + '\n')
            for allele in sorted (proteins[schema], key=int):
                proteins_fh.write(schema + '\t' + str(allele) + '\t' + proteins[schema][allele] + '\n')

    variability_length, gene_length, coding_cds , error_type, direction = summary_schema( schema_info, outputdir, logger) 
    
    
    return True

if __name__ == '__main__' :
    version = 'analyze_schema  version 0.0.1'
    if len(sys.argv) == 1 :
        print( 'Mandatory parameters are missing to execute the program. \n ' ,'Usage: "analyze_schema -help " for more information \n')
        exit (0)
    if sys.argv[1] == '-v' or sys.argv[1] == '--version':
        print( version, '\n')
        exit (0)
    arguments = check_arg(sys.argv[1:])
    start_time = datetime.now()
    print('Start the execution at :', start_time )
    # open log file
    logger = open_log ('analyze_schema.log')
    

    try:
        os.makedirs(arguments.output_dir)
    except:
        print('The output directory is not empty')
        choice_value = input('Enter yes to delete directory. Any other character to exit the program >>  ')
        if choice_value == 'yes' or choice_value == 'YES' :
            logger.info('Deleting the result  directory for a previous execution without cleaning up')
            shutil.rmtree(arguments.output_dir)
            try:
                os.makedirs(arguments.output_dir)
                logger.info ( 'Result folder %s  has been created again', arguments.output_dir)
            except:
                logger.info('Unable to create again the result directory %s', arguments.output_dir)
                print('Cannot create result directory on ', arguments.output_dir)
                exit(0)
        else:
            print('Aborting the execution')
            exit(0)
    if arguments.chosen_option =='evaluate' :
        evaluate_schema (arguments.input_dir, arguments.output_dir, logger)
    else:
        pass # compare 2 schema
    end_time = datetime.now()
    print('completed execution at :', end_time )

