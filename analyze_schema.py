#!/usr/bin/env python3
import argparse
import os
import shutil
import sys
import glob
from datetime import datetime
import statistics
#import matplotlib.pyplot as plt
import plotly.graph_objs as go
import plotly.io as pio
#import numpy as np
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
from progressbar import ProgressBar
from utils.taranis_utils import *
  
def extract_info_schema (schema_files,  alt_codon_start, logger) :
    not_cds_dict = {}
    schema_sequence_dict ={}
    
    schema_info_dict = {}
    reverse_alleles_dict = {}
    protein_dict = {}
    allele_duplicated = {}
    pbar = ProgressBar ()
    for schema_file in pbar (schema_files) :
        schema_fasta_dict ={}
        tmp_gene_name = os.path.basename(schema_file).split('.')
        gene_name = tmp_gene_name[0]
        #print('analyzing : ' ,gene_name)
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
            if alt_codon_start == True and sequence.startswith('GTG') :
                alt_table =2
            else:
                alt_table =1
            try:
                protein = str(sequence.translate(cds=True, table =alt_table))
                #protein = str(sequence.translate(cds=True))
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

def create_bar_graphic (x_data, y_data, x_label, x_prefix ,y_label, title , rotation, file_name) :
    '''
    index = np.arange(len(x_data))
    plt.bar(index, y_data)
    plt.xlabel(x_label, fontsize=5)
    plt.ylabel(y_label, fontsize=5)
    plt.xticks(index, x_data, fontsize= 7, rotation=rotation)

    plt.title(title)
    #plt.show()
    plt.savefig(file_name)
    plt.close()
    '''
    
    
    trace0 = go.Bar(
                #x=['Product A', 'Product B', 'Product C'],
                #y=[20, 14, 23],
                x = x_data,
                y = y_data,
                text = y_data,
                
                #text=['27% market share', '24% market share', '19% market share'],
                textposition = 'auto',
                marker=dict( color='rgb(158,202,225)',
                    line=dict(
                    color='rgb(8,48,107)',
                    width=1.5, )
                ),
                opacity=0.6
                )
    
    data = [trace0]
    #import pdb; pdb.set_trace()
    layout = go.Layout( title=title,
                    xaxis = dict(title = x_label,
                    tickformat = '%' +x_prefix),
                    yaxis = dict(title = y_label),
                    )
    fig = go.Figure(data=data, layout=layout)
    pio.write_image(fig, file_name)
    return True

def find_proteins_in_gene (raw_proteins_per_genes, logger) :
    proteins_sequence_per_gene ={}
    proteins_percent_per_gene ={}
    logger.info('Start handling the raw_proteins to get the unique coding proteins')
    for gene in raw_proteins_per_genes :
        proteins = []

        #num_alleles = len (proteins_per_genes[gene])
        for allele, value in sorted(raw_proteins_per_genes[gene].items()) :
            if value != 'NOT CDS' :
                proteins.append(value)
        proteins_sequence_per_gene[gene] = list(set(proteins))
        if len(proteins) == 0 :
            proteins_percent_per_gene[gene] = '0'
        else:
            proteins_percent_per_gene[gene] = format(len(list(set(proteins))) / len(proteins) , '.2f')
        
    logger.info('Complete the protein handling')
    return proteins_sequence_per_gene, proteins_percent_per_gene


def summary_schema_info ( schema_info,  output_dir , logger) :
    logger.info('Start processing the information in schema info')
    header_variability_length = ['Gene name', 'Length variability']
    header_gene_length = ['Gene name', 'Length']
    header_percent_allele_not_cds =['Gene name', 'Allele Percentage that is not coding CDS']
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

        mode_length=statistics.mode(g_length)
        min_length = min(g_length)
        max_length = max(g_length)
        gene_length[gene] = mode_length
        variability_length[gene]=format(max((mode_length-min_length), (max_length-mode_length))/mode_length, '.2f')
    
    logger.info('Create the summary folder')
    os.makedirs(os.path.join(output_dir, 'summary'))
    
    logger.info('Dumping the variability length from the schema to file')
    variability_length_file =  os.path.join(output_dir, 'summary' , 'variability_length.tsv')
    save_simple_dict_to_file (variability_length,  header_variability_length, variability_length_file, logger)
    '''
    with open (variability_length_file , 'w') as variability_length_fh :
        variability_length_fh.write('\t'.join(header_variability_length) + '\n')
        for gene, value in sorted (variability_length.items()) :
            variability_length_fh.write(gene + '\t' +  value + '\n')
    '''
    logger.info('Dumping completed')
    
    logger.info('Dumping the gene length from the schema to file')
    gene_length_file = os.path.join(output_dir, 'summary' , 'gene_length.tsv')
    save_simple_dict_to_file (gene_length,  header_gene_length, gene_length_file, logger) 

    logger.info('Processing the picture for gene length')
    # Length of the gene will be clustered in 10 groups to be presented in the graphic bar
    x_axis = [150, 250, 500, 1000, 1500, 2000, 2500, 3000, 4000 , 5000]
    gene_length_values = 10 *[0]
    #summary_length = {}    
    #set_of_length = []
    #number_of_set_length = []
    
    for value  in gene_length.values() :
        if value > 5000 :
            # if gene length is bigger than 5000 it will be assigned to 5000
            gene_length_values[len(x_axis)-1] += 1
        else:
            for index in range(len(x_axis)) :
                if value <= x_axis[index] :
                    gene_length_values[index] += 1
                    break
      

    x_axis_label = ['<= {0}'.format(element) for element in x_axis]

    length_graphic_file = os.path.join(output_dir, 'graphic_gene_length.png')
    rotation = 30
    x_prefix = ''
    create_bar_graphic (x_axis_label, gene_length_values, 'Gene length', x_prefix ,'Number of gene with the same length', 'Sequence length for genes defined in the schema ' , rotation,  length_graphic_file) 

    #create_bar_graphic (set_of_length, number_of_set_length, 'length of gene', 'Number of gene with the same length', 'Length of the sequence for each gene defined in the schema ' , rotation,  length_graphic_file) 
    
    logger.info('Processing the picture for variablity length')
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
    
    x_axis_label = ['{0}%'.format(int(float(element)*100)) for element in index_variation]
    varation_length_graphic_file = os.path.join(output_dir, 'graphic_varation_length.png')
    rotation = 30
    x_prefix =''
    create_bar_graphic (x_axis_label, value_varation, 'length variability of gene', x_prefix,  'Numbers of gene variability', 'Variability length of the sequence for each gene defined in the schema ' , rotation,  varation_length_graphic_file) 
    logger.info('Complete picture for variability length')
    
    # combine the number of times that an allele is not protein coding
    summary_coding_cds = {}
    #count_conting_cds = {}
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
        

    logger.info('Dumping the allele percentage that are not codings CDS to file')
    percent_allele_not_coding_file =  os.path.join(output_dir, 'summary' , 'percent_allele_not_coding.tsv')
    save_simple_dict_to_file (summary_coding_cds,  header_percent_allele_not_cds, percent_allele_not_coding_file, logger)
    
  
    
    # create the plot file for the (cdc/non cds) percent relation 
    percent_coding_one_decimal = []
    for per_values in  summary_coding_cds.values() :
        percent_coding_one_decimal.append(str(round(float(per_values), 1)))
    
    percent_number = []
    percent_list = sorted(list(set(percent_coding_one_decimal)))
    for item in percent_list :
        percent_number.append(percent_coding_one_decimal.count(item))
    
    x_axis_label = ['{0}%'.format(int(float(element)*100)) for element in percent_list]
    
    percent_not_contig_graphic_file = os.path.join(output_dir, 'graphic_allele_percent_not_coding.png')
    rotation = 30
    x_prefix = ''
    create_bar_graphic (x_axis_label, percent_number, 'Percent of non coding CDS', x_prefix, 'Number of genes ', 'Alleles that are not coding CDS ( in % ) ' , rotation, percent_not_contig_graphic_file) 
    
    
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
    x_prefix = ''
    create_bar_graphic (error_name, error_value, 'Error type when converting to CDS', x_prefix,  'Number of errors', 'Type of errors that are generated when trying to convert to CDS ' , rotation , error_type_graphic_file) 
    
    logger.info('Schema info has been completed processed ')
    
    
    
    return True


def save_simple_dict_list_to_files (dict_to_save,  heading_text, folder_name ,file_name, logger) :
    logger.info('Saving file %s', file_name)
    for gene , value_list in sorted(dict_to_save.items()):
        f_name = os.path.join(folder_name, str(gene + file_name))
        with open (f_name , 'w') as f_name_fh :
            f_name_fh.write('\t'.join(heading_text) + '\n')
            for item in value_list :
                f_name_fh.write(gene + '\t' + item + '\n')
    logger.info('Saved file  %s', file_name)
    return True

def save_simple_dict_to_file (dict_to_save,  heading_text, file_name, logger) :
    logger.info('Saving file %s', file_name)
    with open (file_name , 'w') as file_name_fh :
        file_name_fh.write('\t'.join(heading_text) + '\n')
        for gene , value in sorted (dict_to_save.items()) :
            file_name_fh.write(gene + '\t' + str(value) + '\n')
    logger.info('Saved file  %s', file_name)
    return True


def summary_proteins (raw_proteins_per_genes, output_dir, logger) :
    logger.info('Start handling protein from the raw information')
    heading_summary_proteins_sequence = ['Gene Name', 'Protein sequence']
    heading_summary_proteins_percent = ['Gene Name', 'Percent of different proteins in the gene']
    proteins_sequence_per_gene, proteins_percent_per_gene = find_proteins_in_gene (raw_proteins_per_genes, logger)
    # Save proteins sequences proteins to file
    os.makedirs(os.path.join(output_dir, 'summary', 'proteins'))
    folder_summary_proteins = os.path.join(output_dir, 'summary', 'proteins')
    proteins_sequence_file = '_summary_protein_sequence.tsv'
    save_simple_dict_list_to_files (proteins_sequence_per_gene, heading_summary_proteins_sequence, folder_summary_proteins, proteins_sequence_file, logger)
    # Save proteins percent to file
    proteins_percent_file = os.path.join(output_dir, 'summary' , 'proteins_percent.tsv')
    save_simple_dict_to_file (proteins_percent_per_gene, heading_summary_proteins_percent, proteins_percent_file ,logger)
    
    # create the diagram to display the percent protoins for each gene
    # round number to 1 decimal to show the graphic
    all_percent = []
    pencent_values = proteins_percent_per_gene.values()
    for percent_value in pencent_values :
        all_percent.append(str(round(float(percent_value), 1)))
    
    #all_percent =  list(proteins_percent_per_gene.values() )
    percent_list = sorted(list(set(all_percent)))
    percent_number = []
    for item in percent_list :
        percent_number.append(all_percent.count(item))
    x_axis_label = ['{0}%'.format(int(float(element)*100)) for element in percent_list]
    protein_percent_graphic_file = os.path.join(output_dir, 'graphic_protein_percent.png')
    rotation = 30
    x_prefix =''
    create_bar_graphic (x_axis_label, percent_number, 'Percent of proteins ', x_prefix ,
                        'Number of genes', 'Percent of Alleles that coding for the same protein (in %)'
                        , rotation, protein_percent_graphic_file)
    
    return True

def evaluate_schema (inputdir, outputdir, alt_codon_start, logger) :

    header_allele_no_cds = ['Gene name', 'Allele id' , 'error description', 'sequence']
    header_reverse_alleles = ['Gene name', 'allele id' , 'sequence']
    header_proteins = ['Gene name', 'allele id' , 'protein']
    header_alleles_duplicated = ['Gene name', 'Duplicated alleles id' ]
    header_schema_info = ['Gene name', 'Allele id' , 'length', 'Coding(Yes/No)' , 'Error description','direction']
    schema_files = get_fasta_file_list(inputdir, logger)
    logger.info('Extract the raw information for each gene in the schema')
    allele_no_cds , reverse_alleles, raw_proteins_per_genes , schema_info , allele_duplicated = extract_info_schema (schema_files, alt_codon_start, logger)
    
    print('saving data to ', outputdir )
    logger.info('Start dumping the raw information to files')
    logger.info('Saving alleles not coding to protein to file..')
    os.makedirs(os.path.join(outputdir, 'raw_info'))
    os.makedirs(os.path.join(outputdir, 'raw_info', 'allele_not_cds'))
    for schema in sorted (allele_no_cds) :
        allele_no_cds_file =  os.path.join(outputdir, 'raw_info' , 'allele_not_cds' ,str(schema + '_allele_no_cds.tsv'))
        with open (allele_no_cds_file , 'w') as allele_no_cds_fh :
            allele_no_cds_fh.write('\t'.join(header_allele_no_cds) + '\n')
            for allele in sorted (allele_no_cds[schema], key=int):
                allele_no_cds_fh.write(schema + '\t' + str(allele) + '\t' + '\t'.join(allele_no_cds[schema][allele]) + '\n')
    
    logger.info('Saving dulicate alleles to file..')
    os.makedirs(os.path.join(outputdir, 'raw_info', 'duplicated_alleles'))
    for gene in sorted (allele_duplicated) :
        allele_duplicated_file =  os.path.join(outputdir, 'raw_info' , 'duplicated_alleles' , str(gene + '_alleles_duplicated.tsv'))
        with open (allele_duplicated_file , 'w') as allele_duplicated_fh :
            allele_duplicated_fh.write('\t'.join(header_allele_duplicated) + '\n')
            for duplication in (allele_duplicated[gene]):
                allele_duplicated_fh.write(gene + '\t'  + '\t'.join(allele_duplicated[gene][duplication]) + '\n')
    
    logger.info('Saving schema info  to file..')
    schema_info_file =  os.path.join(outputdir, 'raw_info', 'schema_information.tsv')
    with open (schema_info_file , 'w') as schema_info_fh :
        schema_info_fh.write('\t'.join(header_schema_info) + '\n')
        for gene in sorted (schema_info) :
            for allele in (schema_info[gene]):
                schema_info_fh.write(gene + '\t' + str(allele) + '\t' + '\t'.join(schema_info[gene][allele]) + '\n')


    logger.info('Saving alleles not coding to protein to file..')
    os.makedirs(os.path.join(outputdir, 'raw_info', 'raw_reverse_alleles'))
    for schema in sorted (reverse_alleles) :
        reverse_alleles_file =  os.path.join(outputdir, 'raw_info', 'raw_reverse_alleles', str(schema + '_reverse_alleles.tsv'))
        with open (reverse_alleles_file , 'w') as reverse_alleles_fh :
            reverse_alleles_fh.write('\t'.join(header_reverse_alleles) + '\n')
            for allele in sorted (reverse_alleles[schema], key=int):
                reverse_alleles_fh.write(schema + '\t' + str(allele) + '\t' + reverse_alleles[schema][allele] + '\n')
                
    logger.info('Saving proteins to file..')
    os.makedirs(os.path.join(outputdir, 'raw_info', 'raw_proteins'))
    for schema in sorted (raw_proteins_per_genes) :
        proteins_file =  os.path.join(outputdir, 'raw_info', 'raw_proteins', str(schema + '_proteins.tsv'))
        with open (proteins_file , 'w') as proteins_fh :
            proteins_fh.write('\t'.join(header_proteins) + '\n')
            for allele in sorted (raw_proteins_per_genes[schema], key=int):
                proteins_fh.write(schema + '\t' + str(allele) + '\t' + raw_proteins_per_genes[schema][allele] + '\n')
    
    logger.info('Completed dumped raw information to files')
    
    logger.info('Analyze the raw proteins to remove the non CDS and duplicated proteins for each gene')
    #proteins_per_gene = find_proteins_in_gene (raw_proteins_per_genes, logger)
    
    logger.info('Dumping proteins to file ')
    
    
    summary_schema_info( schema_info, outputdir, logger) 
    summary_proteins (raw_proteins_per_genes, outputdir, logger) 
    
    return True

'''
if __name__ == '__main__' :
    version = 'analyze_schema  version 0.0.1'
    if len(sys.argv) == 1 :
        print( 'Mandatory parameters are missing to execute the program. \n ' ,'Usage: "analyze_schema -help " for more information \n')
        exit (0)
    if sys.argv[1] == '-v' or sys.argv[1] == '--version':
        print( version, '\n')
        exit (0)
    arguments = check_arg(sys.argv[1:])
'''
def processing_evaluate_schema (arguments) :
    start_time = datetime.now()
    print('Start the execution at :', start_time )
    # open log file
    logger = open_log ('analyze_schema.log')
    

    try:
        os.makedirs(arguments.outputdir)
    except:
        print('The output directory is not empty')
        choice_value = input('Enter yes to delete directory. Any other character to exit the program >>  ')
        if choice_value == 'yes' or choice_value == 'YES' :
            logger.info('Deleting the result  directory for a previous execution without cleaning up')
            shutil.rmtree(arguments.outputdir)
            try:
                os.makedirs(arguments.outputdir)
                logger.info ( 'Result folder %s  has been created again', arguments.outputdir)
            except:
                logger.info('Unable to create again the result directory %s', arguments.outputdir)
                print('Cannot create result directory on ', arguments.outputdir)
                exit(0)
        else:
            print('Aborting the execution')
            exit(0)
    evaluate_schema (arguments.inputdir, arguments.outputdir, arguments.alt,  logger)
    
    end_time = datetime.now()
    print('completed execution at :', end_time )
    return True

