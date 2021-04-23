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


# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * #  
# Extract info from schema: duplicates, subsets, quality, lenght statistics, annotation and general info  #
# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * #

def extract_info_schema (schema_files, outputdir, genus, species, usegenus, logger) :
    schema_info_dict = {}
    protein_dict = {}
    schema_sequence_dict ={} ## auxiliar duplicados
    allele_duplicated = {}
    allele_subsets = {} 
    schema_quality_per_class_ids = {} 
    schema_statistics = {} 
    schema_variability_count = {} ## auxiliar estadística longitud
    annotation_core_dict = {}

    print('Analyzing schema...')
    pbar = ProgressBar ()
    for schema_file in pbar (schema_files) :

        gene_name = os.path.basename(schema_file).split('.')[0]
        
        protein_dict[gene_name] = {}
        schema_info_dict[gene_name] = {}
        schema_sequence_dict[gene_name] = {}
        schema_quality_per_class_ids[gene_name] = {'good_quality': [], 'bad_quality: no_start': [], 'bad_quality: no_stop': [], 'bad_quality: no_start_stop': [], 'bad_quality: multiple_stop': []}
        schema_statistics[gene_name] = []
        schema_variability_count[gene_name] = {} ## auxiliar estadística longitud

        alleles_len = []


        # ········································ #
        # Get schema alleles quality for core gene #
        # ········································ #

        locus_quality = check_core_gene_quality(schema_file, logger)

        for allele in locus_quality:
            schema_quality_per_class_ids[gene_name][locus_quality[allele]].append(allele)


        # ············································· #
        # Get gene and product annotation for core gene #
        # ············································· #

        gene_annot, product_annot = get_gene_annotation (gene_name, outputdir, genus, species, usegenus, logger)
        if gene_name not in annotation_core_dict.keys():
            annotation_core_dict[gene_name] = {}
        annotation_core_dict[gene_name] = [gene_annot, product_annot]


        alleles_in_locus = list(SeqIO.parse(schema_file, "fasta"))

        for allele_1 in alleles_in_locus:
            #check if allele id contain string characters . If yes get only the nummber
            if '_' in allele_1.id: ## X
                tmp_id = allele_1.id.split('_')
                allele_1_id = int(tmp_id[-1])
            else:
                allele_1_id = int(allele_1.id)


            # ··························································· #
            # Get alleles which are subsets of other locus schema alleles #
            # ··························································· #

            for allele_2 in alleles_in_locus :
                
                if str(allele_1.seq) in str(allele_2.seq) or str(allele_2.seq) in str(allele_1.seq) :
                    if len(str(allele_1.seq)) != len(str(allele_2.seq)) :
                        if '_' in allele_2.id: ## X
                            tmp_id = allele_2.id.split('_')
                            allele_2_id = int(tmp_id[-1])
                        else:
                            allele_2_id = int(allele_2.id)

                        if len(str(allele_1.seq)) > len(str(allele_2.seq)) :
                            no_subset = allele_1_id
                            subset = allele_2_id

                        else:
                            no_subset = allele_2_id
                            subset = allele_1_id


                        if not gene_name in allele_subsets :
                            allele_subsets [gene_name] = {}

                        if not no_subset in allele_subsets [gene_name] :
                            allele_subsets [gene_name][no_subset] = []  
                        if not subset in allele_subsets [gene_name][no_subset]:
                            allele_subsets [gene_name][no_subset].append(subset)


            sequence = allele_1.seq            
            sequence_str = str(sequence)


            # ··············································································· #
            # Get protein sequence for each CDS encoding allele sequence in each schema locus #
            # ··············································································· #

            query_direction = check_sequence_order(sequence_str, logger)
            if query_direction == 'reverse' : 
                sequence = sequence.reverse_complement()

            try:
                protein = str(sequence.translate(cds=True, table = 11))
                protein_dict[gene_name][allele_1_id] = protein
                coding_cds = 'Yes' 
            except Exception as error:
                logger.error('Not CDS for gene %s in the allele %s ', gene_name, allele_1_id)
                protein = '-' 
                coding_cds = 'No'
                protein_dict[gene_name][allele_1_id] = 'NOT CDS'


            # ························································································································································································ #
            # Create schema info summary including for each allele in each locus: nucleotide sequence, protein sequence, nucleotide sequence length, CDS encoding, allele quality and allele direction #
            # ························································································································································································ #            
    
            for quality_class in schema_quality_per_class_ids[gene_name]:
                if str(allele_1_id) in schema_quality_per_class_ids[gene_name][quality_class]:
                    allele_quality = quality_class

            schema_info_dict[gene_name][allele_1_id] = [sequence_str, str(len(sequence_str)), coding_cds, allele_quality, query_direction, protein]

            # Get core gene alleles length to keep length variability and statistics info
            alleles_len.append(len(sequence_str))


            # ··········································· #
            # Get duplicated alleles in each locus schema #
            # ··········································· #
            if not sequence_str in schema_sequence_dict[gene_name] :
                schema_sequence_dict [gene_name][sequence_str]= [allele_1_id]
            else:
                schema_sequence_dict [gene_name][sequence_str].append(allele_1_id)

        for allele in schema_sequence_dict [gene_name] :
            if len(schema_sequence_dict [gene_name][allele]) > 1 :
                if not gene_name in allele_duplicated :
                    allele_duplicated[gene_name] = []
                allele_duplicated [gene_name].append(sorted(schema_sequence_dict [gene_name][allele]))


        # ······························································· #
        # Get length variability and statistics for alleles in this locus #
        # ······························································· #

        if len(alleles_len) == 1:
            stdev = 0
        else:
            stdev = statistics.stdev(alleles_len)
        schema_statistics[gene_name]=[statistics.mode(alleles_len), statistics.mean(alleles_len), stdev, min(alleles_len), max(alleles_len)]

        for length in list(set(alleles_len)):
            schema_variability_count[gene_name][str(length)] = str(alleles_len.count(length))
    
    return schema_info_dict, schema_quality_per_class_ids, allele_duplicated, allele_subsets, schema_statistics, schema_variability_count, annotation_core_dict, protein_dict


def create_bar_graphic (x_data, y_data, x_label, x_prefix ,y_label, title , rotation, file_name) : ## X
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


def find_proteins_in_gene (raw_proteins_per_genes, logger) : ## X
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


def summary_schema_info (schema_info_dict,  output_dir, logger) : ## X
    logger.info('Start processing the information in schema info')
    header_variability_length = ['Gene name', 'Length variability'] 
    header_gene_length = ['Gene name', 'Length']
    header_percent_allele_not_cds = ['Gene name', 'Allele Percentage that is not coding CDS']
    summary_info = {}
    variability_length = {}
    coding_cds = {}
    #error_type = {}
    allele_quality = {}
    gene_length = {}
    direction = {}
    
    # join all individual information to one item per gene
    for gene in sorted(schema_info_dict) :

        g_length = [] # longitud
        coding_cds[gene] = {} # coding cds
        allele_quality[gene] = {} # tipo de error --> CALIDAD
        direction[gene] = {} # dirección

        logger.debug('dumping g_length for gene  %s ' ,gene)
        for allele in schema_info_dict[gene] :
            values = schema_info_dict[gene][allele]
            g_length.append(int(values[1]))
            #g_coding.append(values[1])
            if not values[2] in coding_cds[gene] :
                coding_cds[gene][values[2]] = 0
            coding_cds[gene][values[2]] += 1
            
            if not values[3] in allele_quality[gene] :
               allele_quality[gene][values[3]] = 0
            allele_quality[gene][values[3]] += 1
            
            if not values[4] in direction [gene]:
                direction[gene][values[4]] = 0
            direction[gene][values[4]] += 1


        mode_length=statistics.mode(g_length)
        min_length = min(g_length) 
        max_length = max(g_length) 
        gene_length[gene] = mode_length 
        variability_length[gene]=format(max((mode_length-min_length), (max_length-mode_length))/mode_length, '.2f')
    
    logger.info('Create the summary folder')
    os.makedirs(os.path.join(output_dir, 'summary'))
    

    #logger.info('Dumping the variability length from the schema to file')
    #variability_length_file =  os.path.join(output_dir, 'summary' , 'variability_length.tsv')
    #save_simple_dict_to_file (variability_length,  header_variability_length, variability_length_file, logger)
    
    '''
    with open (variability_length_file , 'w') as variability_length_fh :
        variability_length_fh.write('\t'.join(header_variability_length) + '\n')
        for gene, value in sorted (variability_length.items()) :
            variability_length_fh.write(gene + '\t' +  value + '\n')
    '''
    
    logger.info('Dumping completed')
    #logger.info('Dumping the gene length from the schema to file')
    #gene_length_file = os.path.join(output_dir, 'summary' , 'gene_length.tsv')
    #save_simple_dict_to_file (gene_length,  header_gene_length, gene_length_file, logger) 

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

    #logger.info('Dumping the allele percentage that are not codings CDS to file')
    #percent_allele_not_coding_file =  os.path.join(output_dir, 'summary' , 'percent_allele_not_coding.tsv')
    #save_simple_dict_to_file (summary_coding_cds,  header_percent_allele_not_cds, percent_allele_not_coding_file, logger)
    
  
    # Create the plot file for the (cdc/non cds) percent relation 
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
    
    # Combine the number of times that the error codo arise when trying to conver to cds
    summary_allele_quality = {}
    error_name = []
    error_value = []
    for gene, errors in allele_quality.items() :
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
                if not error_code in summary_allele_quality :
                    summary_allele_quality[error_code] = 0
                summary_allele_quality[error_code] += value_error
    for error , value in summary_allele_quality.items():
        error_name.append(error)
        error_value.append(value)
    
    # Create the plot file for error types when trying to convert to cds
    allele_quality_graphic_file = os.path.join(output_dir, 'graphic_allele_quality_cds.png')
    rotation = 0
    x_prefix = ''
    create_bar_graphic (error_name, error_value, 'Error type when converting to CDS', x_prefix,  'Number of errors', 'Type of errors that are generated when trying to convert to CDS ' , rotation , allele_quality_graphic_file) 
    
    logger.info('Schema info has been completed processed ')
    
    return True


def save_simple_dict_list_to_files (dict_to_save,  heading_text, folder_name ,file_name, logger) : ## X
    logger.info('Saving file %s', file_name)
    for gene , value_list in sorted(dict_to_save.items()):
        f_name = os.path.join(folder_name, str(gene + file_name))
        with open (f_name , 'w') as f_name_fh :
            f_name_fh.write('\t'.join(heading_text) + '\n')
            for item in value_list :
                f_name_fh.write(gene + '\t' + item + '\n')
    logger.info('Saved file  %s', file_name)
    return True


def save_simple_dict_to_file (dict_to_save,  heading_text, file_name, logger) : ## X
    logger.info('Saving file %s', file_name)
    with open (file_name , 'w') as file_name_fh :
        file_name_fh.write('\t'.join(heading_text) + '\n')
        for gene , value in sorted (dict_to_save.items()) :
            file_name_fh.write(gene + '\t' + str(value) + '\n')
    logger.info('Saved file  %s', file_name)
    return True


def summary_proteins (raw_proteins_per_genes, output_dir, logger) : ## X
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
    
    # create the diagram to display the percent proteins for each gene
    # round number to 1 decimal to show the graphic
    all_percent = []
    percent_values = proteins_percent_per_gene.values()
    for percent_value in percent_values :
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


# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * ·  #  
# Keep info extracted from schema: duplicates, subsets, quality, lenght statistics, annotation and general info  #
# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * ·  #

#def analyze_schema (inputdir, outputdir, alt_codon_start, logger) :
def analyze_schema (inputdir, outputdir, genus, species, usegenus, logger) :

    header_schema_quality = ['Core Gene', 'Good quality', 'Bad quality: no start', 'Bad quality: no stop', 'Bad quality: no start stop', 'Bad quality: multiple stop', 'Total']
    header_schema_statistics = ['Core Gene', 'Mode', 'Mean', 'Standard Deviation', 'Min Length', 'Max Length', 'Schema Variability', 'Total']
    header_annotation = ['Core Gene', 'Gene Annotation', 'Product Annotation'] ####
    header_alleles_duplicated = ['Core Gene', 'Duplicated alleles group IDs' ]
    header_alleles_subsets = ['Core Gene', 'Allele', 'Subsets alleles group IDs' ] 
    header_schema_info = ['Core Gene', 'Allele', 'Nucleotides Sequence', 'Length', 'Encoding CDS' , 'Allele Quality', 'Direction', 'Protein Sequence']
    
    schema_files = get_fasta_file_list(inputdir, logger)
    
    logger.info('Extract the raw information for each gene in the schema')

    #schema_info_dict, schema_quality_per_class_ids, allele_duplicated, allele_subsets, raw_proteins_per_genes = extract_info_schema (schema_files, logger)
    schema_info_dict, schema_quality_per_class_ids, allele_duplicated, allele_subsets, schema_statistics, schema_variability_count, annotation_core_dict, raw_proteins_per_genes = extract_info_schema (schema_files, outputdir, genus, species, usegenus, logger)

    print('Saving data to ', outputdir )
    logger.info('Start dumping the raw information to files')
    os.makedirs(os.path.join(outputdir, 'raw_info'))

    # Saving schema info to file
    logger.info('Saving schema info to file..')
    os.makedirs(os.path.join(outputdir, 'raw_info', 'raw_schema_information'))
    for core in sorted(schema_info_dict) :
        schema_info_file = os.path.join(outputdir, 'raw_info', 'raw_schema_information', str(core + '_schema_information.tsv'))
        with open (schema_info_file , 'w') as schema_info_fh :
            schema_info_fh.write('\t'.join(header_schema_info) + '\n')
            for allele in (schema_info_dict[core]) :
                schema_info_fh.write(core + '\t' + str(allele) + '\t' + '\t'.join(schema_info_dict[core][allele]) + '\n')

    # Saving duplicated alleles to file
    logger.info('Saving duplicated alleles to file..')
    allele_duplicated_file =  os.path.join(outputdir, 'raw_info' , 'duplicated_alleles.tsv')
    with open (allele_duplicated_file , 'w') as allele_duplicated_fh :
        allele_duplicated_fh.write('\t'.join(header_alleles_duplicated) + '\n')
        for core in sorted(allele_duplicated) :
            for duplication in (allele_duplicated[core]):
                allele_duplicated_fh.write(core + '\t'  + ', '.join(map(str, list(duplication))) + '\n')
 
    # Saving alleles subsets to file
    logger.info('Saving subsets alleles to file..') 
    #os.makedirs(os.path.join(outputdir, 'raw_info', 'subsets_alleles'))
    allele_subsets_file =  os.path.join(outputdir, 'raw_info' , 'alleles_subsets.tsv')    
    with open (allele_subsets_file , 'w') as allele_subsets_fh :
        allele_subsets_fh.write('\t'.join(header_alleles_subsets) + '\n')
        for core in sorted(allele_subsets) :
            for allele_id in allele_subsets[core]: 
                allele_subsets_fh.write(core + '\t' + str(allele_id) + '\t' + ', '.join(map(str, list(allele_subsets[core][allele_id]))) + '\n')

    # Saving schema quality to file
    logger.info('Saving schema quality information to file..')
    quality_file =  os.path.join(outputdir, 'raw_info', 'schema_quality.tsv')
    with open (quality_file , 'w') as quality_fh :
        quality_fh.write('\t'.join(header_schema_quality) + '\n')
        for core in sorted(schema_quality_per_class_ids) :
            len_quality_class_type = [len(schema_quality_per_class_ids[core]['good_quality']), len(schema_quality_per_class_ids[core]['bad_quality: no_start']), \
                                        len(schema_quality_per_class_ids[core]['bad_quality: no_stop']), len(schema_quality_per_class_ids[core]['bad_quality: no_start_stop']), \
                                        len(schema_quality_per_class_ids[core]['bad_quality: multiple_stop'])]

            ### orden alfabético? (['bad_quality: multiple_stop', 'bad_quality: no_start', 'bad_quality: no_start_stop', 'bad_quality: no_stop', 'good_quality'])            
            #len_quality_class_type = []
            #for quality_class in sorted(schema_quality_per_class_ids[core]):
             #   len_quality_class_type.append(len(quality_class))
            #len_quality_class_type = [len(value) for value in list(schema_quality_per_class_ids[core].values())] 

            quality_fh.write(core + '\t' + '\t'.join (map(str, len_quality_class_type)) + '\t' + str(sum(len_quality_class_type)) + '\n')

    # Saving length statistics to file
    logger.info('Saving schema length statistics information to file..')
    statistics_file =  os.path.join(outputdir, 'raw_info', 'length_statistics.tsv')
    with open (statistics_file , 'w') as stat_fh :
        stat_fh.write('\t'.join(header_schema_statistics) + '\n')
        for core in sorted (schema_statistics):
            length_number = []
            total_alleles = 0
            for length in schema_variability_count[core]:
                length_number.append(length + ': ' + schema_variability_count[core][length])
                total_alleles += int(schema_variability_count[core][length])

            stat_fh.write(core + '\t' + '\t'.join (map(str,schema_statistics[core])) + '\t' + ', '.join(length_number) + '\t' + str(total_alleles) + '\n')

    # Saving schema annotation to file
    logger.info('Saving core gene schema annotation to file..') 
    annotation_file =  os.path.join(outputdir, 'raw_info' , 'annotation.tsv')    
    with open (annotation_file , 'w') as annot_fh :
        annot_fh.write('\t'.join(header_annotation) + '\n')
        for core in sorted(annotation_core_dict) :
            annot_fh.write(core + '\t' + '\t'.join(annotation_core_dict[core]) + '\n')

    logger.info('Completed dumped raw information to files')
        
    #summary_schema_info(schema_info_dict, outputdir, logger) 
    #summary_proteins (raw_proteins_per_genes, outputdir, logger) 

    return schema_quality_per_class_ids, allele_duplicated, allele_subsets, schema_files


# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * ·  #  
# Filter schema removing subsets, duplicates and bad quality alleles from each locus #
# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * ·  #

def remove_alleles_from_schema (schema_files, remove_subsets, remove_duplicates, remove_no_cds, new_schema, allele_subsets, allele_duplicated, schema_quality_per_class, schema_dir, outputdir, logger): 
                                ## logger?

    ## Create a copy of core genes schema if updateschema = 'new' / 'New'
    if new_schema == 'true':
      #  print("Ha entrado a true new schema para filtrar el esquema", '\n')
        no_filtered_schema_dir = schema_dir
        schema_dir_name = os.path.basename(no_filtered_schema_dir)
        schema_dir = os.path.join(outputdir, schema_dir_name + '_filtered') 
        shutil.copytree(no_filtered_schema_dir, schema_dir)
        logger.info('Copying core genes fasta files to filter schema')
        
        schema_files = get_fasta_file_list(schema_dir, logger)


    print('\n', 'Filtering schema...')
    pbar = ProgressBar ()
    for schema_file in pbar (schema_files) :
        core_name = os.path.basename (schema_file).split('.') [0]
        
        alleles_to_remove = []

        if remove_subsets == 'true':           
            if core_name in allele_subsets:
                subsets_alleles = sum(list(allele_subsets [core_name].values()), []) 
                alleles_to_remove += [str(x) for x in subsets_alleles]

        if remove_duplicates == 'true':
            if core_name in allele_duplicated:
                for duplicates_group in allele_duplicated[core_name]:
                    for id_index in range(1, len(duplicates_group)):
                        alleles_to_remove += str(duplicates_group[id_index])

        if remove_no_cds == 'true':
            for quality_class in schema_quality_per_class [core_name]:
                if 'bad_quality' in quality_class:
                    alleles_to_remove += schema_quality_per_class [core_name][quality_class]

        alleles_to_remove_unique = list(set(alleles_to_remove))

        alleles_in_locus_dict = {}
        allele_str_id = ''
        
        for allele in SeqIO.parse(schema_file, 'fasta'):
            if '_' in str(allele.id):
                split_id = str(allele.id).split('_')
                allele_id = int(split_id[-1])
                allele_str_id = '_'.join(split_id[0:len(split_id)])
            else:
                allele_id = int(allele.id)
            
            alleles_in_locus_dict[allele_id] = str(allele.seq)


        with open(schema_file, 'w') as schema_fh :
            for allele_id in sorted(alleles_in_locus_dict):
                if str(allele_id) not in alleles_to_remove_unique:
                    if len(allele_str_id) > 0:
                        allele_id_comp = allele_str_id + '_' + str(allele_id)
                    else:
                        allele_id_comp = str(allele_id)
                    schema_fh.write('>' + allele_id_comp + '\n' + alleles_in_locus_dict[allele_id] + '\n' + '\n' )
                                                                                                                                 
    return True


# · * · * · * · * · * · * · * · * · * · * · * · * #  
# Processing schema analysis and schema filtering #
# · * · * · * · * · * · * · * · * · * · * · * · * #

def processing_analyze_schema(arguments) :
    
    start_time = datetime.now()
    print('Start the execution at :', start_time )
    
    # Open log file
    logger = open_log ('analyze_schema.log')

    #############################
    ## Create output directory ##
    #############################
    try:
        os.makedirs(arguments.outputdir)
    except:
        logger.info('Deleting the result  directory for a previous execution without cleaning up')
        shutil.rmtree(arguments.outputdir)
        try:
            os.makedirs(arguments.outputdir)
            logger.info ( 'Results folder %s  has been created again', arguments.outputdir)
        except:
            logger.info('Unable to create again the result directory %s', arguments.outputdir)
            print('Cannot create result directory on ', arguments.outputdir)
            exit(0)

    #########################
    ## Get schema analysis ##
    #########################
    #analyze_schema (arguments.inputdir, arguments.outputdir, arguments.alt,  logger)
    schema_quality_per_class_ids, allele_duplicated, allele_subsets, schema_files = analyze_schema (arguments.inputdir, arguments.outputdir, arguments.genus, arguments.species, arguments.usegenus, logger)        
    if not schema_quality_per_class_ids:
        print('There is an error while processing the schema analysis. Check the log file to get more information \n')
        exit(0)

    ###################################################################################
    ## Remove allele subsets, duplicated alleles and bad quality alleles from schema ##
    ###################################################################################
    if str(arguments.removesubsets).lower() == 'true' or str(arguments.removeduplicates).lower() == 'true' or str(arguments.removenocds).lower() == 'true' :
        if not remove_alleles_from_schema (schema_files, str(arguments.removesubsets).lower(), str(arguments.removeduplicates).lower(), str(arguments.removenocds).lower(), str(arguments.newschema).lower(), allele_subsets, allele_duplicated, schema_quality_per_class_ids, arguments.inputdir, arguments.outputdir, logger):
            print('There is an error while processing the schema allele filtering. Check the log file to get more information \n')
            exit(0)

    end_time = datetime.now()
    print('completed execution at :', end_time )
    
    return True

