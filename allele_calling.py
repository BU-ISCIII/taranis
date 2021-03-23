#!/usr/bin/env python3

import argparse
import sys
import io
import os
import re
import statistics
import logging
from logging.handlers import RotatingFileHandler

from datetime import datetime
import glob
import pickle

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Blast.Applications import NcbiblastnCommandline
from io import StringIO
from Bio.Blast import NCBIXML
from BCBio import GFF
import pandas as pd
import shutil
from progressbar import ProgressBar

from utils.taranis_utils import *

import math ### c/m: import math para obtención length_trhesholds
import numpy as np ### c/m: import numpy para obtener alelo de referencia para cada locus
import csv ### c/m: para obtener anotación gen y producto de cada locus

def check_blast (reference_allele, sample_files, db_name, logger) : # No se está usando
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

def parsing_fasta_file_to_dict (fasta_file, logger): # No se usa si no se obtienen diccionarios de locus y muestras
    fasta_dict = {}
    fasta_dict_ordered = {}
    for contig in SeqIO.parse(fasta_file, "fasta", generic_dna):
        fasta_dict[str(contig.id)] = str(contig.seq.upper())
    logger.debug('file %s parsed to dictionary', fasta_file)

    for key in sorted(list(fasta_dict.keys())):
        fasta_dict_ordered[key] = fasta_dict[key]
    return fasta_dict_ordered

def prepare_core_gene (core_gene_file_list, store_dir, ref_alleles_dir, outputdir, logger): # añadiendo outputdir para guardar report de calidad y stat en el directorio de resultados
    #check if directory exists and files have fasta files
    #valid_core_gene_files = get_fasta_file_list(core_gene_dir, logger)
    #if not valid_core_gene_files :
    #    return False
    #logger.debug('Schema files to be processed are : %s', valid_core_gene_files)
    
    #processing the files in the schema
    schema_quality = {} ### c/m: introducción del diccionario de calidad del esquema incluyendo todos los locus
    annotation_core_dict = {} ### c/m: diccionaro locus-anotación
    schema_variability = {}
    schema_statistics = {}
    file_list = []
  ###  reference_alleles_list = [] ### CAMBIANDO/MODIFICANDO: añadiendo lista para meter el nombre de los archivos de los alelos de referencia
    ### first_alleles_list = [] ### cambiando/modificando: comentado todo lo relacionado con el uso del primer alelo del locus como referencia
    
    schema_variability_count = {} # diccionario para report estadisticas longitud
    schema_quality_per_class = {} ### diccionario para report de calidad

    alleles_in_locus_dict = {} # diccionario para ir guardando, para cada locus, id-secuencia de todos los alelos

    blast_dir = os.path.join(store_dir,'blastdb')
    logger.info('start preparation  of core genes files')
    for fasta_file in core_gene_file_list:

        f_name = os.path.basename(fasta_file).split('.')
        file_list.append(os.path.join(store_dir, f_name[0]))

        # parsing fasta file and get in the dictionary the id and the sequence
        fasta_file_parsed_dict = parsing_fasta_file_to_dict(fasta_file, logger) #--------> Generacion de diccionario de alelos por locus comentado. Vuelvo a generar los diccionarios pero esta vez sin generar archivos para ver si así tarda menos tiempo
        
        if f_name[0] not in alleles_in_locus_dict.keys():
            alleles_in_locus_dict[f_name[0]] = {}
        alleles_in_locus_dict[f_name[0]] = fasta_file_parsed_dict ### guardando el diccionario id-secuencia de este locus en el diccionario general con todos los locus

        # dump fasta file into pickle file
        #with open (file_list[-1],'wb') as f:
         #   pickle.dump(fasta_file_parsed_dict, f) --------> Generacion de diccionario de alelos por locus comentado
        
        # obteniendo calidad del esquema
        #locus_quality = check_core_gene_quality(fasta_file_parsed_dict, logger) ### cambiando/modificando: obteniendo diccionario de calidad de cada locus
        locus_quality = check_core_gene_quality(fasta_file, logger) ### cambiando/modificando: obteniendo diccionario de calidad de cada locus introduciendo path al fasta del locus en cuestión en lugar del diccionario
        ### INTRODUCIR PARSEO DEL ARCHIVO EN LUGAR DE DICCIONARIO
    
        if f_name[0] not in schema_quality.keys(): ### cambiando/modificando: guardando output de check_core_gene_quality para cada locus en dict schema_quality
            schema_quality[f_name[0]] = {}   ### no sé si esto es necesario o si directamente puedo hacer lo de abajo
        schema_quality[f_name[0]] = locus_quality

        ### cambiando/modificando: he comentado todo lo que tiene que ver con el uso del primer alelo del locus como referencia
        # create the first allele for each core gene file
        #### used only for gene annotation
        #### first_alleles_list.append(write_first_allele_seq(fasta_file, store_dir, logger))
        
     #   reference_alleles_list.append(get_reference_allele(locus_quality, fasta_file, store_dir, logger)) ### CAMBIANDO/MODIFICANDO: esta lista sustituye a la anterior (first_alleles_list). Listado de los paths de los archivos que contienen el alelo de referencia

        ###print("reference_alleles_list: ", reference_alleles_list) ### print, borrar

        ### core schema annotation  #########

        # se obtiene anotacion gen/producto de cada locus a partir del alelo de referencia
        ref_allele = os.path.join(ref_alleles_dir, f_name[0] + '.fasta') # path al alelo de referencia del locus en cuestion

        gene_annot, product_annot = get_gene_annotation (ref_allele, store_dir, logger)
        if f_name[0] not in annotation_core_dict.keys(): ### cambiando/modificando: guardando output de check_core_gene_quality para cada locus en dict schema_quality
            annotation_core_dict[f_name[0]] = {}   ### no sé si esto es necesario o si directamente puedo hacer lo de abajo
        annotation_core_dict[f_name[0]] = [gene_annot, product_annot]

        ### Obteniendo longitudes de los alelos del locus
        #alleles_len = []
        #for allele in fasta_file_parsed_dict : ### sacando a partir de diccionario, comentado
         #   alleles_len.append(len(fasta_file_parsed_dict[allele]))

        alleles_len = []
        alleles_in_locus = list (SeqIO.parse(fasta_file, "fasta"))
        for allele in alleles_in_locus :
            alleles_len.append(len(str(allele.seq)))

        ######## schema_variability_count para report estadisticas
        schema_variability_count[f_name[0]] = {} ## añadiendo locus al diccionario
        for length in list(set(alleles_len)): ## por cada longitud
            schema_variability_count[f_name[0]][str(length)] = str(alleles_len.count(length))
        ########

        schema_variability[f_name[0]]=list(set(alleles_len))
        schema_statistics[f_name[0]]=[statistics.mode(alleles_len), statistics.mean(alleles_len), statistics.stdev(alleles_len), min(alleles_len), max(alleles_len)] ### Cambiando/modificando: he añadido mean y stdev de las longitudes de los alelos de cada locus para usarlo para la clasificación de las secuencias encontradas en cada tipo

    ############################
    ###### Report calidad ######
    ############################

    header_schema_quality = ['Core gene', 'Good quality', 'Bad quality: no start', 'Bad quality: no stop', 'Bad quality: no start stop', 'Bad quality: multiple stop', 'Total']

    ### Dos formas de obtener el report de calidad
    # la primera que hice: 
    for core_gene in schema_quality:
        schema_quality_per_class[core_gene] = {'good_quality': 0, 'bad_quality: no_start': 0, 'bad_quality: no_stop': 0, 'bad_quality: no_start_stop': 0, 'bad_quality: multiple_stop': 0}
        for quality_class in schema_quality_per_class[core_gene]:
            #print("quality_class: ", quality_class, '\n')
            for allele in schema_quality[core_gene]:
                if schema_quality[core_gene][allele] == quality_class:
                    schema_quality_per_class[core_gene][quality_class] += 1

    # la segunda que he hecho, que aún no he probado, mirar si funciona:
    # a esta segunda forma le he añadido la creación del diccionario donde guardar los ids de los alelos del locus con cada tipo de calidad, además del diccionario de la contabilidad
    # de hecho, igual se podría sacar únicamente este primer diccionario y para saber la contabilidad pues hacer un len de cada tipo de calidad y ya
    # el dict de ids me hace falta para saber qué alelos borrar del esquema en analyze schema

    """
    schema_quality_per_class_ids = {}
    
    for core_gene in schema_quality:
        schema_quality_per_class[core_gene] = {}
        for allele in schema_quality[core_gene]:
            if schema_quality[core_gene][allele] not in schema_quality_per_class[core_gene].keys():
                schema_quality_per_class[core_gene][schema_quality[core_gene][allele]] = 0 ## para llevar contabilidad de alelos de cada tipo de calidad
                schema_quality_per_class_ids[core_gene][schema_quality[core_gene][allele]] = [] ## para guardar id de alelos de cada tipo de calidad
            schema_quality_per_class[core_gene][schema_quality[core_gene][allele]] += 1
            schema_quality_per_class_ids[core_gene][schema_quality[core_gene][allele]].append(allele)
    """
    



    #print("schema_quality_per_class: ", schema_quality_per_class, '\n')

    # saving schema quality to file
    logger.info('Saving schema quality information to file..')
    quality_file =  os.path.join(outputdir, 'schema_quality.tsv')
    with open (quality_file , 'w') as quality_fh :
        quality_fh.write('\t'.join(header_schema_quality) + '\n')
        for core in sorted (schema_quality_per_class) :

           # allele_number = []
            #for quality_class in schema_quality_per_class[core]:
             #   allele_number.append(schema_quality_per_class[core][quality_class])
            quality_fh.write(core + '\t' + '\t'.join (map(str, list(schema_quality_per_class[core].values()))) + '\t' + str(sum(list(schema_quality_per_class[core].values()))) + '\n')

    ##################################################
    ###### Report estadisticas longitud esquema ######
    ##################################################

    header_schema_statistics = ['Core gene', 'Mode', 'Mean', 'Standard deviation', 'Min length', 'Max length', 'Schema variability', 'Total']

    # saving lengtt statistics to file
    logger.info('Saving schema length statistics information to file..')
    statistics_file =  os.path.join(outputdir, 'length_statistics.tsv')
    with open (statistics_file , 'w') as stat_fh :
        stat_fh.write('\t'.join(header_schema_statistics) + '\n')
        for core in sorted (schema_statistics):
            length_number = []
            total_alleles = 0
            for length in schema_variability_count[core]:
                length_number.append(length + ': ' + schema_variability_count[core][length]) ## disponiendo length : numero de alelos con esa length
                total_alleles += int(schema_variability_count[core][length]) ## numero total de alelos en el locus

            stat_fh.write(core + '\t' + '\t'.join (map(str,schema_statistics[core])) + '\t' + ', '.join(length_number) + '\t' + str(total_alleles) + '\n')

  #  print("schema_quality: ", schema_quality, '\n')

    return file_list, alleles_in_locus_dict, annotation_core_dict, schema_variability, schema_statistics, schema_quality ### cambiando/modificando: he borrado first_alleles_list del output de la función, he añadido reference_alleles_list y el diccionario schema_quality



def prodigal_training(reference_genome_file, prodigal_dir, logger): ### cambiando/modificando: función añadida para generar el archivo de training a partir del genoma de referencia 
    ### No sé si poner el training file para que se meta desde fuera a taranis, de modo que hay que hacer el paso previo de prodigal o si es mejor que genere el training file dentro y se le pasaría a taranis el genoma de referencia para generar el training file (si lo hago de esta última forma cada vez que se ejecute taranis va a generar un archivo de entrenamiento de prodigal y no se puede skipear la generación del archivo) 

    f_name = os.path.basename(reference_genome_file).split('.')[0] ### reference_genome_file_name sería el path del genoma de ref. Se está cogiendo solo el nombre del archivo
    prodigal_train_dir = os.path.join(prodigal_dir, 'training')

    output_prodigal_train_dir = os.path.join(prodigal_train_dir, f_name + '.trn')
    
    if not os.path.exists(prodigal_train_dir):
        try:
            os.makedirs(prodigal_train_dir)
            logger.debug('Created prodigal directory for training file %s', f_name)
        except:
            logger.info('Cannot create prodigal directory for training file %s', f_name)
            print ('Error when creating the directory %s for training file', prodigal_train_dir) 
            exit(0)

        prodigal_command = ['prodigal' , '-i', reference_genome_file, '-t', output_prodigal_train_dir]
        prodigal_result = subprocess.run(prodigal_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    ### He comentado este if porque prodigal devuelve comentarios del proceso que se guardan en stderr y me lo toma como errores y que no se ha podido generar el archivo de training, cuando realmente sí
      #  if prodigal_result.stderr:
       #     logger.error('cannot create training file for %s', f_name)
        #    logger.error('prodigal returning error code %s', prodigal_result.stderr)
         #   return False
    else:
        logger.info('Skeeping prodigal training file creation for %s, as it has already been created', f_name)
    
    return output_prodigal_train_dir


def prodigal_prediction(file_name, prodigal_dir, prodigal_train_dir, logger): ### cambiando/modificando: función añadida para la predicción de genes en las muestras con prodigal

    f_name = os.path.basename(file_name).split('.')[0] ### file_name sería el path de cada muestra. Se coge el nombre de solo una muestra en cuestión
    prodigal_dir_sample = os.path.join(prodigal_dir,f_name)

    output_prodigal_coord = os.path.join(prodigal_dir_sample, f_name + '_coord.gff') # no sería necesario
    output_prodigal_prot = os.path.join(prodigal_dir_sample, f_name + '_prot.faa') # no sería necesario
    output_prodigal_dna = os.path.join(prodigal_dir_sample, f_name + '_dna.faa')

    if not os.path.exists(prodigal_dir_sample):
        try:
            os.makedirs(prodigal_dir_sample)
            logger.debug('Created prodigal directory for Core Gene %s', f_name)
        except:
            logger.info('Cannot create prodigal directory for Core Gene %s' , f_name)
            print ('Error when creating the directory %s for prodigal genes prediction', prodigal_dir_sample)
            exit(0)

        prodigal_command = ['prodigal' , '-i', file_name , '-t', prodigal_train_dir, '-f', 'gff', '-o', output_prodigal_coord, '-a', output_prodigal_prot, '-d', output_prodigal_dna]
        prodigal_result = subprocess.run(prodigal_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        ### He comentado este if porque prodigal devuelve comentarios del proceso que se guardan en stderr y me lo toma como errores y que no se ha podido generar el archivo de training, cuando realmente sí

        # if prodigal_result.stderr:
          #  logger.error('cannot predict genes for %s ', f_name)
           # logger.error('prodigal returning error code %s', prodigal_result.stderr)
            #return False
    else:
        logger.info('Skeeping prodigal genes prediction for %s, as it has already been made', f_name)
    
    return True


# Función get_prodigal_sequence anterior
"""
def get_prodigal_sequence(sstart, send, contig_id, prodigal_directory, sample_name): ### Cambiando/modificando: función para obtener la secuencia completa del gen de interés a partir de la predicción de prodigal
    # sstart = posición inicio sseq BLAST
    # send = posición final sseq BLAST
    # prodigal_directory* = directorio donde se encuentran los archivos de prodigal si se va a parsear directamente desde aquí y no se va a crear un diccionario de cada muestra previamente
    # sample_name* = nombre de la sample en cuestión para generar el path completo para los datos de esta muestra. Solo si se va a parsear directamente desde aquí y no se va a crear un diccionario de cada muestra previamente

    # prodigal_directory = os.path.join(tmp_samples_dir,'prodigal') ### para saber qué es prodigal_directory
    # sample_name = os.path.basename(sample_file) ### para saber qué es sample_name. Esto se le pasaría así, el nombre de la muestra en cuestión, si se parsea directamente desde allele_call, si no habría que pasar el sample_file que es el path completo de cada muestra a donde se cree diccionario o lo que sea de la muestra, y aquí habría que pasar el diccionario en cuestión, habría que pasar el path del diccionario, o el diccionario en sí se se parsea el diccionario al inicio de cada muestra, así no hay que parsear el diccionario para cada locus

    prodigal_directory_sample = os.path.join(prodigal_directory, sample_name)
    #coords_file = os.path.join(prodigal_directory_sample, sample_name + '_coord.gff')
    genes_file = os.path.join(prodigal_directory_sample, sample_name + '_dna.faa')

    #records = GFF.parse(coords_file)

    if int(sstart) < int(send):
        start_blast = int(sstart)
        end_blast = int(send)
    else:
        start_blast = int(send)
        end_blast = int(sstart)
    
    dist_final = 10000000 ### número muuuy grande
    id_final = ''
    new_sseq = ''
    
    
    #for rec in records:
     #   for gene_feature in rec.features:
      #      if (gene_feature.id).startswith(str(contig_id)): ### mira si el id del gen predicho empieza por el número del contig y si no, pasa al siguiente gen
                
       #         start_prodigal = int(gene_feature.location.start)
        #        end_prodigal = int(gene_feature.location.end)
              #  print("end_prodigal: ", end_prodigal, '\n')
              #  print("end_blast: ", end_blast, '\n')
              #  print("sstart: ", send, '\n', '\n')

         #       dist_start = abs(start_prodigal - start_blast)
              #  print('dist_start: ', dist_start, '\n')
          #      dist_end = abs(end_prodigal - end_blast)
              #  print('dist_end: ', dist_end, '\n')
           #     dist_total = dist_start + dist_end
              #  print('dist_total: ', dist_total, '\n', '\n')

            #    if dist_total < dist_final:
             #       dist_final = dist_total
              #      id_final = gene_feature.id

         #   else:
          #      next
    

    predicted_genes = list(SeqIO.parse(genes_file, "fasta"))

    for rec in predicted_genes:
        if (rec.id).startswith(str(contig_id)): ### mira si el id del gen predicho empieza por el número del contig y si no, pasa al siguiente gen

            start_prodigal = int(rec.description.split('#')[1])
            end_prodigal = int(rec.description.split('#')[2])

            dist_start = abs(start_prodigal - start_blast)
            dist_end = abs(end_prodigal - end_blast)
            dist_total = dist_start + dist_end

            if dist_total < dist_final:
                #print("rec.id: ", rec.id)
                #print("dist_start: ", dist_start)
                #print("dist_end: ", dist_end)
                #print("dist_total: ", dist_total)

                dist_final = dist_total
                id_final = rec.id
                new_sseq = rec.seq
    
    print("El id_final de la seq predicha por prodigal es: ", id_final, '\n')

    ### esto es para printear y ver la info de la sec que pilla
    #for rec_2 in predicted_genes:
     #   if rec_2.id == id_final:
      #      start_prodigal = int(rec_2.description.split('#')[1])
       #     end_prodigal = int(rec_2.description.split('#')[2])
        #    description_final = rec_2.description
         #   rec_final = rec_2
          #  dist_start = abs(start_prodigal - start_blast)
           # dist_end = abs(end_prodigal - end_blast)
           # dist_total = dist_start + dist_end

    
    #for rec in predicted_genes:
     #   if rec.id == id_final:
      #      break

    #new_sseq = rec.seq
    
    #print("rec (fasta dna): ", rec_final, '\n')
    #print("rec.description (fasta dna): ", description_final, '\n', '\n')

    #print('contig_id: ', contig_id, '\n', '\n')

    #print("start_prodigal_fasta: ", start_prodigal, '\n')
    #print("start_blast: ", start_blast, '\n')
    #print("sstart: ", sstart, '\n', '\n')

    #print("end_prodigal_fasta: ", end_prodigal, '\n')
    #print("end_blast: ", end_blast, '\n')
    #print("send: ", send, '\n', '\n')

    #print('dist_start: ', dist_start, '\n')
    #print('dist_end: ', dist_end, '\n')
    #print('dist_total: ', dist_total, '\n', '\n')

    #print('id_final: ', id_final, '\n', '\n')

    return new_sseq
"""

### get_prodigal_sequence con make blast db y BLAST
def get_prodigal_sequence(blast_sseq, contig_blast_id, prodigal_directory, sample_name, blast_parameters, logger): ### Cambiando/modificando: función para obtener la secuencia completa del gen de interés a partir de la predicción de prodigal
    # blast_sseq = secuencia encontrada tras BLAST sin gaps
    # contig_blast_id = id del contig de la muestra donde se ha encontrado la secuencia en BLASt --> Para coger los genes de prodigal que se han predicho en este contig para pillar esas 
        # secuencias para comparar con la secuencia encontrada con BLAST para intentar pillar la secuencia predicha con prodigal más similar a la encontrada con BLAST.
    # prodigal_directory
    # sample_name (nombre de la muestra)
    # no haría falta sstart y send pero dejar para generar informe con todas las secuencias que se sacan y sus sstart y end en prodigal junto el sstart y end en blast?
    
    ### argumentos anteriores
    # sstart = posición inicio sseq BLAST
    # send = posición final sseq BLAST
    # prodigal_directory* = directorio donde se encuentran los archivos de prodigal si se va a parsear directamente desde aquí y no se va a crear un diccionario de cada muestra previamente
    # sample_name* = nombre de la sample en cuestión para generar el path completo para los datos de esta muestra. Solo si se va a parsear directamente desde aquí y no se va a crear un diccionario de cada muestra previamente

    # prodigal_directory = os.path.join(tmp_samples_dir,'prodigal') ### para saber qué es prodigal_directory
    # sample_name = os.path.basename(sample_file) ### para saber qué es sample_name. Esto se le pasaría así, el nombre de la muestra en cuestión, si se parsea directamente desde allele_call, si no habría que pasar el sample_file que es el path completo de cada muestra a donde se cree diccionario o lo que sea de la muestra, y aquí habría que pasar el diccionario en cuestión, habría que pasar el path del diccionario, o el diccionario en sí se se parsea el diccionario al inicio de cada muestra, así no hay que parsear el diccionario para cada locus

    #print("contig_blast_id: ", contig_blast_id, '\n')

    prodigal_directory_sample = os.path.join(prodigal_directory, sample_name)
    genes_file = os.path.join(prodigal_directory_sample, sample_name + '_dna.faa')


    ### Creacion del directorio mash en directorio prodigal de la muestra en cuestion
    ### aqui se van a guardar las bases de datos de los fastas con genes por contig
    blastdb_per_contig_directory = 'blastdb_per_contig'
    full_path_blastdb_per_contig = os.path.join(prodigal_directory_sample, blastdb_per_contig_directory)
    if not os.path.exists(full_path_blastdb_per_contig):
        try:
            os.makedirs(full_path_blastdb_per_contig)
            logger.info('Directory %s has been created', full_path_blastdb_per_contig)
        except:
            print ('Cannot create the directory ', full_path_blastdb_per_contig)
            logger.info('Directory %s cannot be created', full_path_blastdb_per_contig)
            exit (0)

    
    ### Creacion del directorio prodigal_genes_per_contig dentro del directorio mash en directorio prodigal de la muestra en cuestion
    ### Aqui se van a guardar los fastas con genes por contig
    prodigal_genes_per_contig_directory = 'prodigal_genes_per_contig'
    full_path_prodigal_genes_per_contig = os.path.join(prodigal_directory_sample, prodigal_genes_per_contig_directory)
    if not os.path.exists(full_path_prodigal_genes_per_contig):
        try:
            os.makedirs(full_path_prodigal_genes_per_contig)
            logger.info('Directory %s has been created', full_path_prodigal_genes_per_contig)
        except:
            print ('Cannot create the directory ', full_path_prodigal_genes_per_contig)
            logger.info('Directory %s cannot be created', full_path_prodigal_genes_per_contig)
            exit (0)

    ### Parseando el archivo de prodigal que contiene todos los genes predichos para la muestra en cuestión
    predicted_genes = SeqIO.parse(genes_file, "fasta")

    #print("predicted_genes: ", predicted_genes, '\n')

    #for rec in SeqIO.parse(genes_file, 'fasta'):

    ### Generando fasta con las secuencias de los genes predichos por prodigal en este contig
    contig_genes_path = os.path.join(full_path_prodigal_genes_per_contig, contig_blast_id + '.fasta')
    
    with open (contig_genes_path, 'w') as out_fh:
       # print("contig_genes_path: ", contig_genes_path, '\n')
        for rec in predicted_genes: ### por cada uno de los genes predichos pro prodigal
            #print("rec in predicted_genes: ", rec, '\n')
            contig_prodigal_id = (rec.id).split("_")[0]
        #    print("for rec in predicted_genes, if contig_prodigal_id == contig_blast_id", '\n')
            if contig_prodigal_id == contig_blast_id: ### mira si el id del gen predicho por prodigal se ha encontrado en el mismo contig que la secuencia encontrada por blast
         #       print("contig_prodigal_id: ", contig_prodigal_id, '\n')
                out_fh.write ('>' + str(rec.description) + '\n' + str(rec.seq) + '\n') ### poniendo descripcion como cabecera del fasta para poder sacar posteriormente start y end para report prodigal

    # cerando base de datos local BLAST para genes prodigal en contig de interes
    if not create_blastdb(contig_genes_path, full_path_blastdb_per_contig, 'nucl', logger):
        print('Error when creating the blastdb for samples files. Check log file for more information. \n ')
        return False
    #blastdb = create_blastdb(contig_genes_path, full_path_blastdb_per_contig, 'nucl', logger) ### creacion de base de datos blast de los genes del contig en cuestion mediante ceate_blastdb

    #print("se ha creado la base de datos blastdb: ", blastdb, '\n') #### blastdb para saber si se genera la base de datos
    blast_db_name = os.path.join(full_path_blastdb_per_contig, contig_blast_id, contig_blast_id) ### path a la base de datos creada par ael contig en cuestion

    #print("blast_sseq: ", blast_sseq, '\n')

    cline = NcbiblastnCommandline(db=blast_db_name, evalue=0.001, perc_identity = 90, outfmt= blast_parameters, max_target_seqs=10, max_hsps=10,num_threads=1)    
    out, err = cline(stdin = blast_sseq) 
    out_lines = out.splitlines()

    #print("len(out_lines) get_prodigal_sequence: ", len(out_lines), '\n')

    #print("out_lines: ", out_lines, '\n')
    bigger_bitscore = 0
    if len (out_lines) > 0 :                  
        for line in out_lines :
            values = line.split('\t')
            if  float(values[8]) > bigger_bitscore: 
                qseqid , sseqid , pident ,  qlen , s_length , mismatch , gapopen , evalue , bitscore , sstart , send , qstart , qend ,sseq , qseq = values
                bigger_bitscore = float(bitscore)

        ### parseo el archivo de prodigal con todas las secuencias de todos los genes o el archivo del contig en cuestion con los genes predichos por prodigal en ese contig
        ### y obtengo la secuencia del gen que ha hecho match tras BLAST mediante el id del gen, sseqid, ademas de la posicion de inicio y final en el contig segun prodigal

        #contig_genes_path = os.path.join(full_path_prodigal_genes_per_contig + contig_blast_id + '.fasta')

        predicted_genes_in_contig = SeqIO.parse(contig_genes_path, "fasta")

        for rec in predicted_genes_in_contig: ### por cada uno de los genes predichos por prodigal en ese contig
            if rec.id == sseqid:
                predicted_gene_sequence = str(rec.seq)
                start_prodigal = str(rec.description.split( '#')[1]) ### para sacar report prodigal
                end_prodigal = str(rec.description.split('#')[2]) ### para sacar report prodigal
                break
   # print("predicted_gene_sequence: ", predicted_gene_sequence, '\n')
   # print("start_prodigal: ", start_prodigal, '\n')
   # print("end_prodigal: ", end_prodigal, '\n')

    if len (out_lines) == 0: ## if añadido para evitar error de código cuando prodigal no ha predicho la secuencia obtenida con blast y el resultado de este BLAST de las secs de prodigal frente a la secuencia de blast es 0
        predicted_gene_sequence = 'Sequence not found by Prodigal' 
        start_prodigal = '-' 
        end_prodigal = '-'

    #print("predicted_gene_sequence: ", predicted_gene_sequence)
    #print("predicted_gene_sequence: ", min_dist_seq, '\n')
    #print("start_prodigal: ", start_prodigal, '\n')
    #print("end_prodigal: ", end_prodigal, '\n')

    return predicted_gene_sequence, start_prodigal, end_prodigal ### añadido start_prodigal, end_prodigal para poder sacar report prodigal, quitar?



### get prodigal MASH
"""
def get_prodigal_sequence(blast_sseq, contig_blast_id, prodigal_directory, sample_name, core_name, logger): ### Cambiando/modificando: función para obtener la secuencia completa del gen de interés a partir de la predicción de prodigal
    # blast_sseq = secuencia encontrada tras BLAST sin gaps
    # contig_blast_id = id del contig de la muestra donde se ha encontrado la secuencia en BLASt --> Para coger los genes de prodigal que se han predicho en este contig para pillar esas 
        # secuencias para comparar con la secuencia encontrada con BLAST para intentar pillar la secuencia predicha con prodigal más similar a la encontrada con BLAST.
    # prodigal_directory
    # sample_name (nombre de la muestra)
    # core_name (nombre del locus, necesario para generar el archivo con la secuencia de blast, ya que va a llevar el nombre del locus con el que se corresponde)
    # no haría falta sstart y send pero dejar para generar informe con todas las secuencias que se sacan y sus sstart y end en prodigal junto el sstart y end en blast?
    
    
    ### argumentos anteriores
    # sstart = posición inicio sseq BLAST
    # send = posición final sseq BLAST
    # prodigal_directory* = directorio donde se encuentran los archivos de prodigal si se va a parsear directamente desde aquí y no se va a crear un diccionario de cada muestra previamente
    # sample_name* = nombre de la sample en cuestión para generar el path completo para los datos de esta muestra. Solo si se va a parsear directamente desde aquí y no se va a crear un diccionario de cada muestra previamente

    # prodigal_directory = os.path.join(tmp_samples_dir,'prodigal') ### para saber qué es prodigal_directory
    # sample_name = os.path.basename(sample_file) ### para saber qué es sample_name. Esto se le pasaría así, el nombre de la muestra en cuestión, si se parsea directamente desde allele_call, si no habría que pasar el sample_file que es el path completo de cada muestra a donde se cree diccionario o lo que sea de la muestra, y aquí habría que pasar el diccionario en cuestión, habría que pasar el path del diccionario, o el diccionario en sí se se parsea el diccionario al inicio de cada muestra, así no hay que parsear el diccionario para cada locus


    #print("core name: ", core_name, '\n', '\n')

    prodigal_directory_sample = os.path.join(prodigal_directory, sample_name)
    genes_file = os.path.join(prodigal_directory_sample, sample_name + '_dna.faa')


    ### Creacion del directorio mash en directorio prodigal de la muestra en cuestion
    mash_directory = 'mash'
    full_path_mash = os.path.join(prodigal_directory_sample, mash_directory)
    if not os.path.exists(full_path_mash):
        try:
            os.makedirs(full_path_mash)
            logger.info('Directory %s has been created', full_path_mash)
        except:
            print ('Cannot create the directory ', full_path_mash)
            logger.info('Directory %s cannot be created', full_path_mash)
            exit (0)


    ### Creacion del directorio prodigal_genes_per_contig dentro del directorio mash en directorio prodigal de la muestra en cuestion
    prodigal_genes_per_contig_directory = 'prodigal_genes_per_contig'
    full_path_prodigal_genes_per_contig = os.path.join(full_path_mash, prodigal_genes_per_contig_directory)
    if not os.path.exists(full_path_prodigal_genes_per_contig):
        try:
            os.makedirs(full_path_prodigal_genes_per_contig)
            logger.info('Directory %s has been created', full_path_prodigal_genes_per_contig)
        except:
            print ('Cannot create the directory ', full_path_prodigal_genes_per_contig)
            logger.info('Directory %s cannot be created', full_path_prodigal_genes_per_contig)
            exit (0)


    ### Creacion del directorio blast_seq_per_locus dentro del directorio mash en directorio prodigal de la muestra en cuestion
    blast_seq_per_locus_directory = 'blast_seq_per_locus'
    full_path_blast_seq_per_locus = os.path.join(full_path_mash, blast_seq_per_locus_directory)
    if not os.path.exists(full_path_blast_seq_per_locus):
        try:
            os.makedirs(full_path_blast_seq_per_locus)
            logger.info('Directory %s has been created', full_path_blast_seq_per_locus)
        except:
            print ('Cannot create the directory ', full_path_blast_seq_per_locus)
            logger.info('Directory %s cannot be created', full_path_blast_seq_per_locus)
            exit (0)


    ### Generando archivo fasta que contiene la secuencia encontrada con BLAST en directorio blast_seq_per_locus
    blast_seq_path = os.path.join(full_path_blast_seq_per_locus, core_name + '_blast.fasta')
    with open (blast_seq_path, 'w') as out_fh:
        out_fh.write ('>' + core_name + '_blast.fasta' + '\n' + blast_sseq)


    ### Creacion del directorio para este contig para guardar los fasta de los genes prodigal de este contig de forma indiviual + sketch de mash
    contig_directory = os.path.join(full_path_prodigal_genes_per_contig, contig_blast_id)
    
    if not os.path.exists(contig_directory):
        try:
            os.makedirs(contig_directory)
            logger.info('Directory %s has been created', contig_directory)
        except:
            print ('Cannot create the directory ', contig_directory)
            logger.info('Directory %s cannot be created', contig_directory)
            exit (0)


    ### Parseando el archivo de prodigal que contiene todos los genes predichos para la muestra en cuestión
    predicted_genes = SeqIO.parse(genes_file, "fasta")

    #print("predicted_genes: ", predicted_genes, '\n')

    #for rec in SeqIO.parse(genes_file, 'fasta'):
    for rec in predicted_genes: ### por cada uno de los genes predichos pro prodigal
        #print("rec in predicted_genes: ", rec, '\n')
        contig_prodigal_id = (rec.id).split("_")[0]

        #print("contig_prodigal_id: ", contig_prodigal_id, '\n')
        
        #print("type(contig_blast_id): ", type(contig_blast_id), '\n')

        #print("type(contig_prodigal_id): ", type(contig_prodigal_id), '\n')

        if contig_prodigal_id == contig_blast_id: ### mira si el id del gen predicho por prodigal se ha encontrado en el mismo contig que la secuencia encontrada por blast

            #print("ha entrado a int(contig_prodigal_id) == int(contig_blast_id)")

            ### AQUÍ CADA UNO DE ESOS GENES LOS VOY GUARDANDO EN ARCHIVOS FASTA INDIVIDUALES PARA HACER SKETCH CON MASH
            ### Habría que generar directorio mash dentro de directorio prodigal de la muestra

            ### Esto se va a generar para cada muestra y cada locus, crear una carpeta dentro de prodigal para cada muestra que sea para cada contig de la muestra, rollo ---> 
            ### nombre: nombremuestra_contig, de modo que si previamente ya se han generado las secuencias fastas individuales para ese contig se comprueba si existe ya la carpeta
            ### y de existir se utiliza esta y no se vuelve a crear y a volver a splitear las secuencias de ese contig y a generar el sketch de mash con estas secuencias, usando este
            ### sketch. path prodigal/muestra/mash/prodigal_genes_per_contig/contigX

            ### Para evitar repetir este proceso para contigs repetidos se borrarían todos los archivos generados para mash al final --> Pero tendría que ser al final de todo el allele_calling
            ### ya que el allele_calling se hace por locus y por muestra, no por muestra y por locus, así que tendría que pasar todos los locus por todas las muestras antes de borrar...
                ### Se estaría almacenando mucha mierda, no sé si al final ocuparía mucho espacio, pero si no lo hago así se borraría y se generarían fastas y sketches del mismo contig 
                ### muchísimas veces ------------------> Sería mejor hacer este proceso de dividir los genes de prodigal en fastas individuales al inicio al preaprar las muestras?

            ### Por otro lado hay que generar el fasta de sseq para comparar con el sketch de prodigal. Esta secuencia se guardaría en la carpeta de mash para cada muestra con el nombre
            ### del locus al cual se corresponde esa secuencia de esa muestra ---> path prodigal/muestra/mash/blast_seq_per_locus/nombrelocus

            ### Generando fasta con la secuencia de este gen predicho por prodigal que se ha encontrado en el mismo contig en el que se encontro la secuencia tras BLAST
            contig_genes_path = os.path.join(contig_directory, str(rec.id) + '.fasta')
            with open (contig_genes_path, 'w') as out_fh:
                #out_fh.write ('>' + str(rec.id) + '\n' + str(rec.seq))
                out_fh.write ('>' + str(rec.description) + '\n' + str(rec.seq)) ### poniendo descripcion como cabecera del fasta para poder sacar posteriormente start y end para report prodigal

    #print("contig_blast_id: ", contig_blast_id, '\n', '\n')

    ### comando mash para generar sketch con todos los fasta files para obtener distancias entre secuencias simultáneamente
    sketch_path = os.path.join(contig_directory, contig_blast_id + '_sketch.msh')
    #print("sketch_path: ", sketch_path, '\n')

    #print("contig_directory de donde se listan los fastas de los genes prodigal del contig en cuestion: ", contig_directory, '\n', '\n')

    ### obtengo la lista de paths de fastas que hay en el directorio contig_directory
    prodigal_gene_files_list = get_fasta_file_list(contig_directory, logger)

    mash_sketch_command = ["mash", "sketch", "-o", sketch_path]

    #print("prodigal_gene_files_list: ", prodigal_gene_files_list, '\n')

    #print("prodigal_gene_files_list: ", prodigal_gene_files_list, '\n', '\n')
    for fasta_path in prodigal_gene_files_list : ### Añadiendo los paths de todos los fastas que se van a incluir en el sketch al comando mash_sketch_command
        #print("fasta_path: ", fasta_path, '\n')
        mash_sketch_command.append(fasta_path)
    #print("mash_sketch_command: ", mash_sketch_command, '\n')

    mash_sketch_result = subprocess.run(mash_sketch_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


    ### comando mash para obtener mash distances entre sencuencias del sketch (secuencias predichas por prodigal que se encuentran en el contig de interes) y la secuencia encontrada con BLAST
    mash_distance_command = ["mash", "dist", sketch_path, blast_seq_path]
    
    mash_distance_result = subprocess.Popen(mash_distance_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    out, err = mash_distance_result.communicate()
    #print("out: ", out)
    out = out.decode('UTF-8').split('\n')
    #print("out_decode: ", out, '\n')

    comp_dist_list = []
    id_genes_list = []
    for n in range(len(out)-1): ### Comprobar que esto pasa también en este caso: restando 1 porque genera una línea final vacía, no se si es algo general que pasara siempre o solo se está dando en el locus que estoy usando como referencia para construir la función
        comp = out[n].split('\t')
        #print("comp: ", comp, '\n')
        comp_dist = float(comp[2]) ### cuidado porque por la cara mete un ultimo elemento en la lista que es un "blanco" y no tiene elemento 2, porque es una lista vaca -.-
        comp_dist_list.append(comp_dist)

        ### probando a hacer lista de ids tal y como compara mash para poder sacar el id por el índice de la min mean
        gene_id = os.path.basename(comp[0]).split('.')[0]
        id_genes_list.append(gene_id)

    #print("comp_dist_list: ", comp_dist_list, '\n')

    min_dist = min(comp_dist_list) ### distancia minima
    min_dist_index = comp_dist_list.index(min_dist) ### sacando indice de la distancia minima en la lista comp_dist_list para obtener el id del gen predicho por prodigal que se aleja menos de la secuencia blast
    min_dist_id = id_genes_list[min_dist_index] ### sacando id del gen para el que se ha obtenido la mínima distancia

    ### se lee el fasta que cuyo nombre se corresponde con el id de la secuencia para la que se he obtenido una menor distancia con respecto a la secuencia blast (min_dist_prodigal_id)
    min_dist_prodigal_path = os.path.join(contig_directory, min_dist_id + '.fasta') ### path del gen para el que se ha obtenido menor distancia

    min_dist_gene = SeqIO.parse(min_dist_prodigal_path, "fasta") ### no sé si puedo sacar directamente la secuencia haciendo min_dist_gene_seq, o si tengo que hacer el for record in min_dist_gene, etc
                                                                ### no se puede. Lo saco haciendo for, pero mirar si hay otra forma de sacar la seq sin tener que hacer for ya que solo hay una secuencia en el fasta
    
    #print("min_dist_id : ", min_dist_id, '\n')


    for gene in min_dist_gene:
        min_dist_seq = str(gene.seq)
        
        start_prodigal = str(gene.description.split('#')[1]) ### para sacar report prodigal
        end_prodigal = str(gene.description.split('#')[2]) ### para sacar report prodigal


    #print("min_dist_seq: ", min_dist_seq, '\n')
    #print("start_prodigal: ", start_prodigal, '\n')
    #print("end_prodigal: ", end_prodigal, '\n')

    return min_dist_seq, start_prodigal, end_prodigal ### añadido start_prodigal, end_prodigal para poder sacar report prodigal, quitar?
"""

# introduciendo directorios de las secuencias d eblast encontradas con alelo de referencia
def prepare_samples(sample_file_list, store_dir, reference_genome_file, logger): ### modificando/cambiando: introducido reference_genome_file, genoma de referencia 
    file_list = []
    contigs_in_sample_dict = {} ### diccionario para guardar para cada muestra los id-contig

    blast_dir = os.path.join(store_dir,'blastdb')
    prodigal_dir = os.path.join(store_dir,'prodigal') ### cambiando/modificando: se introduce el prodigal_dir
    blast_results_seq_dir = os.path.join(store_dir,'blast_results', 'blast_results_seq') #---> a este path hay que ir añadiendo el nombre de la muestra

    # create training file for genes prediction
    output_prodigal_train_dir = prodigal_training(reference_genome_file, prodigal_dir, logger) ### cambiando/modificando: se llama a la función prodigal_training
    if not output_prodigal_train_dir:
        print('Error when creating training file for genes prediction. Check log file for more information. \n ')
        return False
    """
    if not prodigal_training(reference_genome_file, prodigal_dir, logger): ### cambiando/modificando: se llama a la función prodigal_training
        print('Error when creating training file for genes prediction. Check log file for more information. \n ')
        return False
    """

    for fasta_file in sample_file_list:
        f_name = os.path.basename(fasta_file).split('.')
        file_list.append(os.path.join(store_dir, f_name[0]))

        # parsing fasta file and get in the dictionary the id and the sequence 
        fasta_file_parsed_dict = parsing_fasta_file_to_dict(fasta_file, logger) #---------> De momento este diccionario de contigs no se utiliza, comentado
        if f_name[0] not in contigs_in_sample_dict.keys():
            contigs_in_sample_dict[f_name[0]] = {}
        contigs_in_sample_dict[f_name[0]] = fasta_file_parsed_dict ### guardando el diccionario id-contigs de esta muestra en el diccionario general con todos los locus

        # dump fasta file into pickle file
        #with open (file_list[-1],'wb') as f: -------------> Comentado generacion de diccionarios de contigs para cada muestra
         #   pickle.dump(fasta_file_parsed_dict, f)
    
        # creacion directorios por muestra para guardar secuencias encontradas con BLAST mediante alelo de referencia de cada locus
        blast_results_seq_per_sample_dir = os.path.join(blast_results_seq_dir, f_name[0]) #---> a este path hay que ir añadiendo el nombre de la muestra
        
        if not os.path.exists(blast_results_seq_per_sample_dir):
            try:
                os.makedirs(blast_results_seq_per_sample_dir)
                logger.debug('Created blast results directory for sample %s', f_name[0])
            except:
                logger.info('Cannot create blast results directory for sample %s', f_name[0])
                print ('Error when creating the directory for blast results', blast_results_seq_per_sample_dir)
                exit(0)

        # predict genes for each sample fasta file
        if not prodigal_prediction(fasta_file, prodigal_dir, output_prodigal_train_dir, logger): ### cambiando/modificando: se llama a la función prodigal_prediction
            print('Error when predicting genes for samples files. Check log file for more information. \n ')
            return False

        # create local blast db for sample fasta file
        if not create_blastdb(fasta_file, blast_dir, 'nucl', logger):
            print('Error when creating the blastdb for samples files. Check log file for more information. \n ')
            return False

    return file_list, contigs_in_sample_dict

# Prepare_samples anterior antes de intentar meter los directorios de las secuencias de blast encontradas con alelo de referencia
"""
def prepare_samples(sample_file_list, store_dir, reference_genome_file, logger): ### modificando/cambiando: introducido reference_genome_file, genoma de referencia 
    file_list = []
    blast_dir = os.path.join(store_dir,'blastdb')
    prodigal_dir = os.path.join(store_dir,'prodigal') ### cambiando/modificando: se introduce el prodigal_dir

    # create training file for genes prediction
    output_prodigal_train_dir = prodigal_training(reference_genome_file, prodigal_dir, logger) ### cambiando/modificando: se llama a la función prodigal_training
    if not output_prodigal_train_dir:
        print('Error when creating training file for genes prediction. Check log file for more information. \n ')
        return False
    
    #if not prodigal_training(reference_genome_file, prodigal_dir, logger): ### cambiando/modificando: se llama a la función prodigal_training
     #   print('Error when creating training file for genes prediction. Check log file for more information. \n ')
     #   return False
    
    for fasta_file in sample_file_list:
        # parsing fasta file and get in the dictionary the id and the sequence 
        # fasta_file_parsed_dict = parsing_fasta_file_to_dict(fasta_file, logger) ---------> De momento este diccionario de contigs no se utiliza, comentado
        f_name = os.path.basename(fasta_file).split('.')
        file_list.append(os.path.join(store_dir, f_name[0]))
        # dump fasta file into pickle file
        #with open (file_list[-1],'wb') as f: -------------> Comentado generacion de diccionarios de contigs para cada muestra
         #   pickle.dump(fasta_file_parsed_dict, f)
    
        # predict genes for each sample fasta file
        if not prodigal_prediction(fasta_file, prodigal_dir, output_prodigal_train_dir, logger): ### cambiando/modificando: se llama a la función prodigal_prediction
            print('Error when predicting genes for samples files. Check log file for more information. \n ')
            return False

        # create local blast db for sample fasta file
        if not create_blastdb(fasta_file, blast_dir, 'nucl', logger):
            print('Error when creating the blastdb for samples files. Check log file for more information. \n ')
            return False

    return file_list
"""

def length_thresholds(core_name, schema_statistics, percent): ### cambiando/modificando: función añadida para obtener el rango de longitud para clasificación en INF, ASM o ALM
                                                                    ### logger?
    locus_mean = int(schema_statistics[core_name][1])

    if percent != "SD": 
        max_length_threshold = math.ceil(locus_mean + ((locus_mean * float(percent)) / 100))
        min_length_threshold = math.floor(locus_mean - ((locus_mean * float(percent)) / 100))
    else:
        percent = float(schema_statistics[core_name][2])

        max_length_threshold = math.ceil(locus_mean + (locus_mean * percent))
        min_length_threshold = math.floor(locus_mean - (locus_mean * percent))

    return max_length_threshold, min_length_threshold

def convert_to_protein (sequence) :

    seq = Seq.Seq(sequence)
    protein = str(seq.translate())

    return protein

def get_snp (sample, query) :

    ### Tiene sentido mirar SNPs en INF, ASM y ALM si realmente al tener inserciones y deleciones puedes llegar a una base de lasecuencia insertada por lo que la estarías comparando
    ### con otra base en la otra secuencia que se correspondería con al siguiente baase realmente y esa coincidiría con un gap? 
    ### Quitar de Exact Match ya que si realmente hay Exact Match no debería haber SNPs.


    prot_annotation = {'S': 'polar' ,'T': 'polar' ,'Y': 'polar' ,'Q': 'polar' ,'N': 'polar' ,'C': 'polar' ,'S': 'polar' ,
                        'F': 'nonpolar' ,'L': 'nonpolar','I': 'nonpolar','M': 'nonpolar','P': 'nonpolar','V': 'nonpolar','A': 'nonpolar','W': 'nonpolar','G': 'nonpolar',
                        'D' : 'acidic', 'E' :'acidic',
                        'H': 'basic' , 'K': 'basic' , 'R' : 'basic',
                        '-': '-----', '*' : 'Stop codon'}
    snp_list = []
    sample = sample.replace('-','')
    #length = max(len(sample), len(query))
    length = len(query)
    # normalize the length of the sample for the iteration
    if len(sample) < length :
        need_to_add = length - len(sample)
        sample = sample + need_to_add * '-'

    # convert to Seq class to translate to protein
    seq_sample = Seq.Seq(sample)
    seq_query = Seq.Seq(query)
    
    for index in range(length):
        if seq_query[index] != seq_sample[index] :
            triple_index = index - (index % 3)
            codon_seq = seq_sample[triple_index : triple_index + 3]
            codon_que = seq_query[triple_index : triple_index + 3]
            if not '-' in str(codon_seq) :
                prot_seq = str(codon_seq.translate())
                prot_que = str(codon_que.translate())
            else:
                prot_seq = '-'
                prot_que = str(seq_query[triple_index: ].translate())
            if prot_annotation[prot_que[0]] == prot_annotation[prot_seq[0]] :
                missense_synonym = 'Synonymous'
            elif prot_seq == '*' :
                    missense_synonym = 'Nonsense'
            else:
                missense_synonym = 'Missense'
            #snp_list.append([str(index+1),str(seq_sample[index]) + '/' + str(seq_query[index]), str(codon_seq) + '/'+ str(codon_que),
            snp_list.append([str(index+1),str(seq_query[index]) + '/' + str(seq_sample[index]), str(codon_que) + '/'+ str(codon_seq),
                             # when one of the sequence ends but not the other we will translate the remain sequence to proteins
                             # in that case we will only annotate the first protein. Using [0] as key of the dictionary  annotation
                             missense_synonym, prot_que + '/' + prot_seq, prot_annotation[prot_que[0]] + ' / ' + prot_annotation[prot_seq[0]]])
            if '-' in str(codon_seq) :
                break
        
    return snp_list

### cambiando/modificando: 
### creo que tiene sentido comparar las secuencias sin gaps porque así se ve cómo se ha alterado la secuencia? en una misma posición de una secuencia y otra qué es lo que había y qué es lo que hay ahora?
### revisar, puede que no tenga que cambiar esta función, pero sí quitarla de exact_match, aplicarla a INF, ASM, ALM
### al menos sí usar convert to protein en lugar de tener que utilizar Seq y translate
"""
def get_snp (sample, query, sample_prot_align, query_prot_align) :

    #esto sería si se busca a partir de proteína
    #sample = new_sseq
    #querry = matching_allele_seq
    #protein_alignment = get_alignment (sample_seq, query_seq, seq_type = "protein")
    #sample_prot_align = protein_alignment[0][1]
    #query_prot_align = protein_alignment[2][1]

    #esto sería si se busca a partir de dna
    #dna_alignment = get_alignment (sample_seq, query_seq, seq_type = "dna")
    #sample_dna_align = dna_alignment[0][1]
    #query_dna_align = dna_alignment[2][1]





    prot_annotation = {'S': 'polar' ,'T': 'polar' ,'Y': 'polar' ,'Q': 'polar' ,'N': 'polar' ,'C': 'polar' ,'S': 'polar' ,
                        'F': 'nonpolar' ,'L': 'nonpolar','I': 'nonpolar','M': 'nonpolar','P': 'nonpolar','V': 'nonpolar','A': 'nonpolar','W': 'nonpolar','G': 'nonpolar',
                        'D' : 'acidic', 'E' :'acidic',
                        'H': 'basic' , 'K': 'basic' , 'R' : 'basic',
                        '-': '-----', '*' : 'Stop codon'}

    snp_list = []
    #sample = sample.replace('-','')
    #length = max(len(sample), len(query))
    length = len(query_dna_align) ### sample_dna_align y query_dna_align deben tener la misma longitud al contener los gaps tras el alineamiento con pairwise2.align.localms en la función get_alignment
    

    # normalize the length of the sample for the iteration
    #if len(sample) < length :
     #   need_to_add = length - len(sample)
      #  sample = sample + need_to_add * '-'

    # convert to Seq class to translate to protein
    ### utilizar aquí la función convert_to_protein
    #seq_sample = Seq.Seq(sample)
    #seq_query = Seq.Seq(query)

    ### obteniendo prot a partir de sample_dna_align eliminando gaps y empleando la función convert_to_protein
    sample_dna = sample_dna_align.replace('-', '')
    query_dna = query_dna_align.replace('-', '')

    sample_prot = convert_to_protein(sample_dna)
    query_prot = convert_to_protein(query_dna)


    for index_prot in range(length):
        index_sample_dna = index_prot
        index_query_dna = index_prot
        if sample_dna[index_sample_dna] != sample_dna[index_query_dna]:
            aa_sample = sample_protein_align[index_prot]
            aa_query = query_protein_align[index_prot]


        if sample_protein_align[index] != sample_dna_align[index] :
            triple_index = index - (index % 3)
            codon_seq = sample_dna_align[triple_index : triple_index + 3]
            codon_que = query_dna_align[triple_index : triple_index + 3]
            if not '-' in str(codon_seq) :
                prot_seq = str(codon_seq.translate())
                prot_que = str(codon_que.translate())
            else:
                prot_seq = '-'
                prot_que = str(seq_query[triple_index: ].translate())
            if prot_annotation[prot_que[0]] == prot_annotation[prot_seq[0]] :
                missense_synonym = 'Synonymous'
            elif prot_seq == '*' :
                    missense_synonym = 'Nonsense'
            else:
                missense_synonym = 'Missense'
            #snp_list.append([str(index+1),str(seq_sample[index]) + '/' + str(seq_query[index]), str(codon_seq) + '/'+ str(codon_que),
            snp_list.append([str(index+1),str(seq_query[index]) + '/' + str(seq_sample[index]), str(codon_que) + '/'+ str(codon_seq),
                             # when one of the sequence ends but not the other we will translate the remain sequence to proteins
                             # in that case we will only annotate the first protein. Using [0] as key of the dictionary  annotation
                             missense_synonym, prot_que + '/' + prot_seq, prot_annotation[prot_que[0]] + ' / ' + prot_annotation[prot_seq[0]]])
            if '-' in str(codon_seq) :
                break
        
    return snp_list
"""

def nucleotide_to_protein_alignment (sample_seq, query_seq ) : ### Sustituir por get_alignment
    aligment = []
    sample_prot = convert_to_protein(sample_seq)
    query_prot = convert_to_protein(query_seq)
    minimun_length = min(len(sample_prot), len(query_prot))
    for i in range(minimun_length):
        if sample_prot[i] == query_prot[i] :
            aligment.append('|')
        else:
            aligment.append(' ')
    protein_alignment = [['sample', sample_prot],['match', ''.join(aligment)], ['schema', query_prot]]

    return protein_alignment

def get_alignment_for_indels (blast_db_name, qseq) : ### Sustituir por get_alignment
    #match_alignment =[]
    cline = NcbiblastnCommandline(db=blast_db_name, evalue=0.001, perc_identity = 80, outfmt= 5, max_target_seqs=10, max_hsps=10,num_threads=1)
    out, err = cline(stdin = qseq)
    psiblast_xml = StringIO(out)
    blast_records = NCBIXML.parse(psiblast_xml)

    #print("blast_records: ", blast_records, '\n') ### printeando, borar
    for blast_record in blast_records:
     #   print("blast_record in blast_records: ", blast_record, '\n') ### printeando, borar
        for alignment in blast_record.alignments:
      #      print("alignment in blast_record.alignments: ", alignment, '\n') ### printeando, borar
            for match in alignment.hsps:
       #         print(" match in alignment.hsps: ", match, '\n') ### printeando, borar
                match_alignment = [['sample', match.sbjct],['match', match.match], ['schema',match.query]]
        #        print("match_alignment: ", match_alignment, '\n', '\n') 
    return match_alignment


def get_alignment_for_deletions (sample_seq, query_seq): ### Sustituir por get_alignment
    index_found = False
    alignments = pairwise2.align.globalxx(sample_seq, query_seq)
    #print("alignments: ", alignments, '\n') ### printeando, borrar
    for index in range(len(alignments)) :
        if alignments[index][4] == len(query_seq) :  ### Esto es para coger el alineamiento (index) que tiene la misma longitud que la secuencia del alelo query_seq
            #print("alignments[index][4]: ", alignments[index][4], '\n') ### printeando, borrar
            index_found = True
            #print("index_found: ", index_found, '\n') ### printeando, borrar
            break
    if not index_found : ### Y si no se da el if de arriba entonces se coge el primer alineamiento que se obtiene
        index = 0
    values = format_alignment(*alignments[index]).split('\n')
    #print("values: ", values, '\n') ### printeando, borrar
    match_alignment = [['sample', values[0]],['match', values[1]], ['schema',values[2]]]
    #print("match_alignment: ", match_alignment, '\n') ### printeando, borrar

    return match_alignment

### Cambiando/modificando: obteniendo get_alignment unificando get_alignment_for_indels, get_alignment_for_deletions y nucleotide_to_protein_alignment. Aplicar sobre INF, deletions, insertions, equal
def get_alignment (sample_seq, query_seq, seq_type = "dna"):
    
    ### si se condiciona el uso de id 85 e id 90 habría que condicionar los params de pairwise2.align
    #param blast_id introducido por checkeo de coverage y para saber qué params utilizar en pairwise2 en función del BLAST empleado (BLAST 90 o BLAST 85 con dif params gaps)

    if seq_type == "protein":
        sample_seq = convert_to_protein(sample_seq)
        query_seq = convert_to_protein(query_seq)

    # arguments pairwise2.align.globalms: match, mismatch, gap opening, gap extending
    alignments = pairwise2.align.localms(sample_seq, query_seq, 1, -2, -1, -1) ### cambiando/modificando: mismos match, mismatch, gap opening, gap extending especificados en BLAST ID 85

    #print("alignments: ", alignments, '\n') ### printeando, borrar

    values = format_alignment(*alignments[0]).split('\n')
    #print("values: ", values, '\n') ### printeando, borrar
    
    match_alignment = [['sample', values[0]],['match', values[1]], ['schema',values[2]]]
    #print("match_alignment: ", match_alignment, '\n') ### printeando, borrar
    #print("len(sample_seq): ", len(sample_seq), '\n') ### printeando, borrar
    #print("len(query_seq): ", len(query_seq), '\n') ### printeando, borrar

    return match_alignment


def lnf_tag(core_name, sample_name, alleles_in_locus_dict, samples_matrix_dict, lnf_dict, locus_alleles_path, qseqid, pident, s_length, new_sequence_length, annotation_core_dict, logger):

    # meter values igual que en inf_asm_alm_tag en lugar de qseqid, pident, s_length?
    # si meto locus_alleles_path es necesario meter core_name o hago basename de locus_alleles_path?

    # obteniendo anotación gen y producto del locus
    gene_annot, product_annot = annotation_core_dict[core_name]

    # obteniendo tag y longitud alelo match
    if qseqid == '-':
        samples_matrix_dict[sample_name].append('LNF')
        matching_allele_length = '-'
    
    else:
        samples_matrix_dict[sample_name].append('LNF_' + str(qseqid)) ### indicar etiqueta así o solo LNF?

        alleles_in_locus = list (SeqIO.parse(locus_alleles_path, "fasta")) ## comentar parse
        for allele in alleles_in_locus : ## comentar parse
            if allele.id == qseqid : ## comentar parse
                break ## comentar parse
        matching_allele_seq = str(allele.seq) ## comentar parse
        matching_allele_length = len(matching_allele_seq) ## comentar parse

        matching_allele_seq = alleles_in_locus_dict[core_name][qseqid] ## Accediendo a diccionario locus. Alternativa a parser
        matching_allele_length = len(matching_allele_seq) ## Accediendo a diccionario locus. Alternativa a parser

    if pident == '-':
        # los dos BLAST sin resultado:
        coverage_blast = '-'
        coverage_new_sequence = '-'
        add_info = 'Locus not found'
        logger.info('Locus not found at sample %s, for gene %s', sample_name, core_name)

    elif 90 > float(pident):  ### mirar cómo lo indica, por si indica 90.0 o algo así
        # BLAST 90 sin resultado y BLAST 70 con resultado:
        coverage_blast = '-'
        coverage_new_sequence = '-'
        add_info = 'BLAST sequence ID under threshold: {}%'.format(perc_identity_ref) ### En este mensaje poner el porcentaje de ID oportuno, de momento es 90%
        logger.info('BLAST sequence ID %s under threshold at sample %s, for gene  %s', pident, sample_name, core_name)

    elif 90 <= float(pident) and new_sequence_length == '-':
        # BLAST 90 con resultado, bajo coverage BLAST:
        coverage_blast = int(s_length) / matching_allele_length ### aquí debería poner matching_allele_length o la media en función de cómo se está calculando los rangos?
        coverage_new_sequence = '-'
        if coverage_blast < 1:
            add_info = 'BLAST sequence coverage under threshold: {}%'.format(coverage)
        else:  
            add_info = 'BLAST sequence coverage above threshold: {}%'.format(coverage)
        logger.info('BLAST sequence coverage %s under threshold at sample %s, for gene  %s', coverage_blast, sample_name, core_name)

    elif 90 <= float(pident) and new_sequence_length != '-':
        # BLAST 90 con resultado, buen coverage BLAST, bajo coverage new_sseq:
        coverage_blast = int(s_length) / matching_allele_length ### aquí debería poner matching_allele_length o la media en función de cómo se está calculando los rangos?
        coverage_new_sequence = new_sequence_length / matching_allele_length ### aquí debería poner matching_allele_length o la media en función de cómo se está calculando los rangos?
        if coverage_new_sequence < 1:
            add_info = 'New sequence coverage under threshold: {}%'.format(coverage)
        else:  
            add_info = 'New sequence coverage above threshold: {}%'.format(coverage)
        logger.info('New sequence coverage %s under threshold at sample %s, for gene  %s', coverage_new_sequence, sample_name, core_name)

    # creando report lnf
    if not core_name in lnf_dict:
        lnf_dict[core_name] = {}
    if not sample_name in lnf_dict[core_name]:
        lnf_dict[core_name][sample_name] = []

    lnf_dict[core_name][sample_name].append([gene_annot, product_annot, qseqid, pident, str(coverage_blast), str(coverage_new_sequence), str(matching_allele_length), str(s_length), str(new_sequence_length), add_info]) ### Meter secuencias alelo, blast y new_sseq (si las hay)?

    return True

# Funcion paralog_exact_tag para etiquetar paralogos y exact matches

def paralog_exact_tag(sample_name, core_name, tag, schema_quality, matching_genes_dict, samples_matrix_dict, allele_found, tag_dict, prodigal_report, prodigal_directory, blast_parameters, annotation_core_dict, logger):

    # tag --> NIPHEM / NIPH X
    # sample_name X
    # samples_matrix_dict X
    # tag_dict X paralog_dict // exact_dict
    # core_name X
    # allele_found X
    # paralog_found (NIPH) X PROBANDO A ELIMINAR PARALOG FOUND
    # prodigal_directory X
    # blast_parameters X    
    # matching_genes_dict X

    logger.info('Found %s at sample %s for core gene %s ', tag, sample_name, core_name)

    # obteniendo anotación gen y producto del locus
    gene_annot, product_annot = annotation_core_dict[core_name]

    if not sample_name in tag_dict :
        tag_dict[sample_name] = {}
    if not core_name in tag_dict[sample_name] :
        tag_dict[sample_name][core_name]= []

   # all_sequences_matrix = {}
   # if tag == 'NIPH':
        # merging the 2 dictionaries
    ## paralog_matrix[sample_name] = {**allele_found, **paralog_found} ### no le veo sentido a crear el diccionario así poniendo el sample_name como clave si cada vez que haya un NIPH en el código original se genera el diccionario paralog_matrix = [] vacío, paralog matrix no se utiliza para guardar archivos al final
    #    all_sequences_matrix = {**allele_found, **paralog_found} ### intentando ponerlo sin sample_name
     #   print("all_sequences_matrix: ", all_sequences_matrix, '\n')
      #  print("allele_found: ", allele_found, '\n')
       # print("paralog_found: ", paralog_found, '\n')
   # else:
    #    all_sequences_matrix = allele_found ### guardando allele_found en caso NIPHEM como paralog_matrix para poder hacer el siguiente paso conjunto para NIPH y NIPHEM

    if tag == 'EXACT':
        allele = list(allele_found.keys())[0]
        qseqid = allele_found[allele][0]
        tag = qseqid

    samples_matrix_dict[sample_name].append(tag) 


    for sequence in allele_found:
        qseqid, sseqid, pident, qlen, s_length, mismatch, gapopen, evalue, bitscore, sstart, send, qstart, qend, sseq, qseq = allele_found[sequence]
        sseq = sseq.replace('-', '') ### MIRAR CÓMO VOY A GUARDAR LOS DATOS EN ALLELE_FOUND Y PARALOG FOUND PARA SABER SI TENGO QUE SACAR LONGITUDES Y SECUENCIAS SIN GAPS O NO
                                    ## No hace falta quitar - en sseq de NIPHEM, pero lo hago en general

        # Obteniendo calidad del alelo
        allele_quality = schema_quality[core_name][qseqid]

        # Obteniendo secuencia predicha por prodigal si la calidad del alelo es bad_quality
        if 'bad_quality' in allele_quality: 
            #complete_predicted_seq = str(get_prodigal_sequence(sstart, send, sseqid, prodigal_directory, sample_name))
            #complete_predicted_seq, start_prodigal, end_prodigal = get_prodigal_sequence(sseq, sseqid, prodigal_directory, sample_name, core_name, logger)
            complete_predicted_seq, start_prodigal, end_prodigal = get_prodigal_sequence(sseq, sseqid, prodigal_directory, sample_name, blast_parameters, logger)

            ##### PARA SACAR INFORME DE PRODIGAL VS BLAST #####
            prodigal_report.append([core_name, sample_name, gene_annot, product_annot, qseqid, tag, sstart, send, start_prodigal, end_prodigal, sseq, complete_predicted_seq])

        else:
            complete_predicted_seq = '-'
        
        print("complete_predicted_seq Prodigal: ", complete_predicted_seq, '\n', '\n') ### printeando, borrar

        if not sseqid in matching_genes_dict[sample_name] :
            matching_genes_dict[sample_name][sseqid] = []
        if sstart > send :
            matching_genes_dict[sample_name][sseqid].append([core_name, sstart, send,'-', tag])
        else:
            matching_genes_dict[sample_name][sseqid].append([core_name, sstart, send,'+', tag])

        if tag == 'NIPH' or tag == 'NIPHEM':
            tag_dict[sample_name][core_name].append([tag, pident, qseqid, allele_quality, sseqid, bitscore, sstart, send, sseq, complete_predicted_seq]) ### cambiando/modificando: introducido allele_quality, complete_predicted_seq y etiqueta NIPHEM
        else:
            tag_dict[sample_name][core_name] = [qseqid, allele_quality, sseqid, s_length, sstart, send, sseq, complete_predicted_seq] ### cambiando/modificando: introducido allele_quality, complete_predicted_seq y etiqueta NIPHEM

    return True 


"""
# Funcion paralog_exact_tag para etiquetar paralogos y exact matches

def paralog_exact_tag(sample_name, core_name, tag, schema_quality, matching_genes_dict, samples_matrix_dict, allele_found, tag_dict, prodigal_report, prodigal_directory, blast_parameters, logger):

    # tag --> NIPHEM / NIPH X
    # sample_name X
    # samples_matrix_dict X
    # tag_dict X paralog_dict // exact_dict
    # core_name X
    # allele_found X
    # paralog_found (NIPH) X PROBANDO A ELIMINAR PARALOG FOUND
    # prodigal_directory X
    # blast_parameters X    
    # matching_genes_dict X

    logger.info('Found %s at sample %s for core gene %s ', tag, sample_name, core_name)

    if not sample_name in tag_dict :
        tag_dict[sample_name] = {}
    if not core_name in tag_dict[sample_name] :
        tag_dict[sample_name][core_name]= []

   # all_sequences_matrix = {}
   # if tag == 'NIPH':
        # merging the 2 dictionaries
    ## paralog_matrix[sample_name] = {**allele_found, **paralog_found} ### no le veo sentido a crear el diccionario así poniendo el sample_name como clave si cada vez que haya un NIPH en el código original se genera el diccionario paralog_matrix = [] vacío, paralog matrix no se utiliza para guardar archivos al final
    #    all_sequences_matrix = {**allele_found, **paralog_found} ### intentando ponerlo sin sample_name
     #   print("all_sequences_matrix: ", all_sequences_matrix, '\n')
      #  print("allele_found: ", allele_found, '\n')
       # print("paralog_found: ", paralog_found, '\n')
   # else:
    #    all_sequences_matrix = allele_found ### guardando allele_found en caso NIPHEM como paralog_matrix para poder hacer el siguiente paso conjunto para NIPH y NIPHEM

    if tag == 'EXACT':
        allele = list(allele_found.keys())[0]
        qseqid = allele_found[allele][0]
        tag = qseqid

    samples_matrix_dict[sample_name].append(tag) 


    for sequence in allele_found:
        qseqid, sseqid, pident, qlen, s_length, mismatch, gapopen, evalue, bitscore, sstart, send, qstart, qend, sseq, qseq = allele_found[sequence]
        sseq = sseq.replace('-', '') ### MIRAR CÓMO VOY A GUARDAR LOS DATOS EN ALLELE_FOUND Y PARALOG FOUND PARA SABER SI TENGO QUE SACAR LONGITUDES Y SECUENCIAS SIN GAPS O NO
                                    ## No hace falta quitar - en sseq de NIPHEM, pero lo hago en general

        # Obteniendo calidad del alelo
        allele_quality = schema_quality[core_name][qseqid]

        # Obteniendo secuencia predicha por prodigal si la calidad del alelo es bad_quality
        if 'bad_quality' in allele_quality: 
            #complete_predicted_seq = str(get_prodigal_sequence(sstart, send, sseqid, prodigal_directory, sample_name))
            #complete_predicted_seq, start_prodigal, end_prodigal = get_prodigal_sequence(sseq, sseqid, prodigal_directory, sample_name, core_name, logger)
            complete_predicted_seq, start_prodigal, end_prodigal = get_prodigal_sequence(sseq, sseqid, prodigal_directory, sample_name, blast_parameters, logger)

            ##### PARA SACAR INFORME DE PRODIGAL VS BLAST #####
            prodigal_report.append([core_name, sample_name, qseqid, tag, sstart, send, start_prodigal, end_prodigal, sseq, complete_predicted_seq])

        else:
            complete_predicted_seq = '-'
        
        print("complete_predicted_seq Prodigal: ", complete_predicted_seq, '\n', '\n') ### printeando, borrar

        if not sseqid in matching_genes_dict[sample_name] :
            matching_genes_dict[sample_name][sseqid] = []
        if sstart > send :
            matching_genes_dict[sample_name][sseqid].append([core_name, sstart, send,'-', tag])
        else:
            matching_genes_dict[sample_name][sseqid].append([core_name, sstart, send,'+', tag])

        #tag_dict[sample_name][core_name].append([tag, pident, qseqid, allele_quality, sseqid, bitscore, sstart, send, sseq, complete_predicted_seq]) ### cambiando/modificando: introducido allele_quality, complete_predicted_seq y etiqueta NIPHEM
        tag_dict[sample_name][core_name] = [tag, pident, qseqid, allele_quality, sseqid, bitscore, sstart, send, sseq, complete_predicted_seq] ### cambiando/modificando: introducido allele_quality, complete_predicted_seq y etiqueta NIPHEM

    return True 
"""

# Funcion inf_asm_alm_tag incorporando PLOTs y unificando ifs

def inf_asm_alm_tag(core_name, sample_name, tag, blast_values, allele_quality, new_sseq, matching_allele_length, tag_dict, list_tag, samples_matrix_dict, matching_genes_dict, prodigal_report, start_prodigal, end_prodigal, complete_predicted_seq, annotation_core_dict, logger): 
    
    print("Ha entrado a función inf_asm_alm_tag para etiquetar la secuencia encontrada", '\n')

    # blast_values ----> allele_found[allele]
    # matching_allele_length ----> meter secuencia y calcular aquí dentro la longitud o meter directamente la longitud? 
    

    # obtener allele_quality aquí dentro?


    # core_name ---> nombre locus X
    # tag_dict ---> asm_dict X
    # new_sseq ---> secuencia encontrada completa X
    # tag ---> INF, ASM, ALM X
    # qseqid ## values ---> id del alelo
    # samples_matrix_dict X
    # sample_name ---> nombre muestra X
    # sseq ## values --> sin gaps, supongo
    # matching_allele_length X --> longitud alelo match, hace falta fuera que no sea para este tag? si no hace falta para otra cosa introducir la secuencia y calcular la longitud aquí dentro
    # list_tag X --> list_asm


    # METER values mejor y sacar todos estos params dentro de la función?

    # qseqid ## values
    # allele_quality X
    # sseqid ## values
    # bitscore ## values
    # str(s_length) --> longitud de la sseq, se puede obtener dentro de la función a partir de la sseq sin gaps. No sé si hará falta sacar la sseq sin gaps aparte de para la clasificación. sI no hace falta para nada más entonces sacar la sseq sin gaps dentro de esta función
    # str(new_sequence_length) --> longitud de la new_sseq --> meter como argumento o calcular dentro de la función aunque se calcule fuera antes para clasificar en INF, ASM y ALM?
    # mismatch ## values
    # gapopen ## values
    # sstart ## values
    # send ## values

    # samples_matrix_dict X
    # matching_genes_dict X

    # obteniendo anotación gen y producto del locus
    gene_annot, product_annot = annotation_core_dict[core_name]

    qseqid, sseqid, pident,  qlen, s_length, mismatch, gapopen, evalue, bitscore, sstart, send, qstart, qend, sseq, qseq = blast_values

    sseq = sseq.replace('-', '') # se utiliza para report plot, así que dejar
    s_length = len(sseq) ### cambiando/modificando: obteniendo longitud sseq sin gaps, ya que s_length incluye gaps en caso de que los haya en el alineamiento
    new_sequence_length = len(new_sseq)

    logger.info('Found %s at sample %s for core gene %s ', tag, sample_name, core_name)                    

    if tag == 'PLOT':
        tag_allele = tag + '_' + str(qseqid)
    else:
        ### adding ASM allele to the asm_allele_matrix if it is not already include
        if not core_name in tag_dict:
            tag_dict[core_name] = []
        if not new_sseq in tag_dict[core_name] :
            tag_dict[core_name].append(new_sseq)
        ### find the index of ASM  to include it in the sample matrix dict
        index_tag = tag_dict[core_name].index(new_sseq)

        tag_allele = tag + '_' + core_name + '_' + str(qseqid) + '_' + str(index_tag)
    
    print("Alelos inf encontrados para este locus: ", tag_dict, '\n')

    samples_matrix_dict[sample_name].append(tag_allele)
    #if new_sequence_length < query_length : ### deleción
    #if new_sequence_length < min(schema_variability[core_name]) : ### inserción
    #if new_sequence_length < min_length_threshold:

    print("tag, inf_asm_alm_tag: ", tag, '\n')

    ### add the deletion into deletion list
    if not core_name in list_tag : # list_tag ---> plot_dict
        list_tag[core_name] = {}
    if not sample_name in list_tag[core_name] :
        list_tag[core_name][sample_name] = {}   

    # Fusionar este if y el siguiente y meter logger para cada etiqueta?
    if tag == 'INF':
        list_tag[core_name][sample_name][tag_allele] = [gene_annot, product_annot, qseqid, allele_quality, sseqid,  bitscore, str(matching_allele_length), str(s_length), str(new_sequence_length), mismatch , gapopen, sstart, send,  new_sseq]
        ### duda cambiando/modificando: no debería cambiar en list_deletions sstart y send por las sstart y send de la nueva secuencia completa? Ahora mismo se está guardando sstart y send de la secuencia sseq encontrada con BLAST
    
        print("HA ENTRADO A ETIQUETA INF en inf_asm_alm_tag", '\n')

    elif tag == 'PLOT':
        list_tag[core_name][sample_name] = [gene_annot, product_annot, qseqid, allele_quality, sseqid, bitscore, sstart, send, sseq, new_sseq] ### new_sseq equivalente a complete_predicted_seq ### Modificando/cambiando: he añadido allele_quality y la seq predicha por prodigal complete_predicted_seq cuando la calidad del alelo es bad_quality

    else :
        if tag == 'ASM':
            newsseq_vs_blastseq = 'shorter'
        elif tag == 'ALM':
            newsseq_vs_blastseq = 'longer'

        if len(sseq) < matching_allele_length: ### si la secuencia obtenida con BLAST tiene una longitud menor a la del alelo con el que hizo match -> DELETION
            add_info = 'Global effect: DELETION. BLAST sequence length shorter than matching allele sequence length / Net result: ' + tag + '. Final gene sequence length ' + newsseq_vs_blastseq + 'than matching allele sequence length'
            #delete_allele = 'ASM_DELETE_' + core_name + '_' + str(qseqid) + '_' + str(index_delete)

        elif len(sseq) == matching_allele_length:### si la secuencia obtenida con BLAST tiene la misma longitud que el alelo con el que hizo match -> 
            add_info = 'Global effect: BASE SUBSTITUTION. BLAST sequence length equal to matching allele sequence length / Net result: ' + tag + '. Final gene sequence length ' + newsseq_vs_blastseq + 'than matching allele sequence length'

        #if new_sequence_length > query_length : ### deleción
        #if new_sequence_length > max(schema_variability[core_name]) : ### inserción
        #if new_sequence_length > max_length_threshold:
        elif len(sseq) > matching_allele_length: ### si la secuencia obtenida con BLAST tiene una longitud mayor a la del alelo con el que hizo match -> INSERTION
            add_info = 'Global effect: INSERTION. BLAST sequence length longer than matching allele sequence length / Net result: ' + tag + '. Final gene sequence length ' + newsseq_vs_blastseq + 'than matching allele sequence length'
            #delete_allele = 'ALM_DELETE_' + core_name + '_' + str(qseqid) + '_' + str(index_delete)

        list_tag[core_name][sample_name][tag_allele] = [gene_annot, product_annot, qseqid, allele_quality, sseqid,  bitscore, str(matching_allele_length), str(s_length), str(new_sequence_length), mismatch , gapopen, sstart, send,  new_sseq, add_info]
        ### duda cambiando/modificando: no debería cambiar en list_deletions sstart y send por las sstart y send de la nueva secuencia completa? Ahora mismo se está guardando sstart y send de la secuencia sseq encontrada con BLAST

    #  print("add_info: ", add_info, '\n')                       

    if not sseqid in matching_genes_dict[sample_name] :
        matching_genes_dict[sample_name][sseqid] = []
    if sstart > send :
        matching_genes_dict[sample_name][sseqid].append([core_name, str(int(sstart)-new_sequence_length -1), sstart,'-', tag_allele]) ### modificando/cambiando: al hacer predecir con prodigal habría que incluir sstart y send obteniéndolos de los archivos que genera prodigal para saber dónde están las coordenadas de la nueva seq completa
    else:
        matching_genes_dict[sample_name][sseqid].append([core_name, sstart,str(int(sstart)+ new_sequence_length),'+', tag_allele])

    ##### PARA SACAR INFORME DE PRODIGAL VS BLAST #####
    ##### introducido como argumento start_prodigal, end_prodigal y complete_predicted_seq para influir en el report de prodigal, temporal
    prodigal_report.append([core_name, sample_name, qseqid, tag_allele, sstart, send, start_prodigal, end_prodigal, sseq, complete_predicted_seq])

    return True


# Funcion inf_asm_alm_tag incorporando PLOTs anterior a unificar ifs
"""
def inf_asm_alm_tag(core_name, sample_name, tag, blast_values, allele_quality, new_sseq, matching_allele_length, tag_dict, list_tag, samples_matrix_dict, matching_genes_dict, prodigal_report, start_prodigal, end_prodigal, complete_predicted_seq, logger):
    
    # blast_values ----> allele_found[allele]
    # matching_allele_length ----> meter secuencia y calcular aquí dentro la longitud o meter directamente la longitud? 
    

    # obtener allele_quality aquí dentro?


    # core_name ---> nombre locus X
    # tag_dict ---> asm_dict X
    # new_sseq ---> secuencia encontrada completa X
    # tag ---> INF, ASM, ALM X
    # qseqid ## values ---> id del alelo
    # samples_matrix_dict X
    # sample_name ---> nombre muestra X
    # sseq ## values --> sin gaps, supongo
    # matching_allele_length X --> longitud alelo match, hace falta fuera que no sea para este tag? si no hace falta para otra cosa introducir la secuencia y calcular la longitud aquí dentro
    # list_tag X --> list_asm


    # METER values mejor y sacar todos estos params dentro de la función?

    # qseqid ## values
    # allele_quality X
    # sseqid ## values
    # bitscore ## values
    # str(s_length) --> longitud de la sseq, se puede obtener dentro de la función a partir de la sseq sin gaps. No sé si hará falta sacar la sseq sin gaps aparte de para la clasificación. sI no hace falta para nada más entonces sacar la sseq sin gaps dentro de esta función
    # str(new_sequence_length) --> longitud de la new_sseq --> meter como argumento o calcular dentro de la función aunque se calcule fuera antes para clasificar en INF, ASM y ALM?
    # mismatch ## values
    # gapopen ## values
    # sstart ## values
    # send ## values

    # samples_matrix_dict X
    # matching_genes_dict X


    qseqid, sseqid, pident,  qlen, s_length, mismatch, gapopen, evalue, bitscore, sstart, send, qstart, qend, sseq, qseq = blast_values


    sseq = sseq.replace('-', '') # se utiliza para report plot, así que dejar
    s_length = len(sseq) ### cambiando/modificando: obteniendo longitud sseq sin gaps, ya que s_length incluye gaps en caso de que los haya en el alineamiento
    new_sequence_length = len(new_sseq)

    logger.info('Found %s at sample %s for core gene %s ', tag, sample_name, core_name)                    

    if tag == 'PLOT':
        tag_allele = 'PLOT_' + str(qseqid)
    else:
        ### adding ASM allele to the asm_allele_matrix if it is not already include
        if not core_name in tag_dict:
            tag_dict[core_name] = []
        if not new_sseq in tag_dict[core_name] :
            tag_dict[core_name].append(new_sseq)
        ### find the index of ASM  to include it in the sample matrix dict
        index_tag = tag_dict[core_name].index(new_sseq)

        tag_allele = tag + '_' + core_name + '_' + str(qseqid) + '_' + str(index_tag)
    
    samples_matrix_dict[sample_name].append(tag_allele)
    #if new_sequence_length < query_length : ### deleción
    #if new_sequence_length < min(schema_variability[core_name]) : ### inserción
    #if new_sequence_length < min_length_threshold:

    print("tag, inf_asm_alm_tag: ", tag, '\n')


    # Fusionar este if y el siguiente y meter logger para cada etiqueta?

    if tag == 'ASM' or tag == 'ALM':

        if tag == 'ASM':
            newsseq_vs_blastseq = 'shorter'
        elif tag == 'ALM':
            newsseq_vs_blastseq = 'longer'

        if len(sseq) < matching_allele_length: ### si la secuencia obtenida con BLAST tiene una longitud menor a la del alelo con el que hizo match -> DELETION
            add_info = 'Global effect: DELETION. BLAST sequence length shorter than matching allele sequence length / Net result: ' + tag + '. Final gene sequence length ' + newsseq_vs_blastseq + 'than matching allele sequence length'
            #delete_allele = 'ASM_DELETE_' + core_name + '_' + str(qseqid) + '_' + str(index_delete)

        if len(sseq) == matching_allele_length:### si la secuencia obtenida con BLAST tiene la misma longitud que el alelo con el que hizo match -> 
            add_info = 'Global effect: BASE SUBSTITUTION. BLAST sequence length equal to matching allele sequence length / Net result: ' + tag + '. Final gene sequence length ' + newsseq_vs_blastseq + 'than matching allele sequence length'

        #if new_sequence_length > query_length : ### deleción
        #if new_sequence_length > max(schema_variability[core_name]) : ### inserción
        #if new_sequence_length > max_length_threshold:
        if len(sseq) > matching_allele_length: ### si la secuencia obtenida con BLAST tiene una longitud mayor a la del alelo con el que hizo match -> INSERTION
            add_info = 'Global effect: INSERTION. BLAST sequence length longer than matching allele sequence length / Net result: ' + tag + '. Final gene sequence length ' + newsseq_vs_blastseq + 'than matching allele sequence length'
            #delete_allele = 'ALM_DELETE_' + core_name + '_' + str(qseqid) + '_' + str(index_delete)

      #  print("add_info: ", add_info, '\n')

    ### add the deletion into deletion list
    if not core_name in list_tag : # list_tag ---> plot_dict
        list_tag[core_name] = {}
    if not sample_name in list_tag[core_name] :
        list_tag[core_name][sample_name] = {}                          

    ######## FUSIONAR ESTE IF CON IF ANTERIOR?                    
    if tag == 'INF':
        list_tag[core_name][sample_name][tag_allele] = [qseqid, allele_quality, sseqid,  bitscore, str(matching_allele_length), str(s_length), str(new_sequence_length), mismatch , gapopen, sstart, send,  new_sseq]
        ### duda cambiando/modificando: no debería cambiar en list_deletions sstart y send por las sstart y send de la nueva secuencia completa? Ahora mismo se está guardando sstart y send de la secuencia sseq encontrada con BLAST
    elif tag == 'PLOT':
        list_tag[core_name][sample_name].append([qseqid, allele_quality, sseqid, bitscore, sstart, send, sseq, new_sseq]) ### new_sseq equivalente a complete_predicted_seq ### Modificando/cambiando: he añadido allele_quality y la seq predicha por prodigal complete_predicted_seq cuando la calidad del alelo es bad_quality

    else:
        list_tag[core_name][sample_name][tag_allele] = [qseqid, allele_quality, sseqid,  bitscore, str(matching_allele_length), str(s_length), str(new_sequence_length), mismatch , gapopen, sstart, send,  new_sseq, add_info]
        ### duda cambiando/modificando: no debería cambiar en list_deletions sstart y send por las sstart y send de la nueva secuencia completa? Ahora mismo se está guardando sstart y send de la secuencia sseq encontrada con BLAST
   
    if not sseqid in matching_genes_dict[sample_name] :
        matching_genes_dict[sample_name][sseqid] = []
    if sstart > send :
        matching_genes_dict[sample_name][sseqid].append([core_name, str(int(sstart)-new_sequence_length -1), sstart,'-', tag_allele]) ### modificando/cambiando: al hacer predecir con prodigal habría que incluir sstart y send obteniéndolos de los archivos que genera prodigal para saber dónde están las coordenadas de la nueva seq completa
    else:
        matching_genes_dict[sample_name][sseqid].append([core_name, sstart,str(int(sstart)+ new_sequence_length),'+', tag_allele])

    ##### PARA SACAR INFORME DE PRODIGAL VS BLAST #####
    ##### introducido como argumento start_prodigal, end_prodigal y complete_predicted_seq para influir en el report de prodigal, temporal
    prodigal_report.append([core_name, sample_name, qseqid, tag_allele, sstart, send, start_prodigal, end_prodigal, sseq, complete_predicted_seq])

    return True
"""

def get_blast_results (core_name, values, sample_contigs, allele_found, logger) : ## introduciendo core_name para acceder a alelo en diccionario locus
    ## get_sample_sequence? en lugar de get_blast_results?
 
    #values
    #records ---> contigs
    #allele_found
    #paralog_found si no lo quito del código
  
    qseqid, sseqid, pident, qlen, s_length, mismatch, gapopen, evalue, bitscore, sstart, send, qstart, qend, sseq, qseq = values

    sseqid_blast = sseqid.split('_')[1] # obteniendo el id del contig a partir del id de la secuencia

    # igual no sería necesario guardar como descripción de cada uand e las secuencias encontradas tras el primer blast todos los valores obtenidos
    # ya que lo único que se neciestia es el sseqid que es el contig y el sstart y send hay que recalcularlo aquí. Entonces el sseqid se podría sacar
    # directamente haciendo split del id de la secuencia guardada que la había guardado del modo numsec_sseqid. De hecho, entonces, no sería necesario
    # parsear la secuencia del fasta, ya que ya se tiene el sseqid del segundo blast que equivale al id de la secuencia, por lo que se haría split del 
    # sseqid para quedarme con el segundo elemento y ya tendría el sseqid del primer blast aka contig. y tampoco necesito parsear la secuencia ya que
    # lo que tengo que hacer es usar la sseq que se ha obtenido tras el segundo blast quitándole los gaps y luego para obtener el otro índice de la seq
    # tendría que sumarle la longitud de la secuencia sin gaps (vale, no, creo que tiene que ser con gaps? cómo tendría que ser para que se correspondiese
    # la longitud con el start y send que sale del segundo blast? -> creo que tiene que ser sin gaps, sep, definitivamente, si no no tendría sentido)

    sseq_no_gaps = sseq.replace('-', '') # sseq sin gaps obtenida tras blast con todos los alelos

    # se obtiene sstart y send de la secuencia en el contig buscando el indice de la secuencia en el mismo
    
    for record in sample_contigs: ## comentar parse
        if record.id == sseqid_blast : ## comentar parse
            break ## comentar parse
    accession_sequence = str(record.seq) ## comentar parse

    accession_sequence = alleles_in_locus_dict[core_name][sseqid] ## Alternativa a parser
    
    try: # Se intenta obtener el indice de la secuencia encontrada tras blast en el contig. Si se obtiene error porque el sentido de la seq no es el correcto, se vuelve a buscar el indice hallando previamente la seq complementaria reversa
        sseq_index_1 = int(accession_sequence.index(sseq_no_gaps)) + 1

    except:
        sseq_no_gaps = str(Seq.Seq(sseq_no_gaps).reverse_complement())
        sseq_index_1 = int(accession_sequence.index(sseq_no_gaps)) + 1

    sseq_index_2 = int(sseq_index_1) + len(sseq_no_gaps) - 1


#    if sstart < send :
  #       if sseq_index_1 < sseq_index_2:
  #          sstart_new = int(sseq_index_1)
    #         send_new = int(sseq_index_2)
    #    else:
      #       sstart_new = int(sseq_index_2)
      #      send_new = int(sseq_index_1)
#    else:
  #       if sseq_index_1 < sseq_index_2:
  #          sstart_new = int(sseq_index_2)
    #         send_new = int(sseq_index_1)
    #    else:
      #       sstart_new = int(sseq_index_1)
      #      send_new = int(sseq_index_2)

    # se adjudican los indices obtenidos de la secuencia en el contig a sstart_new y send_new en funcion de sstart y send tras el blast con todas los alelos (sstart < send o viceversa)
    if int(sstart) < int(send):
        sstart_new = str(min(sseq_index_1, sseq_index_2))
        send_new = str(max(sseq_index_1, sseq_index_2))
    else:
        sstart_new = str(max(sseq_index_1, sseq_index_2))
        send_new = str(min(sseq_index_1, sseq_index_2))

    allele_is_subset = False
    if len(allele_found) > 0 :

        # check if the new match is a subset of the previous allele found in blast
        for allele_id in allele_found :
#            if allele_found[allele][9] == sstart_new or allele_found[allele][10] == send_new :
            #if (int(sstart_new) or int(send_new)) in range(min(int(allele_found[allele_id][9]), int(allele_found[allele_id][10])),  max(int(allele_found[allele_id][9]), int(allele_found[allele_id][10])) +1):

            # formulando if de este modo para arreglar falsos NIPH (estaba cogiendo subsets, ejemplo locus lmo2810)
            min_index = min(int(allele_found[allele_id][9]), int(allele_found[allele_id][10]))
            max_index = max(int(allele_found[allele_id][9]), int(allele_found[allele_id][10]))
            if int(sstart_new) in range(min_index, max_index + 1) or  int(send_new) in range(min_index, max_index + 1):    
                if sseqid_blast == allele_found[allele_id][1]: # metiendo condicion de igual contig para evitar que se escapen parálogos que se encuentra en otros contigs pero cuyas posiciones solapan con las de las secuencias ya encontradas (en otros contigs) y guardadas
                    logger.info('Found allele %s that starts or ends at the same position as %s ' , qseqid, allele_id)
                    allele_is_subset = True
                    break # ? es necesario?

    if len(allele_found) == 0 or not allele_is_subset :
        contig_id_start = str(sseqid_blast + '_'+ sstart_new)
        ## skip the allele found in the 100% identity and 100% alignment
        if not contig_id_start in allele_found:
            allele_found[contig_id_start] = [qseqid, sseqid_blast, pident, qlen, s_length, mismatch, gapopen, evalue, bitscore, sstart_new, send_new, '-', '-', sseq, qseq]

    return True

def keep_snp_alignment_info(new_sseq, matching_allele_seq, qseqid, query_direction, core_name, sample_name, snp_dict, match_alignment_dict, protein_dict, logger):

    # new_sseq --> secuencia encontrada completa
    # matching_allele_seq --> secuencia de alelo
    # core_name --> nombre locus
    # sample_name --> nombre muestra
    # qseqid --> id alelo
    # snp_dict --> diccionario snp
    # match_alignment_dict --> diccionario alineamiento adn
    # protein_dict --> diccionario alineamiento proteina


    #if check_sequence_order(matching_allele_seq, logger) == 'reverse':
    if query_direction == 'reverse':
        #matching_allele_seq = reverse_complement(matching_allele_seq)   ### (cambiando/modificando: he cambiado allele_sequence por matching_allele_seq)
        matching_allele_seq = str(matching_allele_seq.reverse_complement())
    else:
        matching_allele_seq = str(matching_allele_seq)

    # get the SNP information
    #snp_information = get_snp(sseq, reference_allele_for_snp)
    snp_information = get_snp(new_sseq, matching_allele_seq) ### cambiando/modificando: he cambiado sseq por new_sseq ya que ahora se está considerando la secuencia completa habiendo encontrado el codón de stop
                                                        ### cambiando/modificando: he cambiado reference_allele_for_snp (que es la secuencia del alelo que se cogía antes como referencia, que era el primero del locus) por el alelo con el que ha hecho match esta secuencia
    if len(snp_information) > 0 :
        if not core_name in snp_dict :
            snp_dict[core_name] = {}
        if not sample_name in snp_dict[core_name] :
            snp_dict[core_name][sample_name] = {}
        snp_dict[core_name][sample_name][qseqid]= snp_information

    # get new sequence-allele sequence dna alignment 
    if not core_name in match_alignment_dict :
        match_alignment_dict[core_name] = {}
        if not sample_name in match_alignment_dict[core_name] :
            #match_alignment_dict[core_name][sample_name] = get_alignment_for_deletions (new_sseq,  matching_allele_seq)           ### (cambiando/modificando: he cambiado el nombre de qqseq por matching_allele_seq para unificar (aunque realmente era lo mismo porque había adjudicad allele_sequence a qseq)
            match_alignment_dict[core_name][sample_name] = get_alignment (new_sseq, matching_allele_seq) ### cambiando/modificando: sustituyendo por get alignment

    # get new sequence-allele sequence protein alignment 
    if not core_name in protein_dict :
        protein_dict[core_name] = {}
    if not sample_name in protein_dict[core_name] :
        protein_dict[core_name][sample_name] = []
    #protein_dict[core_name][sample_name] = nucleotide_to_protein_alignment(new_sseq, matching_allele_seq)              ### (cambiando/modificando: he cambiado el nombre de qqseq por matching_allele_seq para unificar (aunque realmente era lo mismo porque había adjudicad allele_sequence a qseq)
    protein_dict[core_name][sample_name] = get_alignment (new_sseq, matching_allele_seq, "protein") ### cambiando/modificando: sustituyendo por get alignment
                        
    return True

def create_summary (samples_matrix_dict, logger) :

    summary_dict = {}
    summary_result_list = []
    summary_heading_list = ['Exact match', 'INF', 'ASM_INSERT', 'ASM_EQUAL', 'ASM_DELETE', 'ALM_INSERT', 'ALM_EQUAL', 'ALM_DELETE', 'LNF', 'NIPH', 'NIPHEM', 'PLOT', 'ERROR']
    summary_result_list.append('File\t' + '\t'.join(summary_heading_list))
    for key in sorted (samples_matrix_dict) :

        summary_dict[key] = {'Exact match':0, 'INF':0, 'ASM_INSERT':0, 'ASM_EQUAL':0, 'ASM_DELETE':0, 'ALM_INSERT':0, 'ALM_EQUAL':0, 'ALM_DELETE':0, 'LNF':0, 'NIPH':0, 'NIPHEM':0, 'PLOT':0, 'ERROR':0}
        for values in samples_matrix_dict[key] :
            if 'INF_' in values :
                summary_dict[key]['INF'] += 1
            elif 'ASM_INSERT' in values :
                summary_dict[key]['ASM_INSERT'] += 1
            elif 'ASM_DELETE' in values :
                summary_dict[key]['ASM_DELETE'] += 1
            elif 'ASM_EQUAL' in values : ### cambiando/modificando: etiqueta ASM_EQUAL añadida
                summary_dict[key]['ASM_EQUAL'] += 1
          #  elif 'AEM_DELETE' in values : ### cambiando/modificando: etiqueta AEM_DELETE eliminada
           #     summary_dict[key]['AEM_DELETE'] += 1
            elif 'ALM_INSERT' in values :
                summary_dict[key]['ALM_INSERT'] += 1
            elif 'ALM_DELETE' in values :
                summary_dict[key]['ALM_DELETE'] += 1
            elif 'ALM_EQUAL' in values : ### cambiando/modificando: etiqueta ALM_EQUAL añadida
                summary_dict[key]['ALM_EQUAL'] += 1
          #  elif 'AEM_INSERT' in values : ### cambiando/modificando: etiqueta AEM_INSERT eliminada
           #     summary_dict[key]['AEM_INSERT'] += 1
            elif 'LNF' in values :
                summary_dict[key]['LNF'] += 1
            elif 'NIPH' == values : 
              #  print("Sumando +1 a NIPH", '\n')
                summary_dict[key]['NIPH'] += 1
            elif 'NIPHEM' == values : ### cambiando/modificando: no cogía NIPHEM y lo contabilizaba como NIPH por 'NIPH' in values
               # print("Sumando +1 a NIPHEM", '\n')
                summary_dict[key]['NIPHEM'] += 1
            elif 'PLOT' in values :
                summary_dict[key]['PLOT'] += 1
            elif 'ERROR' in values :
                summary_dict[key]['ERROR'] += 1
            else:
                try:
                    number = int(values)
                    summary_dict[key]['Exact match'] +=1
                except:
                    if '_' in values :
                        tmp_value = values
                        try:
                            number = int(tmp_value[-1])
                            summary_dict[key]['Exact match'] +=1
                        except:
                            logger.debug('The value %s, was found when collecting summary information for the %s', values, summary_dict[key] )
                    else:
                        logger.debug('The value %s, was found when collecting summary information for the %s', values, summary_dict[key] )
        summary_sample_list = []
        for item in summary_heading_list :
            summary_sample_list.append(str(summary_dict[key][item]))
        summary_result_list.append(key + '\t' +'\t'.join(summary_sample_list))
    return summary_result_list


### cómo añadir opciones de prokka más avanzadas para obtener mejores resultados
### cómo generar solo el archivo que necesito y no todos
### borrar carpeta prokka (annotation) al acabar?
### meter algún logger dentro de get_reference_allele rollo "obteniendo alelo der ef de lcous X"y en get_gene_annotation rollo "obteniendo anotación de locus X".
def get_gene_annotation (annotation_file, annotation_dir, logger) :
    name_file = os.path.basename(annotation_file).split('.')
    annotation_dir = os.path.join (annotation_dir, 'annotation', name_file[0])
    
    #annotation_result = subprocess.run (['prokka', annotation_file , '--outdir' , str(annotation_dir + 'prokka_anotation' + name_file[0]),
    annotation_result = subprocess.run (['prokka', annotation_file , '--outdir' , annotation_dir ,
                                        '--prefix', name_file[0], '--quiet'])
    annot_tsv = []
    tsv_path = os.path.join (annotation_dir, name_file[0] + '.tsv')
    with open(tsv_path) as tsvfile:
        tsvreader = csv.reader(tsvfile, delimiter="\t")
        for line in tsvreader:
            annot_tsv.append(line)

#    print("annot_tsv: ", annot_tsv, '\n')
    if len(annot_tsv) > 1: # si hay anotacion 
        gene_annot = annot_tsv[1][2].split('_')[0]
        product_annot = annot_tsv[1][4]
    else:
        gene_annot = 'Annotation not found by Prokka'
        product_annot = 'Annotation not found by Prokka'

    print("gene_annot: ", gene_annot, '\n')
    print("product_annot: ", product_annot, '\n')

    return gene_annot, product_annot

"""
def get_gene_annotation (annotation_file, annotation_dir, logger) : ### No se está usando de momento
    name_file = os.path.basename(annotation_file).split('.')
    annotation_dir = os.path.join (annotation_dir, 'annotation', name_file[0])
    '''
    if not os.path.exists(annotation_dir):
        try:
            os.makedirs(annotation_dir)
            logger.debug(' Created local annotation directory for  %s', name_file[0])
        except:
            logger.info('Cannot create directory for local annotation %s' , f_name[0])
            print ('Error when creating the directory %s for prokka. ', annotation_dir)
            exit(0)
    '''
    #annotation_result = subprocess.run (['prokka', annotation_file , '--outdir' , str(annotation_dir + 'prokka_anotation' + name_file[0]),
    annotation_result = subprocess.run (['prokka', annotation_file , '--outdir' , annotation_dir ,
                                         '--prefix', name_file[0]])

    return str(annotation_dir + 'prokka_anotation' + name_file[0] + name_file[0] + '.gff')
"""

def analize_annotation_files (in_file, logger) : ### No se está usando de momento
    examiner = GFF.GFFExaminer()
    file_fh = open(in_file)
    datos = examiner.available_limits(in_file)
    file_fh.close()
    return True

def get_inferred_allele_number(core_dict, logger): ### No se está usando ### utiliza core_dict comentado
    #This function will look for the highest locus number and it will return a safe high value
    # that will be added to the schema database
    logger.debug('running get_inferred_allele_number function')
    int_keys = []
    for key in core_dict.keys():
        int_keys.append(key)
    max_value = max(int_keys)
    digit_length = len(str(max_value))
    return  True   #str 1 ( #'1'+ '0'*digit_length + 2)

def save_results (full_gene_list, samples_matrix_dict, exact_dict, paralog_dict, inf_dict, plot_dict, matching_genes_dict, list_asm, list_alm, lnf_dict, snp_dict, match_alignment, protein_dict, prodigal_report, shorter_seq_coverage, longer_seq_coverage, equal_seq_coverage, shorter_blast_seq_coverage, longer_blast_seq_coverage, equal_blast_seq_coverage):
    header_matching_alleles_contig = ['Sample Name', 'Contig', 'Core Gene','start', 'stop', 'direction', 'codification']
    header_exact = ['Core Gene', 'Sample Name', 'Gene Annotation', 'Product Annotation', 'Allele', 'Allele Quality', 'Contig', 'Query length', 'Contig start', 'Contig end', 'Sequence', 'Predicted Sequence'] # c/m: header para diccionario exact_dict añadido
    header_paralogs = ['Core Gene','Sample Name', 'Gene Annotation', 'Product Annotation', 'Paralog Type', 'ID %', 'Allele', 'Allele Quality', 'Contig', 'Bit Score', 'Contig start', 'Contig end', 'Sequence', 'Predicted Sequence'] # c/m: introducido paralog type, Allele Quality y Predicted Sequence
    header_inferred = ['Core Gene','Sample Name', 'Gene Annotation', 'Product Annotation','Inferred Allele name', 'Allele', 'Allele Quality', 'Contig', 'Bitscore', 'Query length', 'Contig length', 'New sequence length' , 'Mismatch' , 'gaps', 'Contig start', 'Contig end',  'New sequence']
    # c/m: añadido header para clasificación asm/alm. quitando clasificación delecion/insercion/equal
    header_asm = [ 'Core Gene', 'Sample Name', 'Gene Annotation', 'Product Annotation', 'ASM item', 'Allele', 'Allele Quality', 'Contig', 'Bitscore', 'Query length', 'Contig length', 'New sequence length' , 'Mismatch' , 'gaps', 'Contig start', 'Contig end',  'New sequence', 'Additional info']
    header_alm = [ 'Core Gene', 'Sample Name', 'Gene Annotation', 'Product Annotation', 'ALM item', 'Allele', 'Allele Quality', 'Contig', 'Bitscore', 'Query length', 'Contig length', 'New sequence length' , 'Mismatch' , 'gaps', 'Contig start', 'Contig end',  'New sequence', 'Additional info']    
    header_plot = ['Core Gene', 'Sample Name', 'Gene Annotation', 'Product Annotation', 'Allele', 'Allele Quality', 'Contig','Bit Score', 'Contig start', 'Contig end', 'Sequence', 'Predicted Sequence']
    # c/m: añadido header para report lnf
    header_lnf = ['Core Gene', 'Sample Name', 'Gene Annotation', 'Product Annotation', 'Allele', 'ID %', 'Blast sequence coverage', 'New sequence coverage', 'Allele length', 'Blast sequence length', 'New sequence length', 'Additional info'] ### meter secuencias alelo y newsseq en caso de que haya?        
    header_snp = ['Core Gene', 'Sample Name', 'Allele number', 'Position', 'Mutation Schema/Sample', 'Codon Schema/Sample','Protein in Schema/Sample', 'Missense/Synonymous','Annotation Sample / Schema']
    header_protein = ['Core Gene','Sample Name', 'Protein in ' , 'Protein sequence']
    header_match_alignment = ['Core Gene','Sample Name','Alignment', 'Sequence']
                
    # Añadido header_prodigal_report para report prodigal
    header_prodigal_report = ['Core gene', 'Sample Name', 'Allele', 'Sequence type', 'BLAST start', 'BLAST end', 'Prodigal start', 'Prodigal end', 'BLAST sequence', 'Prodigal sequence']
    # Añadido header_newsseq_coverage_report para determinar coverage threshold a imponer
    header_newsseq_coverage_report = ['Core gene', 'Sample Name', 'Query length', 'New sequence length', 'Locus mean', 'Coverage (new sequence/allele)', 'Coverage (new sequence/locus mean)']
    # Añadido header_blast_coverage_report para determinar coverage threshold a imponer
    header_blast_coverage_report = ['Core gene', 'Sample Name', 'Query length', 'Blast sequence length', 'Locus mean', 'Coverage (blast sequence/allele)', 'Coverage (blast sequence/locus mean)']

    print ('Saving results to files \n')
    result_file = os.path.join ( outputdir, 'result.tsv')
    # saving the result information to file
    logger.info('Saving result information to file..')
    with open (result_file, 'w') as out_fh:
        out_fh.write ('Sample Name\t'+'\t'.join( full_gene_list) + '\n')
        for key in sorted (samples_matrix_dict):
            out_fh.write (key + '\t' + '\t'.join(samples_matrix_dict[key])+ '\n')

    ###########################################################################################
    # Guardando report de prodigal. temporal

    prodigal_report_file = os.path.join (outputdir, 'prodigal_report.tsv')
    # saving the result information to file
    with open (prodigal_report_file, 'w') as out_fh:
        out_fh.write ('\t'.join(header_prodigal_report)+ '\n')
        for prodigal_result in prodigal_report:
            out_fh.write ('\t'.join(prodigal_result)+ '\n')

    ###########################################################################################

    ###########################################################################################
    # Guardando coverage de new_sseq para estimar el threshold a establecer. temporal

    newsseq_coverage_file = os.path.join (outputdir, 'newsseq_coverage_report.tsv')
    # saving the result information to file
    with open (newsseq_coverage_file, 'w') as out_fh:
        out_fh.write ('\t' + '\t'.join(header_newsseq_coverage_report)+ '\n')
        for coverage in shorter_seq_coverage:
            out_fh.write ('Shorter new sequence' + '\t' + '\t'.join(coverage)+ '\n')
        for coverage in longer_seq_coverage:
            out_fh.write ('Longer new sequence' + '\t' + '\t'.join(coverage)+ '\n')
        for coverage in equal_seq_coverage:
            out_fh.write ('Same length new sequence' + '\t' + '\t'.join(coverage)+ '\n')

    # Guardando coverage de la sseq obtenida tras blast para estimar el threshold a establecer. temporal

    blast_coverage_file = os.path.join (outputdir, 'blast_coverage_report.tsv')
    # saving the result information to file
    with open (blast_coverage_file, 'w') as out_fh:
        out_fh.write ('\t' + '\t'.join(header_blast_coverage_report)+ '\n')
        for coverage in shorter_blast_seq_coverage:
            out_fh.write ('Shorter blast sequence' + '\t' + '\t'.join(coverage)+ '\n')
        for coverage in longer_blast_seq_coverage:
            out_fh.write ('Longer blast sequence' + '\t' + '\t'.join(coverage)+ '\n')
        for coverage in equal_blast_seq_coverage:
            out_fh.write ('Same length blast sequence' + '\t' + '\t'.join(coverage)+ '\n')

    ###########################################################################################

    # saving exact matches to file

    logger.info('Saving exact matches information to file..')
    exact_file =  os.path.join(outputdir, 'exact.tsv')
    with open (exact_file , 'w') as exact_fh :
        exact_fh.write('\t'.join(header_exact)+ '\n')
        for core in sorted(exact_dict): ### añadido sorted aquí y abajo, si da mal quitar
            for sample in sorted(exact_dict[core]):
                exact_fh.write(core + '\t' + sample + '\t' + '\t'.join(exact_dict[core][sample]) + '\n')

    # saving paralog sequence to file
    logger.info('Saving paralog information to file..')
    paralog_file =  os.path.join(outputdir, 'paralog.tsv')
    with open (paralog_file , 'w') as paralog_fh :
        paralog_fh.write('\t'.join(header_paralogs) + '\n')
        for sample in sorted (paralog_dict) :
            for core in sorted (paralog_dict[sample]):
                for paralog in paralog_dict[sample][core] :
                    paralog_fh.write(core + '\t' + sample + '\t' + '\t'.join (paralog) + '\n')

    # saving inferred alleles to file
    logger.info('Saving inferred alleles information to file..')
    inferred_file =  os.path.join(outputdir, 'inferred_alleles.tsv')
    with open (inferred_file , 'w') as infer_fh :
        infer_fh.write('\t'.join(header_inferred) + '\n')
        for core in sorted (inf_dict) :
            for sample in sorted (inf_dict[core]) :
                for inferred in inf_dict[core][sample]: 
                    #   seq_in_inferred_allele = '\t'.join (inf_dict[sample])
                    infer_fh.write(core + '\t' + sample + '\t' + inferred + '\t' + '\t'.join(inf_dict[core][sample][inferred]) + '\n')
    
    # saving PLOT  to file
    logger.info('Saving PLOT information to file..')
    plot_file =  os.path.join(outputdir, 'plot.tsv')
    with open (plot_file , 'w') as plot_fh :
        plot_fh.write('\t'.join(header_plot) + '\n')
        for core in sorted (plot_dict) :
            for sample in sorted (plot_dict[core]):
                plot_fh.write(core + '\t' + sample + '\t' + '\t'.join(plot_dict[core][sample]) + '\n')

    # saving matching contigs to file
    logger.info('Saving matching information to file..')
    matching_file =  os.path.join(outputdir, 'matching_contigs.tsv')
    with open (matching_file , 'w') as matching_fh :
        matching_fh.write('\t'.join(header_matching_alleles_contig ) + '\n')
        for samples in sorted ( matching_genes_dict) :
            for contigs in matching_genes_dict[samples] :
                for contig in matching_genes_dict[samples] [contigs]:
                        matching_alleles = '\t'.join (contig)
                        matching_fh.write(samples + '\t' + contigs +'\t' + matching_alleles + '\n')
 
    # saving ASMs to file
    logger.info('Saving asm information to file..')
    asm_file =  os.path.join(outputdir, 'asm.tsv')
    with open (asm_file , 'w') as asm_fh :
        asm_fh.write('\t'.join(header_asm)+ '\n')
        for core in list_asm :
            for sample in list_asm[core] :
                for asm in list_asm[core][sample]:
                    asm_fh.write(core + '\t' + sample + '\t' + asm + '\t' + '\t'.join(list_asm[core][sample][asm]) + '\n')

    # saving ALMs to file
    logger.info('Saving alm information to file..')
    alm_file =  os.path.join(outputdir, 'alm.tsv')
    with open (alm_file , 'w') as alm_fh :
        alm_fh.write('\t'.join(header_alm)+ '\n')
        for core in list_alm :
            for sample in list_alm[core] :
                for alm in list_alm[core][sample]:
                    alm_fh.write(core + '\t' + sample + '\t' + alm + '\t' + '\t'.join(list_alm[core][sample][alm]) + '\n')

    # saving LNFs to file
    logger.info('Saving lnf information to file..')
    lnf_file =  os.path.join(outputdir, 'lnf.tsv')
    with open (lnf_file , 'w') as lnf_fh :
        lnf_fh.write('\t'.join(header_lnf)+ '\n')
        for core in lnf_dict :
            for sample in lnf_dict[core] :
                for lnf in lnf_dict[core][sample] :
                    lnf_fh.write(core + '\t' + sample + '\t' + '\t'.join(lnf) + '\n')

    # saving SNPs  to file
    logger.info('Saving SNPs information to file..')
    snp_file =  os.path.join(outputdir, 'snp.tsv')
    with open (snp_file , 'w') as snp_fh :
        snp_fh.write('\t'.join(header_snp) + '\n')
        for core in sorted (snp_dict) :
            for sample in sorted (snp_dict[core]):
                for allele_id_snp in snp_dict[core][sample] :
                    for snp in snp_dict[core][sample][allele_id_snp] :
                        snp_fh.write(core + '\t' + sample + '\t' + allele_id_snp + '\t' + '\t'.join (snp) + '\n')


    logger.info('Saving matching alignment information to files..')
    alignment_dir = os.path.join(outputdir,'alignments')
    if os.path.exists(alignment_dir) :
        shutil.rmtree(alignment_dir)
        logger.info('deleting the alignment files from previous execution')
    os.makedirs(alignment_dir)
    for core in sorted(match_alignment_dict) :
        for sample in sorted (match_alignment_dict[core]) :
            match_alignment_file = os.path.join(alignment_dir, str('match_alignment_' + core + '_' + sample + '.txt'))
            with open(match_alignment_file, 'w') as match_alignment_fh :
                match_alignment_fh.write( '\t'.join(header_match_alignment) + '\n')
                for match_align in match_alignment_dict[core][sample] :
                    match_alignment_fh.write(core + '\t'+ sample +'\t'+ '\t'.join(match_align) + '\n')

    # saving protein in a separated file
    logger.info('Saving protein information to files..')
    protein_dir = os.path.join(outputdir,'proteins')
    if os.path.exists(protein_dir) :
        shutil.rmtree(protein_dir)
        logger.info('deleting the proteins files from previous execution')
    os.makedirs(protein_dir)
    for core in sorted(protein_dict) :
        for sample in sorted (protein_dict[core]) :
            protein_file = os.path.join(protein_dir, str('protein_' + core + '_' + sample + '.txt'))
            with open(protein_file, 'w') as protein_fh :
                protein_fh.write( '\t'.join(header_protein) + '\n')
                for protein in protein_dict[core][sample] :
                    protein_fh.write(core + '\t'+ sample +'\t'+ '\t'.join(protein) + '\n')

    ### create summary information
    logger.info('Saving summary information to file..')
    summary_result = create_summary (samples_matrix_dict, logger)
    summary_file = os.path.join( outputdir, 'summary_result.tsv')
    with open (summary_file , 'w') as summ_fh:
        for line in summary_result :
            summ_fh.write(line + '\n')

    ### modify the result file to remove the PLOT_ string for creating the file to use in the tree diagram
    logger.info('Saving result information for tree diagram')
    tree_diagram_file = os.path.join ( outputdir, 'result_for_tree_diagram.tsv')
    with open (result_file, 'r') as result_fh:
        with open(tree_diagram_file, 'w') as td_fh:
            for line in result_fh:
                tree_line = line.replace('PLOT_','')
                td_fh.write(tree_line)

    return True


def allele_call_nucleotides (core_gene_list_files, sample_list_files, alleles_in_locus_dict, contigs_in_sample_dict, query_directory, reference_alleles_directory, blast_db_directory, prodigal_directory, blast_results_seq_directory, blast_results_db_directory, inputdir, outputdir, cpus, percentlength, coverage, evalue, perc_identity_ref, perc_identity_loc, reward, penalty, gapopen, gapextend, max_target_seqs, max_hsps, num_threads, flankingnts, schema_variability, schema_statistics, schema_quality, annotation_core_dict, logger ): ### CAMBIANDO/MODIFICANDO: he añadido schema statistics para poder tomar la mean, stdev, etc. No sé si eliminar schema variability. He añadido prodigal_directory, prodigal_train_db_directory y schema_quality
    
    prodigal_report = [] # TEMPORAL. prodigal_report para checkear las secuencias obtenidas con prodigal vs blast y las posiciones sstart y send
    # listas añadidas para calcular coverage medio de new_sseq con respecto a alelo para establecer coverage mínimo por debajo del cual considerar LNF
    shorter_seq_coverage = [] # TEMPORAL
    longer_seq_coverage = [] # TEMPORAL
    equal_seq_coverage = [] # TEMPORAL
    # listas añadidas para calcular coverage medio de sseq con respecto a alelo tras blast para establecer coverage mínimo por debajo del cual considerar LNF
    shorter_blast_seq_coverage = [] # TEMPORAL
    longer_blast_seq_coverage = [] # TEMPORAL
    equal_blast_seq_coverage = [] # TEMPORAL

    full_gene_list = []
    samples_matrix_dict = {} # to keep allele number
    matching_genes_dict = {} # to keep start and stop positions
    exact_dict = {} # c/m: to keep exact matches found for each sample
    inferred_alleles_dict = {} # to keep track of the new inferred alleles
    inf_dict = {} # to keep inferred alleles found for each sample
    paralog_dict = {} # to keep paralogs found for each sample
    asm_dict = {} # c/m: to keep track of asm
    alm_dict = {} # c/m: to keep track of alm
    list_asm = {} # c/m: to keep asm found for each sample
    list_alm = {} # c/m: to keep alm found for each sample
    lnf_dict = {} # c/m: to keep locus not found for each sample
    plot_dict = {} # c/m: to keep plots for each sample
    snp_dict = {} # c/m: to keep snp information for each sample
    protein_dict = {}
    match_alignment_dict = {}

    blast_parameters = '"6, qseqid, sseqid, pident,  qlen, length, mismatch, gapopen, evalue, bitscore, sstart, send, qstart, qend, sseq, qseq"'

    print('Allele calling starts')
    pbar = ProgressBar ()

    ## # # # # # # # # # # # # # # # # # # # # # # # # ##
    ## Processing the search for each schema core gene ##
    ## # # # # # # # # # # # # # # # # # # # # # # # # ##

    for core_file in pbar(core_gene_list_files) :
        core_name = os.path.basename(core_file)
        full_gene_list.append(core_name)
        logger.info('Processing core gene file %s ', core_file) 

        # Get path to this locus fasta file
        locus_alleles_path = os.path.join(query_directory, str(core_name + '.fasta'))
        
        # Get path to reference allele fasta file for this locus
        core_reference_allele_path = os.path.join(reference_alleles_directory, core_name + '.fasta')
        
        # Get length thresholds for INF, ASM and ALM classification
        max_length_threshold, min_length_threshold = length_thresholds(core_name, schema_statistics, percentlength)

        # Get length thresholds for LNF, ASM and ALM classification
        max_coverage_threshold, min_coverage_threshold = length_thresholds(core_name, schema_statistics, coverage)

        ## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ##
        ## Processing the search for each schema core gene in each sample  ##
        ## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ##
        
        for sample_file in sample_list_files:
        logger.info('Processing sample file %s ', sample_file)

            sample_name = os.path.basename(sample_file)

            # Initialize the sample list to add the number of alleles and the start, stop positions
            if not sample_name in samples_matrix_dict:
                samples_matrix_dict[sample_name] = []
                matching_genes_dict[sample_name] = {}

            # Path to this sample BLAST database created when processing samples 
            blast_db_name = os.path.join(blast_db_directory, sample_name, sample_name)

            # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
            # Sample contigs VS reference allele(s) BLAST for locus detection in sample #
            # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #

            cline = NcbiblastnCommandline(db=blast_db_name, evalue=evalue, perc_identity=perc_identity_ref, reward=reward, penalty=penalty, gapopen=gapopen, gapextend=gapextend, outfmt = blast_parameters, max_target_seqs=max_target_seqs, max_hsps=max_hsps, num_threads=num_threads, query=core_reference_allele_path)
            out, err = cline()
            out_lines = out.splitlines()
            
            bigger_bitscore = 0 

            # ······························································ #
            # LNF if there are no BLAST results for this gene in this sample #
            # ······························································ #
            if len (out_lines) == 0:

                print("Muestra ", sample_name, '\n') ### printeando, borrar
                print("LNF", '\n', '\n') ### printeando, borrar

                # DUDA: tiene sentido hacer este BLAST si se hace frente al alelo de referencia e igualmente no se va a saber cuál es el alelo que hace match? Hacerlo con todos los alelos? no hacerlo?
                
                # Trying to get the allele number to avoid that a bad quality assembly impact on the tree diagram
                cline = NcbiblastnCommandline(db=blast_db_name, evalue=evalue, perc_identity = 70, reward=reward, penalty=penalty, gapopen=gapopen, gapextend=gapextend, outfmt=blast_parameters, max_target_seqs=1, max_hsps=1, num_threads=1, query=core_reference_allele_path)
                out, err = cline()
                out_lines = out.splitlines()

                #print("Número de resultados tras BLAST 1: ", len(out_lines), '\n')

                if len (out_lines) > 0 :
                
                    for line in out_lines :
                        values = line.split('\t')
                        if  float(values[8]) > bigger_bitscore:
                            qseqid , sseqid , pident ,  qlen , s_length , mismatch , gapopen , evalue , bitscore , sstart , send , qstart , qend ,sseq , qseq = values
                            bigger_bitscore = float(bitscore)

                    # Keep LNF info
                    lnf_tag(core_name, sample_name, alleles_in_locus_dict, samples_matrix_dict, lnf_dict, locus_alleles_path, qseqid, pident, '-', '-', logger)

                else:
                    # Keep LNF info
                    lnf_tag(core_name, sample_name, '-', samples_matrix_dict, lnf_dict, locus_alleles_path, '-', '-', '-', '-', logger)

                continue

            ## Continue classification process if the core gene has been detected in sample after BLAST search
            if len (out_lines) > 0:

                # Parse contigs for this sample 
                contig_file = os.path.join(inputdir, str(sample_name + ".fasta")) ## comentar parse
                records = list(SeqIO.parse(contig_file, "fasta")) ## comentar parse

                ## Keep BLAST results after locus detection in sample using reference allele
                
                # Path to BLAST results fasta file
                path_to_blast_seq = os.path.join(blast_results_seq_directory, sample_name, core_name + "_blast.fasta")
                
                with open (path_to_blast_seq, 'w') as outblast_fh:
                    seq_number = 1
                    for line in out_lines :
                        values = line.split('\t')
                        qseqid = values[0]
                        sseqid = values[1]
                        sstart = values[9]
                        send = values[10]

                        # Get flanked found BLAST sequences from contig for correct allele tagging
                        for record in records: ## comentar parse
                            if record.id == sseqid : ## comentar parse
                                break ## comentar parse
                        accession_sequence = str(record.seq) ## comentar parse

                        accession_sequence = contigs_in_sample_dict[core_name][sseqid] ## accediendo a diccionario. Alternativa a parse

                        if int(send) > int(sstart):
                            max_index = int(send)
                            min_index = int(sstart)
                        else:
                            max_index = int(sstart)
                            min_index = int(send)

                        if (flankingnts + 1) <= min_index:
                            if flankingnts <= (len(accession_sequence) - max_index):
                                flanked_sseq = accession_sequence[ min_index -1 -flankingnts : max_index + flankingnts ]
                            else:
                                flanked_sseq = accession_sequence[ min_index -1 -flankingnts : ]
                        else:
                            flanked_sseq = accession_sequence[ : max_index + flankingnts ]

                        outblast_fh.write('>' + str(seq_number) + '_' + sseqid + ' # ' + ' # '.join(values[0:13]) + '\n' + flanked_sseq + '\n' )                                                                                                                                        
                        
                        seq_number += 1

                ## Create local BLAST database for BLAST results after locus detection in sample using reference allele
                db_name = os.path.join(blast_results_db_directory, sample_name)
                if not create_blastdb(path_to_blast_seq, db_name, 'nucl', logger):
                    print('Error when creating the blastdb for blast results file for locus %s at sample %s. Check log file for more information. \n ', core_name, sample_name)
                    return False

                # Path to local BLAST database for BLAST results after locus detection in sample using reference allele
                locus_blast_db_name = os.path.join(blast_results_db_directory, sample_name, os.path.basename(core_name) + '_blast', os.path.basename(core_name) + '_blast')


                # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
                # BLAST result sequences VS ALL alleles in locus BLAST for allele identification detection in sample  #
                # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
                cline = NcbiblastnCommandline(db=locus_blast_db_name, evalue=evalue, perc_identity=perc_identity_loc, reward=reward, penalty=penalty, gapopen=gapopen, gapextend=gapextend, outfmt = blast_parameters, max_target_seqs=max_target_seqs, max_hsps=max_hsps, num_threads=num_threads, query=locus_alleles_path)

                out, err = cline()
                out_lines = out.splitlines()

                allele_found = {} # To keep filtered BLAST results

                print("Comprobando si algún resultado tiene ID 100")
                
                ## Check if there is any BLAST result with ID = 100 ##
                for line in out_lines:

                    values = line.split('\t')
                    pident = values[2]

                    if float(pident) == 100:

                        qseqid, sseqid, pident, qlen, s_length, mismatch, gapopen, evalue, bitscore, sstart, send, qstart, qend, sseq, qseq = values

                        # Parse core gene fasta file to get matching allele sequence and length
                        alleles_in_locus = list (SeqIO.parse(locus_alleles_path, "fasta")) ## comentar parse
                        for allele in alleles_in_locus : ## comentar parse
                            if allele.id == qseqid : ## comentar parse
                                break ## comentar parse
                        matching_allele_seq = str(allele.seq) ## comentar parse
                        matching_allele_length = len(matching_allele_seq) ## comentar parse

                        matching_allele_seq = alleles_in_locus_dict[core_name][qseqid] ## Accediendo a diccionario contigs. Alternativa a parser
                        matching_allele_length = len(matching_allele_seq) ## Accediend a diccionario contigs. Alternativa a parser

                        # Keep BLAST results with ID = 100 and same length as matching allele
                        if int(s_length) == matching_allele_length:
                            get_blast_results (core_name, values, records, allele_found, logger)

                # ·································································································································· #
                # NIPHEM (paralog) if there are multiple BLAST results with ID = 100 and same length as matching allele for this gene in this sample #
                # ·································································································································· #
                if len(allele_found) > 1:
                    print("Se encontraron múltiples resultados con ID 100: NIPHEM")

                    print("Muestra ", sample_name, '\n') ### printeando, borrar
                    print("Locus ", core_name, '\n') ### printeando, borrar
                    print('NIPHEM', '\n') ### printeando, borrar
                    print('allele_found: ', allele_found, '\n') ### printeando, borrar

                    # Keep NIPHEM info
                    paralog_exact_tag(sample_name, core_name, 'NIPHEM', schema_quality, matching_genes_dict, samples_matrix_dict, allele_found, paralog_dict, prodigal_report, prodigal_directory, blast_parameters, logger)

                    continue
                
                ## Check for possible paralogs with ID < 100 if there is only one BLAST result with ID = 100 and same length as matching allele
                elif len(allele_found) == 1 :
                    print("se encontró un resultado con ID 100", '\n')

                    for line in out_lines :
                        values = line.split('\t')

                        sseq_no_gaps = values[13].replace('-', '')
                        s_length_no_gaps = len(sseq_no_gaps)

                        # Keep BLAST result if its coverage is within min and max thresholds
                        if min_length_threshold <= s_length_no_gaps <= max_length_threshold:                            
                            get_blast_results (core_name, values, records, allele_found, logger)

                    print("allele_found tras búsqueda de NIPHs: ", allele_found, '\n')

                    # ································································ #
                    # EXACT MATCH if there is any paralog for this gene in this sample #
                    # ································································ #
                    if len(allele_found) == 1 :

                        paralog_exact_tag(sample_name, core_name, 'EXACT', schema_quality, matching_genes_dict, samples_matrix_dict, allele_found, exact_dict, prodigal_report, prodigal_directory, blast_parameters, logger)

                        print("No se encontraron parálogos del resultado con ID 100: EXACT MATCH", '\n')
                        print("Muestra ", sample_name, '\n') ### printeando, borrar
                        print("Match con alelo ", qseqid, "del locus ", core_name, '\n') ### printeando, borrar
                        print("EXACT MATCH", '\n') ### printeando, borrar
                        print("sseq encontrada en BLAST: ", sseq, '\n') ### printeando, borrar

                        continue
                    
                    # ··········································································· #
                    # NIPH if there there are paralogs with ID < 100 for this gene in this sample #
                    # ··········································································· #
                    else:
                        print("se encontraron parálogos del resultado con ID 100: NIPH", '\n')

                        paralog_exact_tag(sample_name, core_name, 'NIPH', schema_quality, matching_genes_dict, samples_matrix_dict, allele_found, paralog_dict, prodigal_report, prodigal_directory, blast_parameters, logger)

                        continue

                ## Look for the best BLAST result if there are no results with ID = 100 ##
                elif len(allele_found) == 0:
                    print("No se encontraron resultados con ID 100. Se busca el resultado con mejor bitscore cuyo coverage entre dentro del rango establecido", '\n')

                    bigger_bitscore_seq_values = []

                    for line in out_lines :
                        values = line.split('\t')

                        if  float(values[8]) > bigger_bitscore:
                            s_length_no_gaps = len(values[13].replace('-', ''))

                            # Keep BLAST result if its coverage is within min and max thresholds and its bitscore is bigger than the one previously kept
                            if min_coverage_threshold <= s_length_no_gaps <= max_coverage_threshold:
                                bigger_bitscore_seq_values = values
                                bigger_bitscore = float(bigger_bitscore_seq_values[8])
                    
                    print("bigger_bitscore_seq_values: ", bigger_bitscore_seq_values,'\n')


                    ## Check if best BLAST result out of coverage thresholds is a possible PLOT or LNF due to low coverage ##
                    #if len(allele_found) == 0:
                    if len(bigger_bitscore_seq_values) == 0:
                
                        print("No se encontraron resultados con ID menor a 100 con coverage dentro del rango establecido: LNF", '\n')

                        # Look for best bitscore BLAST result out of coverage thresholds to check possible PLOT or reporting LNF due to low coverage
                        for line in out_lines :
                            values = line.split('\t')

                            if  float(values[8]) > bigger_bitscore:
                                qseqid, sseqid, pident,  qlen, s_length, mismatch, gapopen, evalue, bitscore, sstart, send, qstart, qend, sseq, qseq = values
                                bigger_bitscore = float(bitscore)

                            s_length_no_gaps = len(sseq.replace('-', ''))


                            # Get contig sequence and length for best bitscore BLAST result ID 
                            for record in records: ## comentar parse
                                if record.id == sseqid : ## comentar parse
                                    break ## comentar parse
                            accession_sequence = record.seq ## comentar parse
                            length_sseqid = len(accession_sequence) ## comentar parse

                            accession_sequence = contigs_in_sample_dict[core_name][sseqid] ## Accediendo a diccionario contigs. Alternativa a parser
                            length_sseqid = len(accession_sequence) ## Accediendo a diccionario contigs. Alternativa a parser

                            # Check if best BLAST result out of coverage thresholds is a possible PLOT. If so, keep result info for later PLOT classification
                            if int(sstart) == length_sseqid or int(send) == length_sseqid or int(sstart) == 1 or int(send) == 1:
                                bigger_bitscore_seq_values = values

                            # ·············································································································································· #
                            # LNF if there are no BLAST results within coverage thresholds for this gene in this sample and best out threshold result is not a possible PLOT #
                            # ·············································································································································· #
                            else:

                                # Keep LNF info
                                lnf_tag(core_name, sample_name, alleles_in_locus_dict, samples_matrix_dict, lnf_dict, locus_alleles_path, qseqid, pident, s_length_no_gaps, '-', logger)


                    ## Keep result with bigger bitscore in allele_found dict and look for possible paralogs ##
                    if len(bigger_bitscore_seq_values) > 0:

                        print("Se ha encontrado secuencia con ID < 100 y con coverage dentro de rango. Se buscan posibles parálogos", '\n')

                        qseqid, sseqid, pident, qlen, s_length, mismatch, gapopen, evalue, bitscore, sstart, send, qstart, qend, sseq, qseq = bigger_bitscore_seq_values

                        get_blast_results (core_name, bigger_bitscore_seq_values, records, allele_found, logger)

                        # Possible paralogs search
                        for line in out_lines :
                            values = line.split('\t')

                            qseqid, sseqid, pident, qlen, s_length, mismatch, gapopen, evalue, bitscore, sstart, send, qstart, qend, sseq, qseq = values ### OBTENER ESTO DESPUÉS DEL SIGUIENTE IF Y SACAR AQUÍ SOLO SSEQ PARA S_LENGTH?
                            sseq_no_gaps = sseq.replace('-', '')
                            s_length_no_gaps = len(sseq_no_gaps) 

                            if min_length_threshold <= s_length_no_gaps <= max_length_threshold: # comprobando que la longitud de la secuencia encontrada sin gaps se encuentra dentro del rango esperado                            

                                get_blast_results (core_name, values, records, allele_found, logger)

                        # ····························································· #
                        # NIPH if there there are paralogs for this gene in this sample #
                        # ····························································· #
                        if len(allele_found) > 1 :

                            print("Se encontraron parálogos del resultado con ID menor a 100 con coverage dentro del rango establecido: NIPH", '\n')

                            paralog_exact_tag(sample_name, core_name, 'NIPH', schema_quality, matching_genes_dict, samples_matrix_dict, allele_found, paralog_dict, prodigal_report, prodigal_directory, blast_parameters, logger)

                            continue 

                        ## Continue classification if there are no paralogs ##
                        elif len(allele_found) == 1 :

                            print("No se encontraron parálogos del resultado con ID menor a 100 con coverage dentro del rango establecido.", '\n')

                            allele_id = str(list(allele_found.keys())[0])
                            qseqid, sseqid, pident,  qlen, s_length, mismatch, gapopen, evalue, bitscore, sstart, send, qstart, qend, sseq, qseq = allele_found[allele_id]
                            
                            sseq_no_gaps = sseq.replace('-', '')
                            s_length_no_gaps = len(sseq_no_gaps) # --> usando en LNF, mirar dónde sacar para no ser repetitivo

                            # Get matching allele quality
                            allele_quality = schema_quality[core_name][qseqid]

                            # Parse core gene fasta file to get matching allele sequence and length
                            alleles_in_locus = list (SeqIO.parse(locus_alleles_path, "fasta")) ## comentar parse
                            for allele in alleles_in_locus : ## comentar parse
                                if allele.id == qseqid : ## comentar parse
                                    break ## comentar parse
                            matching_allele_seq = allele.seq ## comentar parse
                            matching_allele_length = len(matching_allele_seq) ## comentar parse               

                            matching_allele_seq = alleles_in_locus_dict [core_name][qseqid] ## Accediendo a diccionario alelos. Alternativa a parser
                            matching_allele_length = len(matching_allele_seq) ## Accediendo a diccionario alelos. Alternativa a parser

                            # Get contig sequence and length for ID found in BLAST
                            for record in records: ## comentar parse
                                if record.id == sseqid : ## comentar parse
                                    break ## comentar parse
                            accession_sequence = record.seq ## comentar parse
                            length_sseqid = len(accession_sequence) ## comentar parse

                            accession_sequence = contigs_in_sample_dict ## Accediendo a diccionario contigs. Alternativa a parser
                            length_sseqid = len(accession_sequence) ## Accediendo a diccionario contigs. Alternativa a parser

                            # ········································································································· #
                            # PLOT if found sequence is shorter than matching allele and it is located on the edge of the sample contig #
                            # ········································································································· #
                            if int(sstart) == length_sseqid or int(send) == length_sseqid or int(sstart) == 1 or int(send) == 1:
                                if int(s_length) < matching_allele_length:

                                    ### No sé si sacar la sec de prodigal para PLOT o quitarlo
                                    # Get prodigal predicted sequence if matching allele quality is "bad quality"
                                    if 'bad_quality' in allele_quality: 
                                        #complete_predicted_seq = str(get_prodigal_sequence(sstart, send, sseqid, prodigal_directory, sample_name))
                                        complete_predicted_seq, start_prodigal, end_prodigal = get_prodigal_sequence(sseq_no_gaps, sseqid, prodigal_directory, sample_name, blast_parameters, logger)

                                        # Keep info for prodigal report
                                        prodigal_report.append([core_name, sample_name, qseqid, 'PLOT', sstart, send, start_prodigal, end_prodigal, sseq_no_gaps, complete_predicted_seq])

                                    else:
                                        complete_predicted_seq = '-'
                                        start_prodigal = '-'
                                        end_prodigal = '-'

                                    # Keep PLOT info
                                    inf_asm_alm_tag(core_name, sample_name, 'PLOT', allele_found[allele_id], allele_quality, '-', matching_allele_length, '-', plot_dict, samples_matrix_dict, matching_genes_dict, prodigal_report, start_prodigal, end_prodigal, complete_predicted_seq, logger)

                                    print("Muestra ", sample_name, '\n') ### printeando, borrar
                                    print("Match con alelo ", qseqid, "del locus ", core_name, '\n') ### printeando, borrar
                                    print("PLOT", '\n') ### printeando, borrar
                                    print("plot_dict: ", plot_dict, '\n', '\n') ### printeando, borrar

                                    continue 

                            # * * * * * * * * * * * * * * * * * * * * #
                            # Search for complete final new sequence  # 
                            # * * * * * * * * * * * * * * * * * * * * #

                            print("Muestra ", sample_name, '\n') ### printeando, borrar
                            print("Match con alelo ", qseqid, "del locus ", core_name, '\n') ### printeando, borrar

                            ## Get Prodigal predicted sequence ##
                            #complete_predicted_seq = str(get_prodigal_sequence(sstart, send, sseqid, prodigal_directory, sample_name))
                            complete_predicted_seq, start_prodigal, end_prodigal = get_prodigal_sequence(sseq_no_gaps, sseqid, prodigal_directory, sample_name, blast_parameters, logger)

                            print("secuencia prodigal predicha: ", complete_predicted_seq, '\n', '\n', '\n')


                            ## Search for new codon stop using contig sequence info ##

                            # Check matching allele sequence direction
                            query_direction = check_sequence_order(matching_allele_seq, logger)

                            print("query_direction: ", query_direction, '\n')

                            # Get extended BLAST sequence for stop codon search
                            if query_direction == 'reverse':
                                if int(send) > int (sstart):
                                    sample_gene_sequence = accession_sequence[ : int(send) ]
                                    sample_gene_sequence = sample_gene_sequence.reverse_complement()
                                else:
                                    sample_gene_sequence = accession_sequence[ int(send) -1 : ]

                            else:
                                if int(sstart) > int (send):
                                    sample_gene_sequence = accession_sequence[ :  int(sstart) ]
                                    sample_gene_sequence = sample_gene_sequence.reverse_complement()
                                else:
                                    sample_gene_sequence = accession_sequence[ int(sstart) -1 : ]

                            # Get new stop codon index
                            stop_index = get_stop_codon_index(sample_gene_sequence) 


                            ## Classification of final new sequence if it is found ##
                            if stop_index != False:
                                new_sequence_length = stop_index +3
                                new_sseq = str(sample_gene_sequence[0:new_sequence_length])

                                #########################################################################################################################
                                ### c/m: introducido para determinar qué umbral de coverage poner. TEMPORAL
                                new_sseq_coverage = new_sequence_length/matching_allele_length ### introduciendo coverage new_sseq /// debería ser con respecto a la media?
                                
                                if new_sseq_coverage < 1:
                                    shorter_seq_coverage.append([core_name, sample_name, str(matching_allele_length), str(new_sequence_length), str(schema_statistics[core_name][1]), str(new_sseq_coverage), str(new_sequence_length/schema_statistics[core_name][1])])
                                elif new_sseq_coverage > 1:
                                    longer_seq_coverage.append([core_name, sample_name, str(matching_allele_length), str(new_sequence_length), str(schema_statistics[core_name][1]), str(new_sseq_coverage), str(new_sequence_length/schema_statistics[core_name][1])])
                                elif new_sseq_coverage == 1:
                                    equal_seq_coverage.append([core_name, sample_name, str(matching_allele_length), str(new_sequence_length), str(schema_statistics[core_name][1]), str(new_sseq_coverage), str(new_sequence_length/schema_statistics[core_name][1])])
                                #########################################################################################################################

                                # Get and keep SNP and DNA and protein alignment
                                keep_snp_alignment_info(new_sseq, matching_allele_seq, qseqid, query_direction, core_name, sample_name, snp_dict, match_alignment_dict, protein_dict, logger)

                                print("new_sseq encontrada por Taranis: ", new_sseq, '\n', '\n')
                                print("sseq encontrada por BLAST: ", sseq, '\n', '\n')
                                print("Inicio de clasificación de la secuencia", '\n', '\n')

                                # ····································································································· #
                                # INF if final new sequence length is within min and max length thresholds for this gene in this sample #
                                # ····································································································· #
                                if min_length_threshold <= new_sequence_length <= max_length_threshold:
                                    ### duda: incluir cuando INF == new_sequence_length? --> Si se toma el rango de longitud por defecto, "media + SD" para clasificar puede que un resultado que tiene la
                                    ### misma longitud que el alelo con el que hizo match se clasifique como ASM o ALM porque se trate del alelo más largo o más corto del locus, por lo que su longitud
                                    ### no entraría dentro del margen permitido por la SD y aunque sea de igual tamaño que el alelo que hizo match no se consideraría INF por esto.
                                                                    
                                    print("INF", '\n')
                                    print("allele_found: ", allele_found, '\n')
                                    print("allele: ", allele, '\n')

                                    inf_asm_alm_tag(core_name, sample_name, 'INF', allele_found[allele_id], allele_quality, new_sseq, matching_allele_length, inferred_alleles_dict, inf_dict, samples_matrix_dict, matching_genes_dict, prodigal_report, start_prodigal, end_prodigal, complete_predicted_seq, logger) ### introducido start_prodigal, end_prodigal, complete_predicted_seq, prodigal_report como argumento a inf_asm_alm_tag para report prodigal, temporal

                                # ············································································································································ #
                                # ASM if final new sequence length is under min length threshold but its coverage is above min coverage threshold for this gene in this sample #
                                # ············································································································································ #
                                elif min_coverage_threshold <= new_sequence_length < min_length_threshold:
                                    # math.floor(schema_statistics[core_name][1] - (schema_statistics[core_name][1]*0.5))

                                    print("ASM")

                                    # Keep ASM info
                                    inf_asm_alm_tag(core_name, sample_name, 'ASM', allele_found[allele_id], allele_quality, new_sseq, matching_allele_length, asm_dict, list_asm, samples_matrix_dict, matching_genes_dict, prodigal_report, start_prodigal, end_prodigal, complete_predicted_seq, logger)

                                # ············································································································································ #
                                # ALM if final new sequence length is above max length threshold but its coverage is under max coverage threshold for this gene in this sample #
                                # ············································································································································ #
                                elif max_length_threshold < new_sequence_length <= max_coverage_threshold:
                                    #math.ceil(schema_statistics[core_name][1] + (schema_statistics[core_name][1]*0.5 
                                    
                                    print("ALM", '\n')

                                    # Keep ASM info
                                    inf_asm_alm_tag(core_name, sample_name, 'ALM', allele_found[allele_id], allele_quality, new_sseq, matching_allele_length, alm_dict, list_alm, samples_matrix_dict, matching_genes_dict, prodigal_report, start_prodigal, end_prodigal, complete_predicted_seq, logger) ### introducido start_prodigal, end_prodigal, complete_predicted_seq, prodigal_report como argumento a inf_asm_alm_tag para report prodigal, temporal

                                # ························································································· #
                                # LNF if final new sequence coverage is not within thresholds for this gene in this sample #
                                # ························································································· #
                                else: 
                                    print("LNF por bajo coverage new_sseq", '\n')

                                    # Keep LNF info
                                    lnf_tag(core_name, sample_name, alleles_in_locus_dict, samples_matrix_dict, lnf_dict, locus_alleles_path, qseqid, pident, s_length_no_gaps, new_sequence_length, logger)

                                print("Fin de la clasificación", '\n')
                    
                            # ········································ #
                            # ERROR if final new sequence is not found #
                            # ········································ #
                            else:
                                logger.error('ERROR : Stop codon was not found for the core %s and the sample %s', core_name, sample_name)
                                samples_matrix_dict[sample_name].append('ERROR not stop codon')
                                if not sseqid in matching_genes_dict[sample_name] :
                                    matching_genes_dict[sample_name][sseqid] = []
                                if sstart > send :
                                    matching_genes_dict[sample_name][sseqid].append([core_name, sstart,send,'-', 'ERROR'])
                                else:
                                    matching_genes_dict[sample_name][sseqid].append([core_name, sstart,send,'+', 'ERROR'])
 
    ## Save results and create reports
    save_results (full_gene_list, samples_matrix_dict, exact_dict, paralog_dict, inf_dict, plot_dict, matching_genes_dict, list_asm, list_alm, lnf_dict, snp_dict, match_alignment, protein_dict, prodigal_report, shorter_seq_coverage, longer_seq_coverage, equal_seq_coverage, shorter_blast_seq_coverage, longer_blast_seq_coverage, equal_blast_seq_coverage)

    return True


def processing_allele_calling (arguments) :
    '''
    Description:
        This is the main function for allele calling.
        With the support of additional functions it will create the output files
        with the summary report.
    Input:
        arguments   # Input arguments given on command line 
    Functions:
        
    Variables:
        run_metric_processed # True or False if there are some rows in
                            StatsRunSummary for this run
    Return:
        experiment_name if the run is updated. Empty if not
    '''

    #logger = logging.getLogger(__name__)
    #logger.debug ('Starting function check_run_metrics_processed')
    start_time = datetime.now()
    print('Start the execution at :', start_time )
    
    # Open log file
    logger = open_log ('taranis_wgMLST.log')
    print('Checking the pre-requisites./n')
    
    ############################################################
    ## Check additional programs are installed in your system ##
    ############################################################
    pre_requisites_list = [['blastp', '2.5'], ['makeblastdb', '2.5']]
    if not check_prerequisites (pre_requisites_list, logger):
        print ('your system does not fulfill the pre-requistes to run the script ')
        exit(0)

    ######################################################
    ## Check that given directories contain fasta files ##
    ######################################################
    print('Validating schema fasta files in ' , arguments.coregenedir , '\n')
    valid_core_gene_files = get_fasta_file_list(arguments.coregenedir, logger)
    if not valid_core_gene_files :
        print ('There are not valid fasta files in ',  arguments.coregenedir , ' directory. Check log file for more information ')
        exit(0)

    print('Validating reference alleles fasta files in ' , arguments.refalleles , '\n')
    valid_reference_alleles_files = get_fasta_file_list(arguments.refalleles, logger)
    if not valid_reference_alleles_files :
        print ('There are not valid reference alleles fasta files in ',  arguments.refalleles, ' directory. Check log file for more information ')
        exit(0)

    print('Validating sample fasta files in ' , arguments.inputdir , '\n')
    valid_sample_files = get_fasta_file_list(arguments.inputdir, logger)
    if not valid_sample_files :
        print ('There are not valid fasta files in ',  arguments.inputdir , ' directory. Check log file for more information ')
        exit(0)

    #################################
    ## Prepare the coreMLST schema ##
    #################################
    tmp_core_gene_dir = os.path.join(arguments.outputdir,'tmp','cgMLST')
    try:
        os.makedirs(tmp_core_gene_dir)
    except:
        logger.info('Deleting the temporary directory for a previous execution without cleaning up')
        shutil.rmtree(os.path.join(arguments.outputdir, 'tmp'))
        try:
            os.makedirs(tmp_core_gene_dir)
            logger.info ( 'Temporary folder %s  has been created again', tmp_core_gene_dir)
        except:
            logger.info('Unable to create again the temporary directory %s', tmp_core_gene_dir)
            print('Cannot create temporary directory on ', tmp_core_gene_dir)
            exit(0)

    core_gene_list_files, alleles_in_locus_dict, annotation_core_dict, schema_variability, schema_statistics, schema_quality = prepare_core_gene (valid_core_gene_files, tmp_core_gene_dir, arguments.refalleles, arguments.outputdir, logger) ### cambiando/modificando: he quitado core_first_alleles_files, relacionado con el primer alelo del locus que se tomaba antes como referencia, he añadido schema_quality
    if not core_gene_list_files :
        print('There is an error while processing the schema preparation phase. Check the log file to get more information \n')
        logger.info('Deleting the temporary directory to clean up the temporary files created')
        shutil.rmtree(os.path.join(arguments.outputdir, 'tmp'))
        exit(0)

    ###############################
    ## Prepare the samples files ##
    ###############################
    tmp_samples_dir = os.path.join(arguments.outputdir,'tmp','samples')
    try:
        os.makedirs(tmp_samples_dir)
    except:
        logger.info('Deleting the temporary directory for a previous execution without properly cleaning up')
        shutil.rmtree(tmp_samples_dir)
        try:
            os.makedirs(tmp_samples_dir)
            logger.info('Temporary folder %s  has been created again', tmp_samples_dir)
        except:
            logger.info('Unable to create again the temporary directory %s', tmp_samples_dir)
            shutil.rmtree(os.path.join(arguments.outputdir, 'tmp'))
            logger.info('Cleaned up temporary directory ', )
            print('Cannot create temporary directory on ', tmp_samples_dir, 'Check the log file to get more information \n')
            exit(0)
    
    sample_list_files, contigs_in_sample_dict = prepare_samples(valid_sample_files, tmp_samples_dir, arguments.refgenome, logger)
    if not sample_list_files :
        print('There is an error while processing the saving temporary files. Check the log file to get more information \n')
        logger.info('Deleting the temporary directory to clean up the temporary files created')
        shutil.rmtree(os.path.join(arguments.outputdir, 'tmp'))
        exit(0)

    query_directory = arguments.coregenedir
    reference_alleles_directory = arguments.refalleles
    blast_db_directory = os.path.join(tmp_samples_dir,'blastdb')
    prodigal_directory = os.path.join(tmp_samples_dir,'prodigal') 
    blast_results_seq_directory = os.path.join(tmp_samples_dir,'blast_results', 'blast_results_seq')  ### path a directorio donde guardar secuencias encontradas tras blast con alelo de referencia
    blast_results_db_directory = os.path.join(tmp_samples_dir,'blast_results', 'blast_results_db') ### path a directorio donde guardar db de secuencias encontradas tras blast con alelo de referencia

    if not allele_call_nucleotides(core_gene_list_files, sample_list_files, alleles_in_locus_dict, contigs_in_sample_dict, query_directory, reference_alleles_directory, blast_db_directory, prodigal_directory, blast_results_seq_directory, blast_results_db_directory, arguments.inputdir, arguments.outputdir,  int(arguments.cpus), arguments.percentlength, arguments.coverage, float(arguments.evalue), float(arguments.perc_identity_ref), float(arguments.perc_identity_loc), int(arguments.reward), int(arguments.penalty), int(arguments.gapopen), int(arguments.gapextend), int(arguments.max_target_seqs), int(arguments.max_hsps), int(arguments.num_threads), int(arguments.flankingnts), schema_variability, schema_statistics, schema_quality, annotation_core_dict, logger): ### CAMBIANDO/MODIFICANDO: He añadido schema_statistics, path a prodigal, prodigal training y schema_quality        
        print('There is an error while processing the allele calling. Check the log file to get more information \n')
        exit(0)

    ################################
    ## Create the distance matrix ##
    ################################
    try:        
        print ('Creating matrix distance\n')
        create_distance_matrix(arguments.outputdir, 'result_for_tree_diagram.tsv')
    except:
        print('There was an error when creating distance matrix\n')
    end_time = datetime.now()
    print('completed execution at :', end_time )

    return True



