#!/usr/bin/env python3
# -coregenedir /srv/project_wgmlst/seqSphere_listeria_cgMLST_test/targets/ -inputdir /srv/project_wgmlst/samples_listeria_test  -outputdir /srv/results
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

import math ### cambiando/modificando: import math para obtención length_trhesholds
import numpy as np ### import numpy para obtener alelo de referencia para cada locus
import enchant ### import enchant para obtener alelo de referencia para cada locus
from Levenshtein import distance ### import enchant para obtener alelo de referencia para cada locus

def check_prerequisites (logger): 
    pre_requisite_list = [['blastp', '2.5'], ['makeblastdb' , '2.5']]
    # check if blast is installed and has the minimum version
    for program, version in pre_requisite_list :
        if not check_program_is_exec_version (program , version, logger):
            return False
    return True

"""
def write_first_allele_seq(file_sequence, store_dir, logger):
    #with open (file_name, 'r' ) as fh :
    #seq_record = SeqIO.parse(open(file_name), "genbank").next()
    first_allele_directory = 'first_alleles'
    # split file_sequence into directory and filename
    f_name = os.path.basename(file_sequence)
    full_path_first_allele = os.path.join(store_dir, first_allele_directory)
    if not os.path.exists(full_path_first_allele):
        try:
            os.makedirs(full_path_first_allele)
            logger.info('Directory %s has been created', full_path_first_allele)
        except:
            print ('Cannot create the directory ', full_path_first_allele)
            logger.info('Directory %s cannot be created', full_path_first_allele)
            exit (0)
    first_record = SeqIO.parse(file_sequence, "fasta").__next__()
    # build the fasta file name to store under first_allele_firectory
    fasta_file = os.path.join(full_path_first_allele, f_name)
    SeqIO.write(first_record, fasta_file, "fasta")
    return fasta_file
"""

def create_blastdb (file_name, db_name,db_type, logger ):
    f_name = os.path.basename(file_name).split('.')
    db_dir = os.path.join(db_name,f_name[0])
    output_blast_dir = os.path.join(db_dir, f_name[0])
    if not os.path.exists(db_dir):
        try:
            os.makedirs(db_dir)
            logger.debug(' Created local blast directory for Core Gene %s', f_name[0])
        except:
            logger.info('Cannot create directory for local blast database on Core Gene file %s' , f_name[0])
            print ('Error when creating the directory %s for blastdb. ', db_dir)
            exit(0)

        blast_command = ['makeblastdb' , '-in' , file_name , '-parse_seqids', '-dbtype',  db_type, '-out' , output_blast_dir]
        blast_result = subprocess.run(blast_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if blast_result.stderr:
            logger.error('cannot create blast db for %s ', f_name[0])
            logger.error('makeblastdb returning error code %s', blast_result.stderr)
            return False
    else:
        logger.info('Skeeping the blastdb creation for %s, as it is already exists', f_name[0])
    return True

def check_blast (reference_allele, sample_files, db_name, logger) : ### No se está usando
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

def parsing_fasta_file_to_dict (fasta_file, logger): 
    fasta_dict = {}
    fasta_dict_ordered = {}
    for contig in SeqIO.parse(fasta_file, "fasta", generic_dna):
        #fasta_dict[int(contig.id)] = str(contig.seq.upper()) ### cambiando/modificando: int en contig.id en lugar de str para ordenar?
        fasta_dict[str(contig.id)] = str(contig.seq.upper())
    logger.debug('file %s parsed to dictionary', fasta_file)

    for key in sorted(list(fasta_dict.keys())): ### cambiando/modificando: ordenando diccionario de alelos en función de id
        fasta_dict_ordered[key] = fasta_dict[key]
  #  print("fasta_dict_ordered.keys(): ", fasta_dict_ordered.keys())
    return fasta_dict_ordered

"""
def check_core_gene_quality(alleles_seqs_dict, logger): ### MODIFICANDO/CAMBIANDO: Función check_core_gene_quality añadida para checkear calidad de los alelos de cada locus del esquema
                                                        ### logger?

    ### logger.info('check quality of locus %s', fasta_file)

    locus_alleles_seqs = list(alleles_seqs_dict.values()) # secuencias de los alelos del locus
    locus_alleles_ids = list(alleles_seqs_dict.keys()) # ids secuencias de los alelos del locus

    #print("locus_alleles_seqs: ", locus_alleles_seqs, '\n') ### print, borrar
    #print("locus_alleles_ids: ", locus_alleles_ids, '\n') ### print, borrar

    start_codons_forward = ['ATG', 'ATA', 'ATT', 'GTG', 'TTG', 'CTG'] ### duda: tener en cuenta codones de inico no clásicos? (prodigal no los considera)
    start_codons_reverse = ['CAT', 'TAT', 'AAT', 'CAC', 'CAA', 'CAG']

    stop_codons_forward = ['TAA', 'TAG', 'TGA']
    stop_codons_reverse = ['TTA', 'CTA', 'TCA']

    locus_quality = {}
    locus_quality["good_quality"] = []
    locus_quality["bad_quality"] = {"no_start": [], "no_stop": [], "no_start_stop": [], "multiple_stop": []}
 

    for allele_index in range(len(locus_alleles_seqs)):
        if locus_alleles_seqs[allele_index][0:3] in start_codons_forward or locus_alleles_seqs[allele_index][-3:] in start_codons_reverse:
            if locus_alleles_seqs[allele_index][-3:] in stop_codons_forward or locus_alleles_seqs[allele_index][0:3] in stop_codons_reverse: ### si tiene codón de stop y codón de inicio
            ### duda: Añadir condición de que sea múltiplo de 3 la secuencia para así saber que el stop y start que presenta el alelo entran dentro del marco de lectura de la secuencia?

                ### Buscando codón de stop para checkear que el primer codón de stop de la secuencia se corresponde con el codón final de la secuencia, de modo que no presenta más de un codón stop
                sequence_order = check_sequence_order(locus_alleles_seqs[allele_index], logger)
                if sequence_order == "reverse":
                    allele_sequence = reverse_complement(locus_alleles_seqs[allele_index])
                else:
                    allele_sequence = locus_alleles_seqs[allele_index]
                stop_index = get_stop_codon_index(allele_sequence)

                if stop_index < (len(allele_sequence) - 3): ### si tiene codón start y stop pero tiene más de un codón stop (-3 --> 1 por índice python y 2 por las 2 bases restantes del codón)
                    locus_quality["bad_quality"]["multiple_stop"].append(locus_alleles_ids[allele_index])
                else: ### si tiene codón start y stop y un solo codón stop
                    locus_quality["good_quality"].append(locus_alleles_ids[allele_index])

            else: ### si tiene codón start pero no stop
                locus_quality["bad_quality"]["no_stop"].append(locus_alleles_ids[allele_index])
        else: ### Si no tiene start
            if locus_alleles_seqs[allele_index][-3:] in stop_codons_forward or locus_alleles_seqs[allele_index][0:3] in stop_codons_reverse: ### si no tiene start pero sí stop
                locus_quality["bad_quality"]["no_start"].append(locus_alleles_ids[allele_index])
            else: ### Si no tiene start ni stop
                locus_quality["bad_quality"]["no_start_stop"].append(locus_alleles_ids[allele_index])

    ### En la función prepare_core_gene el resultado de esta función para cada locus se guardará en el diccionario schema_quality y devolverá el diccionario como resultado de llamar a la función prepare_core_gene, al igual que devuelve schema_statistics, etc
    ### Se le introducirá este diccionario como argumento a la función allele_call_nucleotides, y no es necesario guardarlo como archivo, ya que se le introduce directamente cogiéndolo del output de llamar a prepare_core_gene (al igual que schema_statistics, etc)
    return locus_quality
"""

### Reformando diccionario check_core_gene_quality
### introduciendo path al fasta del locus para parsear en lugar de diccionario del locus generado previamente
### Tengo que cambiar alleles_seqs_dict por el path del locus y así parsear con seqIO sacando las seqs y los ids
def check_core_gene_quality(fasta_file_path, logger): ### MODIFICANDO/CAMBIANDO: Función check_core_gene_quality añadida para checkear calidad de los alelos de cada locus del esquema
                                                        ### logger?

    ### logger.info('check quality of locus %s', fasta_file)

    #locus_alleles_seqs = list(alleles_seqs_dict.values()) # secuencias de los alelos del locus
    #locus_alleles_ids = list(alleles_seqs_dict.keys()) # ids secuencias de los alelos del locus

    #print("locus_alleles_seqs: ", locus_alleles_seqs, '\n') ### print, borrar
    #print("locus_alleles_ids: ", locus_alleles_ids, '\n') ### print, borrar

    start_codons_forward = ['ATG', 'ATA', 'ATT', 'GTG', 'TTG', 'CTG'] ### duda: tener en cuenta codones de inico no clásicos? (prodigal no los considera)
    start_codons_reverse = ['CAT', 'TAT', 'AAT', 'CAC', 'CAA', 'CAG']

    stop_codons_forward = ['TAA', 'TAG', 'TGA']
    stop_codons_reverse = ['TTA', 'CTA', 'TCA']

    locus_quality = {}
  #  locus_quality["good_quality"] = []
  #  locus_quality["bad_quality"] = {"no_start": [], "no_stop": [], "no_start_stop": [], "multiple_stop": []}
 


    alleles_in_locus = list (SeqIO.parse(fasta_file_path, "fasta"))
    for allele_item in alleles_in_locus :
        print("type allele_item.seq: ", type(allele_item.seq), '\n')
        print("primeras 3 bases de allele_item.seq: ", allele_item.seq[0:3], '\n')
        print("primeras 3 bases de str(allele_item.seq): ", str(allele_item.seq[0:3]), '\n')
        print("allele_item.seq.reverse_complement(): ", allele_item.seq.reverse_complement(), '\n')

  #  for allele_index in range(len(locus_alleles_seqs)):
        if allele_item.seq[0:3] in start_codons_forward or allele_item.seq[-3:] in start_codons_reverse:
            if allele_item.seq[-3:] in stop_codons_forward or allele_item.seq[0:3] in stop_codons_reverse: ### si tiene codón de stop y codón de inicio
            ### duda: Añadir condición de que sea múltiplo de 3 la secuencia para así saber que el stop y start que presenta el alelo entran dentro del marco de lectura de la secuencia?

                ### Buscando codón de stop para checkear que el primer codón de stop de la secuencia se corresponde con el codón final de la secuencia, de modo que no presenta más de un codón stop
                sequence_order = check_sequence_order(allele_item.seq, logger)
                if sequence_order == "reverse":
                    #allele_sequence = reverse_complement(str(allele_item.seq)) ####### utilizar Seq y reverse.complement y a continuación convertirlo a str
                    allele_sequence = str(allele_item.seq.reverse_complement())
                    print("allele_sequence: ", allele_sequence, '\n')
                else:
                    allele_sequence = str(allele_item.seq)
                stop_index = get_stop_codon_index(allele_sequence)

                if stop_index < (len(allele_sequence) - 3): ### si tiene codón start y stop pero tiene más de un codón stop (-3 --> 1 por índice python y 2 por las 2 bases restantes del codón)
                   # locus_quality["bad_quality"]["multiple_stop"].append(locus_alleles_ids[allele_index])
                    locus_quality[str(allele_item.id)] = 'bad_quality: multiple_stop'
                else: ### si tiene codón start y stop y un solo codón stop
                  #  locus_quality["good_quality"].append(locus_alleles_ids[allele_index])
                    locus_quality[str(allele_item.id)] = 'good_quality'

            else: ### si tiene codón start pero no stop
                #locus_quality["bad_quality"]["no_stop"].append(locus_alleles_ids[allele_index])
                locus_quality[str(allele_item.id)] = 'bad_quality: no_stop'
        else: ### Si no tiene start
            if allele_item.seq[-3:] in stop_codons_forward or allele_item.seq[0:3] in stop_codons_reverse: ### si no tiene start pero sí stop
                #locus_quality["bad_quality"]["no_start"].append(locus_alleles_ids[allele_index])
                locus_quality[str(allele_item.id)] = 'bad_quality: no_start'
            else: ### Si no tiene start ni stop
                #locus_quality["bad_quality"]["no_start_stop"].append(locus_alleles_ids[allele_index])
                locus_quality[str(allele_item.id)] = 'bad_quality: no_start_stop'

    ### En la función prepare_core_gene el resultado de esta función para cada locus se guardará en el diccionario schema_quality y devolverá el diccionario como resultado de llamar a la función prepare_core_gene, al igual que devuelve schema_statistics, etc
    ### Se le introducirá este diccionario como argumento a la función allele_call_nucleotides, y no es necesario guardarlo como archivo, ya que se le introduce directamente cogiéndolo del output de llamar a prepare_core_gene (al igual que schema_statistics, etc)
    return locus_quality




"""
### Reformando diccionario check_core_gene_quality
### Tengo que cambiar alleles_seqs_dict por el path del locus y así parsear con seqIO sacando las seqs y los ids
def check_core_gene_quality(alleles_seqs_dict, logger): ### MODIFICANDO/CAMBIANDO: Función check_core_gene_quality añadida para checkear calidad de los alelos de cada locus del esquema
                                                        ### logger?

    ### logger.info('check quality of locus %s', fasta_file)

    locus_alleles_seqs = list(alleles_seqs_dict.values()) # secuencias de los alelos del locus
    locus_alleles_ids = list(alleles_seqs_dict.keys()) # ids secuencias de los alelos del locus

    #print("locus_alleles_seqs: ", locus_alleles_seqs, '\n') ### print, borrar
    #print("locus_alleles_ids: ", locus_alleles_ids, '\n') ### print, borrar

    start_codons_forward = ['ATG', 'ATA', 'ATT', 'GTG', 'TTG', 'CTG'] ### duda: tener en cuenta codones de inico no clásicos? (prodigal no los considera)
    start_codons_reverse = ['CAT', 'TAT', 'AAT', 'CAC', 'CAA', 'CAG']

    stop_codons_forward = ['TAA', 'TAG', 'TGA']
    stop_codons_reverse = ['TTA', 'CTA', 'TCA']

    locus_quality = {}
  #  locus_quality["good_quality"] = []
  #  locus_quality["bad_quality"] = {"no_start": [], "no_stop": [], "no_start_stop": [], "multiple_stop": []}
 

    for allele_index in range(len(locus_alleles_seqs)):
        if locus_alleles_seqs[allele_index][0:3] in start_codons_forward or locus_alleles_seqs[allele_index][-3:] in start_codons_reverse:
            if locus_alleles_seqs[allele_index][-3:] in stop_codons_forward or locus_alleles_seqs[allele_index][0:3] in stop_codons_reverse: ### si tiene codón de stop y codón de inicio
            ### duda: Añadir condición de que sea múltiplo de 3 la secuencia para así saber que el stop y start que presenta el alelo entran dentro del marco de lectura de la secuencia?

                ### Buscando codón de stop para checkear que el primer codón de stop de la secuencia se corresponde con el codón final de la secuencia, de modo que no presenta más de un codón stop
                sequence_order = check_sequence_order(locus_alleles_seqs[allele_index], logger)
                if sequence_order == "reverse":
                    allele_sequence = reverse_complement(locus_alleles_seqs[allele_index])
                else:
                    allele_sequence = locus_alleles_seqs[allele_index]
                stop_index = get_stop_codon_index(allele_sequence)

                if stop_index < (len(allele_sequence) - 3): ### si tiene codón start y stop pero tiene más de un codón stop (-3 --> 1 por índice python y 2 por las 2 bases restantes del codón)
                   # locus_quality["bad_quality"]["multiple_stop"].append(locus_alleles_ids[allele_index])
                    locus_quality[locus_alleles_ids[allele_index]] = 'bad_quality: multiple_stop'
                else: ### si tiene codón start y stop y un solo codón stop
                  #  locus_quality["good_quality"].append(locus_alleles_ids[allele_index])
                    locus_quality[locus_alleles_ids[allele_index]] = 'good_quality'

            else: ### si tiene codón start pero no stop
                #locus_quality["bad_quality"]["no_stop"].append(locus_alleles_ids[allele_index])
                locus_quality[locus_alleles_ids[allele_index]] = 'bad_quality: no_stop'
        else: ### Si no tiene start
            if locus_alleles_seqs[allele_index][-3:] in stop_codons_forward or locus_alleles_seqs[allele_index][0:3] in stop_codons_reverse: ### si no tiene start pero sí stop
                #locus_quality["bad_quality"]["no_start"].append(locus_alleles_ids[allele_index])
                locus_quality[locus_alleles_ids[allele_index]] = 'bad_quality: no_start'
            else: ### Si no tiene start ni stop
                #locus_quality["bad_quality"]["no_start_stop"].append(locus_alleles_ids[allele_index])
                locus_quality[locus_alleles_ids[allele_index]] = 'bad_quality: no_start_stop'

    ### En la función prepare_core_gene el resultado de esta función para cada locus se guardará en el diccionario schema_quality y devolverá el diccionario como resultado de llamar a la función prepare_core_gene, al igual que devuelve schema_statistics, etc
    ### Se le introducirá este diccionario como argumento a la función allele_call_nucleotides, y no es necesario guardarlo como archivo, ya que se le introduce directamente cogiéndolo del output de llamar a prepare_core_gene (al igual que schema_statistics, etc)
    return locus_quality
"""

"""
def get_reference_allele(alleles_seqs_dict, locus_quality, file_sequence, store_dir, logger): ### MODIFICANDO/CAMBIANDO: Función reference_allele añadida para encontrar el alelo de referencia
                                                                                                ### PROBANDO obtención de alelo de referencia con función distance de paquete Levenshtein
                                                                                                ### logger?

    ### logger.info('search reference allele')
    
    # alleles_seqs_dict --> diccionario generado previamente key:id alelo - value:seq alelo

    locus_alleles_seqs = list(alleles_seqs_dict.values()) # secuencias de los alelos del locus
    locus_alleles_ids = list(alleles_seqs_dict.keys()) # ids de los alelos del locus

    alleles_matrix_distance = np.zeros((len(locus_alleles_seqs), len(locus_alleles_seqs)))

    for allele_index_1 in range(len(locus_alleles_seqs)):
        for allele_index_2 in range(len(locus_alleles_seqs)):
            levenshtein_distance = distance(locus_alleles_seqs[allele_index_1], locus_alleles_seqs[allele_index_2])
            alleles_matrix_distance[allele_index_1][allele_index_2] = levenshtein_distance

    allele_mean_distance = np.mean(alleles_matrix_distance, 0)


    # casos en los que todos los alelos del locus están imcompletos? -> si se ignoran los que están incompletos, se estarían ignorando todos los alelos del locus

    min_mean = max(allele_mean_distance) # inicializando con la mayor mean
    ref_allele_index = int()

    for mean_index in range(len(allele_mean_distance)):
        if allele_mean_distance[mean_index] < min_mean:
            if locus_alleles_ids[mean_index] in locus_quality['good_quality']:
                min_mean = allele_mean_distance[mean_index]
                ref_allele_index = mean_index

   # print("allele_mean_distance[ref_allele_index]: ", allele_mean_distance[ref_allele_index]) ### printeando, borrar

    min_index_allele = min(allele_mean_distance)
    index_alelo_min_mean = list(allele_mean_distance).index(min(allele_mean_distance))
    #print("min_index_allele: ", min_index_allele)
    #print("Alelo que tiene una mean diferente que el alelo que toma, checkeando si está incompleto: ", locus_alleles_seqs[index_alelo_min_mean]) ### printeando, borrar

    reference_allele = locus_alleles_seqs[ref_allele_index] # alelo representante del locus, cuya media de distancias al resto de alelos es la más pequeña (el que menos difiere del resto) y que está completo, con codón de stop y de start
    id_reference_allele = locus_alleles_ids[ref_allele_index] # id del alelo representante

    reference_allele_directory = 'reference_alleles'
 
    full_path_reference_allele = os.path.join(store_dir, reference_allele_directory)
    if not os.path.exists(full_path_reference_allele):
        try:
            os.makedirs(full_path_reference_allele)
            logger.info('Directory %s has been created', full_path_reference_allele)
        except:
            print ('Cannot create the directory ', full_path_reference_allele)
            logger.info('Directory %s cannot be created', full_path_reference_allele)
            exit (0)
 
    # build the fasta file name to store under first_allele_firectory
    f_name = os.path.basename(file_sequence)
    fasta_file = os.path.join(full_path_reference_allele, f_name)

    with open (fasta_file, 'w') as out_fh:
        out_fh.write ('>' + str(id_reference_allele) + '\n' + reference_allele)

    return fasta_file
"""

### FUNCIÓN ANTERIOR A PROBAR CON MASH. SI SE USA ESTA HABRÍA QUE QUITAR alleles_seqs_dict Y METER EL PATH AL FASTA DEL LOCUS EN CUESTIÓN Y PARSEAR DENTRO LAS SECUENCIAS Y SUS IDS.
### Reformando función check_core_gene_quality
"""
def get_reference_allele(alleles_seqs_dict, locus_quality, file_sequence, store_dir, logger): ### MODIFICANDO/CAMBIANDO: Función reference_allele añadida para encontrar el alelo de referencia
                                                                                                ### PROBANDO obtención de alelo de referencia con función distance de paquete Levenshtein
                                                                                                ### logger?

    ### logger.info('search reference allele')
    
    # alleles_seqs_dict --> diccionario generado previamente key:id alelo - value:seq alelo

    locus_alleles_seqs = list(alleles_seqs_dict.values()) # secuencias de los alelos del locus
    locus_alleles_ids = list(alleles_seqs_dict.keys()) # ids de los alelos del locus

    alleles_matrix_distance = np.zeros((len(locus_alleles_seqs), len(locus_alleles_seqs)))

    for allele_index_1 in range(len(locus_alleles_seqs)):
        for allele_index_2 in range(len(locus_alleles_seqs)):
            levenshtein_distance = distance(locus_alleles_seqs[allele_index_1], locus_alleles_seqs[allele_index_2])
            alleles_matrix_distance[allele_index_1][allele_index_2] = levenshtein_distance

    allele_mean_distance = np.mean(alleles_matrix_distance, 0)


    # casos en los que todos los alelos del locus están imcompletos? -> si se ignoran los que están incompletos, se estarían ignorando todos los alelos del locus

    min_mean = max(allele_mean_distance) # inicializando con la mayor mean
    ref_allele_index = int()

    for mean_index in range(len(allele_mean_distance)):
        if allele_mean_distance[mean_index] < min_mean:
            #if locus_alleles_ids[mean_index] in locus_quality['good_quality']:
            if locus_quality[locus_alleles_ids[mean_index]] == 'good_quality':
                min_mean = allele_mean_distance[mean_index]
                ref_allele_index = mean_index

   # print("allele_mean_distance[ref_allele_index]: ", allele_mean_distance[ref_allele_index]) ### printeando, borrar

    min_index_allele = min(allele_mean_distance)
    index_alelo_min_mean = list(allele_mean_distance).index(min(allele_mean_distance))
    #print("min_index_allele: ", min_index_allele)
    #print("Alelo que tiene una mean diferente que el alelo que toma, checkeando si está incompleto: ", locus_alleles_seqs[index_alelo_min_mean]) ### printeando, borrar

    reference_allele = locus_alleles_seqs[ref_allele_index] # alelo representante del locus, cuya media de distancias al resto de alelos es la más pequeña (el que menos difiere del resto) y que está completo, con codón de stop y de start
    id_reference_allele = locus_alleles_ids[ref_allele_index] # id del alelo representante

    reference_allele_directory = 'reference_alleles'
 
    full_path_reference_allele = os.path.join(store_dir, reference_allele_directory)
    if not os.path.exists(full_path_reference_allele):
        try:
            os.makedirs(full_path_reference_allele)
            logger.info('Directory %s has been created', full_path_reference_allele)
        except:
            print ('Cannot create the directory ', full_path_reference_allele)
            logger.info('Directory %s cannot be created', full_path_reference_allele)
            exit (0)
 
    # build the fasta file name to store under first_allele_firectory
    f_name = os.path.basename(file_sequence)
    fasta_file = os.path.join(full_path_reference_allele, f_name)

    with open (fasta_file, 'w') as out_fh:
        out_fh.write ('>' + str(id_reference_allele) + '\n' + reference_allele)

    return fasta_file
"""


### PROBANDO get_reference_allele USANDO MASH
### Reformando función check_core_gene_quality
### Reformando get_reference_allele con mash
### Quitar uso de diccionario de locus con ids-alelos, parsear con seq.IO

####################################################################################################################################################################################
def get_reference_allele(locus_quality, fasta_file, store_dir, logger): ### MODIFICANDO/CAMBIANDO: Función reference_allele añadida para encontrar el alelo de referencia
                                                                                                ### PROBANDO obtención de alelo de referencia con función distance de paquete Levenshtein
                                                                                                ### logger?
    ### logger.info('search reference allele')
    
    # alleles_seqs_dict --> diccionario generado previamente key:id alelo - value:seq alelo
    # fasta_file --> Es el path entero del locus en la carpeta de core_gene que se le mete en la terminal

    ### Me llevo la generación del directorio de alelos de referencia al inicio del código para poder así generar en reference_alleles el directorio mash donde se va a guardar el split del multifasta para cada locus y el sketch de cada locus
    reference_allele_directory = 'reference_alleles'
 
    full_path_reference_allele = os.path.join(store_dir, reference_allele_directory)
    if not os.path.exists(full_path_reference_allele):
        try:
            os.makedirs(full_path_reference_allele)
            logger.info('Directory %s has been created', full_path_reference_allele)
        except:
            print ('Cannot create the directory ', full_path_reference_allele)
            logger.info('Directory %s cannot be created', full_path_reference_allele)
            exit (0)
 

    ### Aquí metería una función para generar los archivos fasta individuales a partir del multifasta, generando previamente el directorio mash dentro del directorio reference_alleles
    ### y dentro del directorio mash crearía un directorio para el locus en cuestión donde se guardarían los fasta individuales
    ### En esa función metería un subprocess de mash sketch para generar el sketch con todos los fastas individuales y haría otro subproces de mash dist para obtener las distancias y las 
    ### guardaría en forma de matriz que va a sacar como output dentro de la función get_reference_allele y no sé si también la lista con todos los archivos que hay en la carpeta que
    ### debería haberse creado para 

    ### Creacion del directorio mash
    mash_directory = 'mash'
    full_path_mash = os.path.join(full_path_reference_allele, mash_directory)
    if not os.path.exists(full_path_mash):
        try:
            os.makedirs(full_path_mash)
            logger.info('Directory %s has been created', full_path_mash)
        except:
            print ('Cannot create the directory ', full_path_mash)
            logger.info('Directory %s cannot be created', full_path_mash)
            exit (0)

    ### Creacion del directory mash para cada locus. Si cada vez que se obtiene la referencia para un locuse se borra el directorio mash igual no haría falta generar un directorio
    ### dentro de mash para cada locus
    
    f_name = os.path.basename(fasta_file).split('.')
    full_path_locus_mash = os.path.join(full_path_mash, f_name[0])
    
    if not os.path.exists(full_path_locus_mash):
        try:
            os.makedirs(full_path_locus_mash)
            logger.info('Directory %s has been created', full_path_locus_mash)
        except:
            print ('Cannot create the directory ', full_path_locus_mash)
            logger.info('Directory %s cannot be created', full_path_locus_mash)
            exit (0)

    ### -
    ### split_multifasta_command y así poder acceder al archivo al cual le quiero hacer split
    #split_fasta_path = os.path.join(full_path_locus_mash, filename) ### path a cada fasta individual que se va a generar con split_multifasta_command
                                                                    ### al hacerlo con python no puedo poner filename, tengo que poner el id de nombre del archivo, id.fasta

    #split_multifasta_command = ["cat", fasta_file, "| awk '{if (substr($0, 1, 1)=='>') {filename=(substr($0,2) '.fa')} print $0 >", split_fasta_path, "}'"]
    #split_multifasta_result = subprocess.run(split_multifasta_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ### -

    ### Divido el multifasta del locus en fastas individuales cada uno con un alelo del locus cuyo nombre del archivo se corresponde con el id del alelo. Se guardan en el
    ### directorio nombrado como el locus en cuestión y que se encuentra dentro de la carpeta mash

    for record in SeqIO.parse(fasta_file, "fasta"):
        split_fasta_path = os.path.join(full_path_locus_mash, str(record.id) + ".fasta")
        with open (split_fasta_path, 'w') as out_fh:
            out_fh.write ('>' + str(record.id) + '\n' + str(record.seq))

    ### obtengo la lista de paths de fastas que hay en el directorio mash del locus
    split_multifasta_files_list = get_fasta_file_list(full_path_locus_mash, logger)

    ### comando mash para generar sketch con todos los fasta files para obtener distancias entre secuencias simultáneamente
    sketch_path = os.path.join(full_path_locus_mash, "reference.msh")
    mash_sketch_command = ["mash", "sketch", "-o", sketch_path]
    for fasta_path in split_multifasta_files_list: ### Añadiendo los paths de todos los fastas que se van a incluir en el sketch al comando mash_sketch_command
        mash_sketch_command.append(fasta_path)

    mash_sketch_result = subprocess.run(mash_sketch_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


    ### comando mash para obtener mash distances entre sencuencias
    mash_distance_command = ["mash", "dist", sketch_path, sketch_path]
    mash_distance_result = subprocess.Popen(mash_distance_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    out, err = mash_distance_result.communicate()
    out = out.decode('UTF-8').split('\n')
    #print("out_decode: ", out, '\n')


    comp_dist_list = []
    for n in range(len(out)-1): ### restando 1 porque genera una línea final vacía, no sé si es algo general que pasará siempre o solo se está dando en el locus que estoy usando como referencia para construir la función
        comp = out[n].split('\t')
        comp_dist = float(comp[2]) ### cuidado porque por la cara mete un último elemento en la lista que es un "blanco" y no tiene elemento 2, porque es una lista vacía -.-
        comp_dist_list.append(comp_dist)


    ### Obteniendo array de distancias y array de means
    comp_dist_list_per_allele = []
    alleles_number = len(split_multifasta_files_list)

    for index_distance in range(0, len(comp_dist_list), alleles_number):  
    
        dist_per_allele = comp_dist_list[index_distance : index_distance + alleles_number] 
        comp_dist_list_per_allele.append(dist_per_allele)

    comp_dist_arr_per_allele = np.asarray(comp_dist_list_per_allele)
    allele_mean_distance = np.mean(comp_dist_arr_per_allele, 0)
    
    #print("allele_mean_distance: ", allele_mean_distance)


    min_mean = max(allele_mean_distance) # inicializando con la mayor mean
    ref_allele_id = str()

    for mean_index in range(len(allele_mean_distance)):
        if allele_mean_distance[mean_index] < min_mean:
            allele_id = os.path.basename(split_multifasta_files_list[mean_index]).split('.')[0]
            #if locus_alleles_ids[mean_index] in locus_quality['good_quality']:
            if locus_quality[allele_id] == 'good_quality':
                min_mean = allele_mean_distance[mean_index]
                ref_allele_id = allele_id

    print("ref_allele_id: ", ref_allele_id, '\n')
    print("min_mean: ", min_mean, '\n')

    ### Moviendo fasta que contiene el alelo de referencia al directorio reference
    #for file_path in split_multifasta_files_list:
     #   f_name = os.path.basename(file_path).split('.')[0]
      #  if f_name == ref_allele_id:
       #     shutil.move(file_path, full_path_reference_allele)
    
    # Probando mover fasta metiendo path obtenido generándolo a partir del path del locus en el directorio mash y el id del alelo de referencia obtenido
    reference_file_path = os.path.join(full_path_locus_mash, ref_allele_id + ".fasta")
    new_reference_file_path = os.path.join(full_path_reference_allele, os.path.basename(fasta_file))
    #shutil.move(reference_file_path, full_path_reference_allele)
    shutil.move(reference_file_path, new_reference_file_path) ### tengo que renombrar el archivo del alelo de referencia poniéndole el nombre del locus
                                                                ### intento de renombrar el archivo al moverlo al nuevo directorio. Si no sirve probar con os.rename

    ### Borrando el directorio mash del locus ya que el archivo de interes se ha movido al directorio reference
    ### Si cada vez que obtengo una referencia para un locus borro mash para el siguiente alelo se volvería a crear
    ### y no seria necesario crear el directorio para el locus dentro del directorio mash
    shutil.rmtree(full_path_locus_mash)

 
    # build the fasta file name to store under reference_allele_directory
    #f_name = os.path.basename(fasta_file)
    #fasta_file = os.path.join(full_path_reference_allele, f_name)

    #with open (fasta_file, 'w') as out_fh:
     #   out_fh.write ('>' + str(id_reference_allele) + '\n' + reference_allele)

    return new_reference_file_path

####################################################################################################################################################################################


def prepare_core_gene(core_gene_file_list, store_dir, outputdir, logger): ### añadiendo outputdir para guardar report de calidad y stat en el directorio de resultados
    #check if directory exists and files have fasta files
    #valid_core_gene_files = get_fasta_file_list(core_gene_dir, logger)
    #if not valid_core_gene_files :
    #    return False
    #logger.debug('Schema files to be processed are : %s', valid_core_gene_files)
    #processing the files in the schema

    schema_quality = {} ### cambiando/modificando: introducción del diccionario de calidad del esquema incluyendo todos los locus
    schema_variability = {}
    schema_statistics = {}
    file_list = []
    ###reference_alleles_list = [] ### CAMBIANDO/MODIFICANDO: añadiendo lista para meter el nombre de los archivos de los alelos de referencia
    ### first_alleles_list = [] ### cambiando/modificando: comentado todo lo relacionado con el uso del primer alelo del locus como referencia
    
    schema_variability_count = {} # diccionario para report estadisticas longitud
    schema_quality_per_class = {} ### diccionario para report de calidad


    blast_dir = os.path.join(store_dir,'blastdb')
    logger.info('start preparation  of core genes files')
    for fasta_file in core_gene_file_list:

        # parsing fasta file and get in the dictionary the id and the sequence
        # fasta_file_parsed_dict = parsing_fasta_file_to_dict(fasta_file, logger) --------> Generacion de diccionario de alelos por locus comentado
        f_name = os.path.basename(fasta_file).split('.')
        file_list.append(os.path.join(store_dir, f_name[0]))
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
        
        ###reference_alleles_list.append(get_reference_allele(locus_quality, fasta_file, store_dir, logger)) ### CAMBIANDO/MODIFICANDO: esta lista sustituye a la anterior (first_alleles_list). Listado de los paths de los archivos que contienen el alelo de referencia

        ###print("reference_alleles_list: ", reference_alleles_list) ### print, borrar

        ### Obteniendo longitudes de los alelos del locus
        #alleles_len = []
        #for allele in fasta_file_parsed_dict : ### sacando a partir de diccionario, comentado
         #   alleles_len.append(len(fasta_file_parsed_dict[allele]))

        alleles_len = []
        alleles_in_locus = list (SeqIO.parse(fasta_file, "fasta"))
        for allele_item in alleles_in_locus :
            alleles_len.append(len(str(allele_item.seq)))

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

    for core_gene in schema_quality:
        schema_quality_per_class[core_gene] = {'good_quality': 0, 'bad_quality: no_start': 0, 'bad_quality: no_stop': 0, 'bad_quality: no_start_stop': 0, 'bad_quality: multiple_stop': 0}
        for quality_class in schema_quality_per_class[core_gene]:
            if schema_quality[core_gene] == quality_class:
                schema_quality_per_class[core_gene][quality_class] += 1

    # saving schema quality to file
    logger.info('Saving schema quality information to file..')
    quality_file =  os.path.join(outputdir, 'schema_quality.tsv')
    with open (quality_file , 'w') as quality_fh :
        quality_fh.write('\t'.join(header_schema_quality) + '\n')
        for core in sorted (schema_quality_per_class) :

           # allele_number = []
            #for quality_class in schema_quality_per_class[core]:
             #   allele_number.append(schema_quality_per_class[core][quality_class])
            quality_fh.write(core + '\t' + '\t'.join (map(str, list(schema_quality_per_class[core].values()))) + '\t' + str(sum(list(schema_quality_per_class[core].values()))) )

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

            stat_fh.write(core + '\t' + '\t'.join (map(str,schema_statistics[core])) + '\t' + ', '.join(length_number) + '\t' + str(total_alleles))

    return file_list, schema_variability, schema_statistics, schema_quality ### cambiando/modificando: he borrado first_alleles_list del output de la función, he añadido reference_alleles_list y el diccionario schema_quality



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
### Reformando codigo get_prodigal_sequence --> Buscando secuencia prodigal por semejanza a secuencia blast (sseq sin gaps)

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

            print("ha entrado a int(contig_prodigal_id) == int(contig_blast_id)")

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


            ### Creacion del directorio para este contig para guardar los fasta de los genes prodigal de este contig en de forma indiviual + sketch de mash
            contig_directory = os.path.join(full_path_prodigal_genes_per_contig, contig_blast_id)
    
            if not os.path.exists(contig_directory):
                try:
                    os.makedirs(contig_directory)
                    logger.info('Directory %s has been created', contig_directory)
                except:
                    print ('Cannot create the directory ', contig_directory)
                    logger.info('Directory %s cannot be created', contig_directory)
                    exit (0)

            ### Generando fasta con la secuencia de este gen predicho por prodigal que se ha encontrado en el mismo contig en el que se encontro la secuencia tras BLAST
            contig_genes_path = os.path.join(contig_directory, str(rec.id) + '.fasta')
            with open (contig_genes_path, 'w') as out_fh:
                #out_fh.write ('>' + str(rec.id) + '\n' + str(rec.seq))
                out_fh.write ('>' + str(rec.description) + '\n' + str(rec.seq)) ### poniendo descripcion como cabecera del fasta para poder sacar posteriormente start y end para report prodigal


    ### comando mash para generar sketch con todos los fasta files para obtener distancias entre secuencias simultáneamente
    sketch_path = os.path.join(contig_directory, contig_blast_id + '_sketch.msh')
    print("sketch_path: ", sketch_path, '\n')

    ### obtengo la lista de paths de fastas que hay en el directorio contig_directory
    prodigal_gene_files_list = get_fasta_file_list(contig_directory, logger)

    mash_sketch_command = ["mash", "sketch", "-o", sketch_path]

    for fasta_path in prodigal_gene_files_list : ### Añadiendo los paths de todos los fastas que se van a incluir en el sketch al comando mash_sketch_command
        mash_sketch_command.append(fasta_path)
    print("mash_sketch_command: ", mash_sketch_command, '\n')

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
        print("comp: ", comp, '\n')
        comp_dist = float(comp[2]) ### cuidado porque por la cara mete un ultimo elemento en la lista que es un "blanco" y no tiene elemento 2, porque es una lista vaca -.-
        comp_dist_list.append(comp_dist)

        ### probando a hacer lista de ids tal y como compara mash para poder sacar el id por el índice de la min mean
        gene_id = os.path.basename(comp[0]).split('.')[0]
        id_genes_list.append(gene_id)

    print("comp_dist_list: ", comp_dist_list, '\n')

    min_dist = min(comp_dist_list) ### distancia minima
    min_dist_index = comp_dist_list.index(min_dist) ### sacando indice de la distancia minima en la lista comp_dist_list para obtener el id del gen predicho por prodigal que se aleja menos de la secuencia blast
    min_dist_id = id_genes_list[min_dist_index] ### sacando id del gen para el que se ha obtenido la mínima distancia

    ### se lee el fasta que cuyo nombre se corresponde con el id de la secuencia para la que se he obtenido una menor distancia con respecto a la secuencia blast (min_dist_prodigal_id)
    min_dist_prodigal_path = os.path.join(contig_directory, min_dist_id + '.fasta') ### path del gen para el que se ha obtenido menor distancia

    min_dist_gene = SeqIO.parse(min_dist_prodigal_path, "fasta") ### no sé si puedo sacar directamente la secuencia haciendo min_dist_gene_seq, o si tengo que hacer el for record in min_dist_gene, etc
                                                                ### no se puede. Lo saco haciendo for, pero mirar si hay otra forma de sacar la seq sin tener que hacer for ya que solo hay una secuencia en el fasta
    
    print("min_dist_id : ", min_dist_id, '\n')


    for gene in min_dist_gene:
        min_dist_seq = str(gene.seq)
        
        start_prodigal = str(gene.description.split('#')[1]) ### para sacar report prodigal
        end_prodigal = str(gene.description.split('#')[2]) ### para sacar report prodigal


    print("min_dist_seq: ", min_dist_seq, '\n')
    print("start_prodigal: ", start_prodigal, '\n')
    print("end_prodigal: ", end_prodigal, '\n')


    return min_dist_seq, start_prodigal, end_prodigal ### añadido start_prodigal, end_prodigal para poder sacar report prodigal, quitar?
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
    """
    if not prodigal_training(reference_genome_file, prodigal_dir, logger): ### cambiando/modificando: se llama a la función prodigal_training
        print('Error when creating training file for genes prediction. Check log file for more information. \n ')
        return False
    """

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


def length_thresholds(core_name, schema_statistics, percentlength): ### cambiando/modificando: función añadida para obtener el rango de longitud para clasificación en INF, ASM o ALM
                                                                    ### logger?
    locus_mean = int(schema_statistics[core_name][1])

    if percentlength != "SD": 
        max_length_threshold = math.ceil(locus_mean + ((locus_mean * int(percentlength)) / 100))
        min_length_threshold = math.floor(locus_mean - ((locus_mean * int(percentlength)) / 100))
    else:
        percentlength = int(schema_statistics[core_name][2])

        max_length_threshold = math.ceil(locus_mean + (locus_mean * percentlength))
        min_length_threshold = math.floor(locus_mean - (locus_mean * percentlength))

    return max_length_threshold, min_length_threshold

### cambiado todo por .reverse_complement()
"""
def reverse_complement(seq): ### CAMBIO TEMPORAL. cambiando/modificando: Función introducida temporalmente para la obtención de seqs reversas complementarias a partir de strings (sacando las secuencias de los diccionarios y no parseando)
    seq_cr = []
    seq = seq[::-1]
    for base in seq:
        if base == "A":
            seq_cr.append("T")
        if base == "T":
            seq_cr.append("A")
        if base == "C":
            seq_cr.append("G")
        if base == "G":
            seq_cr.append("C")
            
    return str(("").join(seq_cr))
"""
"""
def get_stop_codon_index(seq, tga_stop_codon, indel_position) :
    stop_codons = ['TAA', 'TAG','TGA']
    seq_len = len(seq)
    index = 0
    for index in range (0, seq_len -2, 3) :
    #while index < seq_len - 2:
        codon = seq[index : index + 3]
        # ignore posible stop codon before the indel position
        if index + 2 < indel_position :
            continue
        if codon in stop_codons :
            if codon == 'TGA' :
                if tga_stop_codon:
                    return index
                else :
                    continue

            return index
        #index +=3
    # Stop condon not foun tn the sequence
    return False
"""

def get_stop_codon_index(seq) : ### CAMBIANDO: He eliminado la parte de checkeo de TGA, he eliminado indel_position porque se le metía algo distinto para deletion e insertion y estoy aunando el proceso de buscar codón de stop
    stop_codons = ['TAA', 'TAG','TGA']
    seq_len = len(seq)
    index = 0
    for index in range (0, seq_len -2, 3) :
    #while index < seq_len - 2:
        codon = seq[index : index + 3]
        if codon in stop_codons :
            return index
        #index +=3
    # Stop condon not found inn the sequence
    return False

"""
def get_stop_codon_index(seq, allele_sequence, indel_position) : ### CAMBIANDO: He eliminado la parte de checkeo de TGA, he añadido allele sequence y checkeo para determinar marco de lectura para búsqueda de codón stop en función de la secuencia del alelo
    
    start_codon_forward = ['ATG','ATA','ATT','GTG', 'TTG']
    start_codon_reverse = ['CAT', 'TAT','AAT','CAC','CAA']

    stop_codons_forward = ['TAA', 'TAG','TGA']
    stop_codons_reverse = ['ATT', 'ATC','ACT']
    
    #seq_len = len(seq)
    #index = 0
    
    # guess reading frame for stop codon searching

    # Esto no serviría para casos de alelos sin stop ni inicio, ya que no se tiene referencia para averiguar el marco de lectura
    #       No podría saberse en estos casos nunca cuál es el marco de lectura? --> Qué hacer en estos casos?
    ###     check_sequence_order estaría gestionando mal el caso de alelos incompletos?
    ###     solo tiene en cuenta el codón de inicio pero podría ser que el alelo no tuviese codón de inicio pero sí de final y se pudiese mirar su orientación a través del codón final
    ###     En los caso en los que los alelos no se tiene ni codón de inicio ni de final, qué orientación devolvería en estos casos?

    if not allele_sequence[0:3] in start_codon_forward and not allele_sequence[-3:] in start_codon_reverse:
        if allele_sequence[0:3] in stop_codons_forward or allele_sequence[-3:] in stop_codons_reverse:
            if len(allele_sequence) % 3 == 0:
                index = 0
            if len(allele_sequence) % 3 == 1:
                index = 1
            if len(allele_sequence) % 3 == 2:
                index = 2
        else:
            index = 0  # Cuando el alelo no tiene ni stop codon ni codón de inicio, de momento, se toma index 0. Aañdir warning de marco de lectura sin referencia y posiblemente erróneo?
    else:
        index = 0 # Cuando la secuencia tiene codón de inicio (aunque falte codón de stop) el índice para la lectura es 0

    for index in range (0, len(seq), 3) :
        #while index < seq_len - 2:
        if len(seq[index:]) >= 3 :
            codon = seq[index : index + 3]
            # ignore posible stop codon before the indel position
            if index + 2 < indel_position :## Tiene sentido 
                continue
            if codon in stop_codons_forward :
                return index
            #index +=3
    # Stop condon not found in the sequence
    return False
"""

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

    print("blast_records: ", blast_records, '\n') ### printeando, borar
    for blast_record in blast_records:
        print("blast_record in blast_records: ", blast_record, '\n') ### printeando, borar
        for alignment in blast_record.alignments:
            print("alignment in blast_record.alignments: ", alignment, '\n') ### printeando, borar
            for match in alignment.hsps:
                print(" match in alignment.hsps: ", match, '\n') ### printeando, borar
                match_alignment = [['sample', match.sbjct],['match', match.match], ['schema',match.query]]
                print("match_alignment: ", match_alignment, '\n', '\n') 
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
def get_alignment (sample_seq, query_seq, blast_id, seq_type = "dna"):
    
    ### si se condiciona el uso de id 85 e id 90 habría que condicionar los params de pairwise2.align
    #param blast_id introducido por checkeo de coverage y para saber qué params utilizar en pairwise2 en función del BLAST empleado (BLAST 90 o BLAST 85 con dif params gaps)

    if seq_type == "protein":
        sample_seq = convert_to_protein(sample_seq)
        query_seq = convert_to_protein(query_seq)

    if not blast_id:
        # arguments pairwise2.align.globalms: match, mismatch, gap opening, gap extending
        alignments = pairwise2.align.localms(sample_seq, query_seq, 1, -2, -1, -1) ### cambiando/modificando: mismos match, mismatch, gap opening, gap extending especificados en BLAST ID 85
                                                                                    ### cambiando/modificando: local alignment en lugar de global para evitar mal alineamiento de secuencias encontradas muy cortas
                                                                                    
    else:
        alignments = pairwise2.align.localxx(sample_seq, query_seq) ### mirar qué parámetros por defecto utiliza blastn para gaps, reward, etc

    #print("alignments: ", alignments, '\n') ### printeando, borrar

    values = format_alignment(*alignments[0]).split('\n')
    #print("values: ", values, '\n') ### printeando, borrar
    
    match_alignment = [['sample', values[0]],['match', values[1]], ['schema',values[2]]]
    #print("match_alignment: ", match_alignment, '\n') ### printeando, borrar
    #print("len(sample_seq): ", len(sample_seq), '\n') ### printeando, borrar
    #print("len(query_seq): ", len(query_seq), '\n') ### printeando, borrar

    return match_alignment


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
                print("Sumando +1 a NIPH", '\n')
                summary_dict[key]['NIPH'] += 1
            elif 'NIPHEM' == values : ### cambiando/modificando: no cogía NIPHEM y lo contabilizaba como NIPH por 'NIPH' in values
                print("Sumando +1 a NIPHEM", '\n')
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

def allele_call_nucleotides (core_gene_dict_files, query_directory, sample_dict_files, blast_db_directory, prodigal_directory, inputdir, outputdir, cpus , percentlength, schema_variability, schema_statistics, schema_quality, logger ): ### CAMBIANDO/MODIFICANDO: he añadido schema statistics para poder tomar la mean, stdev, etc. No sé si eliminar schema variability. He añadido prodigal_directory, prodigal_train_db_directory y schema_quality
                                                                                                                ### cambiando/modificando: REVISAR argumentos que puedo quitarle a la función
    
    prodigal_report = [] # prodigal_report TEMPORAL para checkear las secuencias obtenidas con prodigal vs blast y las posiciones sstart y send


    shorter_seq_coverage = []
    longer_seq_coverage = [] # listas añadidas para calcular coverage medio de new_sseq con respecto a alelo para establecer coverage mínimo por debajo del cual considerar LNF
    equal_seq_coverage = []

    shorter_blast_seq_coverage = []
    longer_blast_seq_coverage = [] # listas añadidas para calcular coverage medio de sseq con respecto a alelo tras blast para establecer coverage mínimo por debajo del cual considerar LNF
    equal_blast_seq_coverage = []

    full_gene_list = []
    samples_matrix_dict = {} # to keep allele number
    matching_genes_dict = {} # to keep start and stop positions
    exact_dict = {} ### cambiando/modificando: diccionario para guardar los exact matches obtenidos junto con su calidad
    inferred_counter = 0
    inferred_alleles_dict = {} # to keep track of the new inferred alleles
    inf_dict = {} # Store the inferred alleles found for each sample
    paralog_dict = {}
    insertions_dict = {}
    deletions_dict = {}
    equal_dict = {} ### cambiando/modificando: añadido

    asm_dict = {} ### cambiando/modificando: añadido para clasificación asm/alm. quitando clasificacion delecion/insercion/equal
    alm_dict = {}


    list_insertions = {} # list all insertions together with Sample file and core gene
    list_deletions = {} # list all deletions together with Sample file and core gene
    list_equal = {} ### cambiando/modificando: añadido
    
    list_asm = {} ### cambiando/modificando: añadido para clasificación asm/alm. quitando clasificacion delecion/insercion/equal
    list_alm = {}

    plot_dict = {}
    snp_dict = {}
    protein_dict = {}
    match_alignment_dict = {}
    blast_parameters = '"6 , qseqid , sseqid , pident ,  qlen , length , mismatch , gapopen , evalue , bitscore , sstart , send , qstart , qend , sseq , qseq"'
    header_macthing_alleles_conting = ['Sample Name', 'Contig', 'Core Gene','start', 'stop', 'direction', 'codification']
    header_exact = ['Core Gene', 'Sample Name', 'Allele', 'Allele Quality', 'Contig', 'Query length', 'Contig start', 'Contig end', 'Sequence', 'Predicted Sequence'] ### cambiando/modificando: header para diccionario exact_dict añadido
    header_paralogs = ['Core Gene','Sample Name', 'Paralog Type', 'ID %', 'Allele', 'Allele Quality', 'Contig', 'Bit Score', 'Contig start', 'Contig end', 'Sequence', 'Predicted Sequence'] ### cambiando/modificando: introducido paralog type, Allele Quality y Predicted Sequence
    header_inferred = ['Core Gene','Sample Name', 'Inferred Allele name', 'Allele', 'Allele Quality', 'Contig', 'Bitscore', 'Query length', 'Contig length', 'New sequence length' , 'Mismatch' , 'gaps', 'Contig start', 'Contig end',  'New sequence']
    
    ### cambiando/modificando: añadido para clasificación asm/alm. quitando clasificacion delecion/insercion/equal
    header_insertions = [ 'Core Gene', 'Sample Name', 'Insertion item', 'Allele', 'Allele Quality', 'Contig', 'Bitscore', 'Query length', 'Contig length', 'New sequence length' , 'Mismatch' , 'gaps', 'Contig start', 'Contig end',  'New sequence']
    header_deletions = [ 'Core Gene', 'Sample Name', 'Deletion item', 'Allele', 'Allele Quality', 'Contig', 'Bitscore', 'Query length', 'Contig length', 'New sequence length' , 'Mismatch' , 'gaps', 'Contig start', 'Contig end',  'New sequence']
    

    header_asm = [ 'Core Gene', 'Sample Name', 'ASM item', 'Allele', 'Allele Quality', 'Contig', 'Bitscore', 'Query length', 'Contig length', 'New sequence length' , 'Mismatch' , 'gaps', 'Contig start', 'Contig end',  'New sequence', 'Additional info']
    header_alm = [ 'Core Gene', 'Sample Name', 'ALM item', 'Allele', 'Allele Quality', 'Contig', 'Bitscore', 'Query length', 'Contig length', 'New sequence length' , 'Mismatch' , 'gaps', 'Contig start', 'Contig end',  'New sequence', 'Additional info']
    
    
    header_plot = ['Core Gene', 'Sample Name' , 'Allele', 'Allele Quality', 'Contig','Bit Score', 'Contig start', 'Contig end', 'Sequence', 'Predicted Sequence']
    header_snp = ['Core Gene', 'Sample Name', 'Allele number', 'Position', 'Mutation Schema/Sample', 'Codon Schema/Sample','Protein in Schema/Sample', 'Missense/Synonymous','Annotation Sample / Schema']
    header_protein = ['Core Gene','Sample Name', 'Protein in ' , 'Protein sequence']
    header_match_alignment = ['Core Gene','Sample Name','Alignment', 'Sequence']
    

    lnf_dict = {} ### cambiando/modificando: introducido lnf_dict y header_lnf para crear report de lnf
    header_lnf = ['Core Gene', 'Sample Name', 'Allele', 'ID %', 'Allele length', 'New sequence length', 'Coverage', 'Additional info'] ### meter secuencias alelo y newsseq en caso de que haya?
                                                                                                                                    ### additional info -> Coverage under threshold: "threshold" 
                                                                                                                                    ###                    ID under threshold: "threshold"
                            
    ### Añadido header_prodigal_report para report prodigal
    header_prodigal_report = ['Core gene', 'Sample Name', 'Allele', 'Sequence type', 'BLAST start', 'BLAST end', 'Prodigal start', 'Prodigal end', 'BLAST sequence', 'Prodigal sequence']
    
    ### Añadido header_newsseq_coverage_report para determinar coverage threshold a imponer
    header_newsseq_coverage_report = ['Core gene', 'Sample Name', 'Query length', 'New sequence length', 'Locus mean', 'Coverage (new sequence/allele)', 'Coverage (new sequence/locus mean)']

    ### Añadido header_blast_coverage_report para determinar coverage threshold a imponer
    header_blast_coverage_report = ['Core gene', 'Sample Name', 'Query length', 'Blast sequence length', 'Locus mean', 'Coverage (blast sequence/allele)', 'Coverage (blast sequence/locus mean)']



    print('Allele calling starts')
    pbar = ProgressBar ()
    for core_file in pbar(core_gene_dict_files) :
    #for core_file in core_gene_dict_files:

        print("schema_quality: ", schema_quality, '\n', '\n') ### printeando, borrar

        full_gene_list.append(os.path.basename(core_file))
        logger.info('Processing core gene file %s ', core_file)
        core_name = os.path.basename(core_file)
        reference_query = os.path.join(query_directory, str( core_name + '.fasta'))
        
        #with open (core_file, 'rb') as core_f: ### comentando lectura de diccionarios de locus 
         #   core_dict = pickle.load(core_f)
        #logger.debug('load in memory the core file %s ', core_file)

        ### cambiando/modificando. Esto es del primer alelo que toma como referencia, comentado ya que no se usa a partir de ahora.
        # get the reference allele to be used to find the SNP 
        ### core_first_allele_file = os.path.join(outputdir, 'tmp', 'cgMLST', 'first_alleles',core_name + '.fasta')
        ### reference_allele_for_snp = str(SeqIO.parse(core_first_allele_file, 'fasta').__next__().seq)
        ### query_length = len(reference_allele_for_snp)
        
        ### cambiando/modificando: obteniendo rango de longitud para clasificación INF, ASM o ALM
        max_length_threshold, min_length_threshold = length_thresholds(core_name, schema_statistics, percentlength)

        for sample_file in sample_dict_files:
            #with open (sample_file,'rb') as sample_f :
            #    sample_dict = pickle.load(sample_f)
            #logger.debug('loaded in memory the sample file %s' , sample_file

            sample_name = os.path.basename(sample_file)
            if not sample_name in samples_matrix_dict:
                # initialize the sample list to add the number of alleles and the start, stop positions
                samples_matrix_dict[sample_name] = []
                matching_genes_dict[sample_name] = {}

            blast_db_name = os.path.join(blast_db_directory, sample_name, sample_name) ### path a la base de datos creada para esta muestra previamente con makeblastdb
            #blast_db_name = '/srv/project_wgmlst/samples_listeria/RA-L2073/blastdb'
            #reference_query = '/srv/project_wgmlst/lmo_test.fasta'

            cline = NcbiblastnCommandline(db=blast_db_name, evalue=0.001, perc_identity = 100, outfmt = blast_parameters , max_target_seqs=2, max_hsps=1, num_threads=1, query=reference_query)
            #cline = NcbiblastnCommandline(db=Gene_Blast_DB_name, evalue=0.001, outfmt=5, max_target_seqs=10, max_hsps=10,num_threads=1, query='/srv/project_wgmlst/seqSphere_listeria_cgMLST_test/targets/lmo0001.fasta')
            out, err = cline()
            out_lines = out.splitlines( )

            if len (out_lines) > 0 :

            ###    bigger_bitscore = 0
                allele_found = {}

                for line in out_lines :
                    values = line.split('\t')

                    s_length = values[4] ### cambiando/modificando: no sé si debería sacar la longitud de la seq encontrada pero quitándole los gaps por si acaso...
                    qseqid = values[0] ### cambiando/modificando: sacando id del alelo que ha hecho match para sacar su secuencia parseando, para coger su longitud y comparar la len de la sec encontrada en la muestra con la len del alelo que ha hecho match, ya que si es exact match esta sec encontrada tiene que tener la misma longtud

                    #matching_allele_seq = core_dict[int(qseqid)] ### obteniendo secuencia alelo a partir de diccionario
                    #matching_allele_seq = core_dict[qseqid]
                    
                    ### Obteniendo la secuencia del alelo que ha hecho match con la muestra parseando
                    alleles_in_locus = list (SeqIO.parse(reference_query, "fasta"))
                    for allele_item in alleles_in_locus :
                       if allele_item.id == qseqid :
                          break
                    matching_allele_seq = str(allele_item.seq)

                    matching_allele_length = len(matching_allele_seq)

                    #if int(s_length) in schema_variability[core_name] :
                    if int(s_length) == matching_allele_length: ### cambiando/modificando:
                                                             ### no sé si tendría sentido eliminar los subsets si se está imponiendo que la longitud de la secuencia encontrada tiene que ser igual a la longitud del alelo que ha hecho match, por lo que no podrían existir subsets (secuencias más cortas), no? Lo que sí se podría obtener es la misma secuencia que se ha encontrado en otras coordenadas del contig
                        contig_id = values[1]
                        gene_start = values[9]
                        gene_end = values[10]
                    ###    sseq = values[13] 
                    ###    qseq = values[14] 
                        allele_is_subset = False
                        if len(allele_found) > 0 :
                            # check if the new match is a subset of the previous allele found in blast
                            for allele in allele_found :
                                if allele_found[allele][9] == gene_start or allele_found[allele][10] == gene_end :
                                    logger.info('Found allele %s that starts or ends as the same position as %s ' , values[0], allele_found[allele][0])
                                    allele_is_subset = True
                                    break

                        if len(allele_found) == 0 or not allele_is_subset :
                            contig_id_start = str(contig_id + '_'+ gene_start)
                            allele_found[contig_id_start] = values
                        #import pdb; pdb.set_trace()
                        ### if  float(values[8]) > bigger_bitscore : ### cambiando/modificando: no se está utilizando el bigger_bitscore para nada, no sé si borrarlo o si puede ser útil
                            #qseqid , sseqid , pident ,  qlen , length , mismatch , gapopen , evalue , bitscore , sstart , send , qstart , qend ,sseq , qseq= values
                            #bigger_bitscore = int(bitscore)
                            ### bigger_bitscore = float(values[8])

                if len(allele_found) > 1:

                    print("Muestra ", sample_name, '\n') ### printeando, borrar
                    print("Match con alelo ", qseqid, "del locus ", core_name, '\n') ### printeando, borrar
                    print('NIPHEM', '\n') ### printeando, borrar
                    print('allele_found: ', allele_found, '\n') ### printeando, borrar

                    # found paralogs in the sample for the core gene
                    samples_matrix_dict[sample_name].append('NIPHEM')
                    print('samples_matrix_dict: ', samples_matrix_dict, '\n')

                    if not sample_name in paralog_dict :
                        paralog_dict[sample_name] = {}
                    if not core_name in paralog_dict[sample_name] :
                        paralog_dict[sample_name] [core_name]= []

                    for allele_item in allele_found : ### Por cada resultado en allele_found (cada match en la muestra con ID 100%)
                        
                        qseqid = allele_found[allele_item][0]
                        sseqid = allele_found[allele_item][1]
                        pident = allele_found[allele_item][2]
                        bitscore = allele_found[allele_item][8]
                        sstart = allele_found[allele_item][9]
                        send = allele_found[allele_item][10]
                        sseq = allele_found[allele_item][13]
                        qseq = allele_found[allele_item][14] ### para sacar informe prodigal

                        #####################################
                        ### Obteniendo calidad del alelo: ### cambiando/modificando
                        #####################################

                        allele_quality = schema_quality[core_name][qseqid]

                        #if int(qseqid) in schema_quality[core_name]["good_quality"]:
                        #if qseqid in schema_quality[core_name]["good_quality"]:
                         #   allele_quality = "good_quality"
                        #elif int(qseqid) in schema_quality[core_name]["bad_quality"]["no_start"]:
                        #elif qseqid in schema_quality[core_name]["bad_quality"]["no_start"]:
                         #   allele_quality = "bad_quality: no_start"
                        #elif int(qseqid) in schema_quality[core_name]["bad_quality"]["no_stop"]:
                        #elif qseqid in schema_quality[core_name]["bad_quality"]["no_stop"]:   
                         #   allele_quality = "bad_quality: no_stop"
                        #elif int(qseqid) in schema_quality[core_name]["bad_quality"]["no_start_stop"]:
                        #elif qseqid in schema_quality[core_name]["bad_quality"]["no_start_stop"]:
                         #   allele_quality = "bad_quality: no_start_stop"
                        #elif int(qseqid) in schema_quality[core_name]["bad_quality"]["multiple_stop"]:
                        #elif qseqid in schema_quality[core_name]["bad_quality"]["multiple_stop"]:
                         #   allele_quality = "bad_quality: multiple_stop"

                        ### Obteniendo secuencia predicha por prodigal si la calidad del alelo es bad_quality
                        if 'bad_quality' in allele_quality: 
                            #complete_predicted_seq = str(get_prodigal_sequence(sstart, send, sseqid, prodigal_directory, sample_name))
                            complete_predicted_seq, start_prodigal, end_prodigal = get_prodigal_sequence(sseq, sseqid, prodigal_directory, sample_name, core_name, logger)

                            ##### PARA SACAR INFORME DE PRODIGAL VS BLAST #####
                            prodigal_report.append([core_name, sample_name, qseqid, 'NIPHEM', sstart, send, start_prodigal, end_prodigal, sseq, complete_predicted_seq])

                        else:
                            complete_predicted_seq = "-"

                        print("complete_predicted_seq Prodigal: ", complete_predicted_seq, '\n', '\n') ### printeando, borrar

                        paralog_dict[sample_name][core_name].append(['NIPHEM', pident, qseqid, allele_quality, sseqid, bitscore, sstart, send, sseq, complete_predicted_seq]) ### cambiando/modificando: introducido allele_quality, complete_predicted_seq y etiqueta NIPHEM
                        if not sseqid in matching_genes_dict[sample_name] :
                            matching_genes_dict[sample_name][sseqid] = []
                        if sstart > send :
                            matching_genes_dict[sample_name][sseqid].append([core_name, sstart,send,'-','NIPHEM'])
                        else:
                            matching_genes_dict[sample_name][sseqid].append([core_name, sstart,send,'+', 'NIPHEM'])
                    continue

                elif len(allele_found) == 1 :

                    ## look for possible paralogs by finding other alleles that identity is  equal to  90%
                    paralog_found ={}
                    ### allele_sequence = allele_found[contig_id_start][14] # cambiando/modificando: Esto no sería necesario sacarlo, porque está el diccionario de core_genes y se puede sacar de ahí la seq como he hecho antes, además de que en caso de que el match haya sido más corto que el alelo del fasta se estaría cogiendo esta secuencia más corta... aunque la secuencia qseq en teoría tal y como he puesto antes que debe tener la misma longitud s_length que el alelo que hace match, no debería ser más corta, pero bueno, a ver cómo lo acabo poniendo al final...
                    
                    print("matching_allele_seq, secuencia del alelo con el que se ha obtenido resultado con BLAST ID 100 y se emplea como query en BLAST ID 90: ", matching_allele_seq) ### printeando, borrar
                    allele_sequence = allele_found[contig_id_start][14] ### sacando la secuencia
                    print("allele_sequence: ", allele_sequence) ### print, borrar


                    #cline = NcbiblastnCommandline(db=blast_db_name, evalue=0.001, perc_identity=85, reward=1, penalty=-2, gapopen=1, gapextend=1, outfmt= blast_parameters, max_target_seqs=10, max_hsps=10,num_threads=1)
                    cline = NcbiblastnCommandline(db=blast_db_name, evalue=0.001, perc_identity=90, outfmt= blast_parameters, max_target_seqs=1, max_hsps=1,num_threads=1, query=reference_query)
                    out, err = cline(stdin = allele_sequence) 
                    out_lines = out.splitlines()
                    for line in out_lines :
                        values = line.split('\t')
                        sseq = values[13]
                        s_length = len(sseq.replace('-', '' )) ### cambiando/modificando: sacando la longitud de la secuencia encontrada sin gaps para comparar longitud


                        #if int(s_length) in schema_variability[core_name] :
                        if min_length_threshold <= s_length <= max_length_threshold: ### cambiando/modificando: he cambiado in schema_variability por el rango permitido, ya que al ser un parálogo puede haber más margen de variablidad de longitud, al contrario que los NIPHEM
                            
                            ### Cambiando/modificando: Check del nuevo match para descartarlo si se trata de una subsecuencia del alelo encontrado previamente con ID 100 y guardado en allele_found
                            gene_start = values[9]
                            gene_end = values[10]

                            allele_is_subset = False
                            allele = str(list(allele_found.keys())[0]) ### cambiando/modificando: obteniendo el id del alelo encontrado con BLAST ID 100 y guardado en el dict allele_found

                            if allele_found[allele][9] == gene_start or allele_found[allele][10] == gene_end :
                                logger.info('Found allele %s that starts or ends as the same position as %s ' , values[0], allele_found[allele][0])
                                allele_is_subset = True

                            if not allele_is_subset:
                                contig_id = values[1]
                                #gene_start = values[9]
                                #gene_end = values[10]
                                #sseq = allele_found[allele_item][13]
                                #qseq = allele_found[allele_item][14]
                                contig_id_start = str(contig_id + '_' + gene_start)
                                ## skip the allele found in the 100% identity and 100% alignment
                                if not contig_id_start in allele_found:
                                    paralog_found[contig_id_start] = values

                            """
                            ### Comentando código anterior a añadir subset para intentar arreglar problema de NIPHs subsets
                            contig_id = values[1]
                            gene_start = values[9]
                            #gene_end = values[10]
                            #sseq = allele_found[allele_item][13]
                            #qseq = allele_found[allele_item][14]
                            contig_id_start = str(contig_id + '_' + gene_start)
                            ## skip the allele found in the 100% identity and 100% alignment
                            if not contig_id_start in allele_found:
                                paralog_found[contig_id_start] = values
                            """
                    
                    if len(paralog_found) == 0 :

                    # exact match found
                        qseqid = allele_found[contig_id_start][0]
                        sseqid = allele_found[contig_id_start][1]
                        s_length = allele_found[contig_id_start][4] ### cambiando/modificando: creo que no es neceasrio sacar s_length porque si hay 0 parálogos solo ha podido sacar justamente antes la s_length de esta secuencia, pero lo dejo de momento por si acaso
                        sstart = allele_found[contig_id_start][9]
                        send = allele_found[contig_id_start][10]
                        sseq = allele_found[contig_id_start][13]
                        qseq = allele_found[contig_id_start][14] ### para sacar report prodigal
                        
                        samples_matrix_dict[sample_name].append(qseqid)

                        #####################################
                        ### Obteniendo calidad del alelo: ### cambiando/modificando
                        #####################################
                        
                        allele_quality = schema_quality[core_name][qseqid]

                        #if int(qseqid) in schema_quality[core_name]["good_quality"]:
                        #if qseqid in schema_quality[core_name]["good_quality"]:
                         #   allele_quality = "good_quality"
                        #elif int(qseqid) in schema_quality[core_name]["bad_quality"]["no_start"]:
                        #elif qseqid in schema_quality[core_name]["bad_quality"]["no_start"]:
                         #   allele_quality = "bad_quality: no_start"
                        #elif int(qseqid) in schema_quality[core_name]["bad_quality"]["no_stop"]:
                        #elif qseqid in schema_quality[core_name]["bad_quality"]["no_stop"]:   
                         #   allele_quality = "bad_quality: no_stop"
                        #elif int(qseqid) in schema_quality[core_name]["bad_quality"]["no_start_stop"]:
                        #elif qseqid in schema_quality[core_name]["bad_quality"]["no_start_stop"]:
                         #   allele_quality = "bad_quality: no_start_stop"
                        #elif int(qseqid) in schema_quality[core_name]["bad_quality"]["multiple_stop"]:
                        #elif qseqid in schema_quality[core_name]["bad_quality"]["multiple_stop"]:
                         #   allele_quality = "bad_quality: multiple_stop"

                        ### Obteniendo secuencia predicha por prodigal si la calidad del alelo es bad_quality
                        if 'bad_quality' in allele_quality: 
                            #complete_predicted_seq = str(get_prodigal_sequence(sstart, send, sseqid, prodigal_directory, sample_name))
                            complete_predicted_seq, start_prodigal, end_prodigal = get_prodigal_sequence(sseq, sseqid, prodigal_directory, sample_name, core_name, logger)

                            ##### PARA SACAR INFORME DE PRODIGAL VS BLAST #####
                            prodigal_report.append([core_name, sample_name, qseqid, 'EXACT', sstart, send, start_prodigal, end_prodigal, sseq, complete_predicted_seq])
                        else:
                            complete_predicted_seq = '-'

                        if not core_name in exact_dict:
                            exact_dict[core_name] = {}
                        if not sample_name in exact_dict[core_name]:
                            exact_dict[core_name][sample_name] = []
                        exact_dict[core_name][sample_name] = [qseqid, allele_quality, sseqid, s_length, sstart, send, sseq, complete_predicted_seq]

                        if not sseqid in matching_genes_dict[sample_name] :
                            matching_genes_dict[sample_name][sseqid] = []
                        # store the matching genes in forward order
                        if sstart > send :
                            matching_genes_dict[sample_name][sseqid].append([core_name, sstart, send,'-','EXACT'])
                        else:
                            matching_genes_dict[sample_name][sseqid].append([core_name, sstart, send,'+','EXACT'])

                        print("Muestra ", sample_name, '\n') ### printeando, borrar
                        print("Match con alelo ", qseqid, "del locus ", core_name, '\n') ### printeando, borrar
                        print("EXACT MATCH", '\n') ### printeando, borrar
                        print("sseq encontrada en BLAST: ", sseq, '\n') ### printeando, borrar

                        # get the snp for the alleles that exact match
                        #alleles_in_locus = list (SeqIO.parse(reference_query, "fasta"))
                        #reference_allele = str(alleles_in_locus[1].seq)
                        
                        ###snp_information = get_snp(sseq, reference_allele_for_snp)
                        snp_information = get_snp(sseq, matching_allele_seq) ### duda cambiando/modificando: tiene sentido sacar información de snp en caso de exact match? no debería obtenerse nada si realmente son EXACT, no? 
                                                                             ### cambiando/modificando: he cambiado reference_allele_for_snp (que es la secuencia del alelo que se cogía antes como referencia, que era el primero del locus) por el alelo con el que ha hecho match esta secuencia habiándolo obtenido de core_dict
                        if len(snp_information) > 0 :
                            if not core_name in snp_dict :
                                snp_dict[core_name] = {}
                            if not sample_name in snp_dict[core_name] :
                                snp_dict[core_name][sample_name] = {}
                            snp_dict[core_name][sample_name][qseqid]= snp_information
                        continue
                    else:
                        # paralog has been found
                        paralog_matrix = {}
                        samples_matrix_dict[sample_name].append('NIPH')
                        if not sample_name in paralog_dict :
                            paralog_dict[sample_name] = {}
                        if not core_name in paralog_dict[sample_name] :
                            paralog_dict[sample_name] [core_name]= []

                        print("allele_found: ", allele_found)

                        allele_found_name = list(allele_found.keys())[0]
                        qseqid = allele_found[allele_found_name] [0] ### cambiando/modificando: sacando el id del alelo que hizo match en el BLAST con ID 100 para poder sacar la calidad
                        
                        #####################################
                        ### Obteniendo calidad del alelo que hizo match con todos los parálogos : ### cambiando/modificando
                        #####################################

                        allele_quality = schema_quality[core_name][qseqid]

                        #if int(qseqid) in schema_quality[core_name]["good_quality"]:
                        #if qseqid in schema_quality[core_name]["good_quality"]:
                         #   allele_quality = "good_quality"
                        #elif int(qseqid) in schema_quality[core_name]["bad_quality"]["no_start"]:
                        #elif qseqid in schema_quality[core_name]["bad_quality"]["no_start"]:
                         #   allele_quality = "bad_quality: no_start"
                        #elif int(qseqid) in schema_quality[core_name]["bad_quality"]["no_stop"]:
                        #elif qseqid in schema_quality[core_name]["bad_quality"]["no_stop"]:   
                         #   allele_quality = "bad_quality: no_stop"
                        #elif int(qseqid) in schema_quality[core_name]["bad_quality"]["no_start_stop"]:
                        #elif qseqid in schema_quality[core_name]["bad_quality"]["no_start_stop"]:
                         #   allele_quality = "bad_quality: no_start_stop"
                        #elif int(qseqid) in schema_quality[core_name]["bad_quality"]["multiple_stop"]:
                        #elif qseqid in schema_quality[core_name]["bad_quality"]["multiple_stop"]:
                         #   allele_quality = "bad_quality: multiple_stop"

                        # merging the 2 dictionary
                        paralog_matrix[sample_name] = {**allele_found, **paralog_found}

                        for paralog in paralog_matrix[sample_name] :
                          #  qseqid = paralog_matrix[sample_name][paralog] [0] ### cambiando/modificando: no saco el qseqid del parálogo para que coja el qseqid del alelo con el que hizo match el parálogo y que en el report no aparezca Query_1
                            sseqid = paralog_matrix[sample_name][paralog] [1]
                            pident = paralog_matrix[sample_name][paralog] [2]
                            bitscore = paralog_matrix[sample_name][paralog] [8]
                            sstart = paralog_matrix[sample_name][paralog][9]
                            send = paralog_matrix[sample_name][paralog] [10]
                            sseq = paralog_matrix[sample_name][paralog] [13]
                            #qseq = paralog_matrix[sample_name][paralog] [14] ### para sacar report prodigal

                            #qseq_no_gaps = qseq.replace('-', '') ### para sacar report prodigal

                            sseq_no_gaps = sseq.replace('-', '')

                            ### Obteniendo secuencia predicha por prodigal para cada secuencia paráloga si la calidad del alelo es bad_quality
                            if 'bad_quality' in allele_quality: 
                                #complete_predicted_seq = str(get_prodigal_sequence(sstart, send, sseqid, prodigal_directory, sample_name))
                                complete_predicted_seq, start_prodigal, end_prodigal= get_prodigal_sequence(sseq_no_gaps, sseqid, prodigal_directory, sample_name, core_name, logger)

                                ##### PARA SACAR INFORME DE PRODIGAL VS BLAST #####
                                prodigal_report.append([core_name, sample_name, qseqid, 'NIPH', sstart, send, start_prodigal, end_prodigal, sseq_no_gaps, complete_predicted_seq])

                            else:
                                complete_predicted_seq = '-'

                            paralog_dict[sample_name][core_name].append(['NIPH', pident, qseqid, allele_quality, sseqid, bitscore, sstart, send, sseq_no_gaps, complete_predicted_seq]) ### cambiando/modificando: incluida ID, allele_quality y complete_predicted_seq
                            if not sseqid in matching_genes_dict[sample_name] :
                                matching_genes_dict[sample_name][sseqid] = []
                            if sstart > send :
                                matching_genes_dict[sample_name][sseqid].append([core_name, sstart, send, '-', 'NIPH'])
                            else:
                                matching_genes_dict[sample_name][sseqid].append([core_name, sstart, send, '+', 'NIPH'])
                       
                

                        """
                        # merging the 2 dictionary
                        paralog_matrix[sample_name] = {**allele_found, **paralog_found}

                        for paralog in paralog_matrix[sample_name] :
                            
                            qseqid = paralog_matrix[sample_name][paralog] [0]
                            sseqid = paralog_matrix[sample_name][paralog] [1]
                            bitscore = paralog_matrix[sample_name][paralog] [8]
                            sstart = paralog_matrix[sample_name][paralog][9]
                            send = paralog_matrix[sample_name][paralog] [10]
                            sseq = paralog_matrix[sample_name][paralog] [13]

                            sseq_no_gaps = sseq.replace('-', '')

                            #####################################
                            ### Obteniendo calidad del alelo: ### cambiando/modificando
                            #####################################

                            if int(qseqid) in schema_quality[core_name]["good_quality"]:
                                allele_quality = "good_quality"
                            elif int(qseqid) in schema_quality[core_name]["bad_quality"]["no_start"]:
                                allele_quality = "bad_quality: no_start"
                            elif int(qseqid) in schema_quality[core_name]["bad_quality"]["no_stop"]:
                                allele_quality = "bad_quality: no_stop"
                            elif int(qseqid) in schema_quality[core_name]["bad_quality"]["no_start_stop"]:
                                allele_quality = "bad_quality: no_start_stop"
                            elif int(qseqid) in schema_quality[core_name]["bad_quality"]["multiple_stop"]:
                                allele_quality = "bad_quality: multiple_stop"

                            ### Obteniendo secuencia predicha por prodigal si la calidad del alelo es bad_quality
                            if "bad_quality" in allele_quality: 
                                #complete_predicted_seq = str(get_prodigal_sequence(sstart, send, sseqid, prodigal_directory, sample_name))
                                complete_predicted_seq = str(get_prodigal_sequence(sseq_no_gaps, sseqid, prodigal_directory, sample_name, core_name, logger))
                            else:
                                complete_predicted_seq = "-"

                            paralog_dict[sample_name][core_name].append(['NIPH', qseqid, allele_quality, sseqid, bitscore, sstart, send, sseq_no_gaps, complete_predicted_seq]) ### cambiando/modificando: incluida etiqueta NIPH para diferenciar de NIPHEM en archivo de resultados, allele_quality y complete_predicted_seq
                            if not sseqid in matching_genes_dict[sample_name] :
                                matching_genes_dict[sample_name][sseqid] = []
                            if sstart > send :
                                matching_genes_dict[sample_name][sseqid].append([core_name, sstart, send, '-', 'NIPH'])
                            else:
                                matching_genes_dict[sample_name][sseqid].append([core_name, sstart, send, '+', 'NIPH'])
                        """
                        print("Muestra ", sample_name, '\n') ### printeando, borrar
                        print("Match con alelo ", qseqid, "del locus ", core_name, '\n') ### printeando, borrar
                        print("NIPH", '\n') ### printeando, borrar
                        print("paralog_found: ", paralog_found, '\n') ### printeando, borrar
                        print("allele_found: ", allele_found, '\n') ### printeando, borrar
                        print("paralog_matrix ---> paralog_found + allele_found", '\n', '\n') ### printeando, borrar
                        
                        continue

            ### duda Cambiando/modificando: aquí no habría que poner un if len(out_lines) == 0: entonces se hace el segundo blast?
            #cline = NcbiblastnCommandline(db=blast_db_name, evalue=0.001, perc_identity=85, reward=1, penalty=-2, gapopen=1, gapextend=1, outfmt= blast_parameters, max_target_seqs=1, max_hsps=1,num_threads=1, query=reference_query)
            cline = NcbiblastnCommandline(db=blast_db_name, evalue=0.001, perc_identity=90, outfmt= blast_parameters, max_target_seqs=1, max_hsps=1,num_threads=1, query=reference_query)
            
            blast_90 = True ### introducido por checkeo coverage e indicar params a utilizar en get_alignments en función del BLAST que se haya utlizado (BLAST 90 o BLAST 85 con dif params gaps)
            
            out, err = cline()
            out_lines = out.splitlines( )
            bigger_bitscore = 0
            if len (out_lines) == 0:

                print("Muestra ", sample_name, '\n') ### printeando, borrar
                print("Match con alelo ", qseqid, "del locus ", core_name, '\n') ### printeando, borrar
                print("LNF", '\n', '\n') ### printeando, borrar

                # trying to get the allele number to avoid that a bad quality assembly impact on the tree diagram
                cline = NcbiblastnCommandline(db=blast_db_name, evalue=0.001, perc_identity = 70, outfmt= blast_parameters, max_target_seqs=1, max_hsps=1,num_threads=1, query=reference_query)
                out, err = cline()
                out_lines = out.splitlines()

                if len (out_lines) > 0 :
                    
                    for line in out_lines :
                        values = line.split('\t')
                        if  float(values[8]) > bigger_bitscore: ### duda cambiando/modificnado: no sé qué sentdo tiene utilizar bigger_bitscore aquí, la verdad, y no en el resto de casos al mirar los resultados del blast...
                           ### qseqid , sseqid , pident ,  qlen , s_length , mismatch , gapopen , evalue , bitscore , sstart , send , qstart , qend ,sseq , qseq = values
                            qseqid = values[0]
                            #s_length = values[4]
                            pident = values[2] # obteniendo ID para report LNF
                            bitscore = values[8]
                            sseq = values[13] # obteniendo sseq para sacar s_length sin gaps

                            #print('q len seq is : ', len(qseq), ' s len seq is : ', len(sseq))
                            bigger_bitscore = float(bitscore)
                    #import pdb; pdb.set_trace()
                    #percent = int(s_length)/ query_length
                    #percent = int(s_length)/ len(core_dict[int(qseqid)]) ### cambiando/modificando: he cambiado query length por la longitud del alelo que hace match
                    ##percent = int(s_length)/ len(core_dict[qseqid]) # mensaje id bajo, no coverage

                    samples_matrix_dict[sample_name].append('LNF_' + str(qseqid)) ### indicar etiqueta así o solo LNF?
                
                    # Obteniendo alelo que ha hecho match a partir de diccionario core_dict
                    #matching_allele_seq = core_dict[qseqid]

                    ### Obteniendo la secuencia del alelo que ha hecho match con la muestra parseando
                    alleles_in_locus = list (SeqIO.parse(reference_query, "fasta"))
                    for allele_item in alleles_in_locus :
                       if allele_item.id == qseqid :
                          break
                    matching_allele_seq = str(allele_item.seq)
                    matching_allele_length = len(matching_allele_seq)

                    s_length = len(sseq.replace('-', '')) ### cambiando/modificando: sacando sseq para obtener longitud de la secuencia con mejor bitscore obtenido sin gaps

                    # obteniendo coverage del resultado con mejor bitscore (sacar coverage con respecto a la media o con respecto al alelo en sí?)
                    #coverage_blast = int(s_length)/matching_allele_length ### aquí debería poner matching_allele_length o la media en función de cómo se está calculando los rangos?
                    coverage_blast = '-'
                
                    coverage_new_sequence = '-'

                    new_sequence_length = '-'

                    add_info = 'BLAST sequence ID under threshold: 90%' ### En este mensaje poner el porcentaje de ID oportuno, de momento es 90%

                    ### creando report lnf
                    if not core_name in lnf_dict:
                        lnf_dict[core_name] = {}
                    if not sample_name in lnf_dict[core_name]:
                        lnf_dict[core_name][sample_name] = []

                    lnf_dict[core_name][sample_name].append([qseqid, pident, coverage_blast, coverage_new_sequence,  str(matching_allele_length), str(s_length), new_sequence_length, add_info]) ### Meter secuencias alelo, blast y new_sseq (si las hay)?

                    logger.info('BLAST sequence ID %s under threshold at sample %s, for gene  %s', pident, sample_name, core_name)


                else:
                    samples_matrix_dict[sample_name].append('LNF')
                    logger.info('Locus not found at sample %s, for gene %s', sample_name, core_name)

                    ### de momento dejo estas variables así hasta que cree alguna función para etiquetado LNF
                    qseqid = '-'
                    pident = '-'    
                    matching_allele_length = '-' ### introduciendo - en longitud de alelo match
                    s_length = '-' ### introduciendo - en longitud de seq encontrada en BLAST
                    coverage_blast = '-'
                    coverage_new_sequence = '-'
                    new_sequence_length = '-'
                    add_info = 'Locus not found'

                    ### creando report lnf
                    if not core_name in lnf_dict:
                        lnf_dict[core_name] = {}
                    if not sample_name in lnf_dict[core_name]:
                        lnf_dict[core_name][sample_name] = []

                    lnf_dict[core_name][sample_name].append([qseqid, pident, coverage_blast, coverage_new_sequence,  matching_allele_length, s_length, new_sequence_length, add_info]) ### Meter secuencias alelo, blast y new_sseq (si las hay)?

                continue

            ### Duda, cambiando/modificando: Creo que aquí se debería poner un if len(out_lines) > 0 y ya lo siguiente


            """
            #################################################################################################################################################
            # Probando a comprobar el coverage de la mejor secuencia obtenida con respecto al alelo que hizo match para descartar o no y correr BLAST 85 o no
            #################################################################################################################################################

            if len(out_lines) != 0: ## sería coomo meter dos veces len(out_lines), esta vez y después de comprobar LNFs pero si no se mete este if hay que hacer este checkeo después
                                    ## comprobar LNFs para BLAST ID 90 e igual no sirve para nada porque se mete en BLAST 85 por el coverage y hay que volver a calcular LNFs para ID 85
                                    ## Si meto este if este código iría antes de if len (out_lines) == 0:

                # Obteniendo los valores del resultado de BLAST ID 90 con mejor bitscore
                for line in out_lines :
                    values = line.split('\t')
            
                    if  float(values[8]) > bigger_bitscore:
                        qseqid , sseqid , pident ,  qlen , s_length , mismatch , gapopen , evalue , bitscore , sstart , send , qstart , qend ,sseq , qseq = values
                        bigger_bitscore = float(bitscore)

                # Obteniendo alelo que ha hecho match a partir de diccionario core_dict
                matching_allele_seq = core_dict[qseqid]
                matching_allele_length = len(matching_allele_seq)

                coverage_percent = (s_length/matching_allele_length)*100 ### duda: utilizar longitud seq quitando gaps a sseq o longitud con posibles gaps (s_length)?

                if coverage_percent < 70: ### duda: qué umbral de coverage poner?
                
                    cline = NcbiblastnCommandline(db=blast_db_name, evalue=0.001, perc_identity=85, reward=1, penalty=-2, gapopen=1, gapextend=1, outfmt= blast_parameters, max_target_seqs=1, max_hsps=1,num_threads=1, query=reference_query)
                    out, err = cline()
                    out_lines = out.splitlines()

                    ### Habría que volver a meter aquí el if len (out_lines) == 0: para obtener LNF en caso que proceda... 
                    ### Meter tras BLAST ID 90 una condición de rapideo para hacer este check que sea if len(out_lines) != 0? y ya meter este paso entero
                    ### y tras esta parte ya continuar normal checkeando LNF, parálogos, INF, etc con las out_lines de este BLAST si se da el caso de que se acaba ejecutando.

            ##############################################################################################################################################
            """


            ### cambiando/modificando: revisando posibles parálogos
            allele_found = {}

            for line in out_lines :
                values = line.split('\t')

                sseq = values[13] ### cambiando/modificando: sacando sseq para obtener longitud de la secuencia sin gaps
                s_length = len(sseq.replace('-', '')) ### sacando longitud de la secuencia encontrada sin gaps
                
                ### Sacando los valores de BLAST del mejor resultado para, si no hay resultados en allele_found tras el siguiente filtro, poder sacar el coverage del mejor resultado obtenido con este BLAST
                ### y poder incluir información en el report de LNFs
                #############################################################################################################################################################################################
                if  float(values[8]) > bigger_bitscore:
                    qseqid_bb , sseqid_bb , pident_bb ,  qlen_bb , s_length_bb , mismatch_bb , gapopen_bb , evalue_bb , bitscore_bb , sstart_bb , send_bb , qstart_bb , qend_bb ,sseq_bb , qseq_bb = values
                    #print('q len seq is : ', len(qseq), ' s len seq is : ', len(sseq))
                    bigger_bitscore = float(bitscore_bb)
                #############################################################################################################################################################################################


                #if min_length_threshold <= s_length <= max_length_threshold: ### Al mirar NIPHs e imponer esta condición, de que la longitud de la secuencia tiene que encontrarse en este rango, ya se estarían evitando las secuencias MUY Largas y MUY cortas tras BLAST
                                                                            ### Si no se mira NIPHs aquí entonces sí tendría sentido hacer un paso de coverage de las secuencias obtenidas tras BLAST para no pillar fragmentos de secuencias, o se cuencias muy largas o muy cortas, etc
                                                                            ### y tendría sentido que ese coverage fuese básicamente este rango de longitud, no? o no
                    
                if (schema_statistics[core_name][1] - schema_statistics[core_name][1]*0.5) <= s_length <= (schema_statistics[core_name][1] + schema_statistics[core_name][1]*0.5):
                                                                                            ### Cambiando rango anterior equivalente a clasificación en INF, ASM y ALM por rango umbral para clasificación de ASM y ALM a LNF
                                                                                            ### De momento pongo 50% de coverage de la seq encontrada tras blast con respecto a LA MEDIA (pongo media en lugar de len del alelo 
                                                                                            ### que hizo match ya que también se usa la media para el rango de clasificación de INF, etc... por determinar qué es mejor utilizar), 
                                                                                            ### aunque aún está por determinar. Hacerlo así o sacar coverage de la seq encontrada y comparar con 0.5
                    contig_id = values[1]
                    bitscore = values[8]
                    gene_start = values[9]
                    gene_end = values[10]
                    ###    sseq = values[13]  
                    ###    qseq = values[14] 
                    allele_is_subset = False
                    if len(allele_found) > 0 :
                        # check if the new match is a subset of the previous allele found in blast
                        for allele in allele_found :
                            if allele_found[allele][9] == gene_start or allele_found[allele][10] == gene_end : ### esto es suficiente? O tendría que poner como condición que no se solapen las posiciones porque un alelo puede hacer match una base arriba o abajo, no sería el mismo start ni end pero estaría solapando, estaría en la misma localización del genoma
                                                                                                                #### si la diferencia de start de allele_found - start del encontrad no es mayor a la longitud del encontrado? (si es menor sería que está solapando...)
                                if float(bitscore) > float(allele_found[allele][8]): ### cambiando/modificando: si se encuentra un resultado que empieza o termina donde otro pero su bitscore es mayor, se sustituye el resultado anterior por el nuevo con mayor bitscore
                                    allele_found[allele] = values
                               #     allele_is_subset = True ### cambiando/modificando: cambiaría a True subset aunque ya haya sobreescrito el alelo precisamente para que no vuelva a añadirlo al diccionario otra vez en el siguiente if
                               # else:
                                logger.info('Found allele %s that starts or ends as the same position as %s ' , values[0], allele_found[allele][0])
                                allele_is_subset = True
                                break 

                    if len(allele_found) == 0 or not allele_is_subset :
                        contig_id_start = str(contig_id + '_'+ gene_start)
                        allele_found[contig_id_start] = values

  
            print("allele_found: ", allele_found, '\n') ### printeando, checkeando el resultado tras mirar si hay NIPHs, borrar

            if len(allele_found) == 0:
                
                samples_matrix_dict[sample_name].append('LNF_' + str(qseqid_bb)) ### indicar etiqueta así o solo LNF?
                
                # Obteniendo alelo que ha hecho match a partir de diccionario core_dict
                #matching_allele_seq = core_dict[qseqid_bb]

                ### Obteniendo la secuencia del alelo que ha hecho match con la muestra parseando
                alleles_in_locus = list (SeqIO.parse(reference_query, "fasta"))
                for allele_item in alleles_in_locus :
                    if allele_item.id == qseqid :
                        break
                matching_allele_seq = str(allele_item.seq)
                matching_allele_length = len(matching_allele_seq)

                s_length_bb = len(sseq_bb.replace('-', '')) ### cambiando/modificando: sacando sseq para obtener longitud de la secuencia con mejor bitscore obtenido sin gaps

                # obteniendo coverage del resultado con mejor bitscore (sacar coverage con respecto a la media o con respecto al alelo en sí?)
                coverage_blast = int(s_length_bb)/ matching_allele_length ### aquí debería poner matching_allele_length o la media en función de cómo se está calculando los rangos?

                if coverage_blast < 1:
                    add_info = 'BLAST sequence coverage under threshold: 50%' ### En este mensaje poner el porcentaje de coverage oportuno, de momento pongo 50%

                else:  
                    add_info = 'BLAST sequence coverage above threshold: 50%' ### En este mensaje poner el porcentaje de coverage oportuno, de momento pongo 50%
                
                coverage_new_sequence = '-'

                new_sequence_length = '-'

                ### creando report lnf
                if not core_name in lnf_dict:
                    lnf_dict[core_name] = {}
                if not sample_name in lnf_dict[core_name]:
                    lnf_dict[core_name][sample_name] = []

                lnf_dict[core_name][sample_name].append([qseqid_bb, pident_bb, str(coverage_blast), str(coverage_new_sequence),  str(matching_allele_length), str(s_length_bb), new_sequence_length, add_info]) ### Meter secuencias alelo, blast y new_sseq (si las hay)?

                logger.info('BLAST sequence coverage %s under threshold at sample %s, for gene  %s', coverage_blast, sample_name, core_name)



            if len(allele_found) > 1: ### si hay más de una sec en allele_found etiqueta como parálogos       
                #print("Muestra ", sample_name, '\n') ### printeando, borrar
                #print("Match con alelo ", qseqid, "del locus ", core_name, '\n') ### printeando, borrar
                #print('NIPH', '\n') ### printeando, borrar
                #print('allele_found: ', allele_found, '\n') ### printeando, borrar

                # found paralogs in the sample for the core gene
                samples_matrix_dict[sample_name].append('NIPH')
                    #print('samples_matrix_dict: ', samples_matrix_dict, '\n')

                if not sample_name in paralog_dict :
                    paralog_dict[sample_name] = {}
                if not core_name in paralog_dict[sample_name] :
                    paralog_dict[sample_name] [core_name]= []
                
                for allele_item in allele_found : ### Por cada resultado en allele_found (cada match en la muestra con ID 100%)
                        
                    qseqid = allele_found[allele_item][0]
                    sseqid = allele_found[allele_item][1]
                    pident = allele_found[allele_item][2]
                    bitscore = allele_found[allele_item][8]
                    sstart = allele_found[allele_item][9]
                    send = allele_found[allele_item][10]
                    sseq = allele_found[allele_item][13]
                    #qseq = allele_found[allele_item][14] ### para sacar report prodigal

                    #qseq_no_gaps = qseq.replace('-', '') ### para sacar report prodigal

                    sseq_no_gaps = sseq.replace('-', '')
                    #####################################
                    ### Obteniendo calidad del alelo: ### cambiando/modificando
                    #####################################

                    allele_quality = schema_quality[core_name][qseqid]

                    #if int(qseqid) in schema_quality[core_name]["good_quality"]:
                    #if qseqid in schema_quality[core_name]["good_quality"]:
                     #   allele_quality = "good_quality"
                    #elif int(qseqid) in schema_quality[core_name]["bad_quality"]["no_start"]:
                    #elif qseqid in schema_quality[core_name]["bad_quality"]["no_start"]:
                     #   allele_quality = "bad_quality: no_start"
                    #elif int(qseqid) in schema_quality[core_name]["bad_quality"]["no_stop"]:
                    #elif qseqid in schema_quality[core_name]["bad_quality"]["no_stop"]:   
                     #   allele_quality = "bad_quality: no_stop"
                    #elif int(qseqid) in schema_quality[core_name]["bad_quality"]["no_start_stop"]:
                    #elif qseqid in schema_quality[core_name]["bad_quality"]["no_start_stop"]:
                     #   allele_quality = "bad_quality: no_start_stop"
                    #elif int(qseqid) in schema_quality[core_name]["bad_quality"]["multiple_stop"]:
                    #elif qseqid in schema_quality[core_name]["bad_quality"]["multiple_stop"]:
                     #   allele_quality = "bad_quality: multiple_stop"

                    ### Obteniendo secuencia predicha por prodigal si la calidad del alelo es bad_quality
                    if 'bad_quality' in allele_quality: 
                        #complete_predicted_seq = str(get_prodigal_sequence(sstart, send, sseqid, prodigal_directory, sample_name))
                        complete_predicted_seq, start_prodigal, end_prodigal = get_prodigal_sequence(sseq_no_gaps, sseqid, prodigal_directory, sample_name, core_name, logger)

                        ##### PARA SACAR INFORME DE PRODIGAL VS BLAST #####
                        prodigal_report.append([core_name, sample_name, qseqid, 'NIPH', sstart, send, start_prodigal, end_prodigal, sseq_no_gaps, complete_predicted_seq])

                    else:
                        complete_predicted_seq = "-"

                    paralog_dict[sample_name][core_name].append(['NIPH', pident, qseqid, allele_quality, sseqid, bitscore, sstart, send, sseq_no_gaps, complete_predicted_seq]) ### cambiando/modificando: introducido allele_quality, complete_predicted_seq y etiqueta NIPHEM
                    if not sseqid in matching_genes_dict[sample_name] :
                        matching_genes_dict[sample_name][sseqid] = []
                    if sstart > send :
                        matching_genes_dict[sample_name][sseqid].append([core_name, sstart,send,'-','NIPH'])
                    else:
                        matching_genes_dict[sample_name][sseqid].append([core_name, sstart,send,'+', 'NIPH'])

            if len(allele_found) == 1: ### Si no hay parálogos continúa el proceso de clasificación
                allele = str(list(allele_found.keys())[0]) ### cambiando/modificando: obteniendo el id del alelo que ha hecho match tras el BLAST con id 100 para continuar con la clasificación
                qseqid , sseqid , pident ,  qlen , s_length , mismatch , gapopen , evalue , bitscore , sstart , send , qstart , qend ,sseq , qseq = allele_found[allele]

                sseq_no_gaps = sseq.replace('-', '')

                """
                ### Este if de abajo iba dentro del for line in out_lines, pero estoy introduciendo testeo de NIPHs
                ### lo comento para probar la nueva versión intentando pillar NIPHs...
                # indento el resto del código en el if len(allele_found) == 1

                if  float(values[8]) > bigger_bitscore:
                    qseqid , sseqid , pident ,  qlen , s_length , mismatch , gapopen , evalue , bitscore , sstart , send , qstart , qend ,sseq , qseq = values
                    final_values = values ### cambiando/modificando: para comprobar resultados, borrar
                    #print('q len seq is : ', len(qseq), ' s len seq is : ', len(sseq))
                    bigger_bitscore = float(bitscore)
                """                

                ### Cambiando/modificando: igual no es necesario sacar todas las variables de arriba, cuando haya mirado cuáles se usan y cuáles no, borrar estas últimas
                ### qseqid hace falta para sacar alelo de core_dict y sseqid para sacar el contig parseando
                ### sstart y send hacen falta para comprobar si se trata de PLOT
                ### bitscore y sseq los guarda en el dict matching_genes_dict

                #####################################
                ### Obteniendo calidad del alelo: ### cambiando/modificando
                #####################################

                allele_quality = schema_quality[core_name][qseqid]

                #if int(qseqid) in schema_quality[core_name]["good_quality"]:
                #if qseqid in schema_quality[core_name]["good_quality"]:
                 #   allele_quality = "good_quality"
                #elif int(qseqid) in schema_quality[core_name]["bad_quality"]["no_start"]:
                #elif qseqid in schema_quality[core_name]["bad_quality"]["no_start"]:
                 #   allele_quality = "bad_quality: no_start"
                #elif int(qseqid) in schema_quality[core_name]["bad_quality"]["no_stop"]:
                #elif qseqid in schema_quality[core_name]["bad_quality"]["no_stop"]:   
                 #   allele_quality = "bad_quality: no_stop"
                #elif int(qseqid) in schema_quality[core_name]["bad_quality"]["no_start_stop"]:
                #elif qseqid in schema_quality[core_name]["bad_quality"]["no_start_stop"]:
                 #   allele_quality = "bad_quality: no_start_stop"
                #elif int(qseqid) in schema_quality[core_name]["bad_quality"]["multiple_stop"]:
                #elif qseqid in schema_quality[core_name]["bad_quality"]["multiple_stop"]:
                 #   allele_quality = "bad_quality: multiple_stop"


                ### Cambiando/modificando: Obteniendo alelo que ha hecho match con la muestra a partir de diccionario core_dict y contig donde se ha obtenido match en la muestra parseando
                #matching_allele_seq = core_dict[int(qseqid)]
                #matching_allele_seq = core_dict[qseqid]

                ### Obteniendo la secuencia del alelo que ha hecho match con la muestra parseando
                alleles_in_locus = list (SeqIO.parse(reference_query, "fasta"))
                for allele_item in alleles_in_locus :
                    if allele_item.id == qseqid :
                        break
                matching_allele_seq = str(allele_item.seq)
                matching_allele_length = len(matching_allele_seq)                


                # Retrieve the contig file for getting the contig sequence for the id found in Blast 
                contig_file = os.path.join(inputdir, str(sample_name + '.fasta'))
                records = list(SeqIO.parse(contig_file, "fasta"))
                for record in records:
                    if record.id == sseqid :
                        break
                accession_sequence = record.seq
                length_sseqid = len(accession_sequence)

                ### Revisando PLOTs
                ## check if the blast alignment could be classified as PLOT

                if int(sstart) == length_sseqid or int(send) == length_sseqid or int(sstart) == 1 or int(send) == 1:
                    #if int(s_length) < int(query_length) :
                    if int(s_length) < int(matching_allele_length): ### Cambiando/modificando/aclaración: Aquí no le estoy quitando los posibles gaps para comparar con la longitud del alelo query ya que si se quitan los gaps la secuencia podría ser más corta simplemente porque tenga deleciones (y además se haya colado en este if porque termina justo en el borde del contig) y no porque esté cortada por encontrarse en el borde del contig
                                                                ### Cambiando/modificando: cambio query_length por la longitud del alelo que ha hecho match
                        samples_matrix_dict[sample_name].append('PLOT_' + str(qseqid))
                        logger.info('PLOT found at sample %s, for gene  %s', sample_name, core_name)

                        ### Cambiando/modificando: Si la calidad del alelo es bad_quality, buscar la seq predicha por prodigal. Puede tratarse de la misma al tratarse de un PLOT y faltar justamente la parte que le faltaba al alelo porque coincide con el borde del contig
                        ### No sé si sacar la sec de prodigal para PLOT o quitarlo
                        if 'bad_quality' in allele_quality: 
                            #complete_predicted_seq = str(get_prodigal_sequence(sstart, send, sseqid, prodigal_directory, sample_name))
                            complete_predicted_seq, start_prodigal, end_prodigal = get_prodigal_sequence(sseq_no_gaps, sseqid, prodigal_directory, sample_name, core_name, logger)

                            ##### PARA SACAR INFORME DE PRODIGAL VS BLAST #####
                            prodigal_report.append([core_name, sample_name, qseqid, 'PLOT', sstart, send, start_prodigal, end_prodigal, sseq_no_gaps, complete_predicted_seq])

                        else:
                            complete_predicted_seq = "-"

                        if core_name not in plot_dict :
                            plot_dict[core_name] = {}
                        if not sample_name in plot_dict[core_name] :
                            plot_dict[core_name][sample_name] = []
                        plot_dict[core_name][sample_name].append([qseqid, allele_quality, sseqid, bitscore, sstart, send, sseq, complete_predicted_seq]) ### Modificando/cambiando: he añadido allele_quality y la seq predicha por prodigal complete_predicted_seq cuando la calidad del alelo es bad_quality
                    
                        if not sseqid in matching_genes_dict[sample_name] :
                            matching_genes_dict[sample_name][sseqid] = []
                        if sstart > send :
                            matching_genes_dict[sample_name][sseqid].append([core_name, sstart,send,'-', 'PLOT'])
                        else:
                            matching_genes_dict[sample_name][sseqid].append([core_name, sstart,send,'+', 'PLOT'])

                        print("Muestra ", sample_name, '\n') ### printeando, borrar
                        print("Match con alelo ", qseqid, "del locus ", core_name, '\n') ### printeando, borrar
                        print("PLOT", '\n') ### printeando, borrar
                        print("plot_dict: ", plot_dict, '\n', '\n') ### printeando, borrar

                        continue

                ###############################################
                ### BÚSQUEDA DE LA SECUENCIA COMPLETA FINAL ### 
                ###############################################

                print("Muestra ", sample_name, '\n') ### printeando, borrar
                print("Match con alelo ", qseqid, "del locus ", core_name, '\n') ### printeando, borrar

                ### probando resultado de prodigal
                #complete_predicted_seq = str(get_prodigal_sequence(sstart, send, sseqid, prodigal_directory, sample_name))
                complete_predicted_seq, start_prodigal, end_prodigal = get_prodigal_sequence(sseq_no_gaps, sseqid, prodigal_directory, sample_name, core_name, logger)
                print('\n', '\n', '\n', "secuencia prodigal predicha: ", complete_predicted_seq, '\n', '\n', '\n')


                query_direction = check_sequence_order(matching_allele_seq, logger)

                ### contig_file = os.path.join(inputdir,str(sample_name + '.fasta')) ### cambiando/modificando: esto no se utiliza para nada, esto sería para sacar el contig donde se ha encontrado la secuencia, pero ya se ha sacado antes de checkear los plots
                ### records = list (SeqIO.parse(contig_file, "fasta"))

                ### Cambiando/modificando: obviando el checkeo de que la secuencia del alelo tenga como codón de stop TGA ya que no es una buena forma de prevenir los casos en los que TGA sea codificado en mitad de la secuencia como serina o triptófano en códigos genéticos diferentes al usual
                ### ya que al hacer esto se está sesgando que los alelos encontrados solo puedan tener TGA como nuevo codón stop si el codón stop del alelo que ha hecho match con la secuencia era TGA, lo cual no tiene sentido
               
                ### if allele_sequence.endswith ('TGA') or  allele_sequence.startswith ('TCA') :
                    ### tga_stop_codon = True
                ### else:
                    ### tga_stop_codon = False

                ### Cambiando/modificando (temporalmente): obteniendo secuencia alargada por el final para buscar el nuevo codón de stop. No se está considerando que el alelo esté incompleto y haya que buscar el codónd e inicio también. 
                ### no se alarga X nts, si no que se coge todo lo que queda de contig hasta el final para evitar casosen los que el número de nts que se coja no sea suficiente para encontrar el codón de stop o para evitar errores en los
                ### casos en los que el número de nts que queden hasta llegar al final del contig sea menor al número de nts que se hubiese establecido que se quieren coger para alargar la secuencia

                if query_direction == 'reverse':
                    if int(send) > int (sstart): ## increasing the number of nucleotides to check if getting  longer protein
                        sample_gene_sequence = accession_sequence[ : int(send) ]
                        sample_gene_sequence = sample_gene_sequence.reverse_complement()

                    else:
                        sample_gene_sequence = accession_sequence[ int(send) -1 : ]
                        # import pdb; pdb.set_trace()
                else:
                    if int(sstart) > int (send):
                        sample_gene_sequence = accession_sequence[ :  int(sstart) ]
                        sample_gene_sequence = sample_gene_sequence.reverse_complement()

                    else:
                        sample_gene_sequence = accession_sequence[ int(sstart) -1 : ]


                #sseq = sseq.replace('-','') ### MODIFICANDO: comentado.
                #stop_index = get_stop_codon_index(sample_gene_sequence, tga_stop_codon, int(s_length)- int(qstart)) ### deleción
                #stop_index = get_stop_codon_index(sseq, tga_stop_codon, qseq.find('-')) ### inserción
                stop_index = get_stop_codon_index(sample_gene_sequence) 
                ### cambiando/MODIFICANDO: No uso sseq para buscar el codón de stop, como estaba puesto antes, si no que meto sample_gene_sequence en la función get_stop_codon_index, que es la secuencia encontrada ampliada para la búsqueda del codón de stop
                ### Además, borro el parámetro tga_stop_codon y no lo introduzco a la función get_stop_codon_index modificada

                if stop_index != False:
                    new_sequence_length = stop_index +3
                    new_sseq = str(sample_gene_sequence[0:new_sequence_length]) ### MODIFICANDO: cambio de sseq por sample_gene_sequence

                    new_sseq_coverage = new_sequence_length/matching_allele_length ### introduciendo coverage new_sseq /// debería ser con respecto a la media?

                    #########################################################################################################################
                    ### Cambiando/modificando: introducido para determinar qué umbral de coverage poner. TEMPORAL
                    if new_sseq_coverage < 1:
                        shorter_seq_coverage.append([core_name, sample_name, str(matching_allele_length), str(new_sequence_length), str(schema_statistics[core_name][1]), str(new_sseq_coverage), str(new_sequence_length/schema_statistics[core_name][1])])
                    elif new_sseq_coverage > 1:
                        longer_seq_coverage.append([core_name, sample_name, str(matching_allele_length), str(new_sequence_length), str(schema_statistics[core_name][1]), str(new_sseq_coverage), str(new_sequence_length/schema_statistics[core_name][1])])
                    elif new_sseq_coverage == 1:
                        equal_seq_coverage.append([core_name, sample_name, str(matching_allele_length), str(new_sequence_length), str(schema_statistics[core_name][1]), str(new_sseq_coverage), str(new_sequence_length/schema_statistics[core_name][1])])
                    #########################################################################################################################

                    print("new_sseq encontrada por Taranis: ", new_sseq, '\n', '\n')
                    print("sseq encontrada por BLAST: ", sseq, '\n', '\n')
                    print("Inicio de clasificación de la secuencia", '\n', '\n')
                    print("Aquí vienen INF", '\n')

                    ### INF
                    if min_length_threshold <= new_sequence_length <= max_length_threshold:
                        ### duda: incluir cuando INF == new_sequence_length? --> Si se toma el rango de longitud por defecto, "media + SD" para clasificar puede que un resultado que tiene la
                        ### misma longitud que el alelo con el que hizo match se clsifique como ASM o ALM porque se trate del alelo más largo o más corto del locus, por lo que su longitud
                        ### no entraría dentro del margen permitido por la SD y aunque sea de igual tamaño que el alelo que hizo match no se consideraría INF por esto.

                        print("Ha entrado a INF", '\n')

                        # print("values tras BLAST: ", values, '\n') ### printeando, borrar
                        #  print("% ID tras BLAST: ", values[2], '\n', '\n') ### printeando, borrar

                        logger.info('Found new allele for core gene %s ', core_name)

                        ### adding new allele to the  inferred allele list if it is not already included
                        if not core_name in inferred_alleles_dict :
                            inferred_alleles_dict[core_name] = []
                        if not new_sseq in inferred_alleles_dict[core_name] : ### cambiando/modificando: he cambiado sseq por new_sseq ya que ahora se está considerando la secuencia completa habiendo encontrado el codón de stop
                            inferred_alleles_dict[core_name].append(new_sseq)
                        ### find the index to include in the sample matrix dict
                        index_inferred = inferred_alleles_dict[core_name].index(new_sseq)
                        inferred_allele = 'INF_' + core_name + '_' + str(index_inferred)
                        samples_matrix_dict[sample_name].append(inferred_allele)
                        #if not sample_name in inf_dict : ### cambiando/modificando: comentado
                        #   inf_dict[sample_name] = {} ### cambiando/modificando: comentado
                        #inf_dict[sample_name][core_name] = inferred_allele ### cambiando/modificando: comentado

                        ### cambiando/modificando: incluyendo para INF misma información que para deleciones e inserciones en los resultados
                        if not core_name in inf_dict:
                            inf_dict[core_name] = {}
                        if not sample_name in inf_dict[core_name]:
                            inf_dict[core_name][sample_name] = {}

                        len_sseq = len(sseq.replace('-', '')) ### cambiando/modificando: obtengo la longitud de la sseq sin los gaps, ya que s_length tiene en cuenta gaps en caso de que los haya en el alineamiento
                        inf_dict[core_name][sample_name][inferred_allele] = [qseqid, allele_quality, sseqid, bitscore, str(matching_allele_length), str(len_sseq), str(new_sequence_length), mismatch , gapopen, sstart, send,  new_sseq]

                        print("inf_dict: ", inf_dict, '\n', '\n', '\n')  ### printeando, borrar
                        print("INF: ", inferred_allele, '\n') ### printeando, borrar

                        # Get the SNP for the new allele inferred

                        #snp_information = get_snp(sseq, reference_allele_for_snp)
                        snp_information = get_snp(new_sseq, matching_allele_seq) ### cambiando/modificando: he cambiado sseq por new_sseq ya que ahora se está considerando la secuencia completa habiendo encontrado el codón de stop
                                                                             ### cambiando/modificando: he cambiado reference_allele_for_snp (que es la secuencia del alelo que se cogía antes como referencia, que era el primero del locus) por el alelo con el que ha hecho match esta secuencia
                        if len(snp_information) > 0 :
                            if not core_name in snp_dict :
                                snp_dict[core_name] = {}
                            if not sample_name in snp_dict[core_name] :
                                snp_dict[core_name][sample_name] = {}
                            snp_dict[core_name][sample_name][qseqid]= snp_information

                        ### Intrduciendo obtención de alineamientos adn y proteína para INFs al igual que en ASM y ALM

                        """
                        # execute again blast with the reference query the previous query found to get the aligment format to get the SNPs
                        if not core_name in match_alignment_dict :
                            match_alignment_dict[core_name] = {}
                            if not sample_name in match_alignment_dict[core_name] :
                                #match_alignment_dict[core_name][sample_name] = get_alignment_for_deletions (new_sseq,  matching_allele_seq)           ### (cambiando/modificando: he cambiado el nombre de qqseq por matching_allele_seq para unificar (aunque realmente era lo mismo porque había adjudicad allele_sequence a qseq)
                                match_alignment_dict[core_name][sample_name] = get_alignment (new_sseq, matching_allele_seq, blast_90) ### cambiando/modificando: sustituyendo por get alignment

                        # convert the sequence to protein
                        if not core_name in protein_dict :
                            protein_dict[core_name] = {}
                        if not sample_name in protein_dict[core_name] :
                            protein_dict[core_name][sample_name] = []
                        #protein_dict[core_name][sample_name] = nucleotide_to_protein_alignment(new_sseq, matching_allele_seq)              ### (cambiando/modificando: he cambiado el nombre de qqseq por matching_allele_seq para unificar (aunque realmente era lo mismo porque había adjudicad allele_sequence a qseq)
                        protein_dict[core_name][sample_name] = get_alignment (new_sseq, matching_allele_seq, blast_90, "protein") ### cambiando/modificando: sustituyendo por get alignment
                        """

                        if not sseqid in matching_genes_dict[sample_name] :
                            matching_genes_dict[sample_name][sseqid] = []
                        if sstart > send :
                            matching_genes_dict[sample_name][sseqid].append([core_name, sstart,send,'-', inferred_allele])
                        else:
                            matching_genes_dict[sample_name][sseqid].append([core_name, sstart,send,'+', inferred_allele])
                        continue


                        ##### PARA SACAR INFORME DE PRODIGAL VS BLAST #####
                        sseq_no_gaps = sseq.replace('-', '') ### para sacar report prodigal
                        prodigal_report.append([core_name, sample_name, qseqid, 'INF', sstart, send, start_prodigal, end_prodigal, sseq_no_gaps, complete_predicted_seq])


                    print("Aquí vienen ASM y ALM")

                    ### ASM/ALM deletion y ASM/ALM insertion (y ASLM/ALM equal)*** 
               
                    if not min_length_threshold <= new_sequence_length <= max_length_threshold: ### Cambiando/modificando: si la longtidud de la new seq NO se encuentra dentro del rango de longitud
                    
                        s_length = len(sseq.replace('-', '')) ### cambiando/modificando: sacando sseq para obtener longitud de la secuencia sin gaps
                                                            ### no sé si esta longitud sin gaps se saca antes
                                                            ### mirar dónde meter esto para hallarlo solo una vez para todos los lugares donde se necesita usarlo

                        ###############################################################################################################################
                        ## Introduciendo filtro coverage new_sseq. si la longitud de new_sseq no se encuentra dentro del rango establecido, se impone un
                        ## coverage mínimo/máximo de new_sseq para clasificar en ASM y ALM o considerar LNF.
                        if (schema_statistics[core_name][1] - schema_statistics[core_name][1]*0.5) <= new_sequence_length <= (schema_statistics[core_name][1] + schema_statistics[core_name][1]*0.5):

                            samples_matrix_dict[sample_name].append('LNF_' + str(qseqid)) ### indicar etiqueta así o solo LNF?
                
                            # Obteniendo alelo que ha hecho match a partir de diccionario core_dict
                            #matching_allele_seq = core_dict[qseqid] ### Esto se ha obtenido antes, pero dejar aquí de momento por si creo función para etiqueta de lnf ver cómo generalizar los distintos casos donde se etiqueta

                            ### Obteniendo la secuencia del alelo que ha hecho match con la muestra parseando
                            alleles_in_locus = list (SeqIO.parse(reference_query, "fasta"))
                            for allele_item in alleles_in_locus :
                                if allele_item.id == qseqid :
                                    break
                            #matching_allele_seq = str(allele_item.seq)
                            matching_allele_seq = allele_item.seq ### probando a no convertir en este paso a str ya que para SNP hay que utilizar reverse_complement()
                            matching_allele_length = len(matching_allele_seq)


                            s_length = len(sseq.replace('-', '')) ### cambiando/modificando: sacando sseq para obtener longitud de la secuencia sin gaps
                                                                ### no sé si esta longitud sin gaps se saca antes

                            # obteniendo coverage de new_sseq (sacar coverage con respecto a la media o con respecto al alelo en sí?)
                            coverage_blast = int(s_length)/ matching_allele_length ### aquí debería poner matching_allele_length o la media en función de cómo se está calculando los rangos?
                        
                            coverage_new_sequence = int(new_sequence_length)/ matching_allele_length ### aquí debería poner matching_allele_length o la media en función de cómo se está calculando los rangos?

                            if coverage_new_sequence < 1:
                                add_info = 'New sequence coverage under threshold: 50%' ### En este mensaje poner el porcentaje de coverage oportuno, de momento pongo 50%
                            else:  
                                add_info = 'New sequence coverage above threshold: 50%' ### En este mensaje poner el porcentaje de coverage oportuno, de momento pongo 50%

                            #new_sequence_length   ### ya se calcula antes, no hace falta calcular aquí

                            ### creando report lnf
                            if not core_name in lnf_dict:
                                lnf_dict[core_name] = {}
                            if not sample_name in lnf_dict[core_name]:
                                lnf_dict[core_name][sample_name] = []

                            lnf_dict[core_name][sample_name].append([qseqid, pident, str(coverage_blast), str(coverage_new_sequence),  str(matching_allele_length), str(s_length), new_sequence_length, add_info]) ### Meter secuencias alelo, blast y new_sseq (si las hay)?

                            logger.info('New sequence coverage %s under threshold at sample %s, for gene  %s', coverage_blast, sample_name, core_name)
                        ###############################################################################################################################

                        else: ## Si el coverage de new_sseq entra dentro del rango establecido, se procede a etiquetar en ASM/ALM

                            print("Ha entrado a ASM y ALM", '\n')
                            print("len(sseq) con gaps: ", len(sseq), '\n', '\n')

                            sseq = sseq.replace('-', '')
                        
                            ### ASM/ALM deletions
                    
                            print("Aquí vienen las deleciones", '\n')
                            print("len(sseq) sin gaps: ", len(sseq), '\n', '\n')
                            print("new_sequence_length: ", new_sequence_length, '\n', '\n')
                            print("matching_allele_length: ", matching_allele_length, '\n', '\n')
                            print("min_length_threshold: ", min_length_threshold, '\n', '\n')
                            print("max_length_threshold: ", max_length_threshold, '\n', '\n')


                            ######################################################
                            ### QUITANDO CLASIFICAIÓN DELETION/INSERTION/EQUAL ###
                            ######################################################

                            ### ASM
                            if new_sequence_length < min_length_threshold:  ### cambiando/modificando: si la nueva secuencia se encuentra por debajo del umbral mínimo de longitud para ser clasificada como INF                        
                                                                            ### modificando/cambiando: realmente la clasificación en deletion e insertion puede ser ambigua porque se puede estar clasificando una secuencia como deletion pero realmente haya ocurrido inserciones y deleciones, solo que la deleción ha sido más grande y por tanto la sseq sería finalmente más corta que la longitud del alelo con el que hizo match
                        
                                print("Ha entrado a ASM", '\n')
                                ### adding ASM allele to the asm_allele_matrix if it is not already include
                                if not core_name in asm_dict:
                                    asm_dict[core_name] = []
                                if not new_sseq in asm_dict[core_name] :
                                    asm_dict[core_name].append(new_sseq)
                                ### find the index of ASM  to include it in the sample matrix dict
                                index_asm = asm_dict[core_name].index(new_sseq)
                                asm_allele = 'ASM_' + core_name + '_' + str(qseqid) + '_' + str(index_asm)
                                
                                #if new_sequence_length < query_length : ### deleción
                                #if new_sequence_length < min(schema_variability[core_name]) : ### inserción
                                #if new_sequence_length < min_length_threshold:
                                if len(sseq) < matching_allele_length: ### si la secuencia obtenida con BLAST tiene una longitud menor a la del alelo con el que hizo match -> DELETION
                                    add_info = 'Global effect: DELETION. BLAST sequence length shorter than matching allele sequence length / Net result: ASM. Final gene sequence length shorter than matching allele sequence length'
                                    #delete_allele = 'ASM_DELETE_' + core_name + '_' + str(qseqid) + '_' + str(index_delete)

                                if len(sseq) == matching_allele_length:### si la secuencia obtenida con BLAST tiene la misma longitud que el alelo con el que hizo match -> 
                                    add_info = 'Global effect: BASE SUBSTITUTION. BLAST sequence length equal to matching allele sequence length / Net result: ASM. Final gene sequence length shorter than matching allele sequence length'

                                #if new_sequence_length > query_length : ### deleción
                                #if new_sequence_length > max(schema_variability[core_name]) : ### inserción
                                #if new_sequence_length > max_length_threshold:
                                if len(sseq) < matching_allele_length: ### si la secuencia obtenida con BLAST tiene una longitud mayor a la del alelo con el que hizo match -> INSERTION
                                    add_info = 'Global effect: INSERTION. BLAST sequence length longer than matching allele sequence length / Net result: ASM. Final gene sequence length shorter than matching allele sequence length'
                                    #delete_allele = 'ALM_DELETE_' + core_name + '_' + str(qseqid) + '_' + str(index_delete)

                                samples_matrix_dict[sample_name].append(asm_allele)

                                if not sseqid in matching_genes_dict[sample_name] :
                                    matching_genes_dict[sample_name][sseqid] = []
                                if sstart > send :
                                    matching_genes_dict[sample_name][sseqid].append([core_name, str(int(sstart)-new_sequence_length -1), sstart,'-', asm_allele]) ### modificando/cambiando: al hacer predecir con prodigal habría que incluir sstart y send obteniéndolos de los archivos que genera prodigal para saber dónde están las coordenadas de la nueva seq completa
                                else:
                                    matching_genes_dict[sample_name][sseqid].append([core_name, sstart,str(int(sstart)+ new_sequence_length),'+', asm_allele])

                                print("ASM: ", asm_allele, '\n') ### printeando, borrar

                                ### add the deletion into deletion list
                                if not core_name in list_asm :
                                    list_asm [core_name] = {}
                                if not sample_name in list_asm[core_name] :
                                    list_asm[core_name][sample_name] = {}

                                len_sseq = len(sseq) ### cambiando/modificando: obteniendo longitud sseq sin gaps, ya que s_length incluye gaps en caso de que los haya en el alineamiento
                                list_asm[core_name][sample_name][asm_allele] = [qseqid, allele_quality, sseqid,  bitscore, str(matching_allele_length), str(len_sseq), str(new_sequence_length), mismatch , gapopen, sstart, send,  new_sseq, add_info]
                                ### duda cambiando/modificando: no debería cambiar en list_deletions sstart y send por las sstart y send de la nueva secuencia completa? Ahora mismo se está guardando sstart y send de la secuencia sseq encontrada con BLAST

                                #########################
                                ### OBTENCIÓN DE SNPs ### Añadir función get_SNPs como en INF, además de sustituir get_alignment_for_deletions y nucleotide_to_protein_alignment por get_alignment
                                #########################

                                if check_sequence_order(matching_allele_seq, logger) == 'reverse':
                                    #matching_allele_seq = reverse_complement(matching_allele_seq)   ### (cambiando/modificando: he cambiado allele_sequence por matching_allele_seq)
                                    matching_allele_seq = str(matching_allele_seq.reverse_complement())
                                
                                # get the SNP for the  delection
                                #if not core_name in snp_dict :
                                #    snp_dict[core_name] = {}
                                #if not sample_name in snp_dict[core_name] :
                                #    snp_dict[core_name][sample_name] = []
                                #snp_dict[core_name][sample_name] = get_snp(new_sseq, matching_allele_seq)

                                # execute again blast with the reference query the previous query found to get the aligment format to get the SNPs
                                if not core_name in match_alignment_dict :
                                    match_alignment_dict[core_name] = {}
                                    if not sample_name in match_alignment_dict[core_name] :
                                        #match_alignment_dict[core_name][sample_name] = get_alignment_for_deletions (new_sseq,  matching_allele_seq)        ### (cambiando/modificando: he cambiado el nombre de qqseq por matching_allele_seq para unificar (aunque realmente era lo mismo porque había adjudicad allele_sequence a qseq)
                                        match_alignment_dict[core_name][sample_name] = get_alignment (new_sseq, matching_allele_seq, blast_90) ### cambiando/modificando: sustituyendo por get alignment

                                # convert the sequence to protein
                                if not core_name in protein_dict :
                                    protein_dict[core_name] = {}
                                if not sample_name in protein_dict[core_name] :
                                    protein_dict[core_name][sample_name] = []
                                #protein_dict[core_name][sample_name] = nucleotide_to_protein_alignment(new_sseq, matching_allele_seq)          ### (cambiando/modificando: he cambiado el nombre de qqseq por matching_allele_seq para unificar (aunque realmente era lo mismo porque había adjudicad allele_sequence a qseq)
                                protein_dict[core_name][sample_name] = get_alignment(new_sseq, matching_allele_seq, blast_90, "protein") ### cambiando/modificando: sustituyendo por get alignment                     

                                ##### PARA SACAR INFORME DE PRODIGAL VS BLAST #####
                                sseq_no_gaps = sseq.replace('-', '') ### para sacar report prodigal
                                prodigal_report.append([core_name, sample_name, qseqid, asm_allele, sstart, send, start_prodigal, end_prodigal, sseq_no_gaps, complete_predicted_seq])


                            ### ALM
                            if new_sequence_length > min_length_threshold:  ### cambiando/modificando: si la nueva secuencia se encuentra por debajo del umbral mínimo de longitud para ser clasificada como INF                        
                                                                            ### modificando/cambiando: realmente la clasificación en deletion e insertion puede ser ambigua porque se puede estar clasificando una secuencia como deletion pero realmente haya ocurrido inserciones y deleciones, solo que la deleción ha sido más grande y por tanto la sseq sería finalmente más corta que la longitud del alelo con el que hizo match
                        
                                print("Ha entrado a ALM", '\n')
                                ### adding ASM allele to the asm_allele_matrix if it is not already include
                                if not core_name in alm_dict:
                                    alm_dict[core_name] = []
                                if not new_sseq in asm_dict[core_name] :
                                    alm_dict[core_name].append(new_sseq)
                                ### find the index of ASM  to include it in the sample matrix dict
                                index_alm = alm_dict[core_name].index(new_sseq)
                                alm_allele = 'ALM_' + core_name + '_' + str(qseqid) + '_' + str(index_alm)
                                
                                #if new_sequence_length < query_length : ### deleción
                                #if new_sequence_length < min(schema_variability[core_name]) : ### inserción
                                #if new_sequence_length < min_length_threshold:
                                if len(sseq) < matching_allele_length: ### si la secuencia obtenida con BLAST tiene una longitud menor a la del alelo con el que hizo match -> DELETION
                                    add_info = 'Global effect: DELETION. BLAST sequence length shorter than matching allele sequence length / Net result: ALM. Final gene sequence length longer than matching allele sequence length'
                                    #delete_allele = 'ASM_INSERT_' + core_name + '_' + str(qseqid) + '_' + str(index_delete)

                                if len(sseq) == matching_allele_length:### si la secuencia obtenida con BLAST tiene la misma longitud que el alelo con el que hizo match -> 
                                    add_info = 'Global effect: BASE SUBSTITUTION. BLAST sequence length equal to matching allele sequence length / Net result: ALM. Final gene sequence length longer than matching allele sequence length'

                                #if new_sequence_length > query_length : ### deleción
                                #if new_sequence_length > max(schema_variability[core_name]) : ### inserción
                                #if new_sequence_length > max_length_threshold:
                                if len(sseq) < matching_allele_length: ### si la secuencia obtenida con BLAST tiene una longitud mayor a la del alelo con el que hizo match -> INSERTION
                                    add_info = 'Global effect: INSERTION. BLAST sequence length longer than matching allele sequence lengtH / Net result: ALM. Final gene sequence length longer than matching allele sequence length'
                                    #delete_allele = 'ALM_INSERT_' + core_name + '_' + str(qseqid) + '_' + str(index_delete)
                                
                                samples_matrix_dict[sample_name].append(alm_allele)

                                if not sseqid in matching_genes_dict[sample_name] :
                                    matching_genes_dict[sample_name][sseqid] = []
                                if sstart > send :
                                    matching_genes_dict[sample_name][sseqid].append([core_name, str(int(sstart)-new_sequence_length -1), sstart,'-', alm_allele]) ### modificando/cambiando: al hacer predecir con prodigal habría que incluir sstart y send obteniéndolos de los archivos que genera prodigal para saber dónde están las coordenadas de la nueva seq completa
                                else:
                                    matching_genes_dict[sample_name][sseqid].append([core_name, sstart,str(int(sstart)+ new_sequence_length),'+', alm_allele])

                                print("ALM: ", alm_allele, '\n') ### printeando, borrar


                                ### add the deletion into deletion list
                                if not core_name in list_alm :
                                    list_alm[core_name] = {}
                                if not sample_name in list_alm[core_name] :
                                    list_alm[core_name][sample_name] = {}

                                len_sseq = len(sseq) ### cambiando/modificando: obteniendo longitud sseq sin gaps, ya que s_length incluye gaps en caso de que los haya en el alineamiento
                                list_alm[core_name][sample_name][alm_allele] = [qseqid, allele_quality, sseqid,  bitscore, str(matching_allele_length), str(len_sseq), str(new_sequence_length), mismatch , gapopen, sstart, send,  new_sseq, add_info]
                                ### duda cambiando/modificando: no debería cambiar en list_deletions sstart y send por las sstart y send de la nueva secuencia completa? Ahora mismo se está guardando sstart y send de la secuencia sseq encontrada con BLAST

                                #########################
                                ### OBTENCIÓN DE SNPs ### Añadir función get_SNPs como en INF, además de sustituir get_alignment_for_deletions y nucleotide_to_protein_alignment por get_alignment
                                #########################

                                if check_sequence_order(matching_allele_seq, logger) == 'reverse':
                                    #matching_allele_seq = reverse_complement(matching_allele_seq)   ### (cambiando/modificando: he cambiado allele_sequence por matching_allele_seq)
                                    matching_allele_seq = str(matching_allele_seq.reverse_complement())    ### (cambiando/modificando: he cambiado str(matching_allele_seq).reverse_complement() por reverse_complement(matching_allele_seq), ya que no puedo emplear .reverse_complement() sobre una string, sino que tiene que ser sobre un objeto y al obtener la seq del diccionario core_dict se trata de una seq. Para que sea un objeto tendría que parsear. De momento he incluido la función reverse_complement para que funcione con string. Ver cómo lo dejo finalmente, como string (entonces tendría que extraer también el contig en caso de necesitarlo en la versión final a partir del ddiccionario de contigs) o parseando (como está ahora mismo la obtención de accession_sequence y en ese caso sería 0 necesario la generación de los archivos con los diccionarios para cada muestra y cada locus, al menos tal y como se generan ahora)
                                
                                # get the SNP for the  delection
                                #if not core_name in snp_dict :
                                #    snp_dict[core_name] = {}
                                #if not sample_name in snp_dict[core_name] :
                                #    snp_dict[core_name][sample_name] = []
                                #snp_dict[core_name][sample_name] = get_snp(new_sseq, matching_allele_seq)


                                # execute again blast with the reference query the previous query found to get the aligment format to get the SNPs
                                if not core_name in match_alignment_dict :
                                    match_alignment_dict[core_name] = {}
                                    if not sample_name in match_alignment_dict[core_name] :
                                        #match_alignment_dict[core_name][sample_name] = get_alignment_for_deletions (new_sseq,  matching_allele_seq)        ### (cambiando/modificando: he cambiado el nombre de qqseq por matching_allele_seq para unificar (aunque realmente era lo mismo porque había adjudicad allele_sequence a qseq)
                                        match_alignment_dict[core_name][sample_name] = get_alignment (new_sseq, matching_allele_seq, blast_90) ### cambiando/modificando: sustituyendo por get alignment

                                # convert the sequence to protein
                                if not core_name in protein_dict :
                                    protein_dict[core_name] = {}
                                if not sample_name in protein_dict[core_name] :
                                    protein_dict[core_name][sample_name] = []
                                #protein_dict[core_name][sample_name] = nucleotide_to_protein_alignment(new_sseq, matching_allele_seq)          ### (cambiando/modificando: he cambiado el nombre de qqseq por matching_allele_seq para unificar (aunque realmente era lo mismo porque había adjudicad allele_sequence a qseq)
                                protein_dict[core_name][sample_name] = get_alignment(new_sseq, matching_allele_seq, blast_90, "protein") ### cambiando/modificando: sustituyendo por get alignment                     

                                ##### PARA SACAR INFORME DE PRODIGAL VS BLAST #####
                                sseq_no_gaps = sseq.replace('-', '') ### para sacar report prodigal
                                prodigal_report.append([core_name, sample_name, qseqid, alm_allele, sstart, send, start_prodigal, end_prodigal, sseq_no_gaps, complete_predicted_seq])




                            """
                            ### CLASIFICACION DELECION/INSERCION/EQUAL ###
                            if len(sseq) < matching_allele_length: ### modificando/cambiando: Estoy añadiendo este if para clasificar en deletion e insertion en función de la longitud de la secuencia encontrada con BLAST (eliminando gaps, si los hubiere) con respecto a la longitud del alelo que hizo match
                                                                ### modificando/cambiando: realmente la clasificación en deletion e insertion puede ser ambigua porque se puede estar clasificando una secuencia como deletion pero realmente haya ocurrido inserciones y deleciones, solo que la deleción ha sido más grande y por tanto la sseq sería finalmente más corta que la longitud del alelo con el que hizo match
                        
                                print("Ha entrado a deleciones", '\n')

                                ### adding ASM allele to the asm_allele_matrix if it is not already include
                                if not core_name in deletions_dict:
                                    deletions_dict[core_name] = []
                                if not new_sseq in deletions_dict[core_name] :
                                    deletions_dict[core_name].append(new_sseq)
                                ### find the index of ASM  to include it in the sample matrix dict
                                index_delete = deletions_dict[core_name].index(new_sseq)

                                ### Cambiando/modificando: clasificando en ASM o ALM dentro de "deletions" en función de si la longitud de new_sseq es menor al umbral mínimo de longitud o mayor al umbral máximo de longitud
                                #if new_sequence_length < query_length : ### deleción
                                #if new_sequence_length < min(schema_variability[core_name]) : ### inserción
                                if new_sequence_length < min_length_threshold:
                                    delete_allele = 'ASM_DELETE_' + core_name + '_' + str(qseqid) + '_' + str(index_delete)
                                #if new_sequence_length > query_length : ### deleción
                                #if new_sequence_length > max(schema_variability[core_name]) : ### inserción
                                if new_sequence_length > max_length_threshold:
                                    delete_allele = 'ALM_DELETE_' + core_name + '_' + str(qseqid) + '_' + str(index_delete)
                                samples_matrix_dict[sample_name].append(delete_allele)

                                if not sseqid in matching_genes_dict[sample_name] :
                                    matching_genes_dict[sample_name][sseqid] = []
                                if sstart > send :
                                    matching_genes_dict[sample_name][sseqid].append([core_name, str(int(sstart)-new_sequence_length -1), sstart,'-', delete_allele]) ### modificando/cambiando: al hacer predecir con prodigal habría que incluir sstart y send obteniéndolos de los archivos que genera prodigal para saber dónde están las coordenadas de la nueva seq completa
                                else:
                                    matching_genes_dict[sample_name][sseqid].append([core_name, sstart,str(int(sstart)+ new_sequence_length),'+', delete_allele])

                                print("DELETION: ", delete_allele, '\n') ### printeando, borrar

                                ### add the deletion into deletion list
                                if not core_name in list_deletions :
                                    list_deletions [core_name] = {}
                                if not sample_name in list_deletions[core_name] :
                                    list_deletions[core_name][sample_name] = {}

                                len_sseq = len(sseq) ### cambiando/modificando: obteniendo longitud sseq sin gaps, ya que s_length incluye gaps en caso de que los haya en el alineamiento
                                list_deletions[core_name][sample_name][delete_allele] = [qseqid, allele_quality, sseqid,  bitscore, str(matching_allele_length), str(len_sseq), str(new_sequence_length), mismatch , gapopen, sstart, send,  new_sseq ]
                                ### duda cambiando/modificando: no debería cambiar en list_deletions sstart y send por las sstart y send de la nueva secuencia completa? Ahora mismo se está guardando sstart y send de la secuencia sseq encontrada con BLAST

                                #########################
                                ### OBTENCIÓN DE SNPs ### Añadir función get_SNPs como en INF, además de sustituir get_alignment_for_deletions y nucleotide_to_protein_alignment por get_alignment
                                #########################

                                if check_sequence_order(matching_allele_seq, logger) == 'reverse':
                                    matching_allele_seq = reverse_complement(matching_allele_seq)           ### (cambiando/modificando: he cambiado allele_sequence por matching_allele_seq)
                                                                                                        ### (cambiando/modificando: he cambiado str(matching_allele_seq).reverse_complement() por reverse_complement(matching_allele_seq), ya que no puedo emplear .reverse_complement() sobre una string, sino que tiene que ser sobre un objeto y al obtener la seq del diccionario core_dict se trata de una seq. Para que sea un objeto tendría que parsear. De momento he incluido la función reverse_complement para que funcione con string. Ver cómo lo dejo finalmente, como string (entonces tendría que extraer también el contig en caso de necesitarlo en la versión final a partir del ddiccionario de contigs) o parseando (como está ahora mismo la obtención de accession_sequence y en ese caso sería 0 necesario la generación de los archivos con los diccionarios para cada muestra y cada locus, al menos tal y como se generan ahora)
                                # get the SNP for the  delection
                                #if not core_name in snp_dict :
                                #    snp_dict[core_name] = {}
                                #if not sample_name in snp_dict[core_name] :
                                #    snp_dict[core_name][sample_name] = []
                                #snp_dict[core_name][sample_name] = get_snp(new_sseq, matching_allele_seq)


                                # execute again blast with the reference query the previous query found to get the aligment format to get the SNPs
                                if not core_name in match_alignment_dict :
                                    match_alignment_dict[core_name] = {}
                                    if not sample_name in match_alignment_dict[core_name] :
                                        #match_alignment_dict[core_name][sample_name] = get_alignment_for_deletions (new_sseq,  matching_allele_seq)        ### (cambiando/modificando: he cambiado el nombre de qqseq por matching_allele_seq para unificar (aunque realmente era lo mismo porque había adjudicad allele_sequence a qseq)
                                        match_alignment_dict[core_name][sample_name] = get_alignment (new_sseq, matching_allele_seq, blast_90) ### cambiando/modificando: sustituyendo por get alignment

                                # convert the sequence to protein
                                if not core_name in protein_dict :
                                    protein_dict[core_name] = {}
                                if not sample_name in protein_dict[core_name] :
                                    protein_dict[core_name][sample_name] = []
                                #protein_dict[core_name][sample_name] = nucleotide_to_protein_alignment(new_sseq, matching_allele_seq)          ### (cambiando/modificando: he cambiado el nombre de qqseq por matching_allele_seq para unificar (aunque realmente era lo mismo porque había adjudicad allele_sequence a qseq)
                                protein_dict[core_name][sample_name] = get_alignment(new_sseq, matching_allele_seq, blast_90, "protein") ### cambiando/modificando: sustituyendo por get alignment                     

                                ##### PARA SACAR INFORME DE PRODIGAL VS BLAST #####
                                sseq_no_gaps = sseq.replace('-', '') ### para sacar report prodigal
                                prodigal_report.append([core_name, sample_name, qseqid, delete_allele, sstart, send, start_prodigal, end_prodigal, sseq_no_gaps, complete_predicted_seq])

                            ### ASM/ALM equal

                            ##############################################################################################################################################
                            #### Probando 3ª opción: EQUAL para casos en los que la secuencia encontrada tras blast tiene la misma longitud que el alelo que hizo match                   
                            #### De este modo se solucionan problemas como el del locus lmo1500 en la muestra 4169, que no podía clasificarse ni en deleción ni en inserción 
                            #### porque su longitud no era ni menor ni mayor que la del alelo query
                            ##############################################################################################################################################

                            if len(sseq) == matching_allele_length: ### modificando/cambiando: Estoy añadiendo este if para clasificar en deletion e insertion en función de la longitud de la secuencia encontrada con BLAST (eliminando gaps, si los hubiere) con respecto a la longitud del alelo que hizo match

                                print("Ha entrado a equal", '\n')

                                ### adding ASM allele to the asm_allele_matrix if it is not already include
                                if not core_name in equal_dict:
                                    equal_dict[core_name] = []
                                if not new_sseq in equal_dict[core_name] :
                                    equal_dict[core_name].append(new_sseq)
                                ### find the index of ASM  to include it in the sample matrix dict
                                index_equal = equal_dict[core_name].index(new_sseq)

                                ### Cambiando/modificando: clasificando en ASM o ALM dentro de "deletions" en función de si la longitud de new_sseq es menor al umbral mínimo de longitud o mayor al umbral máximo de longitud
                                if new_sequence_length < min_length_threshold:
                                    equal_allele = 'ASM_EQUAL_' + core_name + '_' + str(qseqid) + '_' + str(index_equal)
                                if new_sequence_length > max_length_threshold:
                                    equal_allele = 'ALM_EQUAL_' + core_name + '_' + str(qseqid) + '_' + str(index_equal)
                                samples_matrix_dict[sample_name].append(equal_allele)

                                if not sseqid in matching_genes_dict[sample_name] :
                                    matching_genes_dict[sample_name][sseqid] = []
                                if sstart > send :
                                    ### modificando/cambiando: al hacer predecir con prodigal habría que incluir sstart y send obteniéndolos de los archivos que genera prodigal para saber dónde están las coordenadas de la nueva seq completa
                                    matching_genes_dict[sample_name][sseqid].append([core_name, str(int(sstart)-new_sequence_length -1), sstart,'-', equal_allele])
                                else:
                                    matching_genes_dict[sample_name][sseqid].append([core_name, sstart,str(int(sstart)+ new_sequence_length),'+', equal_allele])

                                print("EQUAL: ", equal_allele, '\n') ### printeando, borrar

                                ### add the deletion into deletion list
                                if not core_name in list_equal :
                                    list_equal [core_name] = {}
                                if not sample_name in list_equal[core_name] :
                                    list_equal[core_name][sample_name] = {}

                                len_sseq = len(sseq) ### cambiando/modificando: obteniendo longitud sseq sin gaps, ya que s_length incluye gaps en caso de que los haya en el alineamiento
                                list_equal[core_name][sample_name][equal_allele] = [qseqid, allele_quality, sseqid,  bitscore, str(matching_allele_length), str(len_sseq), str(new_sequence_length), mismatch , gapopen, sstart, send,  new_sseq ]
                                ### duda cambiando/modificando: no debería cambiar en list_deletions sstart y send por las sstart y send de la nueva secuencia completa? Ahora mismo se está guardando sstart y send de la secuencia sseq encontrada con BLAST

                                #########################
                                ### OBTENCIÓN DE SNPs ### Añadir función get_SNPs como en INF, además de sustituir get_alignment_for_deletions y nucleotide_to_protein_alignment por get_alignment
                                #########################

                                if check_sequence_order(matching_allele_seq, logger) == 'reverse':
                                    matching_allele_seq = reverse_complement(matching_allele_seq) ### cambiando/modificando: he cambiado allele_sequence por matching_allele_seq
                                                                                          ### cambiando/modificando: he cambiado str(matching_allele_seq).reverse_complement() por reverse_complement(matching_allele_seq), ya que no puedo emplear .reverse_complement() sobre una string, sino que tiene que ser sobre un objeto y al obtener la seq del diccionario core_dict se trata de una seq. Para que sea un objeto tendría que parsear. De momento he incluido la función reverse_complement para que funcione con string. Ver cómo lo dejo finalmente, como string (entonces tendría que extraer también el contig en caso de necesitarlo en la versión final a partir del ddiccionario de contigs) o parseando (como está ahora mismo la obtención de accession_sequence y en ese caso sería 0 necesario la generación de los archivos con los diccionarios para cada muestra y cada locus, al menos tal y como se generan ahora)
                                # get the SNP for the  delection
                                #if not core_name in snp_dict :
                                #    snp_dict[core_name] = {}
                                #if not sample_name in snp_dict[core_name] :
                                #    snp_dict[core_name][sample_name] = []
                                #snp_dict[core_name][sample_name] = get_snp(new_sseq, matching_allele_seq)


                                # execute again blast with the reference query the previous query found to get the aligment format to get the SNPs
                                if not core_name in match_alignment_dict :
                                    match_alignment_dict[core_name] = {}
                                    if not sample_name in match_alignment_dict[core_name] :
                                        #match_alignment_dict[core_name][sample_name] = get_alignment_for_deletions (new_sseq,  matching_allele_seq)           ### (cambiando/modificando: he cambiado el nombre de qqseq por matching_allele_seq para unificar (aunque realmente era lo mismo porque había adjudicad allele_sequence a qseq)
                                        match_alignment_dict[core_name][sample_name] = get_alignment (new_sseq, matching_allele_seq, blast_90) ### cambiando/modificando: sustituyendo por get alignment

                                # convert the sequence to protein
                                if not core_name in protein_dict :
                                    protein_dict[core_name] = {}
                                if not sample_name in protein_dict[core_name] :
                                    protein_dict[core_name][sample_name] = []
                                #protein_dict[core_name][sample_name] = nucleotide_to_protein_alignment(new_sseq, matching_allele_seq)              ### (cambiando/modificando: he cambiado el nombre de qqseq por matching_allele_seq para unificar (aunque realmente era lo mismo porque había adjudicad allele_sequence a qseq)
                                protein_dict[core_name][sample_name] = get_alignment (new_sseq, matching_allele_seq, blast_90, "protein") ### cambiando/modificando: sustituyendo por get alignment

                            
                            
                                ##### PARA SACAR INFORME DE PRODIGAL VS BLAST #####
                                sseq_no_gaps = sseq.replace('-', '') ### para sacar report prodigal
                                prodigal_report.append([core_name, sample_name, qseqid, equal_allele, sstart, send, start_prodigal, end_prodigal, sseq_no_gaps, complete_predicted_seq])


                            ### ASM/ALM insertion

                            print("Aquí vienen las inserciones", '\n')
                    
                            if len(sseq) > matching_allele_length: ### modificando/cambiando: Estoy añadiendo este if para clasificar en deletion e insertion en función de la longitud de la secuencia encontrada con BLAST (eliminando gaps, si los hubiere) con respecto a la longitud del alelo que hizo match
                                               
                                print("Ha entrado a inserción", '\n')
                        
                                if not core_name in insertions_dict:
                                    insertions_dict[core_name] = []
                                if not new_sseq in insertions_dict[core_name]:
                                    insertions_dict[core_name].append(new_sseq)
                    
                                ### find the index of ASM  to include it in the sample matrix dict
                                index_insert = insertions_dict[core_name].index(new_sseq)

                                ### Cambiando/modificando: clasificando en ASM o ALM dentro de "insertions" en función de si la longitud de new_sseq es menor al umbral mínimo de longitud o mayor al umbral máximo de longitud
                                if new_sequence_length < min_length_threshold:
                                    print("Se trata de un ASM", '\n') ### print, borrar

                                    insert_allele = 'ASM_INSERT_' + core_name + '_' + str(qseqid) + '_' + str(index_insert)
                                if new_sequence_length > max_length_threshold:
                                    insert_allele = 'ALM_INSERT_' + core_name + '_' + str(qseqid) + '_' + str(index_insert)
                                samples_matrix_dict[sample_name].append(insert_allele)

                                if not sseqid in matching_genes_dict[sample_name]:
                                    matching_genes_dict[sample_name][sseqid] = []
                                    if sstart > send :
                                        matching_genes_dict[sample_name][sseqid].append([core_name, str(int(sstart) - new_sequence_length - 1), sstart, '-', insert_allele])
                                    else:
                                        matching_genes_dict[sample_name][sseqid].append([core_name, sstart, str(int(sstart) + new_sequence_length),'+', insert_allele])
                        
                                ### add the insertion into insertion list
                                if not core_name in list_insertions:
                                    list_insertions [core_name] = {}
                                if not sample_name in list_insertions[core_name]:
                                    list_insertions[core_name][sample_name] = {}

                                len_sseq = len(sseq) ### cambiando/modificando: obteniendo longitud sseq sin gaps, ya que s_length incluye gaps en caso de que los haya en el alineamiento
                                list_insertions[core_name][sample_name][insert_allele] = [qseqid, allele_quality, sseqid,  bitscore, str(matching_allele_length), str(len_sseq), str(new_sequence_length), mismatch , gapopen, sstart, send,  new_sseq ]
                                ### duda cambiando/modificando: no debería cambiar en list_deletions sstart y send por las sstart y send de la nueva secuencia completa? Ahora mismo se está guardando sstart y send de la secuencia sseq encontrada con BLAST

                                print("INSERTION: ", insert_allele, '\n') ### printeando, borrar
                        
                                #########################
                                ### OBTENCIÓN DE SNPs ### Añadir función get_SNPs como en INF, además de sustituir get_alignment_for_indels y nucleotide_to_protein_alignment por get_alignment
                                #########################

                                if check_sequence_order(matching_allele_seq, logger) == 'reverse': ### modificando/cambiando: matching_allele_seq
                                    matching_allele_seq = reverse_complement(matching_allele_seq) ### cambiando/modificando: he cambiado str(matching_allele_seq).reverse_complement() por reverse_complement(matching_allele_seq), ya que no puedo emplear .reverse_complement() sobre una string, sino que tiene que ser sobre un objeto y al obtener la seq del diccionario core_dict se trata de una seq. Para que sea un objeto tendría que parsear. De momento he incluido la función reverse_complement para que funcione con string. Ver cómo lo dejo finalmente, como string (entonces tendría que extraer también el contig en caso de necesitarlo en la versión final a partir del ddiccionario de contigs) o parseando (como está ahora mismo la obtención de accession_sequence y en ese caso sería 0 necesario la generación de los archivos con los diccionarios para cada muestra y cada locus, al menos tal y como se generan ahora)

                                # get the SNP for the  delection
                                #if not core_name in snp_dict :
                                #    snp_dict[core_name] = {}
                                #if not sample_name in snp_dict[core_name] :
                                #    snp_dict[core_name][sample_name] = []
                                #snp_dict[core_name][sample_name] = get_snp(new_sseq, qseq)

                                if not core_name in match_alignment_dict :
                                    match_alignment_dict[core_name] = {}
                                if not sample_name in match_alignment_dict[core_name] :
                                    #match_alignment_dict[core_name][sample_name] = get_alignment_for_indels (blast_db_name, matching_allele_seq) ### modificando/cambiado: he cambiado qseq por matching_allele_seq igual que antes
                                    #match_alignment_dict[core_name][sample_name] = get_alignment_for_indels (new_sseq, matching_allele_seq)
                                    # index_not_match = [m.start() for m in re.finditer(' ', match.match)]
                                    match_alignment_dict[core_name][sample_name] = get_alignment (new_sseq, matching_allele_seq, blast_90) ### cambiando/modificando: sustituyendo por get alignment
                            
                                # convert the sequence to protein
                                if not core_name in protein_dict :
                                    protein_dict[core_name] = {}
                                if not sample_name in protein_dict[core_name] :
                                    #protein_dict[core_name][sample_name] = []
                                    #protein_dict[core_name][sample_name] = nucleotide_to_protein_alignment(new_sseq, matching_allele_seq) ### modificando/cambiado: he cambiado qseq por matching_allele_seq igual que antes
                                    protein_dict[core_name][sample_name] = get_alignment (new_sseq, matching_allele_seq, blast_90, "protein") ### cambiando/modificando: sustituyendo por get alignment


                                # print("El stop codon de esta secuencia es TGA: ", tga_stop_codon) ### printeando, borrar
                                # print("La posición del stop codon en la secuencia encontrada en la muestra ampliada es: ", stop_index, '\n', '\n', '\n') ### printeando, borrar
                                # print("La nueva secuencia obtenida habiendo encontrado el stop codon (new_sseq) es: ", new_sseq, '\n', '\n') ### printeando, borrar
                                # print("La longitud de la nueva secuencia obtenida habiendo encontrado el stop codon (new_sseq) es: ", new_sequence_length, '\n', '\n', '\n') ### printeando, borrar
                                # print("La lista de secuencias con inserciones encontradas para este locus contenidas en el diccionario insertions_dict es: ", insertions_dict,  '\n', '\n') ### printeando, borrar
                                # print("El nombre otorgado al alelo encontrado es: ", insert_allele, '\n', '\n', '\n') ### printeando, borrar
                                # print("El alineamiento entre la secuencia con deleción encontrada en la muestra habiendo encontrado su stop codon y el alelo con el que ha hecho match es (match_alignment_dict): ", match_alignment_dict) ### printeando, borrar
                                # print("Esta es la proteína obtenida a partir de la secuencia con inserción encontrada en la muestra (protein_dict): ", protein_dict) ### printeando, borrar
                                # print("list_insertions: ", list_insertions)
                

                                ##### PARA SACAR INFORME DE PRODIGAL VS BLAST #####
                                sseq_no_gaps = sseq.replace('-', '') ### para sacar report prodigal
                                prodigal_report.append([core_name, sample_name, qseqid, insert_allele, sstart, send, start_prodigal, end_prodigal, sseq_no_gaps, complete_predicted_seq])

                            """
                        print("Fin de la clasificación", '\n')
            
                    #################################################################################################################
                    #### Este else va con el if stop_codon_index != False, es decir, por si no se encuentra el codón de stop se guarda como error
                    #### Va después de INF, ASM y ALM
                    else:
                        # print("No se ha encontrado stop codon: ERROR not stop codon when insertion", '\n') ### printeando, BORRAR
                        logger.error('ERROR : Stop codon was not found for the core %s and the sample %s', core_name, sample_name)
                        samples_matrix_dict[sample_name].append('ERROR not stop codon')
                        if not sseqid in matching_genes_dict[sample_name] :
                            matching_genes_dict[sample_name][sseqid] = []
                        if sstart > send :
                            matching_genes_dict[sample_name][sseqid].append([core_name, sstart,send,'-', 'ERROR'])
                        else:
                            matching_genes_dict[sample_name][sseqid].append([core_name, sstart,send,'+', 'ERROR'])
                            #print ('ERROR when looking the allele match for core gene ', core_name, 'at sample ', sample_name )  ### printeando, BORRAR
                    #################################################################################################################


            ######################################################################################################
            ############# COMENTANDO ESTA PARTE DEL CÓDIGO, QUE ES EL QUE HABÍA ANTES DE REORGANIZAR #############
            ######################################################################################################
        """
            if int(s_length) in schema_variability[core_name] :

                logger.info('Found new allele for core gene %s ', core_name)
                if not sample_name in inf_dict :
                    inf_dict[sample_name] = {}

                ### adding new allele to the  inferred allele list if it is not already included
                if not core_name in inferred_alleles_dict :
                    inferred_alleles_dict[core_name] = []
                if not sseq in inferred_alleles_dict[core_name] :
                    inferred_alleles_dict[core_name].append(sseq)
                ### find the index to include in the sample matrix dict
                index_inferred = inferred_alleles_dict[core_name].index(sseq)
                inferred_allele = 'INF_' + core_name + '_' + str(index_inferred)
                samples_matrix_dict[sample_name].append(inferred_allele)
                inf_dict[sample_name][core_name] = inferred_allele

                # Get the SNP for the new allele inferred

                snp_information = get_snp(sseq, reference_allele_for_snp)
                if len(snp_information) > 0 :
                    if not core_name in snp_dict :
                        snp_dict[core_name] = {}
                    if not sample_name in snp_dict[core_name] :
                        snp_dict[core_name][sample_name] = {}
                    snp_dict[core_name][sample_name][qseqid]= snp_information
                
                if not sseqid in matching_genes_dict[sample_name] :
                    matching_genes_dict[sample_name][sseqid] = []
                if sstart > send :
                    matching_genes_dict[sample_name][sseqid].append([core_name, sstart,send,'-',inferred_allele])
                else:
                    matching_genes_dict[sample_name][sseqid].append([core_name, sstart,send,'+',inferred_allele])
                continue

            alleles_in_locus = list (SeqIO.parse(reference_query, "fasta"))
            for allele_item in alleles_in_locus :
                if allele_item.id == qseqid :
                    break
            allele_sequence = allele_item.seq

            # Retrieve the contig file for getting the contig sequence for the id found in Blast 
            contig_file = os.path.join(inputdir,str(sample_name + '.fasta'))
            records = list (SeqIO.parse(contig_file, "fasta"))

            for record in records:
                if record.id == sseqid :
                    break
            accession_sequence = record.seq

            if int(s_length) < int(query_length) :
                ## check if the blast alignment could be classified as PLOT
                length_sseqid = len(accession_sequence)
                if int(sstart) == length_sseqid or int(send) == length_sseqid or int(sstart) == 1 or int(send) == 1:
                    samples_matrix_dict[sample_name].append('PLOT_' + str(qseqid))
                    logger.info('PLOT found at sample %s, for gene  %s', sample_name, core_name)
                    if sample_name not in plot_dict :
                        plot_dict[sample_name] = {}
                    if not core_name in plot_dict[sample_name] :
                        plot_dict[sample_name][core_name] = []
                    plot_dict[sample_name][core_name].append([qseqid,sseqid,bitscore,sstart, send, sseq])
                    
                    if not sseqid in matching_genes_dict[sample_name] :
                        matching_genes_dict[sample_name][sseqid] = []
                    if sstart > send :
                        matching_genes_dict[sample_name][sseqid].append([core_name, sstart,send,'-', 'PLOT'])
                    else:
                        matching_genes_dict[sample_name][sseqid].append([core_name, sstart,send,'+', 'PLOT'])

                    continue
                else:

                    query_direction = check_sequence_order(allele_sequence, logger)
                    contig_file = os.path.join(inputdir,str(sample_name + '.fasta'))
                    records = list (SeqIO.parse(contig_file, "fasta"))

                    if allele_sequence.endswith ('TGA') or  allele_sequence.startswith ('TCA') :
                        tga_stop_codon = True
                    else:
                        tga_stop_codon = False

                    if int(query_length) > int(s_length):
                        difference_q_s_length = int(query_length)- int(s_length)
                    else:
                        difference_q_s_length = 0
                        
                    if query_direction == 'reverse' :
                        if int(send) > int (sstart): ## increasing the number of nucleotides to check if getting  longer protein
                            #sample_gene_sequence = accession_sequence[int(sstart) - 51 :  int(send)  ]
                            sample_gene_sequence = accession_sequence[int(sstart) - difference_q_s_length - 81 :  int(send)  ]
                            sample_gene_sequence = sample_gene_sequence.reverse_complement()
                        else:
                            #sample_gene_sequence = accession_sequence[int(send) -1 : int(sstart)  + 51]
                            sample_gene_sequence = accession_sequence[int(send) -1 : int(sstart) + difference_q_s_length + 81]
                            # import pdb; pdb.set_trace()
                    else:
                        if int(sstart) > int (send):
                            
                            #sample_gene_sequence = accession_sequence[int(send) - 51 :  int(sstart)  ]
                            sample_gene_sequence = accession_sequence[int(send) - difference_q_s_length - 81 :  int(sstart)  ]
                            sample_gene_sequence = sample_gene_sequence.reverse_complement()
                        else:
                            
                            #sample_gene_sequence = accession_sequence[int(sstart) -1 : int(send)  + 51]
                            sample_gene_sequence = accession_sequence[int(sstart) -1 : int(send)  + difference_q_s_length + 81]
                    
                    #stop_index = get_stop_codon_index(sample_gene_sequence, tga_stop_codon, int(qlen)- int(qstart))
                    stop_index = get_stop_codon_index(sample_gene_sequence, tga_stop_codon, int(s_length)- int(qstart))
                    if stop_index != False:
                        new_sequence_length = stop_index +3
                        new_sseq = str(sample_gene_sequence[0:new_sequence_length])

                        ### adding ASM allele to the asm_allele_matrix if it is not already include
                        if not core_name in deletions_dict :
                            deletions_dict[core_name] = []
                        if not new_sseq in deletions_dict[core_name] :
                            deletions_dict[core_name].append(new_sseq)
                        ### find the index of ASM  to include it in the sample matrix dict
                        index_delete = deletions_dict[core_name].index(new_sseq)
                        if new_sequence_length < query_length :
                            delete_allele = 'ASM_DELETE_' + core_name + '_' + str(qseqid) + '_' + str(index_delete)
                        elif new_sequence_length == query_length:
                            delete_allele = 'AEM_DELETE_' + core_name + '_' + str(qseqid) + '_' + str(index_delete)
                        else:
                            delete_allele = 'ALM_DELETE_' + core_name + '_' + str(qseqid) + '_' + str(index_delete)
                        samples_matrix_dict[sample_name].append(delete_allele)

                        if not sseqid in matching_genes_dict[sample_name] :
                            matching_genes_dict[sample_name][sseqid] = []
                        if sstart > send :
                            matching_genes_dict[sample_name][sseqid].append([core_name, str(int(sstart)-new_sequence_length -1), sstart,'-', delete_allele])
                        else:
                            matching_genes_dict[sample_name][sseqid].append([core_name, sstart,str(int(sstart)+ new_sequence_length),'+', delete_allele])

                        ### add the deletion into deletion list
                        if not core_name in list_deletions :
                            list_deletions [core_name] = {}
                        if not sample_name in list_deletions[core_name] :
                            list_deletions[core_name][sample_name] = {}
                        list_deletions[core_name][sample_name][delete_allele] = [qseqid, sseqid,  bitscore, str(query_length) , s_length, str(new_sequence_length), mismatch , gapopen, sstart, send,  new_sseq ]

                        if check_sequence_order(qseq, logger) == 'reverse' :
                            qseq = str(allele_sequence.reverse_complement())
                        else:
                            qseq = str(allele_sequence)
                        # get the SNP for the  delection
                        #if not core_name in snp_dict :
                        #    snp_dict[core_name] = {}
                        #if not sample_name in snp_dict[core_name] :
                        #    snp_dict[core_name][sample_name] = []
                        #snp_dict[core_name][sample_name] = get_snp(new_sseq, qseq)


                        # execute again blast with the reference query the previous query found to get the aligment format to get the SNPs
                        if not core_name in match_alignment_dict :
                            match_alignment_dict[core_name] = {}
                            if not sample_name in match_alignment_dict[core_name] :
                                match_alignment_dict[core_name][sample_name] = get_aligments_for_deletions (new_sseq,  str(qseq))

                        # convert the sequence to protein
                        if not core_name in protein_dict :
                            protein_dict[core_name] = {}
                        if not sample_name in protein_dict[core_name] :
                            protein_dict[core_name][sample_name] = []
                        protein_dict[core_name][sample_name] = nucleotide_to_protein_aligment(new_sseq, qseq )
                    else:
                        # import pdb; pdb.set_trace()
                        #print(sample_gene_sequence)
                        #print('\n')
                        #print(accession_sequence)
                        #print('\n\n')
                        logger.error('ERROR : Stop codon was not found for the core %s and the sample %s', core_name, sample_name)
                        samples_matrix_dict[sample_name].append('ERROR not stop codon when deletion')
                        if not sseqid in matching_genes_dict[sample_name] :
                            matching_genes_dict[sample_name][sseqid] = []
                        if sstart > send :
                            matching_genes_dict[sample_name][sseqid].append([core_name, sstart,send,'-', 'ERROR'])
                        else:
                            matching_genes_dict[sample_name][sseqid].append([core_name, sstart,send,'+', 'ERROR'])

            elif int(s_length) > int(query_length) :
                print("Muestra ", sample_name, '\n') ### printeando, borrar

                ### print("reference_allele_for_snp: ", reference_allele_for_snp) ### printeando, borrar
                ### print("query_length, longitud de reference_allele_for_snp: ", query_length) ### printeando, borrar

                print("Match con alelo ", qseqid, "del locus ", core_name, '\n') ### printeando, borrar

                print("schema variability: ", schema_variability)

                print("La secuencia del alelo query con id ", qseqid, " es (qseq): ", qseq, '\n', '\n') ### printeando, borrar
                print("Longitud de la secuencia del alelo query (qseq): ", len(qseq), qlen, query_length, '\n', '\n', '\n') ### printeando, borrar
                print("La secuencia encontrada en la muestra con BLAST con id ", sseqid, " es (sseq): ", sseq, '\n', '\n') ### printeando, borrar
                print("Longitud de la secuencia encontrada en la muestra con BLAST (sseq): ", len(sseq), s_length, '\n', '\n', '\n') ### printeando, borrar
                print("La posición inicial y final de la secuencia encontrada con BLAST son: ", sstart, "y", send, ", respectivamente", '\n')

                print("La secuencia se ha encontrado en el contig: ", sseqid) ### printeando, borrar
                print("La secuencia del contig donde se ha encontrado la secuencia es: MIRAR EN FASTA DE LA MUESTRA", '\n', '\n' ) ### printeando, borrar
                print("Longitud de la secuencia del contig donde se ha encontrado la secuencia: ", len(accession_sequence),  '\n', '\n', '\n') ### printeando, borrar


               # tga_stop_codon = qseq.endswith('TGA') ### MODIFICANDO: Silenciado.
                query_direction = check_sequence_order(allele_sequence, logger)
               # difference_s_q_length = int(s_length) - int(len(allele_sequence))

                ### MODIFICANDO: Añadido. inicio
                ###              Cambio de signo de difference_s_q_length a continuación en relación con DELECIONES
                ###              Para compensar el número de bases insertadas ya que puede ser un número no múltiplo de 3. De esta forma me aseguro de que la cantidad de bases que se alarga la secuencia contando con la inserción es múltiplo de 3
                if query_direction == 'reverse' :
                    if int(send) > int (sstart): ## increasing the number of nucleotides to check if getting  longer protein
                        #sample_gene_sequence = accession_sequence[int(sstart) - 51 :  int(send)  ]
                        sample_gene_sequence = accession_sequence[ : int(send) ]
                      #  print("La secuencia encontrada con BLAST habiendo incrementando el número de nts es: ", sample_gene_sequence, '\n', '\n') #### printeando, borrar
                      #  print("Longitud secuencia encontrada con BLAST habiendo incrementado nts: ", len(sample_gene_sequence),  '\n', '\n', '\n') ### printeando, borrar
                        sample_gene_sequence = sample_gene_sequence.reverse_complement()
                      #  print("La secuencia encontrada con BLAST habiendo incrementando el número de nts y habiendo hayado la sec complementaria si se trataba de reverse es: ", sample_gene_sequence, '\n', '\n') ### printeando, borrar
                      #  print("Longitud secuencia encontrada con BLAST habiendo incrementado nts: ", len(sample_gene_sequence), '\n', '\n', '\n') ### printeando, borrar
                      #  inicio = int(sstart) + difference_s_q_length - 81 ### printeando, borrar
                      #  final = int(send) ### printeando, borrar
                    else:
                        #sample_gene_sequence = accession_sequence[int(send) -1 : int(sstart)  + 51]
                        sample_gene_sequence = accession_sequence[ int(send) -1 : ]
                      #  print("La secuencia encontrada con BLAST habiendo incrementando el número de nts es: ", sample_gene_sequence, '\n', '\n') #### printeando, borrar
                      #  print("Longitud secuencia encontrada con BLAST habiendo incrementado nts: ", len(sample_gene_sequence), '\n', '\n', '\n') ### printeando, borrar
                        # import pdb; pdb.set_trace()
                      #  inicio = int(send) -1 ### printeando, borrar
                      #  final = int(sstart) + difference_s_q_length + 81 ### printeando, borrar
                      #  print("La diferencia en número de bases entre alelo query y secuencia encontrada con BLAST es de: ", difference_s_q_length, '\n', '\n', '\n') ### printeando, borrar

                else:
                    if int(sstart) > int (send):
                        #sample_gene_sequence = accession_sequence[int(send) - 51 :  int(sstart)  ]
                        sample_gene_sequence = accession_sequence[ :  int(sstart) ]
                      #  print("La secuencia encontrada con BLAST habiendo incrementando el número de nts es: ", sample_gene_sequence, '\n', '\n') #### printeando, borrar
                      #  print("Longitud secuencia encontrada con BLAST habiendo incrementado nts: ", len(sample_gene_sequence),  '\n', '\n', '\n') ### printeando, borrar
                        sample_gene_sequence = sample_gene_sequence.reverse_complement()
                      #  print("La secuencia encontrada con BLAST habiendo incrementando el número de nts y habiendo hayado la sec complementaria si se trataba de reverse es: ", sample_gene_sequence, '\n', '\n') ### printeando, borrar
                      #  print("Longitud secuencia encontrada con BLAST habiendo incrementado nts: ", len(sample_gene_sequence),  '\n', '\n', '\n') ### printeando, borrar
                      #  inicio = int(send) + difference_s_q_length - 81  ### printeando, borrar
                      #  final = int(sstart) ### printeando, borrar

                    else:
                        #sample_gene_sequence = accession_sequence[int(sstart) -1 : int(send)  + 51]
                        sample_gene_sequence = accession_sequence[ int(sstart) -1 : ]
                      #  print("La secuencia encontrada con BLAST habiendo incrementando el número de nts es: ", sample_gene_sequence, '\n', '\n') #### printeando, borrar
                      #  print("Longitud secuencia encontrada con BLAST habiendo incrementado nts: ", len(sample_gene_sequence), '\n', '\n', '\n') ### printeando, borrar
                      #  inicio = int(sstart) -1  ### printeando, borrar
                      #  final = int(send) - difference_s_q_length + 81 ### printeando, borrar
                ### MODIFICANDO: Añadido. final

                # sseq = sseq.replace('-','') ### MODIFICANDO: Silenciado.
                stop_index = get_stop_codon_index(sample_gene_sequence, allele_sequence, qseq.find('-'))
                ### MODIFICANDO: No uso sseq, que de todos modos no entiendo por qué quitarle los gaps si sseq no debería tener gaps al contar con inserción, si no más bien qseq, meto sample_gene_sequence en la función get_stop_codon_index, que es la secuencia encontrada ampliada para la búsqueda del codón
                ###             Además, borro el parámetro tga_stop_codon y no lo introduzco a la función get_stop_codon_index que la había modificado para quitarle este parámetro (al menos de momento)

                if stop_index != False:
                    new_sequence_length = stop_index +3
                    ### adding ASM allele to the asm_allele_matrix if it is not already include
                    new_sseq = str(sample_gene_sequence[0:new_sequence_length]) ### MODIFICANDO: cambio de sseq por sample_gene_sequence

                    if not core_name in insertions_dict :
                        insertions_dict[core_name] = []
                    if not new_sseq in insertions_dict[core_name] :
                        insertions_dict[core_name].append(new_sseq)
                    ### find the index of ASM  to include it in the sample matrix dict
                    index_insert = insertions_dict[core_name].index(new_sseq)
                    #if new_sequence_length < query_length :
                    if new_sequence_length < len(allele_sequence) :
                        insert_allele = 'ASM_INSERT_' + core_name + '_' + str(qseqid) + '_' + str(index_insert)
                    elif new_sequence_length == len(allele_sequence) :
                        insert_allele = 'AEM_INSERT_' + core_name + '_' + str(qseqid) + '_' + str(index_insert)
                    else:
                        insert_allele = 'ALM_INSERT_' + core_name + '_' + str(qseqid) + '_' + str(index_insert)
                    samples_matrix_dict[sample_name].append(insert_allele)
                
                    if not sseqid in matching_genes_dict[sample_name] :
                        matching_genes_dict[sample_name][sseqid] = []
                        ### MODIFICANDO: a continuación he cambiado el sstart y send y lo he puesto igual que en deleción, pero no sé si está bien, pensar
                        if sstart > send :
                            matching_genes_dict[sample_name][sseqid].append([core_name, str(int(sstart) - new_sequence_length - 1), sstart, '-', insert_allele])
                        else:
                            matching_genes_dict[sample_name][sseqid].append([core_name, sstart, str(int(sstart) + new_sequence_length),'+', insert_allele])
                    ### add the insertion into insertion list
                    if not core_name in list_insertions :
                        list_insertions [core_name] = {}
                    if not sample_name in list_insertions[core_name] :
                        list_insertions[core_name][sample_name] = {}
                    list_insertions[core_name][sample_name][insert_allele] = [qseqid, sseqid,  bitscore, str(query_length) , s_length, str(new_sequence_length), mismatch , gapopen, sstart, send,  new_sseq ]

                    if check_sequence_order(qseq, logger) == 'reverse' :
                        qseq = str(allele_sequence.reverse_complement())
                    else:
                        qseq = str(allele_sequence)
                    # get the SNP for the  delection
                    #if not core_name in snp_dict :
                    #    snp_dict[core_name] = {}
                    #if not sample_name in snp_dict[core_name] :
                    #    snp_dict[core_name][sample_name] = []
                    #snp_dict[core_name][sample_name] = get_snp(new_sseq, qseq)

                    if not core_name in match_alignment_dict :
                        match_alignment_dict[core_name] = {}
                    if not sample_name in match_alignment_dict[core_name] :
                        match_alignment_dict[core_name][sample_name] = get_alignment_for_indels (blast_db_name, qseq)
                    # index_not_match = [m.start() for m in re.finditer(' ', match.match)]

                    # convert the sequence to protein
                    if not core_name in protein_dict :
                        protein_dict[core_name] = {}
                    if not sample_name in protein_dict[core_name] :
                        #protein_dict[core_name][sample_name] = []
                        protein_dict[core_name][sample_name] = nucleotide_to_protein_aligment(new_sseq, qseq )

                    # get the SNP from the alignment
                    # print("El stop codon de esta secuencia es TGA: ", tga_stop_codon) ### printeando, borrar
                    print("La posición del stop codon en la secuencia encontrada en la muestra ampliada es: ", stop_index, '\n', '\n', '\n') ### printeando, borrar
                    print("La nueva secuencia obtenida habiendo encontrado el stop codon (new_sseq) es: ", new_sseq, '\n', '\n') ### printeando, borrar
                    print("La longitud de la nueva secuencia obtenida habiendo encontrado el stop codon (new_sseq) es: ", new_sequence_length, '\n', '\n', '\n') ### printeando, borrar
                    print("La lista de secuencias con inserciones encontradas para este locus contenidas en el diccionario insertions_dict es: ", insertions_dict,  '\n', '\n') ### printeando, borrar
                    print("El nombre otorgado al alelo encontrado es: ", insert_allele, '\n', '\n', '\n') ### printeando, borrar
                    print("El alineamiento entre la secuencia con deleción encontrada en la muestra habiendo encontrado su stop codon y el alelo con el que ha hecho match es (match_alignment_dict): ", match_alignment_dict) ### printeando, borrar
                    print("Esta es la proteína obtenida a partir de la secuencia con inserción encontrada en la muestra (protein_dict): ", protein_dict) ### printeando, borrar
                    print("list_insertions: ", list_insertions)

                else:
                    print("No se ha encontrado stop codon: ERROR not stop codon when insertion", '\n') ### printeando, BORRAR
                    logger.error('ERROR : Stop codon was not found for the core %s and the sample %s', core_name, sample_name)
                    samples_matrix_dict[sample_name].append('ERROR not stop codon when insertion')
                    if not sseqid in matching_genes_dict[sample_name] :
                        matching_genes_dict[sample_name][sseqid] = []
                    if sstart > send :
                        matching_genes_dict[sample_name][sseqid].append([core_name, sstart,send,'-', 'ERROR'])
                    else:
                        matching_genes_dict[sample_name][sseqid].append([core_name, sstart,send,'+', 'ERROR'])

                        print ('ERROR when looking the allele match for core gene ', core_name, 'at sample ', sample_name )
        """


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
    logger.info('Saving insert bases information to file..')
    exact_file =  os.path.join(outputdir, 'exact.tsv')
    with open (exact_file , 'w') as exact_fh :
        exact_fh.write('\t'.join(header_exact)+ '\n')
        for core in exact_dict:
            for sample in exact_dict[core]:
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
    '''
    inf_file = os.path.join(outputdir, 'infe_l.tsv')
    with open (inf_file , 'w') as inf_fh :
        for key in sorted (inf_dict) :
            inf_in_sample = '\t'.join (inf_dict[key])
            infer_fh.write(key + '\t' + inf_in_sample + '\n')
    '''
    # saving matching contigs to file
    logger.info('Saving matching information to file..')
    matching_file =  os.path.join(outputdir, 'matching_contigs.tsv')
    with open (matching_file , 'w') as matching_fh :
        matching_fh.write('\t'.join(header_macthing_alleles_conting ) + '\n')
        for samples in sorted ( matching_genes_dict) :
            for contigs in matching_genes_dict[samples] :
                for contig in matching_genes_dict[samples] [contigs]:
                #matching_alleles = '\t'.join (matching_genes_dict[samples][contig])
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
                lnf_fh.write(core + '\t' + sample + '\t' + '\t'.join(lnf_dict[core][sample]) + '\n')


    ### QUITAR al quitar clasificacion en delecion/insercion/equal
    # saving insertions bases in contigs to file
    logger.info('Saving insert bases  information to file..')
    insertions_file =  os.path.join(outputdir, 'insertions.tsv')
    with open (insertions_file , 'w') as insertions_fh :
        insertions_fh.write('\t'.join(header_insertions )+ '\n')
        for core in list_insertions :
            for sample in list_insertions[core] :
                for insert in list_insertions[core][sample]:
                    insertions_fh.write(core + '\t' + sample + '\t' + insert + '\t' + '\t'.join(list_insertions[core][sample][insert]) + '\n')

    ### QUITAR al quitar clasificacion en delecion/insercion/equal
    # saving deletions bases in contigs to file
    logger.info('Saving deleted bases information to file..')
    deletions_file =  os.path.join(outputdir, 'deletions.tsv')
    with open (deletions_file , 'w') as deletions_fh :
        deletions_fh.write('\t'.join(header_deletions) + '\n')
        for core in list_deletions :
            for sample in list_deletions[core] :
                for delete in list_deletions[core][sample]:
                    deletions_fh.write(core + '\t' + sample + '\t' + delete + '\t' + '\t'.join(list_deletions[core][sample][delete]) + '\n')

   # saving PLOT  to file
    logger.info('Saving PLOT information to file..')
    plot_file =  os.path.join(outputdir, 'plot.tsv')
    with open (plot_file , 'w') as plot_fh :
        plot_fh.write('\t'.join(header_plot) + '\n')
        for sample in sorted (plot_dict) :
            for core in sorted (plot_dict[sample]):
                for plot in plot_dict[sample][core] :
                    plot_fh.write(sample + '\t' + core + '\t' + '\t'.join (plot) + '\n')

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
    # open log file
    logger = open_log ('taranis_wgMLST.log')
    print('Checking the pre-requisites./n')
    # check additional programs are installed in your system
    if not check_prerequisites (logger):
        print ('your system does not fulfill the pre-requistes to run the script ')
        exit(0)
    ##############################################
    # Check that given directories contain fasta files
    ##############################################
    print('Validating schema fasta files in ' , arguments.coregenedir , '\n')
    valid_core_gene_files = get_fasta_file_list(arguments.coregenedir, logger)
    if not valid_core_gene_files :
        print ('There are not valid fasta files in ',  arguments.coregenedir , ' directory. Check log file for more information ')
        exit(0)
    print('Validating sample fasta files in ' , arguments.inputdir , '\n')
    valid_sample_files = get_fasta_file_list(arguments.inputdir, logger)
    if not valid_sample_files :
        print ('There are not valid fasta files in ',  arguments.inputdir , ' directory. Check log file for more information ')
        exit(0)
    ###############################
    # Prepare the coreMLST schema .
    ###############################
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

    core_gene_dict_files, schema_variability, schema_statistics, schema_quality = prepare_core_gene (valid_core_gene_files , tmp_core_gene_dir , arguments.outputdir, logger) ### cambiando/modificando: he quitado core_first_alleles_files, relacionado con el primer alelo del locus que se tomaba antes como referencia, he añadido schema_quality
    if not core_gene_dict_files :
        print('There is an error while processing the schema preparation phase. Check the log file to get more information \n')
        logger.info('Deleting the temporary directory to clean up the temporary files created')
        shutil.rmtree(os.path.join(arguments.outputdir, 'tmp'))
        exit(0)

    #######################################################
    # Prepare the samples files
    #######################################################
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
    sample_dict_files = prepare_samples(valid_sample_files, tmp_samples_dir, arguments.refgenome, logger) ### modificando/cambiando: He introducido arguments.refgenome a prepare_samples para que se lo pase a prodigal_training y genere el training file para la predicción de genes
    if not sample_dict_files :
        print('There is an error while processing the saving temporary files. Check the log file to get more information \n')
        logger.info('Deleting the temporary directory to clean up the temporary files created')
        shutil.rmtree(os.path.join(arguments.outputdir, 'tmp'))
        exit(0)

    ### core schema annotation  #########
    '''


    annotation_core_list =[]
    for annotation_core_file in core_first_alleles_files :
        annotation_core_list.append(get_gene_annotation (annotation_core_file, tmp_core_gene_dir, logger))
    ### sample  annotation    ##############
    annotation_sample_list =[]
    for annotation_sample_file in valid_sample_files :
        annotation_sample_list.append(get_gene_annotation (annotation_sample_file, tmp_samples_dir, logger))


    '''
    #query_directory = os.path.join(tmp_core_gene_dir,'first_alleles')
    ##########   Modified to get all alleles instead of the first one  #############
    query_directory = arguments.coregenedir
    blast_db_directory = os.path.join(tmp_samples_dir,'blastdb')
    prodigal_directory = os.path.join(tmp_samples_dir,'prodigal') ### cambiando/modificando: he incluido el path a prodigal
    if not allele_call_nucleotides( core_gene_dict_files, query_directory, sample_dict_files, blast_db_directory, prodigal_directory, arguments.inputdir, arguments.outputdir,  int(arguments.cpus), arguments.percentlength, schema_variability, schema_statistics, schema_quality, logger): ### CAMBIANDO/MODIFICANDO: He añadido schema_statistics, path a prodigal, prodigal training y schema_quality
        print('There is an error while processing the allele calling. Check the log file to get more information \n')
        exit(0)
    # Create the distance matrix
    try:
        
        print ('Creating matrix distance\n')
        create_distance_matrix(arguments.outputdir, 'result_for_tree_diagram.tsv')
    except:
        print('There was an error when creating distance matrix\n')
    end_time = datetime.now()
    print('completed execution at :', end_time )
    return True



