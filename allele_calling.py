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
#from BCBio import GFF
import pandas as pd
import shutil
from progressbar import ProgressBar

from utils.taranis_utils import *

import math 
import csv 

def check_blast (reference_allele, sample_files, db_name, logger) : ## N
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

# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · #  
# Parse samples and core genes schema fasta files to dictionary #
# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · #

def parsing_fasta_file_to_dict (fasta_file, logger):
    fasta_dict = {}
    fasta_dict_ordered = {}
    for contig in SeqIO.parse(fasta_file, "fasta", generic_dna):
        fasta_dict[str(contig.id)] = str(contig.seq.upper())
    logger.debug('file %s parsed to dictionary', fasta_file)

    for key in sorted(list(fasta_dict.keys())):
        fasta_dict_ordered[key] = fasta_dict[key]
    return fasta_dict_ordered


# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · #  
# Get core genes schema info before allele calling analysis #
# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · #

#def prepare_core_gene (core_gene_file_list, store_dir, ref_alleles_dir, logger):
def prepare_core_gene (core_gene_file_list, store_dir, ref_alleles_dir, genus, species, usegenus, logger):

    ## Initialize dict for keeping id-allele, quality, length variability, length statistics and annotation info for each schema core gene
    alleles_in_locus_dict = {}
    schema_quality = {} 
    annotation_core_dict = {}
    schema_variability = {}
    schema_statistics = {}


    ## Process each schema core gene
    blast_dir = os.path.join(store_dir,'blastdb')
    logger.info('start preparation of core genes files')
    for fasta_file in core_gene_file_list:

        f_name = os.path.basename(fasta_file).split('.')

        # Parse core gene fasta file and keep id-sequence info in dictionary
        fasta_file_parsed_dict = parsing_fasta_file_to_dict(fasta_file, logger)
        if f_name[0] not in alleles_in_locus_dict.keys():
            alleles_in_locus_dict[f_name[0]] = {}
        alleles_in_locus_dict[f_name[0]] = fasta_file_parsed_dict

        # dump fasta file into pickle file
        #with open (file_list[-1],'wb') as f:
         #   pickle.dump(fasta_file_parsed_dict, f)
        
        # Get core gene alleles quality 
        locus_quality = check_core_gene_quality(fasta_file, logger)  
        if f_name[0] not in schema_quality.keys(): 
            schema_quality[f_name[0]] = {}
        schema_quality[f_name[0]] = locus_quality

        # Get gene and product annotation for core gene using reference allele(s)
        ref_allele = os.path.join(ref_alleles_dir, f_name[0] + '.fasta')

        gene_annot, product_annot = get_gene_annotation (ref_allele, store_dir, genus, species, usegenus, logger)
        #gene_annot, product_annot = get_gene_annotation (ref_allele, store_dir, logger)
        if f_name[0] not in annotation_core_dict.keys():
            annotation_core_dict[f_name[0]] = {}
        annotation_core_dict[f_name[0]] = [gene_annot, product_annot]

        # Get core gene alleles length to keep length variability and statistics info
        alleles_len = []
        for allele in fasta_file_parsed_dict :
            alleles_len.append(len(fasta_file_parsed_dict[allele]))

        #alleles_in_locus = list (SeqIO.parse(fasta_file, "fasta")) ## parse
        #for allele in alleles_in_locus : ## parse
            #alleles_len.append(len(str(allele.seq))) ## parse

        schema_variability[f_name[0]]=list(set(alleles_len))
        
        if len(alleles_len) == 1:
            stdev = 0
        else:
            stdev = statistics.stdev(alleles_len)
        schema_statistics[f_name[0]]=[statistics.mode(alleles_len), statistics.mean(alleles_len), stdev, min(alleles_len), max(alleles_len)]

    return alleles_in_locus_dict, annotation_core_dict, schema_variability, schema_statistics, schema_quality


# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · #  
# Get Prodigal training file from reference genome for samples gene prediction  #
# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · #

def prodigal_training(reference_genome_file, prodigal_dir, logger):

    f_name = os.path.basename(reference_genome_file).split('.')[0]
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
    
      #  if prodigal_result.stderr:
       #     logger.error('cannot create training file for %s', f_name)
        #    logger.error('prodigal returning error code %s', prodigal_result.stderr)
         #   return False
    else:
        logger.info('Skeeping prodigal training file creation for %s, as it has already been created', f_name)
    
    return output_prodigal_train_dir


# · * · * · * · * · * · * · * · * · * #  
# Get Prodigal sample gene prediction #
# · * · * · * · * · * · * · * · * · * #

def prodigal_prediction(file_name, prodigal_dir, prodigal_train_dir, logger):

    f_name = os.path.basename(file_name).split('.')[0]
    prodigal_dir_sample = os.path.join(prodigal_dir,f_name)

    #output_prodigal_coord = os.path.join(prodigal_dir_sample, f_name + '_coord.gff')
    #output_prodigal_prot = os.path.join(prodigal_dir_sample, f_name + '_prot.faa')
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
        
        # if prodigal_result.stderr:
          #  logger.error('cannot predict genes for %s ', f_name)
           # logger.error('prodigal returning error code %s', prodigal_result.stderr)
            #return False
    else:
        logger.info('Skeeping prodigal genes prediction for %s, as it has already been made', f_name)
    
    return True


# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · *  #  
# Get Prodigal predicted gene sequence equivalent to BLAST result matching bad quality allele or to no Exact Match BLAST result in allele calling analysis #
# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · *  #

def get_prodigal_sequence(blast_sseq, contig_blast_id, prodigal_directory, sample_name, blast_parameters, logger): 

    prodigal_directory_sample = os.path.join(prodigal_directory, sample_name)
    genes_file = os.path.join(prodigal_directory_sample, sample_name + '_dna.faa')

    ## Create directory for storing prodigal genes prediction per contig BLAST databases 
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

    ## Create directory for storing prodigal genes prediction sequences per contig
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

    ## Parse prodigal genes prediction fasta file
    predicted_genes = SeqIO.parse(genes_file, "fasta")

    ## Create fasta file containing Prodigal predicted genes sequences for X contig in sample
    contig_genes_path = os.path.join(full_path_prodigal_genes_per_contig, contig_blast_id + '.fasta')
    with open (contig_genes_path, 'w') as out_fh:
        for rec in predicted_genes:
            contig_prodigal_id = (rec.id).split("_")[0]
            if contig_prodigal_id == contig_blast_id: 
                out_fh.write ('>' + str(rec.description) + '\n' + str(rec.seq) + '\n')

    ## Create local BLAST database for Prodigal predicted genes sequences for X contig in sample 
    if not create_blastdb(contig_genes_path, full_path_blastdb_per_contig, 'nucl', logger):
        print('Error when creating the blastdb for samples files. Check log file for more information. \n ')
        return False

    ## Local BLAST Prodigal predicted genes sequences database VS BLAST sequence obtained in sample in allele calling analysis
    blast_db_name = os.path.join(full_path_blastdb_per_contig, contig_blast_id, contig_blast_id)

    cline = NcbiblastnCommandline(db=blast_db_name, evalue=0.001, perc_identity = 90, outfmt= blast_parameters, max_target_seqs=10, max_hsps=10, num_threads=1)    
    out, err = cline(stdin = blast_sseq) 
    out_lines = out.splitlines()

    bigger_bitscore = 0
    if len (out_lines) > 0 :                  
        for line in out_lines :
            values = line.split('\t')
            if  float(values[8]) > bigger_bitscore: 
                qseqid , sseqid , pident ,  qlen , s_length , mismatch , r_gapopen , r_evalue , bitscore , sstart , send , qstart , qend ,sseq , qseq = values
                bigger_bitscore = float(bitscore)

        ## Get complete Prodigal sequence matching allele calling BLAST sequence using ID
        predicted_genes_in_contig = SeqIO.parse(contig_genes_path, "fasta")

        for rec in predicted_genes_in_contig: 
            if rec.id == sseqid:
                predicted_gene_sequence = str(rec.seq)
                start_prodigal = str(rec.description.split( '#')[1]) 
                end_prodigal = str(rec.description.split('#')[2]) 
                break

    ## Sequence not found by Prodigal when there are no BLAST results matching allele calling BLAST sequence
    if len (out_lines) == 0:
        predicted_gene_sequence = 'Sequence not found by Prodigal' 
        start_prodigal = '-' 
        end_prodigal = '-'

    return predicted_gene_sequence, start_prodigal, end_prodigal ### start_prodigal y end_prodigal para report prodigal


# · * · * · * · * · * · * · * · * · * · * · * · * #  
# Get samples info before allele calling analysis #
# · * · * · * · * · * · * · * · * · * · * · * · * #

def prepare_samples(sample_file_list, store_dir, reference_genome_file, logger):
    
    ## Initialize dictionary for keeping id-contig
    contigs_in_sample_dict = {} 

    ## Paths for samples blastdb, Prodigal genes prediction and BLAST results
    blast_dir = os.path.join(store_dir,'blastdb')
    prodigal_dir = os.path.join(store_dir,'prodigal')
    blast_results_seq_dir = os.path.join(store_dir,'blast_results', 'blast_results_seq')

    ## Get training file for Prodigal genes prediction
    output_prodigal_train_dir = prodigal_training(reference_genome_file, prodigal_dir, logger)
    if not output_prodigal_train_dir:
        print('Error when creating training file for genes prediction. Check log file for more information. \n ')
        return False

    for fasta_file in sample_file_list:
        f_name = os.path.basename(fasta_file).split('.')

        # Get samples id-contig dictionary
        fasta_file_parsed_dict = parsing_fasta_file_to_dict(fasta_file, logger)
        if f_name[0] not in contigs_in_sample_dict.keys():
            contigs_in_sample_dict[f_name[0]] = {}
        contigs_in_sample_dict[f_name[0]] = fasta_file_parsed_dict 

        # dump fasta file into pickle file
        #with open (file_list[-1],'wb') as f: # generación de diccionarios de contigs para cada muestra
         #   pickle.dump(fasta_file_parsed_dict, f)
    
        # Create directory for storing BLAST results using reference allele(s)
        blast_results_seq_per_sample_dir = os.path.join(blast_results_seq_dir, f_name[0])
        
        if not os.path.exists(blast_results_seq_per_sample_dir):
            try:
                os.makedirs(blast_results_seq_per_sample_dir)
                logger.debug('Created blast results directory for sample %s', f_name[0])
            except:
                logger.info('Cannot create blast results directory for sample %s', f_name[0])
                print ('Error when creating the directory for blast results', blast_results_seq_per_sample_dir)
                exit(0)

        # Prodigal genes prediction for each sample 
        if not prodigal_prediction(fasta_file, prodigal_dir, output_prodigal_train_dir, logger):
            print('Error when predicting genes for samples files. Check log file for more information. \n ')
            return False

        # Create local BLAST db for each sample fasta file
        if not create_blastdb(fasta_file, blast_dir, 'nucl', logger):
            print('Error when creating the blastdb for samples files. Check log file for more information. \n ')
            return False

    return contigs_in_sample_dict

# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * #  
# Get established length thresholds for allele tagging in allele calling analysis #
# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * #

def length_thresholds(core_name, schema_statistics, percent): ### logger
                                                                    
    locus_mean = int(schema_statistics[core_name][1])

    if percent != "SD": 
        max_length_threshold = math.ceil(locus_mean + ((locus_mean * float(percent)) / 100))
        min_length_threshold = math.floor(locus_mean - ((locus_mean * float(percent)) / 100))
    else:
        percent = float(schema_statistics[core_name][2])

        max_length_threshold = math.ceil(locus_mean + (locus_mean * percent))
        min_length_threshold = math.floor(locus_mean - (locus_mean * percent))

    return max_length_threshold, min_length_threshold


# · * · * · * · * · * · * · * · * · * · * · #  
# Convert dna sequence to protein sequence  #
# · * · * · * · * · * · * · * · * · * · * · #

def convert_to_protein (sequence) :

    seq = Seq.Seq(sequence)
    protein = str(seq.translate())

    return protein

# · * · * · * · * · * · * · * · * · * · * · * · * · * #  
# Get SNPs between BLAST sequence and matching allele #
# · * · * · * · * · * · * · * · * · * · * · * · * · * #

def get_snp (sample, query) :

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
    

def nucleotide_to_protein_alignment (sample_seq, query_seq ) : ### Sustituido por get_alignment
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


def get_alignment_for_indels (blast_db_name, qseq) : ### Sustituido por get_alignment
    #match_alignment =[]
    cline = NcbiblastnCommandline(db=blast_db_name, evalue=0.001, perc_identity = 80, outfmt= 5, max_target_seqs=10, max_hsps=10,num_threads=1)
    out, err = cline(stdin = qseq)
    psiblast_xml = StringIO(out)
    blast_records = NCBIXML.parse(psiblast_xml)
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for match in alignment.hsps:
                match_alignment = [['sample', match.sbjct],['match', match.match], ['schema',match.query]]
    return match_alignment


def get_alignment_for_deletions (sample_seq, query_seq): ### Sustituido por get_alignment
    index_found = False
    alignments = pairwise2.align.globalxx(sample_seq, query_seq)
    for index in range(len(alignments)) :
        if alignments[index][4] == len(query_seq) : 
            index_found = True
            break
    if not index_found : 
        index = 0
    values = format_alignment(*alignments[index]).split('\n')
    match_alignment = [['sample', values[0]],['match', values[1]], ['schema',values[2]]]
    return match_alignment


# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · *  #  
# Get DNA and protein alignment between the final sequence found in the sample and the matching allele #
# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · *  #

def get_alignment (sample_seq, query_seq, reward, penalty, gapopen, gapextend, seq_type = "dna"):
    
    ## If sequences alignment type desired is "protein" convert dna sequences to protein
    if seq_type == "protein":
        sample_seq = convert_to_protein(sample_seq)
        query_seq = convert_to_protein(query_seq)

    ## Get dna/protein alignment between final sequence found and matching allele
    # arguments pairwise2.align.globalms: match, mismatch, gap opening, gap extending
    alignments = pairwise2.align.localms(sample_seq, query_seq, reward, penalty, -gapopen, -gapextend)
    values = format_alignment(*alignments[0]).split('\n')
    match_alignment = [['sample', values[0]],['match', values[1]], ['schema',values[2]]]

    return match_alignment


# · * · * · * · * · * · * · * · * #  
# Tag LNF cases and keep LNF info #
# · * · * · * · * · * · * · * · * #

def lnf_tpr_tag(core_name, sample_name, alleles_in_locus_dict, samples_matrix_dict, lnf_tpr_dict, schema_statistics, locus_alleles_path, qseqid, pident, s_length, new_sequence_length, perc_identity_ref, coverage, annotation_core_dict, logger):

    gene_annot, product_annot = annotation_core_dict[core_name]

    if qseqid == '-':
        samples_matrix_dict[sample_name].append('LNF')
        tag_report = 'LNF'
        matching_allele_length = '-'
    
    else:
        if new_sequence_length == '-':
            samples_matrix_dict[sample_name].append('LNF_' + str(qseqid))
            tag_report = 'LNF'
        else:
            samples_matrix_dict[sample_name].append('TPR_' + str(qseqid))
            tag_report = 'TPR'

        matching_allele_seq = alleles_in_locus_dict[core_name][qseqid]
        matching_allele_length = len(matching_allele_seq) 

        #alleles_in_locus = list (SeqIO.parse(locus_alleles_path, "fasta")) ## parse
        #for allele in alleles_in_locus : ## parse
            #if allele.id == qseqid : ## parse
                #break ## parse
        #matching_allele_seq = str(allele.seq) ## parse
        #matching_allele_length = len(matching_allele_seq) ## parse

    if pident == '-':
        # (los dos BLAST sin resultado)
        coverage_blast = '-'
        coverage_new_sequence = '-'
        add_info = 'Locus not found'
        logger.info('Locus not found at sample %s, for gene %s', sample_name, core_name)

    elif 90 > float(pident): 
        # (BLAST 90 sin resultado y BLAST 70 con resultado)
        coverage_blast = '-'
        coverage_new_sequence = '-'
        add_info = 'BLAST sequence ID under threshold: {}%'.format(perc_identity_ref)
        logger.info('BLAST sequence ID %s under threshold at sample %s, for gene  %s', pident, sample_name, core_name)

    elif 90 <= float(pident) and new_sequence_length == '-':
        # (BLAST 90 con resultado, bajo coverage BLAST)
        locus_mean = int(schema_statistics[core_name][1]) 
        coverage_blast = int(s_length) / locus_mean 
        #coverage_blast = int(s_length) / matching_allele_length
        coverage_new_sequence = '-'
        if coverage_blast < 1:
            add_info = 'BLAST sequence coverage under threshold: {}%'.format(coverage)
        else:  
            add_info = 'BLAST sequence coverage above threshold: {}%'.format(coverage)
        logger.info('BLAST sequence coverage %s under threshold at sample %s, for gene  %s', coverage_blast, sample_name, core_name)

    elif 90 <= float(pident) and new_sequence_length != '-':
        # (BLAST 90 con resultado, buen coverage BLAST, bajo coverage new_sseq)
        locus_mean = int(schema_statistics[core_name][1]) 
        coverage_blast = int(s_length) / locus_mean
        #coverage_blast = int(s_length) / matching_allele_length  
        coverage_new_sequence = new_sequence_length / matching_allele_length 
        if coverage_new_sequence < 1:
            add_info = 'New sequence coverage under threshold: {}%'.format(coverage)
        else:  
            add_info = 'New sequence coverage above threshold: {}%'.format(coverage)
        logger.info('New sequence coverage %s under threshold at sample %s, for gene  %s', coverage_new_sequence, sample_name, core_name)

    ## Keeping LNF and TPR report info
    if not core_name in lnf_tpr_dict:
        lnf_tpr_dict[core_name] = {}
    if not sample_name in lnf_tpr_dict[core_name]:
        lnf_tpr_dict[core_name][sample_name] = []

    lnf_tpr_dict[core_name][sample_name].append([gene_annot, product_annot, tag_report, qseqid, pident, str(coverage_blast), str(coverage_new_sequence), str(matching_allele_length), str(s_length), str(new_sequence_length), add_info]) ### Meter secuencias alelo, blast y new_sseq (si las hay)?

    return True


# · * · * · * · * · * · * · * · * · * · * · * · * #  
# Tag paralog and exact match cases and keep info #
# · * · * · * · * · * · * · * · * · * · * · * · * #

def paralog_exact_tag(sample_name, core_name, tag, schema_quality, matching_genes_dict, samples_matrix_dict, allele_found, tag_dict, prodigal_report, prodigal_directory, blast_parameters, annotation_core_dict, logger):

    logger.info('Found %s at sample %s for core gene %s ', tag, sample_name, core_name)

    gene_annot, product_annot = annotation_core_dict[core_name]

    if not sample_name in tag_dict :
        tag_dict[sample_name] = {}
    if not core_name in tag_dict[sample_name] :
        tag_dict[sample_name][core_name]= []

    if tag == 'EXACT':
        allele = list(allele_found.keys())[0]
        qseqid = allele_found[allele][0]
        tag = qseqid

    samples_matrix_dict[sample_name].append(tag) 

    for sequence in allele_found:
        qseqid, sseqid, pident, qlen, s_length, mismatch, r_gapopen, r_evalue, bitscore, sstart, send, qstart, qend, sseq, qseq = allele_found[sequence]
        sseq = sseq.replace('-', '')

        # Get allele quality 
        allele_quality = schema_quality[core_name][qseqid]

        # Get prodigal gene prediction if allele quality is 'bad_quality'
        if 'bad_quality' in allele_quality: 
            complete_predicted_seq, start_prodigal, end_prodigal = get_prodigal_sequence(sseq, sseqid, prodigal_directory, sample_name, blast_parameters, logger)

            ##### informe prodigal #####
            prodigal_report.append([core_name, sample_name, qseqid, tag, sstart, send, start_prodigal, end_prodigal, sseq, complete_predicted_seq])

        else:
            complete_predicted_seq = '-'
        
        if not sseqid in matching_genes_dict[sample_name] :
            matching_genes_dict[sample_name][sseqid] = []
        if sstart > send :
            matching_genes_dict[sample_name][sseqid].append([core_name, sstart, send,'-', tag])
        else:
            matching_genes_dict[sample_name][sseqid].append([core_name, sstart, send,'+', tag])

        ## Keeping paralog NIPH/NIPHEM report info
        if tag == 'NIPH' or tag == 'NIPHEM':
            tag_dict[sample_name][core_name].append([gene_annot, product_annot, tag, pident, qseqid, allele_quality, sseqid, bitscore, sstart, send, sseq, complete_predicted_seq])
        else:
            tag_dict[sample_name][core_name] = [gene_annot, product_annot, qseqid, allele_quality, sseqid, s_length, sstart, send, sseq, complete_predicted_seq]

    return True 


# · * · * · * · * · * · * · * · * · * · *  #  
# Tag INF/ASM/ALM/PLOT cases and keep info #
# · * · * · * · * · * · * · * · * · * · *  #

def inf_asm_alm_tag(core_name, sample_name, tag, blast_values, allele_quality, new_sseq, matching_allele_length, tag_dict, list_tag, samples_matrix_dict, matching_genes_dict, prodigal_report, start_prodigal, end_prodigal, complete_predicted_seq, annotation_core_dict, logger): 

    gene_annot, product_annot = annotation_core_dict[core_name]

    qseqid, sseqid, pident,  qlen, s_length, mismatch, r_gapopen, r_evalue, bitscore, sstart, send, qstart, qend, sseq, qseq = blast_values

    sseq = sseq.replace('-', '')
    s_length = len(sseq)
    new_sequence_length = len(new_sseq)

    logger.info('Found %s at sample %s for core gene %s ', tag, sample_name, core_name)                    

    if tag == 'PLOT':
        tag_allele = tag + '_' + str(qseqid)
    else:
        # Adding ASM/ALM/INF allele to the allele_matrix if it is not already include
        if not core_name in tag_dict:
            tag_dict[core_name] = []
        if not new_sseq in tag_dict[core_name] :
            tag_dict[core_name].append(new_sseq)
        # Find the index of ASM/ALM/INF to include it in the sample matrix dict
        index_tag = tag_dict[core_name].index(new_sseq)

        tag_allele = tag + '_' + core_name + '_' + str(qseqid) + '_' + str(index_tag)
    
    samples_matrix_dict[sample_name].append(tag_allele)

    # Keeping INF/ASM/ALM/PLOT report info
    if not core_name in list_tag :
        list_tag[core_name] = {}
    if not sample_name in list_tag[core_name] :
        list_tag[core_name][sample_name] = {}   

    if tag == 'INF':
        list_tag[core_name][sample_name][tag_allele] = [gene_annot, product_annot, qseqid, allele_quality, sseqid,  bitscore, str(matching_allele_length), str(s_length), str(new_sequence_length), mismatch , r_gapopen, sstart, send,  new_sseq, complete_predicted_seq]
    
    elif tag == 'PLOT':
        list_tag[core_name][sample_name] = [gene_annot, product_annot, qseqid, allele_quality, sseqid, bitscore, sstart, send, sseq, new_sseq]

    else :
        if tag == 'ASM':
            newsseq_vs_blastseq = 'shorter'
        elif tag == 'ALM':
            newsseq_vs_blastseq = 'longer'

        if len(sseq) < matching_allele_length:
            add_info = 'Global effect: DELETION. BLAST sequence length shorter than matching allele sequence length / Net result: ' + tag + '. Final gene sequence length ' + newsseq_vs_blastseq + ' than matching allele sequence length'

        elif len(sseq) == matching_allele_length: 
            add_info = 'Global effect: BASE SUBSTITUTION. BLAST sequence length equal to matching allele sequence length / Net result: ' + tag + '. Final gene sequence length ' + newsseq_vs_blastseq + ' than matching allele sequence length'

        elif len(sseq) > matching_allele_length:
            add_info = 'Global effect: INSERTION. BLAST sequence length longer than matching allele sequence length / Net result: ' + tag + '. Final gene sequence length ' + newsseq_vs_blastseq + ' than matching allele sequence length'

        list_tag[core_name][sample_name][tag_allele] = [gene_annot, product_annot, qseqid, allele_quality, sseqid,  bitscore, str(matching_allele_length), str(s_length), str(new_sequence_length), mismatch , r_gapopen, sstart, send,  new_sseq, add_info, complete_predicted_seq]

    if not sseqid in matching_genes_dict[sample_name] :
        matching_genes_dict[sample_name][sseqid] = []
    if sstart > send :
        matching_genes_dict[sample_name][sseqid].append([core_name, str(int(sstart)-new_sequence_length -1), sstart,'-', tag_allele])
    else:
        matching_genes_dict[sample_name][sseqid].append([core_name, sstart,str(int(sstart)+ new_sequence_length),'+', tag_allele])

    ##### informe prodigal #####
    prodigal_report.append([core_name, sample_name, qseqid, tag_allele, sstart, send, start_prodigal, end_prodigal, sseq, complete_predicted_seq])

    return True


# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · #  
# Keep best results info after BLAST using results from previous reference allele BLAST as database VS ALL alleles in locus as query in allele calling analysis #
# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · #

def get_blast_results (sample_name, values, contigs_in_sample_dict, allele_found, logger) :
  
    qseqid, sseqid, pident, qlen, s_length, mismatch, r_gapopen, r_evalue, bitscore, sstart, send, qstart, qend, sseq, qseq = values

    ## Get contig ID and BLAST sequence
    sseqid_blast = sseqid.split('_')[1]
    sseq_no_gaps = sseq.replace('-', '') 

    
    ## Get start and end positions in contig searching for BLAST sequence index in contig sequence
    
    # Get contig sequence
    accession_sequence = contigs_in_sample_dict[sample_name][sseqid_blast]
    
    #for record in sample_contigs: ## parse
        #if record.id == sseqid_blast : ## parse
            #break ## parse
    #accession_sequence = str(record.seq) ## parse

    # Try to get BLAST sequence index in contig. If index -> error because different contig sequence and BLAST sequence 
    # direction, obtain reverse complement BLAST sequence and try again.
    try: 
        sseq_index_1 = int(accession_sequence.index(sseq_no_gaps)) + 1

    except:
        sseq_no_gaps = str(Seq.Seq(sseq_no_gaps).reverse_complement())
        sseq_index_1 = int(accession_sequence.index(sseq_no_gaps)) + 1

    sseq_index_2 = int(sseq_index_1) + len(sseq_no_gaps) - 1

    # Assign found indexes to start and end possitions depending on BLAST sequence and allele sequence direction
    if int(sstart) < int(send):
        sstart_new = str(min(sseq_index_1, sseq_index_2))
        send_new = str(max(sseq_index_1, sseq_index_2))
    else:
        sstart_new = str(max(sseq_index_1, sseq_index_2))
        send_new = str(min(sseq_index_1, sseq_index_2))


    ## Keep BLAST results info discarding subsets
    allele_is_subset = False

    if len(allele_found) > 0 :
        for allele_id in allele_found :
            min_index = min(int(allele_found[allele_id][9]), int(allele_found[allele_id][10]))
            max_index = max(int(allele_found[allele_id][9]), int(allele_found[allele_id][10]))
            if int(sstart_new) in range(min_index, max_index + 1) or  int(send_new) in range(min_index, max_index + 1): # if both genome locations overlap
                if sseqid_blast == allele_found[allele_id][1]: # if both sequences are in the same contig
                    logger.info('Found allele %s that starts or ends at the same position as %s ' , qseqid, allele_id)
                    allele_is_subset = True
                    break 

    if len(allele_found) == 0 or not allele_is_subset :
        contig_id_start = str(sseqid_blast + '_'+ sstart_new)
        
        # Skip the allele found in the 100% identity and 100% alignment
        if not contig_id_start in allele_found:
            allele_found[contig_id_start] = [qseqid, sseqid_blast, pident, qlen, s_length, mismatch, r_gapopen, r_evalue, bitscore, sstart_new, send_new, '-', '-', sseq, qseq]

    return True


# · * · * · * · * · * · * · * · * · * ·  #  
# Get SNPs and ADN and protein alignment #
# · * · * · * · * · * · * · * · * · * ·  #

def keep_snp_alignment_info(sseq, new_sseq, matching_allele_seq, qseqid, query_direction, core_name, sample_name, reward, penalty, gapopen, gapextend, snp_dict, match_alignment_dict, protein_dict, logger):

    ## Check allele sequence direction
    if query_direction == 'reverse':
        matching_allele_seq = str(Seq.Seq(matching_allele_seq).reverse_complement())
    else:
        matching_allele_seq = str(matching_allele_seq)

    ## Get the SNP information
    snp_information = get_snp(sseq, matching_allele_seq) 
    if len(snp_information) > 0 :
        if not core_name in snp_dict :
            snp_dict[core_name] = {}
        if not sample_name in snp_dict[core_name] :
            snp_dict[core_name][sample_name] = {}
        snp_dict[core_name][sample_name][qseqid]= snp_information

    ## Get new sequence-allele sequence dna alignment 
    if not core_name in match_alignment_dict :
        match_alignment_dict[core_name] = {}
        if not sample_name in match_alignment_dict[core_name] :
            match_alignment_dict[core_name][sample_name] = get_alignment (new_sseq, matching_allele_seq, reward, penalty, gapopen, gapextend) 

    ## Get new sequence-allele sequence protein alignment 
    if not core_name in protein_dict :
        protein_dict[core_name] = {}
    if not sample_name in protein_dict[core_name] :
        protein_dict[core_name][sample_name] = []
    protein_dict[core_name][sample_name] = get_alignment (new_sseq, matching_allele_seq, reward, penalty, gapopen, gapextend, "protein")
                        
    return True


# · * · * · * · * · * · * · * · * · * · * ·  #  
# Create allele tag summary for each sample  #
# · * · * · * · * · * · * · * · * · * · * ·  #

def create_summary (samples_matrix_dict, logger) :

    summary_dict = {}
    summary_result_list = []
    summary_heading_list = ['Exact match', 'INF', 'ASM', 'ALM', 'LNF', 'TPR', 'NIPH', 'NIPHEM', 'PLOT', 'ERROR']
    summary_result_list.append('File\t' + '\t'.join(summary_heading_list))
    for key in sorted (samples_matrix_dict) :

        summary_dict[key] = {'Exact match':0, 'INF':0, 'ASM':0, 'ALM':0, 'LNF':0, 'TPR':0,'NIPH':0, 'NIPHEM':0, 'PLOT':0, 'ERROR':0}
        for values in samples_matrix_dict[key] :
            if 'INF_' in values :
                summary_dict[key]['INF'] += 1
            elif 'ASM_' in values :
                summary_dict[key]['ASM'] += 1
            elif 'ALM_' in values :
                summary_dict[key]['ALM'] += 1
            elif 'LNF' in values :
                summary_dict[key]['LNF'] += 1
            elif 'TPR' in values :
                summary_dict[key]['TPR'] += 1
            elif 'NIPH' == values : 
                summary_dict[key]['NIPH'] += 1
            elif 'NIPHEM' == values :
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


# · * · * · * · * · * · * · * · * · * · * · * · * · * · * ·  #  
# Get gene and product annotation for core gene using Prokka #
# · * · * · * · * · * · * · * · * · * · * · * · * · * · * ·  #

### (tsv para algunos locus? Utils para analyze schema?)
def get_gene_annotation (annotation_file, annotation_dir, genus, species, usegenus, logger) :
    
    name_file = os.path.basename(annotation_file).split('.')
    annotation_dir = os.path.join (annotation_dir, 'annotation', name_file[0])
    
    if usegenus == 'true':
        annotation_result = subprocess.run (['prokka', annotation_file, '--outdir', annotation_dir,
                                            '--genus', genus, '--species', species, '--usegenus', 
                                            '--gcode', '11', '--prefix', name_file[0], '--quiet'])

    elif usegenus == 'false':
        annotation_result = subprocess.run (['prokka', annotation_file, '--outdir', annotation_dir,
                                            '--genus', genus, '--species', species, 
                                            '--gcode', '11', '--prefix', name_file[0], '--quiet'])
    
    annot_tsv = []
    tsv_path = os.path.join (annotation_dir, name_file[0] + '.tsv')

    try:
        with open(tsv_path) as tsvfile:
            tsvreader = csv.reader(tsvfile, delimiter="\t")
            for line in tsvreader:
                annot_tsv.append(line)

        if len(annot_tsv) > 1:
            try:
                if '_' in annot_tsv[1][2]:
                    gene_annot = annot_tsv[1][2].split('_')[0]
                else:
                    gene_annot = annot_tsv[1][2]
            except:
                gene_annot = 'Not found by Prokka'
            
            try: 
                product_annot = annot_tsv[1][4]
            except:
                product_annot = 'Not found by Prokka'
        else:
            gene_annot = 'Not found by Prokka'
            product_annot = 'Not found by Prokka'
    except:
        gene_annot = 'Not found by Prokka'
        product_annot = 'Not found by Prokka'

    return gene_annot, product_annot

"""
def get_gene_annotation (annotation_file, annotation_dir, logger) :
    name_file = os.path.basename(annotation_file).split('.')
    annotation_dir = os.path.join (annotation_dir, 'annotation', name_file[0])
    
    annotation_result = subprocess.run (['prokka', annotation_file, '--outdir', annotation_dir,
                                        '--prefix', name_file[0], '--quiet'])
    
    annot_tsv = []
    tsv_path = os.path.join (annotation_dir, name_file[0] + '.tsv')

    try:
        with open(tsv_path) as tsvfile:
            tsvreader = csv.reader(tsvfile, delimiter="\t")
            for line in tsvreader:
                annot_tsv.append(line)

        if len(annot_tsv) > 1:
            try:
                if '_' in annot_tsv[1][2]:
                    gene_annot = annot_tsv[1][2].split('_')[0]
                else:
                    gene_annot = annot_tsv[1][2]
            except:
                gene_annot = 'Not found by Prokka'
            
            try: 
                product_annot = annot_tsv[1][4]
            except:
                product_annot = 'Not found by Prokka'
        else:
            gene_annot = 'Not found by Prokka'
            product_annot = 'Not found by Prokka'
    except:
        gene_annot = 'Not found by Prokka'
        product_annot = 'Not found by Prokka'

    return gene_annot, product_annot
"""

"""
def get_gene_annotation (annotation_file, annotation_dir, logger) :
    name_file = os.path.basename(annotation_file).split('.')
    annotation_dir = os.path.join (annotation_dir, 'annotation', name_file[0])
    
    annotation_result = subprocess.run (['prokka', annotation_file, '--outdir', annotation_dir,
                                        '--prefix', name_file[0], '--quiet'])
    
    annot_tsv = []
    tsv_path = os.path.join (annotation_dir, name_file[0] + '.tsv')
    with open(tsv_path) as tsvfile:
        tsvreader = csv.reader(tsvfile, delimiter="\t")
        for line in tsvreader:
            annot_tsv.append(line)

    if len(annot_tsv) > 1:
        try:
            if '_' in annot_tsv[1][2]:
                gene_annot = annot_tsv[1][2].split('_')[0]
            else:
                gene_annot = annot_tsv[1][2]
        except:
            gene_annot = 'Not found by Prokka'
        
        try: 
            product_annot = annot_tsv[1][4]
        except:
            product_annot = 'Not found by Prokka'
    else:
        gene_annot = 'Not found by Prokka'
        product_annot = 'Not found by Prokka'

    return gene_annot, product_annot
"""


def analize_annotation_files (in_file, logger) : ## N
    examiner = GFF.GFFExaminer()
    file_fh = open(in_file)
    datos = examiner.available_limits(in_file)
    file_fh.close()
    return True


def get_inferred_allele_number(core_dict, logger): ## N
    #This function will look for the highest locus number and it will return a safe high value
    # that will be added to the schema database
    logger.debug('running get_inferred_allele_number function')
    int_keys = []
    for key in core_dict.keys():
        int_keys.append(key)
    max_value = max(int_keys)
    digit_length = len(str(max_value))
    return  True   #str 1 ( #'1'+ '0'*digit_length + 2)


# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * #  
# Get ST profile for each samples based on allele calling results #
# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * #

## (Revisar)
def get_ST_profile(profile_csv_path, exact_dict, core_gene_list_files):
    
    csv_read = []
    ST_profiles_dict = {}
    samples_profiles_dict = {}

    with open(profile_csv_path) as csvfile:
        csvreader = csv.reader(csvfile, delimiter="\t")
        for line in csvreader:
            csv_read.append(line)

    profile_header = csv_read[0][1:len(core_gene_list_files) + 1]

    for ST_index in range(1, len(csv_read)):
       # if ST_index == 0:
        #    profile_header = csv_read[ST_index][1:len(core_gene_list_files) + 1]
       # else:
        ST_profiles_dict[csv_read[ST_index][0]] = {}  
        for core_index in range(len(profile_header)):
            ST_profiles_dict[csv_read[ST_index][0]][profile_header[core_index]] = csv_read[ST_index][core_index]

    for sample_name in exact_dict:
        for ST in ST_profiles_dict:
            core_counter = 0
            for core_name in profile_header: # for core_name in ST_profiles_dict[ST]:
                if core_name in exact_dict[sample_name]:

                    allele_in_sample = exact_dict[sample_name][core_name][2]
                    allele_in_ST = ST_profiles_dict[ST][core_name]

                    if not '_' in allele_in_ST:
                        if '_' in allele_in_sample:
                            allele_in_sample = allele_in_sample.split('_')[1]

                    #if exact_dict[sample_name][core_name][2] == ST_profiles_dict[ST][core_name]: ## si el locus en cuestión se encuentra entre los exact matches y coincide con el alelo del profile en cuestión, se continua comparando los siguientes locus
                    if allele_in_sample == allele_in_ST:
                        core_counter += 1
                       # next
                    else: 
                        break
                else:
                    #if ST_profiles_dict[ST][core_name] == 'N':
                    if allele_in_ST == 'N':
                        core_counter += 1
                        #next
                    else:
                        break
            
            if core_counter == len(profile_header):
                samples_profiles_dict[sample_name] = ST
                break

        if sample_name not in samples_profiles_dict:
            samples_profiles_dict[sample_name] = '-'    

    return samples_profiles_dict


# · * · * · * · * · * · * · * · * · * · #  
# Create allele calling results reports #
# · * · * · * · * · * · * · * · * · * · #

def save_results (outputdir, full_gene_list, samples_matrix_dict, exact_dict, paralog_dict, inf_dict, plot_dict, matching_genes_dict, list_asm, list_alm, lnf_tpr_dict, snp_dict, match_alignment_dict, protein_dict, prodigal_report, shorter_seq_coverage, longer_seq_coverage, equal_seq_coverage, shorter_blast_seq_coverage, longer_blast_seq_coverage, equal_blast_seq_coverage, samples_profiles_dict, logger):
    header_matching_alleles_contig = ['Sample Name', 'Contig', 'Core Gene', 'Start', 'Stop', 'Direction', 'Codification']
    header_exact = ['Core Gene', 'Sample Name', 'Gene Annotation', 'Product Annotation', 'Allele', 'Allele Quality', 'Contig', 'Query length', 'Contig start', 'Contig end', 'Sequence', 'Predicted Sequence']
    header_paralogs = ['Core Gene','Sample Name', 'Gene Annotation', 'Product Annotation', 'Paralog Tag', 'ID %', 'Allele', 'Allele Quality', 'Contig', 'Bit Score', 'Contig start', 'Contig end', 'Sequence', 'Predicted Sequence']
    header_inferred = ['Core Gene','Sample Name', 'INF tag', 'Gene Annotation', 'Product Annotation', 'Allele', 'Allele Quality', 'Contig', 'Bitscore', 'Query length', 'Contig length', 'New sequence length' , 'Mismatch' , 'gaps', 'Contig start', 'Contig end',  'New sequence', 'Predicted Sequence']
    header_asm = ['Core Gene', 'Sample Name', 'ASM tag', 'Gene Annotation', 'Product Annotation', 'Allele', 'Allele Quality', 'Contig', 'Bitscore', 'Query length', 'Contig length', 'New sequence length' , 'Mismatch' , 'gaps', 'Contig start', 'Contig end',  'New sequence', 'Additional info', 'Predicted Sequence']
    header_alm = ['Core Gene', 'Sample Name', 'ALM tag', 'Gene Annotation', 'Product Annotation', 'Allele', 'Allele Quality', 'Contig', 'Bitscore', 'Query length', 'Contig length', 'New sequence length' , 'Mismatch' , 'gaps', 'Contig start', 'Contig end',  'New sequence', 'Additional info', 'Predicted Sequence']    
    header_plot = ['Core Gene', 'Sample Name', 'Gene Annotation', 'Product Annotation', 'Allele', 'Allele Quality', 'Contig','Bit Score', 'Contig start', 'Contig end', 'Sequence', 'Predicted Sequence']
    header_lnf_tpr = ['Core Gene', 'Sample Name', 'Gene Annotation', 'Product Annotation', 'Tag', 'Allele', 'ID %', 'Blast sequence coverage', 'New sequence coverage', 'Allele length', 'Blast sequence length', 'New sequence length', 'Additional info']        
    header_snp = ['Core Gene', 'Sample Name', 'Allele number', 'Position', 'Mutation Schema/Sample', 'Codon Schema/Sample','Protein in Schema/Sample', 'Missense/Synonymous','Annotation Sample / Schema']
    header_protein = ['Core Gene','Sample Name', 'Protein in ' , 'Protein sequence']
    header_match_alignment = ['Core Gene','Sample Name','Alignment', 'Sequence']
    header_stprofile = ['Sample Name', 'ST']


    # Añadido header_prodigal_report para report prodigal
    header_prodigal_report = ['Core gene', 'Sample Name', 'Allele', 'Sequence type', 'BLAST start', 'BLAST end', 'Prodigal start', 'Prodigal end', 'BLAST sequence', 'Prodigal sequence']
    # Añadido header_newsseq_coverage_report para determinar coverage threshold a imponer
    header_newsseq_coverage_report = ['Core gene', 'Sample Name', 'Query length', 'New sequence length', 'Locus mean', 'Coverage (new sequence/allele)', 'Coverage (new sequence/locus mean)']
    # Añadido header_blast_coverage_report para determinar coverage threshold a imponer
    header_blast_coverage_report = ['Core gene', 'Sample Name', 'Query length', 'Blast sequence length', 'Locus mean', 'Coverage (blast sequence/allele)', 'Coverage (blast sequence/locus mean)']

    ## Saving the result information to file
    print ('Saving results to files \n')
    result_file = os.path.join ( outputdir, 'result.tsv')
    logger.info('Saving result information to file..')
    with open (result_file, 'w') as out_fh:
        out_fh.write ('Sample Name\t'+'\t'.join( full_gene_list) + '\n')
        for key in sorted (samples_matrix_dict):
            out_fh.write (key + '\t' + '\t'.join(samples_matrix_dict[key])+ '\n')

    ## Saving exact matches to file
    logger.info('Saving exact matches information to file..')
    exact_file =  os.path.join(outputdir, 'exact.tsv')
    with open (exact_file , 'w') as exact_fh :
        exact_fh.write('\t'.join(header_exact)+ '\n')
        for core in sorted(exact_dict): 
            for sample in sorted(exact_dict[core]):
                exact_fh.write(core + '\t' + sample + '\t' + '\t'.join(exact_dict[core][sample]) + '\n')

    ## Saving paralog alleles to file
    logger.info('Saving paralog information to file..')
    paralog_file =  os.path.join(outputdir, 'paralog.tsv')
    with open (paralog_file , 'w') as paralog_fh :
        paralog_fh.write('\t'.join(header_paralogs) + '\n')
        for sample in sorted (paralog_dict) :
            for core in sorted (paralog_dict[sample]):
                for paralog in paralog_dict[sample][core] :
                    paralog_fh.write(core + '\t' + sample + '\t' + '\t'.join (paralog) + '\n')

    ## Saving inferred alleles to file
    logger.info('Saving inferred alleles information to file..')
    inferred_file =  os.path.join(outputdir, 'inferred_alleles.tsv')
    with open (inferred_file , 'w') as infer_fh :
        infer_fh.write('\t'.join(header_inferred) + '\n')
        for core in sorted (inf_dict) :
            for sample in sorted (inf_dict[core]) :
                for inferred in inf_dict[core][sample]: 
                    #   seq_in_inferred_allele = '\t'.join (inf_dict[sample])
                    infer_fh.write(core + '\t' + sample + '\t' + inferred + '\t' + '\t'.join(inf_dict[core][sample][inferred]) + '\n')
    
    ## Saving PLOTs to file
    logger.info('Saving PLOT information to file..')
    plot_file =  os.path.join(outputdir, 'plot.tsv')
    with open (plot_file , 'w') as plot_fh :
        plot_fh.write('\t'.join(header_plot) + '\n')
        for core in sorted (plot_dict) :
            for sample in sorted (plot_dict[core]):
                plot_fh.write(core + '\t' + sample + '\t' + '\t'.join(plot_dict[core][sample]) + '\n')

    ## Saving matching contigs to file
    logger.info('Saving matching information to file..')
    matching_file =  os.path.join(outputdir, 'matching_contigs.tsv')
    with open (matching_file , 'w') as matching_fh :
        matching_fh.write('\t'.join(header_matching_alleles_contig ) + '\n')
        for samples in sorted ( matching_genes_dict) :
            for contigs in matching_genes_dict[samples] :
                for contig in matching_genes_dict[samples] [contigs]:
                        matching_alleles = '\t'.join (contig)
                        matching_fh.write(samples + '\t' + contigs +'\t' + matching_alleles + '\n')
 
    ## Saving ASMs to file
    logger.info('Saving asm information to file..')
    asm_file =  os.path.join(outputdir, 'asm.tsv')
    with open (asm_file , 'w') as asm_fh :
        asm_fh.write('\t'.join(header_asm)+ '\n')
        for core in list_asm :
            for sample in list_asm[core] :
                for asm in list_asm[core][sample]:
                    asm_fh.write(core + '\t' + sample + '\t' + asm + '\t' + '\t'.join(list_asm[core][sample][asm]) + '\n')

    ## Saving ALMs to file
    logger.info('Saving alm information to file..')
    alm_file =  os.path.join(outputdir, 'alm.tsv')
    with open (alm_file , 'w') as alm_fh :
        alm_fh.write('\t'.join(header_alm)+ '\n')
        for core in list_alm :
            for sample in list_alm[core] :
                for alm in list_alm[core][sample]:
                    alm_fh.write(core + '\t' + sample + '\t' + alm + '\t' + '\t'.join(list_alm[core][sample][alm]) + '\n')

    ## Saving LNFs to file
    logger.info('Saving lnf information to file..')
    lnf_file =  os.path.join(outputdir, 'lnf.tsv')
    with open (lnf_file , 'w') as lnf_fh :
        lnf_fh.write('\t'.join(header_lnf_tpr)+ '\n')
        for core in lnf_tpr_dict :
            for sample in lnf_tpr_dict[core] :
                for lnf in lnf_tpr_dict[core][sample] :
                    lnf_fh.write(core + '\t' + sample + '\t' + '\t'.join(lnf) + '\n')

    ## Saving SNPs information to file
    logger.info('Saving SNPs information to file..')
    snp_file =  os.path.join(outputdir, 'snp.tsv')
    with open (snp_file , 'w') as snp_fh :
        snp_fh.write('\t'.join(header_snp) + '\n')
        for core in sorted (snp_dict) :
            for sample in sorted (snp_dict[core]):
                for allele_id_snp in snp_dict[core][sample] :
                    for snp in snp_dict[core][sample][allele_id_snp] :
                        snp_fh.write(core + '\t' + sample + '\t' + allele_id_snp + '\t' + '\t'.join (snp) + '\n')

    ## Saving DNA sequences alignments to file
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

    ## Saving protein sequences alignments to file
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

    ## Saving summary information to file 
    logger.info('Saving summary information to file..')
    summary_result = create_summary (samples_matrix_dict, logger)
    summary_file = os.path.join( outputdir, 'summary_result.tsv')
    with open (summary_file , 'w') as summ_fh:
        for line in summary_result :
            summ_fh.write(line + '\n')

    ## Modify the result file to remove the PLOT_ string for creating the file to use in the tree diagram
    logger.info('Saving result information for tree diagram')
    tree_diagram_file = os.path.join ( outputdir, 'result_for_tree_diagram.tsv')
    with open (result_file, 'r') as result_fh:
        with open(tree_diagram_file, 'w') as td_fh:
            for line in result_fh:
                tree_line = line.replace('PLOT_','')
                td_fh.write(tree_line)

    if samples_profiles_dict != '':
        ## Saving ST profile to file
        logger.info('Saving ST profile information to file..')
        stprofile_file =  os.path.join(outputdir, 'stprofile.tsv')
        with open (stprofile_file , 'w') as st_fh :
            st_fh.write('\t'.join(header_stprofile)+ '\n')
            for sample in sorted(samples_profiles_dict): 
                st_fh.write(sample + '\t' + samples_profiles_dict[sample] + '\n')

    ###########################################################################################
    # Guardando report de prodigal. Temporal
    prodigal_report_file = os.path.join (outputdir, 'prodigal_report.tsv')
    # saving prodigal predictions to file
    with open (prodigal_report_file, 'w') as out_fh:
        out_fh.write ('\t'.join(header_prodigal_report)+ '\n')
        for prodigal_result in prodigal_report:
            out_fh.write ('\t'.join(prodigal_result)+ '\n')

    # Guardando coverage de new_sseq para estimar el threshold a establecer. Temporal
    newsseq_coverage_file = os.path.join (outputdir, 'newsseq_coverage_report.tsv')
    # saving the coverage information to file
    with open (newsseq_coverage_file, 'w') as out_fh:
        out_fh.write ('\t' + '\t'.join(header_newsseq_coverage_report)+ '\n')
        for coverage in shorter_seq_coverage:
            out_fh.write ('Shorter new sequence' + '\t' + '\t'.join(coverage)+ '\n')
        for coverage in longer_seq_coverage:
            out_fh.write ('Longer new sequence' + '\t' + '\t'.join(coverage)+ '\n')
        for coverage in equal_seq_coverage:
            out_fh.write ('Same length new sequence' + '\t' + '\t'.join(coverage)+ '\n')

    # Guardando coverage de la sseq obtenida tras blast para estimar el threshold a establecer. Temporal
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

    return True


# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · *  #  
# Update core genes schema adding new inferred alleles found for each locus in allele calling analysis #
# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · *  #

## intentando tener en cuenta los esquemas cuyos alelos contengan ids alfanuméricos
def update_schema (updateschema, schemadir, storedir, core_gene_list_files, inferred_alleles_dict, alleles_in_locus_dict, logger):
    
    ## Create a copy of core genes schema if updateschema = 'new' / 'New'
    if updateschema == 'new' or updateschema == 'New':
        no_updated_schemadir = schemadir
        schemadir_name = os.path.basename(no_updated_schemadir)
        schemadir = os.path.join(storedir, schemadir_name + '_updated') 
        shutil.copytree(no_updated_schemadir, schemadir)
        logger.info('Copying core genes fasta files to update schema')
        
    ## Get INF alleles for each core gene and update each locus fasta file
    for core_file in core_gene_list_files:
        core_name = os.path.basename(core_file).split('.')[0]
        if core_name in inferred_alleles_dict:
            logger.info('Updating core gene file %s adding new INF alleles', core_file)
            
            inf_list = inferred_alleles_dict[core_name]

            try:
                alleles_ids = [int(allele) for allele in alleles_in_locus_dict[core_name]]
                allele_number = max(alleles_ids)
                
                locus_schema_file = os.path.join(schemadir, core_name + '.fasta') 
                with open (locus_schema_file, 'a') as core_fh:
                    for inf in inf_list:
                        allele_number += 1
                        complete_inf_id = core_name + '_' + str(allele_number)
                        core_fh.write('\n' + '>' + str(allele_number) + ' # ' + 'INF by Taranis' + '\n' + inf + '\n')
            except:
                alleles_ids = [int(allele.split('_')[-1]) for allele in alleles_in_locus_dict[core_name]]
                allele_number = max(alleles_ids)

                locus_schema_file = os.path.join(schemadir, core_name + '.fasta') 
                with open (locus_schema_file, 'a') as core_fh:
                    for inf in inf_list:
                        allele_number += 1
                        complete_inf_id = core_name + '_' + str(allele_number)
                        core_fh.write('\n' + '>' + complete_inf_id + ' # ' + 'INF by Taranis' + '\n' + inf + '\n')

    return True

"""
def update_schema (updateschema, schemadir, storedir, core_gene_list_files, inferred_alleles_dict, alleles_in_locus_dict, logger):
    
    ## Create a copy of core genes schema if updateschema = 'new' / 'New'
    if updateschema == 'new' or updateschema == 'New':
        no_updated_schemadir = schemadir
        schemadir_name = os.path.basename(no_updated_schemadir)
        schemadir = os.path.join(storedir, schemadir_name + '_updated') 
        shutil.copytree(no_updated_schemadir, schemadir)
        logger.info('Copying core genes fasta files to update schema')
        
    ## Get INF alleles for each core gene and update each locus fasta file
    for core_file in core_gene_list_files:
        core_name = os.path.basename(core_file).split('.')[0]
        if core_name in inferred_alleles_dict:
            logger.info('Updating core gene file %s adding new INF alleles', core_file)
            
            inf_list = inferred_alleles_dict[core_name]
        
            alleles_ids = [int(allele) for allele in alleles_in_locus_dict[core_name]]
            allele_number = max(alleles_ids)

            locus_schema_file = os.path.join(schemadir, core_name + '.fasta') 
            with open (locus_schema_file, 'a') as core_fh:
                for inf in inf_list:
                    allele_number += 1
                    core_fh.write('\n' + '>' + str(allele_number) + ' # ' + 'INF by Taranis' + '\n' + inf + '\n')

    return True
"""


# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · *  #  
# Allele calling analysis to find each core gene in schema and its variants in samples #
# · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · * · *  #

def allele_call_nucleotides (core_gene_list_files, sample_list_files, alleles_in_locus_dict, contigs_in_sample_dict, query_directory, reference_alleles_directory, blast_db_directory, prodigal_directory, blast_results_seq_directory, blast_results_db_directory, inputdir, outputdir, cpus, percentlength, coverage, evalue, perc_identity_ref, perc_identity_loc, reward, penalty, gapopen, gapextend, max_target_seqs, max_hsps, num_threads, flankingnts, schema_variability, schema_statistics, schema_quality, annotation_core_dict, profile_csv_path, logger ):
    
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
    lnf_tpr_dict = {} # c/m: to keep locus not found for each sample
    plot_dict = {} # c/m: to keep plots for each sample
    snp_dict = {} # c/m: to keep snp information for each sample
    protein_dict = {}
    match_alignment_dict = {}

    blast_parameters = '"6 , qseqid , sseqid , pident ,  qlen , length , mismatch , gapopen , evalue , bitscore , sstart , send , qstart , qend , sseq , qseq"'

    print('Allele calling starts')
    pbar = ProgressBar ()


    ## # # # # # # # # # # # # # # # # # # # # # # # # ##
    ## Processing the search for each schema core gene ##
    ## # # # # # # # # # # # # # # # # # # # # # # # # ##

    for core_file in pbar(core_gene_list_files) :
        core_name = os.path.basename(core_file).split('.')[0]
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

            sample_name = os.path.basename(sample_file).split('.')[0]

            # Initialize the sample list to add the number of alleles and the start, stop positions
            if not sample_name in samples_matrix_dict:
                samples_matrix_dict[sample_name] = []
                matching_genes_dict[sample_name] = {}

            # Path to this sample BLAST database created when processing samples 
            blast_db_name = os.path.join(blast_db_directory, sample_name, sample_name)

            # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
            # Sample contigs VS reference allele(s) BLAST for locus detection in sample #
            # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #

            cline = NcbiblastnCommandline(db=blast_db_name, evalue=evalue, perc_identity=perc_identity_ref, reward=reward, penalty=penalty, gapopen=gapopen, gapextend=gapextend, outfmt=blast_parameters, max_target_seqs=max_target_seqs, max_hsps=max_hsps, num_threads=num_threads, query=core_reference_allele_path)
            out, err = cline()
            out_lines = out.splitlines()
            
            bigger_bitscore = 0 

            # ······························································ #
            # LNF if there are no BLAST results for this gene in this sample #
            # ······························································ #
            if len (out_lines) == 0:

                # Trying to get the allele number to avoid that a bad quality assembly impact on the tree diagram
                cline = NcbiblastnCommandline(db=blast_db_name, evalue=evalue, perc_identity = 70, reward=reward, penalty=penalty, gapopen=gapopen, gapextend=gapextend, outfmt=blast_parameters, max_target_seqs=1, max_hsps=1, num_threads=1, query=core_reference_allele_path)
                out, err = cline()
                out_lines = out.splitlines()

                if len (out_lines) > 0 :
                
                    for line in out_lines :
                        values = line.split('\t')
                        if  float(values[8]) > bigger_bitscore:
                            qseqid , sseqid , pident ,  qlen , s_length , mismatch , r_gapopen , r_evalue , bitscore , sstart , send , qstart , qend ,sseq , qseq = values
                            bigger_bitscore = float(bitscore)

                    # Keep LNF info
                    lnf_tpr_tag(core_name, sample_name, alleles_in_locus_dict, samples_matrix_dict, lnf_tpr_dict, schema_statistics, locus_alleles_path, qseqid, pident, '-', '-', perc_identity_ref, '-', annotation_core_dict, logger)

                else:
                    # Keep LNF info
                    lnf_tpr_tag(core_name, sample_name, '-', samples_matrix_dict, lnf_tpr_dict, schema_statistics, locus_alleles_path, '-', '-', '-', '-', '-', '-', annotation_core_dict, logger)

                continue

            ## Continue classification process if the core gene has been detected in sample after BLAST search
            if len (out_lines) > 0:

                # Parse contigs for this sample 
                #contig_file = os.path.join(inputdir, sample_name + ".fasta") ## parse
                #records = list(SeqIO.parse(contig_file, "fasta")) ## parse

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

                        # Get flanked BLAST sequences from contig for correct allele tagging

                        accession_sequence = contigs_in_sample_dict[sample_name][sseqid]
                        #for record in records: ## parse
                            #if record.id == sseqid : ## parse
                                #break ## parse
                        #accession_sequence = str(record.seq) ## parse

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

                        seq_id = str(seq_number) + '_' + sseqid
                        outblast_fh.write('>' + seq_id + ' # ' + ' # '.join(values[0:13]) + '\n' + flanked_sseq + '\n' + '\n' )                                                                                                                                        
                    
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
                
                ## Check if there is any BLAST result with ID = 100 ##
                for line in out_lines:

                    values = line.split('\t')
                    pident = values[2]

                    if float(pident) == 100:

                        qseqid, sseqid, pident, qlen, s_length, mismatch, r_gapopen, r_evalue, bitscore, sstart, send, qstart, qend, sseq, qseq = values

                        # Parse core gene fasta file to get matching allele sequence and length
                        #alleles_in_locus = list (SeqIO.parse(locus_alleles_path, "fasta")) ## parse
                        #for allele in alleles_in_locus : ## parse
                            #if allele.id == qseqid : ## parse
                                #break ## comentar parse
                        #matching_allele_seq = str(allele.seq) ## parse
                        #matching_allele_length = len(matching_allele_seq) ## parse

                        matching_allele_seq = alleles_in_locus_dict[core_name][qseqid]
                        matching_allele_length = len(matching_allele_seq) 

                        # Keep BLAST results with ID = 100 and same length as matching allele
                        if int(s_length) == matching_allele_length:
                            #get_blast_results (values, records, allele_found, logger)
                            get_blast_results (sample_name, values, contigs_in_sample_dict, allele_found, logger)

                # ·································································································································· #
                # NIPHEM (paralog) if there are multiple BLAST results with ID = 100 and same length as matching allele for this gene in this sample #
                # ·································································································································· #
                if len(allele_found) > 1:

                    # Keep NIPHEM info
                    paralog_exact_tag(sample_name, core_name, 'NIPHEM', schema_quality, matching_genes_dict, samples_matrix_dict, allele_found, paralog_dict, prodigal_report, prodigal_directory, blast_parameters, annotation_core_dict, logger)

                    continue
                
                ## Check for possible paralogs with ID < 100 if there is only one BLAST result with ID = 100 and same length as matching allele
                elif len(allele_found) == 1 :

                    for line in out_lines :
                        values = line.split('\t')

                        sseq_no_gaps = values[13].replace('-', '')
                        s_length_no_gaps = len(sseq_no_gaps)

                        # Keep BLAST result if its coverage is within min and max thresholds
                        if min_length_threshold <= s_length_no_gaps <= max_length_threshold:                            
                            #get_blast_results (values, records, allele_found, logger) 
                            get_blast_results (sample_name, values, contigs_in_sample_dict, allele_found, logger)

                    # ································································ #
                    # EXACT MATCH if there is any paralog for this gene in this sample #
                    # ································································ #
                    if len(allele_found) == 1 :

                        paralog_exact_tag(sample_name, core_name, 'EXACT', schema_quality, matching_genes_dict, samples_matrix_dict, allele_found, exact_dict, prodigal_report, prodigal_directory, blast_parameters, annotation_core_dict, logger)

                        continue
                    
                    # ··········································································· #
                    # NIPH if there there are paralogs with ID < 100 for this gene in this sample #
                    # ··········································································· #
                    else:

                        paralog_exact_tag(sample_name, core_name, 'NIPH', schema_quality, matching_genes_dict, samples_matrix_dict, allele_found, paralog_dict, prodigal_report, prodigal_directory, blast_parameters, annotation_core_dict, logger)

                        continue

                ## Look for the best BLAST result if there are no results with ID = 100 ##
                elif len(allele_found) == 0:

                    bigger_bitscore_seq_values = []

                    for line in out_lines :
                        values = line.split('\t')

                        if  float(values[8]) > bigger_bitscore:
                            s_length_no_gaps = len(values[13].replace('-', ''))

                            # Keep BLAST result if its coverage is within min and max thresholds and its bitscore is bigger than the one previously kept
                            if min_coverage_threshold <= s_length_no_gaps <= max_coverage_threshold:
                                bigger_bitscore_seq_values = values
                                bigger_bitscore = float(bigger_bitscore_seq_values[8])

                    ## Check if best BLAST result out of coverage thresholds is a possible PLOT or LNF due to low coverage ##
                    #if len(allele_found) == 0:
                    if len(bigger_bitscore_seq_values) == 0:
                
                        # Look for best bitscore BLAST result out of coverage thresholds to check possible PLOT or reporting LNF due to low coverage
                        for line in out_lines :
                            values = line.split('\t')

                            if  float(values[8]) > bigger_bitscore:
                                qseqid, sseqid, pident,  qlen, s_length, mismatch, r_gapopen, r_evalue, bitscore, sstart, send, qstart, qend, sseq, qseq = values
                                bigger_bitscore = float(bitscore)

                            s_length_no_gaps = len(sseq.replace('-', ''))

                            # Get contig sequence and length for best bitscore BLAST result ID 
                            
                            #for record in records: ## parse
                                #if record.id == sseqid : ## parse
                                    #break ## parse
                            #accession_sequence = record.seq ## parse
                            #length_sseqid = len(accession_sequence) ## parse

                            accession_sequence = contigs_in_sample_dict[sample_name][sseqid]
                            length_sseqid = len(accession_sequence)

                            # Check if best BLAST result out of coverage thresholds is a possible PLOT. If so, keep result info for later PLOT classification
                            if int(sstart) == length_sseqid or int(send) == length_sseqid or int(sstart) == 1 or int(send) == 1:
                                bigger_bitscore_seq_values = values

                            # ·············································································································································· #
                            # LNF if there are no BLAST results within coverage thresholds for this gene in this sample and best out threshold result is not a possible PLOT #
                            # ·············································································································································· #
                            else:

                                # Keep LNF info
                                lnf_tpr_tag(core_name, sample_name, alleles_in_locus_dict, samples_matrix_dict, lnf_tpr_dict, schema_statistics, locus_alleles_path, qseqid, pident, s_length_no_gaps, '-', '-', coverage, annotation_core_dict, logger)


                    ## Keep result with bigger bitscore in allele_found dict and look for possible paralogs ##
                    if len(bigger_bitscore_seq_values) > 0:

                        qseqid, sseqid, pident, qlen, s_length, mismatch, r_gapopen, r_evalue, bitscore, sstart, send, qstart, qend, sseq, qseq = bigger_bitscore_seq_values

                        #get_blast_results (bigger_bitscore_seq_values, records, allele_found, logger)
                        get_blast_results (sample_name, values, contigs_in_sample_dict, allele_found, logger)

                        # Possible paralogs search
                        for line in out_lines :
                            values = line.split('\t')

                            qseqid, sseqid, pident, qlen, s_length, mismatch, r_gapopen, r_evalue, bitscore, sstart, send, qstart, qend, sseq, qseq = values
                            sseq_no_gaps = sseq.replace('-', '')
                            s_length_no_gaps = len(sseq_no_gaps) 

                            if min_length_threshold <= s_length_no_gaps <= max_length_threshold:                            

                                #get_blast_results (values, records, allele_found, logger)
                                get_blast_results (sample_name, values, contigs_in_sample_dict, allele_found, logger)

                        # ····························································· #
                        # NIPH if there there are paralogs for this gene in this sample #
                        # ····························································· #
                        if len(allele_found) > 1 :

                            paralog_exact_tag(sample_name, core_name, 'NIPH', schema_quality, matching_genes_dict, samples_matrix_dict, allele_found, paralog_dict, prodigal_report, prodigal_directory, blast_parameters, annotation_core_dict, logger)

                            continue 

                        ## Continue classification if there are no paralogs ##
                        elif len(allele_found) == 1 :

                            allele_id = str(list(allele_found.keys())[0])
                            qseqid, sseqid, pident, qlen, s_length, mismatch, r_gapopen, r_evalue, bitscore, sstart, send, qstart, qend, sseq, qseq = allele_found[allele_id]
                            
                            sseq_no_gaps = sseq.replace('-', '')
                            s_length_no_gaps = len(sseq_no_gaps)

                            # Get matching allele quality
                            allele_quality = schema_quality[core_name][qseqid]

                            # Get matching allele sequence and length
                            
                            #alleles_in_locus = list (SeqIO.parse(locus_alleles_path, "fasta")) ## parse
                            #for allele in alleles_in_locus : ## parse
                                #if allele.id == qseqid : ## parse
                                    #break ## parse
                            #matching_allele_seq = allele.seq ## parse
                            #matching_allele_length = len(matching_allele_seq) ## parse               

                            matching_allele_seq = alleles_in_locus_dict [core_name][qseqid]
                            matching_allele_length = len(matching_allele_seq)

                            # Get contig sequence and length for ID found in BLAST
                            
                            #for record in records: ## parse
                                #if record.id == sseqid : ## parse
                                    #break ## parse
                            #accession_sequence = record.seq ## parse
                            #length_sseqid = len(accession_sequence) ## parse

                            accession_sequence = contigs_in_sample_dict[sample_name][sseqid]
                            length_sseqid = len(accession_sequence)

                            # ········································································································· #
                            # PLOT if found sequence is shorter than matching allele and it is located on the edge of the sample contig #
                            # ········································································································· #
                            if int(sstart) == length_sseqid or int(send) == length_sseqid or int(sstart) == 1 or int(send) == 1:
                                if int(s_length) < matching_allele_length:

                                    ### sacar sec prodigal para PLOT?
                                    # Get prodigal predicted sequence if matching allele quality is "bad quality"
                                    if 'bad_quality' in allele_quality: 
                                        complete_predicted_seq, start_prodigal, end_prodigal = get_prodigal_sequence(sseq_no_gaps, sseqid, prodigal_directory, sample_name, blast_parameters, logger)

                                        # Keep info for prodigal report
                                        prodigal_report.append([core_name, sample_name, qseqid, 'PLOT', sstart, send, start_prodigal, end_prodigal, sseq_no_gaps, complete_predicted_seq])

                                    else:
                                        complete_predicted_seq = '-'
                                        start_prodigal = '-'
                                        end_prodigal = '-'

                                    # Keep PLOT info
                                    inf_asm_alm_tag(core_name, sample_name, 'PLOT', allele_found[allele_id], allele_quality, '-', matching_allele_length, '-', plot_dict, samples_matrix_dict, matching_genes_dict, prodigal_report, start_prodigal, end_prodigal, complete_predicted_seq, annotation_core_dict, logger)

                                    continue 

                            # * * * * * * * * * * * * * * * * * * * * #
                            # Search for complete final new sequence  # 
                            # * * * * * * * * * * * * * * * * * * * * #

                            ## Get Prodigal predicted sequence ##
                            complete_predicted_seq, start_prodigal, end_prodigal = get_prodigal_sequence(sseq_no_gaps, sseqid, prodigal_directory, sample_name, blast_parameters, logger)

                            ## Search for new codon stop using contig sequence info ##

                            # Check matching allele sequence direction
                            query_direction = check_sequence_order(matching_allele_seq, logger)

                            # Get extended BLAST sequence for stop codon search
                            if query_direction == 'reverse':
                                if int(send) > int (sstart):
                                    sample_gene_sequence = accession_sequence[ : int(send) ]
                                    sample_gene_sequence = str(Seq.Seq(sample_gene_sequence).reverse_complement())
                                else:
                                    sample_gene_sequence = accession_sequence[ int(send) -1 : ]

                            else:
                                if int(sstart) > int (send):
                                    sample_gene_sequence = accession_sequence[ :  int(sstart) ]
                                    sample_gene_sequence = str(Seq.Seq(sample_gene_sequence).reverse_complement())
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
                                keep_snp_alignment_info(sseq, new_sseq, matching_allele_seq, qseqid, query_direction, core_name, sample_name, reward, penalty, gapopen, gapextend, snp_dict, match_alignment_dict, protein_dict, logger)

                                # ····································································································· #
                                # INF if final new sequence length is within min and max length thresholds for this gene in this sample #
                                # ····································································································· #
                                if min_length_threshold <= new_sequence_length <= max_length_threshold:

                                    # Keep INF info
                                    inf_asm_alm_tag(core_name, sample_name, 'INF', allele_found[allele_id], allele_quality, new_sseq, matching_allele_length, inferred_alleles_dict, inf_dict, samples_matrix_dict, matching_genes_dict, prodigal_report, start_prodigal, end_prodigal, complete_predicted_seq, annotation_core_dict, logger) ### introducido start_prodigal, end_prodigal, complete_predicted_seq, prodigal_report como argumento a inf_asm_alm_tag para report prodigal, temporal

                                # ············································································································································ #
                                # ASM if final new sequence length is under min length threshold but its coverage is above min coverage threshold for this gene in this sample #
                                # ············································································································································ #
                                elif min_coverage_threshold <= new_sequence_length < min_length_threshold:

                                    # Keep ASM info
                                    inf_asm_alm_tag(core_name, sample_name, 'ASM', allele_found[allele_id], allele_quality, new_sseq, matching_allele_length, asm_dict, list_asm, samples_matrix_dict, matching_genes_dict, prodigal_report, start_prodigal, end_prodigal, complete_predicted_seq, annotation_core_dict, logger)

                                # ············································································································································ #
                                # ALM if final new sequence length is above max length threshold but its coverage is under max coverage threshold for this gene in this sample #
                                # ············································································································································ #
                                elif max_length_threshold < new_sequence_length <= max_coverage_threshold:

                                    # Keep ASM info
                                    inf_asm_alm_tag(core_name, sample_name, 'ALM', allele_found[allele_id], allele_quality, new_sseq, matching_allele_length, alm_dict, list_alm, samples_matrix_dict, matching_genes_dict, prodigal_report, start_prodigal, end_prodigal, complete_predicted_seq, annotation_core_dict, logger) ### introducido start_prodigal, end_prodigal, complete_predicted_seq, prodigal_report como argumento a inf_asm_alm_tag para report prodigal, temporal

                                # ························································································· #
                                # LNF if final new sequence coverage is not within thresholds for this gene in this sample  #
                                # ························································································· #
                                else: 

                                    # Keep LNF info
                                    lnf_tpr_tag(core_name, sample_name, alleles_in_locus_dict, samples_matrix_dict, lnf_tpr_dict, schema_statistics, locus_alleles_path, qseqid, pident, s_length_no_gaps, new_sequence_length, '-', coverage, annotation_core_dict, logger)

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
 

    ## Get ST profile for each sample
    
    if profile_csv_path != '':
        samples_profiles_dict = get_ST_profile(profile_csv_path, exact_dict, core_gene_list_files)
    else:
        samples_profiles_dict = ''


    ## Save results and create reports

    if not save_results (outputdir, full_gene_list, samples_matrix_dict, exact_dict, paralog_dict, inf_dict, plot_dict, matching_genes_dict, list_asm, list_alm, lnf_tpr_dict, snp_dict, match_alignment_dict, protein_dict, prodigal_report, shorter_seq_coverage, longer_seq_coverage, equal_seq_coverage, shorter_blast_seq_coverage, longer_blast_seq_coverage, equal_blast_seq_coverage, samples_profiles_dict, logger):
        print('There is an error while saving the allele calling results. Check the log file to get more information \n')
       # exit(0)

    return True, inferred_alleles_dict


# * * * * * * * * * * * * * * * * * * *  #
# Processing gene by gene allele calling #
# * * * * * * * * * * * * * * * * * * *  #

def processing_allele_calling (arguments) :
    '''
    Description:
        This is the main function for allele calling.
        With the support of additional functions it will create the output files
        with the summary report.
    Input:
        arguments   # Input arguments given on command line 
    Functions:
        ????
    Variables: 
        ????
    Return:
        ????
    '''

    start_time = datetime.now()
    print('Start the execution at :', start_time )
    
    # Open log file
    logger = open_log ('taranis_wgMLST.log')
    print('Checking the pre-requisites.')
    
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
            logger.info ('Temporary folder %s  has been created again', tmp_core_gene_dir)
        except:
            logger.info('Unable to create again the temporary directory %s', tmp_core_gene_dir)
            print('Cannot create temporary directory on ', tmp_core_gene_dir)
            exit(0)

    alleles_in_locus_dict, annotation_core_dict, schema_variability, schema_statistics, schema_quality = prepare_core_gene (valid_core_gene_files, tmp_core_gene_dir, arguments.refalleles, arguments.outputdir, arguments.genus, arguments.species, str(arguments.usegenus).lower(), logger)
    #alleles_in_locus_dict, annotation_core_dict, schema_variability, schema_statistics, schema_quality = prepare_core_gene (valid_core_gene_files, tmp_core_gene_dir, arguments.refalleles, arguments.outputdir, logger)
    if not alleles_in_locus_dict:
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
    
    contigs_in_sample_dict = prepare_samples(valid_sample_files, tmp_samples_dir, arguments.refgenome, logger)
    if not contigs_in_sample_dict :
        print('There is an error while processing the saving temporary files. Check the log file to get more information \n')
        logger.info('Deleting the temporary directory to clean up the temporary files created')
        shutil.rmtree(os.path.join(arguments.outputdir, 'tmp'))
        exit(0)

    ##################################
    ## Run allele callling analysis ##
    ##################################
    query_directory = arguments.coregenedir
    reference_alleles_directory = arguments.refalleles
    blast_db_directory = os.path.join(tmp_samples_dir,'blastdb')
    prodigal_directory = os.path.join(tmp_samples_dir,'prodigal') 
    blast_results_seq_directory = os.path.join(tmp_samples_dir,'blast_results', 'blast_results_seq')  ### path a directorio donde guardar secuencias encontradas tras blast con alelo de referencia
    blast_results_db_directory = os.path.join(tmp_samples_dir,'blast_results', 'blast_results_db') ### path a directorio donde guardar db de secuencias encontradas tras blast con alelo de referencia

    complete_allele_call, inferred_alleles_dict = allele_call_nucleotides(valid_core_gene_files, valid_sample_files, alleles_in_locus_dict, contigs_in_sample_dict, query_directory, reference_alleles_directory, blast_db_directory, prodigal_directory, blast_results_seq_directory, blast_results_db_directory, arguments.inputdir, arguments.outputdir,  int(arguments.cpus), arguments.percentlength, arguments.coverage, float(arguments.evalue), int(arguments.perc_identity_ref), int(arguments.perc_identity_loc), int(arguments.reward), int(arguments.penalty), int(arguments.gapopen), int(arguments.gapextend), int(arguments.max_target_seqs), int(arguments.max_hsps), int(arguments.num_threads), int(arguments.flankingnts), schema_variability, schema_statistics, schema_quality, annotation_core_dict, arguments.profile, logger) ### CAMBIANDO/MODIFICANDO: He añadido schema_statistics, path a prodigal, prodigal training y schema_quality        
    if not complete_allele_call:
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

    #########################################################
    ## Update core gene schema adding new inferred alleles ##
    #########################################################
    if inferred_alleles_dict:
        if str(arguments.updateschema).lower() == 'true' or str(arguments.updateschema).lower() == 'new':
            if not update_schema (str(arguments.updateschema).lower(), arguments.coregenedir, tmp_core_gene_dir, valid_core_gene_files, inferred_alleles_dict, alleles_in_locus_dict, logger):        
                print('There is an error adding new inferred alleles found to the core genes schema. Check the log file to get more information \n')
                exit(0)

    end_time = datetime.now()
    print('completed execution at :', end_time )

    return True



