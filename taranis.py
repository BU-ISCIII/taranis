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
import tempfile
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



import subprocess
#from subprocess import check_output
import shutil 

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

def check_program_is_exec_version (program, version, logger):
    # The function will check if the program is installed in your system and if the version
    # installed matched with the pre-requisites
    if shutil.which(program) is not None :
        # check version
        version_str= str(subprocess.check_output([program , '-version']))
        if version_str == "b''" :
            version_str = subprocess.getoutput( str (program + ' -version'))
        if not re.search(version, version_str):
            logger.info('%s program does not have the right version ', program)
            print ('Exiting script \n, Version of ' , program, 'does not fulfill the requirements')
            return False
        return True
    else:
        logger.info('Cannot find %s installed on your system', program)
        return False
        

def check_prerequisites (logger):
    pre_requisite_list = [['blastp', '2.6'], ['makeblastdb' , '2.6'], ['prokka', '1.12']]
    # check if blast is installed and has the minimum version 
    for program, version in pre_requisite_list :
        if not check_program_is_exec_version (program , version, logger):
            return False
    return True

def check_arg(args=None):
    
    parser = argparse.ArgumentParser(prog = 'get_subset', description="This program will make the Allele Calling using a predefined core Schema.")
    
    parser.add_argument('-coregenedir', help = 'Directory where the core gene files are located ')
    parser.add_argument('-inputdir', help ='Directory where are located the sample fasta files')
    parser.add_argument('-outputdir', help = 'Directory where the result files will be stored')
    parser.add_argument('-cpus', required= False, help = 'Number of CPUS to be used in the program. Default is 3.', default = 3)
    parser.add_argument('-updateschema' , required=False, help = 'Create a new schema with the new locus found. Default is True.', default = True)
    parser.add_argument('-percentlength', required=False, help = 'Allowed length percentage to be considered as ASM or ALM. Outside of this limit it is considered as LNF Default is 20.', default = 20)
    return parser.parse_args()

def is_fasta_file (file_name):
    with open (file_name, 'r') as fh:
        fasta = SeqIO.parse(fh, 'fasta')
        return any(fasta)

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

def check_blast (reference_allele, sample_files, db_name, logger) :
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

def parsing_fasta_file_to_dict (fasta_file, logger):
    fasta_dict = {}
    for contig in SeqIO.parse(fasta_file, "fasta", generic_dna):
            fasta_dict[contig.id] = str(contig.seq.upper())
    logger.debug('file %s parsed to dictionary', fasta_file)
    return fasta_dict

def prepare_core_gene(core_gene_file_list, store_dir, logger):
    #check if directory exists and files have fasta files
    #valid_core_gene_files = get_fasta_file_list(core_gene_dir, logger)
    #if not valid_core_gene_files :
    #    return False
    #logger.debug('Schema files to be processed are : %s', valid_core_gene_files)
    #processing the files in the schema
    schema_variability = {}
    schema_statistics = {}
    file_list = []
    first_alleles_list =[]
    blast_dir = os.path.join(store_dir,'blastdb')
    logger.info('start preparation  of core genes files')
    for fasta_file in core_gene_file_list:
        
        # parsing fasta file and get in the dictionary the id and the sequence
        fasta_file_parsed_dict = parsing_fasta_file_to_dict(fasta_file, logger)
        f_name = os.path.basename(fasta_file).split('.')
        file_list.append(os.path.join(store_dir, f_name[0]))
        # dump fasta file into pickle file
        with open (file_list[-1],'wb') as f:
            pickle.dump(fasta_file_parsed_dict, f)
        # create the first allele for each core gene file    
        #### used only for gene annotation
        first_alleles_list.append(write_first_allele_seq(fasta_file, store_dir, logger))
        alleles_len = []
        for allele in fasta_file_parsed_dict :
            alleles_len.append(len(fasta_file_parsed_dict[allele]))
        
        schema_variability[f_name[0]]=list(set(alleles_len))
        schema_statistics[f_name[0]]=[statistics.mode(alleles_len), min(alleles_len), max(alleles_len)]

    return file_list , first_alleles_list , schema_variability, schema_statistics
    
def prepare_samples( sample_file_list, store_dir, logger):
    file_list = []
    blast_dir = os.path.join(store_dir,'blastdb')
    
    for fasta_file in sample_file_list:
        # parsing fasta file and get in the dictionary the id and the sequence
        fasta_file_parsed_dict = parsing_fasta_file_to_dict(fasta_file, logger)
        f_name = os.path.basename(fasta_file).split('.')
        file_list.append(os.path.join(store_dir, f_name[0]))
        # dump fasta file into pickle file
        with open (file_list[-1],'wb') as f:
            pickle.dump(fasta_file_parsed_dict, f)
        
        # create local blast db for each core gene fasta file
        if not create_blastdb(fasta_file, blast_dir, 'nucl' ,logger):
            print('Error when creating the blastdb for core gene files. Check log file for more information. \n ')
            return False

    return file_list


def get_gene_annotation (annotation_file, annotation_dir, logger) :
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

def analize_annotation_files (in_file, logger) :
    examiner = GFF.GFFExaminer()
    file_fh = open(in_file) 
    datos = examiner.available_limits(in_file)
    file_fh.close()
    return True


    
def get_inferred_allele_number(core_dict, logger):
    #This function will look for the highest locus number and it will return a safe high value
    # that will be added to the schema database 
    logger.debug('running get_inferred_allele_number function')
    int_keys = []
    for key in core_dict.keys():
        int_keys.append(key)
    max_value = max(int_keys)
    digit_length = len(str(max_value))
    # return any of the values ( 10000, 100000, 1000000 and so on ) according to the bigest allele number used
    return  True   #str 1 ( #'1'+ '0'*digit_length + 2)
    
    

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


def check_sequence_order(allele_sequence, logger) :
    start_codon_forward= ['ATG','ATA','ATT','GTG', 'TTG']
    start_codon_reverse= ['CAT', 'TAT','AAT','CAC','CAA']
    # check forward direction
    if allele_sequence[0:3] in start_codon_forward :
        return 'forward'
    if allele_sequence[len(allele_sequence) -3: len(allele_sequence)] in start_codon_reverse :
        return 'reverse'
    return False

def get_snp_2(sample, query) :
    snp_list = []
    for i in range(len(sample)):
        try:
            if sample[i] != query[i] :
                snp_list.append([str(i+1), query[i], sample[i]])
        except:
            snp_list.append([str(i+1), '-', sample[i]])
    return snp_list

def get_snp (sample, query) :
    prot_annotation = {'S': 'polar' ,'T': 'polar' ,'Y': 'polar' ,'Q': 'polar' ,'N': 'polar' ,'C': 'polar' ,'S': 'polar' ,
                        'F': 'nonpolar' ,'L': 'nonpolar','I': 'nonpolar','M': 'nonpolar','P': 'nonpolar','V': 'nonpolar','A': 'nonpolar','W': 'nonpolar','G': 'nonpolar',
                        'D' : 'acidic', 'E' :'acidic',
                        'H': 'basic' , 'K': 'basic' , 'R' : 'basic',
                        '-': '-----', '*' : 'Stop codon'}
    
    snp_list = []
    length = max(len(sample), len(query))
    # normalize the lenght of the sample for the iteration
    if len(sample) < length :
        need_to_add = length - len(sample)
        sample = sample + need_to_add * '-'
    if len(query) < length :
        need_to_add = length - len(query)
        query = query + need_to_add * '-'
    # convert to Seq class to translate to protein
    seq_sample = Seq.Seq(sample)
    seq_query = Seq.Seq(query)
    
    
    for index in range (0, length -2, 3) :
        codon_seq = seq_sample[index : index + 3]
        codon_que = seq_query[index : index + 3]
        if codon_seq != codon_que :
            if str(codon_seq) != '---' :
                prot_seq = str(codon_seq.translate())
            else:
                prot_seq = '-'
            if str(codon_que) != '---' :
                prot_que = str(codon_que.translate())
            else:
                prot_que = '-'
            snp_list.append([str(index+1),str(codon_seq) + '/'+ str(codon_que), prot_seq + '/' + prot_que, prot_annotation[prot_seq] + ' / ' + prot_annotation[prot_que]])
        
            
    
    
    return snp_list

def convert_to_protein (sequence) :

    seq = Seq.Seq(sequence)
    protein = str(seq.translate())
    
    return protein

def nucleotide_to_protein_aligment (sample_seq, query_seq ) :
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

def get_alignment_for_indels (blast_db_name, qseq) :
    #match_alignment =[]
    cline = NcbiblastnCommandline(db=blast_db_name, evalue=0.001, perc_identity = 80, outfmt= 5, max_target_seqs=10, max_hsps=10,num_threads=3)
    out, err = cline(stdin = qseq)
    psiblast_xml = StringIO(out)
    blast_records = NCBIXML.parse(psiblast_xml)   
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for match in alignment.hsps:
                match_alignment = [['sample', match.sbjct],['match', match.match], ['schema',match.query]]
    return match_alignment


def get_aligments_for_deletions (sample_seq, query_seq):
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

def create_summary (samples_matrix_dict, logger) :
    summary_dict = {}
    summary_result_list = []
    summary_heading_list = ['Exact match', 'INF', 'ASM_INSERT', 'ASM_DELETE','ALM_INSERT' ,'ALM_DELETE', 'LNF','NIPH','NIPHEM','PLOT']
    summary_result_list.append('File\t' + '\t'.join(summary_heading_list))
    for key in sorted (samples_matrix_dict) :
        summary_dict[key] = {'Exact match':0, 'INF':0, 'ASM_INSERT':0, 'ASM_DELETE':0, 'ALM_INSERT':0, 'ALM_DELETE':0, 'LNF':0, 'NIPH':0, 'NIPHEM':0, 'PLOT':0}
        for values in samples_matrix_dict[key] :
            if 'INF_' in values :
                summary_dict[key]['INF'] += 1
            elif 'ASM_INSERT' in values :
                summary_dict[key]['ASM_INSERT'] += 1
            elif 'ASM_DELETE' in values :
                summary_dict[key]['ASM_DELETE'] += 1
            elif 'ALM_INSERT' in values :
                summary_dict[key]['ALM_INSERT'] += 1
            elif 'ALM_DELETE' in values :
                summary_dict[key]['ALM_DELETE'] += 1
            elif 'LNF' in values :
                summary_dict[key]['LNF'] += 1
            elif 'NIPH' in values :
                summary_dict[key]['NIPH'] += 1
            elif 'NIPHEM' in values :
                summary_dict[key]['NIPHEM'] += 1    
            elif 'PLOT' in values :
                summary_dict[key]['PLOT'] += 1
            else:
                try:
                    number =int(values)
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
    



def allele_call_nucleotides ( core_gene_dict_files, reference_query_directory,  sample_dict_files, blast_db_directory, inputdir, outputdir, cpus , percentlength, schema_variability, logger ):
    full_gene_list = []
    samples_matrix_dict = {} # to keep allele number
    matching_genes_dict = {} # to keep start and stop positions
    inferred_counter = 0
    inferred_alleles_dict = {} # to keep track of the new inferred alleles
    inf_dict = {} # Store the inferred alleles found for each sample
    paralog_dict = {}
    insertions_dict = {}
    deletions_dict = {}
    list_insertions = {} # list all insertions together with Sample file and core gene
    list_deletions = {} # list all deletions together with Sample file and core gene
    plot_dict = {}
    snp_dict = {}
    protein_dict = {}
    match_alignment_dict = {}
    blast_parameters = '"6 , qseqid , sseqid , pident ,  qlen , length , mismatch , gapopen , evalue , bitscore , sstart , send , qstart , qend , sseq , qseq"'
    header_macthing_alleles_conting = ['Sample Name', 'Contig', 'Core Gene','start', 'stop', 'direction', 'codification']
    header_paralogs = ['Sample Name','Core Gene', 'Allele','Contig','Bit Score', 'Start Seq', 'End Seq','Sequence']
    header_inferred = ['Sample Name','Core Gene', 'Allele','Contig','Bit Score', 'Start Seq', 'End Seq','Sequence']
    header_insertions = [ 'Core Gene', 'Sample Name' , 'Insertion item' ,'Allele', 'Contig', 'Bitscore', 'Query length' , 'Contig length', 'New sequence length' , 'Mismatch' , 'gaps', 'Contig start', 'Contig end',  'New sequence' ]
    header_deletions = [ 'Core Gene', 'Sample Name' , 'Deletion item' ,'Allele', 'Contig', 'Bitscore', 'Query length' , 'Contig length', 'New sequence length' , 'Mismatch' , 'gaps', 'Contig start', 'Contig end',  'New sequence' ]
    header_plot = ['Core Gene', 'Sample Name' , 'Allele','Contig','Bit Score', 'Start Seq', 'End Seq','Sequence']
    header_snp = ['Sample Name','Core Gene', 'Position','Sequence Sample/Schema','Protein in Sample/Schema', 'Annotation Sample / Schema']
    header_protein = ['Sample Name','Core Gene', 'Protein in ' , 'Protein sequence']
    header_match_alignment = ['Sample Name','Core Gene','Alignment', 'Sequence']
    
    for core_file in core_gene_dict_files:
        print ( 'Analyzing core file : ', core_file)
        full_gene_list.append(os.path.basename(core_file))
        logger.info('Processing core gene file %s ', core_file)
        core_name = os.path.basename(core_file)
        reference_query = os.path.join(reference_query_directory, str( core_name + '.fasta'))
        with open (core_file, 'rb') as core_f:
            core_dict = pickle.load(core_f)
        logger.debug('load in memory the core file %s ', core_file)      
        ref_query_parse = list (SeqIO.parse(reference_query, "fasta"))
        query_length = len(ref_query_parse[0].seq)
        #query_length_list =[]
        '''
        for allele in ref_query_parse :
            allele_length =  len(allele.seq)
            if not allele_length in query_length_list :
                query_length_list.append(allele_length)
        '''
        
        #create new_allele_dict to infer
        new_allele_dict ={}
        samples_inferred = []
        #allele_list_per_sample = []
        for sample_file in sample_dict_files:
            #print('sample file is: ', sample_file)
            #with open (sample_file,'rb') as sample_f :
            #    sample_dict = pickle.load(sample_f)
            #logger.debug('loaded in memory the sample file %s' , sample_file)
            
            sample_value = os.path.basename(sample_file)
            if not sample_value in samples_matrix_dict:
                # initialize the sample list to add the number of alleles and the start, stop positions 
                samples_matrix_dict[sample_value] = []
                matching_genes_dict[sample_value] = {}
            #intersection = set(core_dict.values()).intersection(gene_dict.values())
            blast_db_name = os.path.join(blast_db_directory, os.path.basename(sample_file),os.path.basename(sample_file))
            #blast_db_name = '/srv/project_wgmlst/samples_listeria/RA-L2073/blastdb'
            #reference_query = '/srv/project_wgmlst/lmo_test.fasta'
            
            cline = NcbiblastnCommandline(db=blast_db_name, evalue=0.001, perc_identity = 100, outfmt = blast_parameters , max_target_seqs=2, max_hsps=1,num_threads=1, query=reference_query)
            #cline = NcbiblastnCommandline(db=Gene_Blast_DB_name, evalue=0.001, outfmt=5, max_target_seqs=10, max_hsps=10,num_threads=1, query='/srv/project_wgmlst/seqSphere_listeria_cgMLST_test/targets/lmo0001.fasta')
            out, err = cline()
            out_lines = out.splitlines( )
            
            if len (out_lines) > 0 :
                
                bigger_bitscore = 0
                allele_found = {}
                
                for line in out_lines :
                    values = line.split('\t')
                    
                    s_length = values[4]
                    
                    #if int(s_length) in query_length_list :
                    if int(s_length) in schema_variability[core_name] :
                    #if int(s_length) == int(query_length) :
                        contig_id = values[1]
                        gene_start = values[9]
                        gene_end = values[10]
                        
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
                        if  int(values[8]) > bigger_bitscore :
                            #qseqid , sseqid , pident ,  qlen , length , mismatch , gapopen , evalue , bitscore , sstart , send , qstart , qend ,sseq , qseq= values
                            #bigger_bitscore = int(bitscore)
                            bigger_bitscore = int(values[8])
                
                if len(allele_found) > 1:
                    # found paralogs in the sample for the core gene
                    samples_matrix_dict[sample_value].append('NIPHEM')
                    if not sample_value in paralog_dict :
                        paralog_dict[sample_value] = {}
                    if not core_name in paralog_dict[sample_value] :
                        paralog_dict[sample_value] [core_name]= []
                    for allele_item in allele_found :
                        sstart = allele_found[allele_item][9]
                        send = allele_found[allele_item][10]
                        sseqid = allele_found[allele_item][1]
                        qseqid = allele_found[allele_item][0]
                        bitscore = allele_found[allele_item][8]
                        sseq = allele_found[allele_item][13]
                        paralog_dict[sample_value][core_name].append([qseqid,sseqid,bitscore,sstart, send, sseq])
                        if not sseqid in matching_genes_dict[sample_value] :
                            matching_genes_dict[sample_value][sseqid] = []
                        if sstart > send :
                            matching_genes_dict[sample_value][sseqid].append([core_name, sstart,send,'-','NIPHEM'])
                        else:
                            matching_genes_dict[sample_value][sseqid].append([core_name, sstart,send,'+', 'NIPHEM'])
                    continue
                    
                elif len(allele_found) == 1 :
                    ## look for possible paralogos by finding other alleles that identity is  > 90%
                    paralog_found ={}
                    allele_sequence = allele_found[contig_id_start][14]
                    cline = NcbiblastnCommandline(db=blast_db_name, evalue=0.001, perc_identity = 90, outfmt= blast_parameters, max_target_seqs=10, max_hsps=10,num_threads=3)
                    out, err = cline(stdin = allele_sequence)
                    out_lines = out.splitlines( )
                    for line in out_lines :
                        values = line.split('\t')
                        s_length = values[4]
                    
                        #if int(s_length) == int(query_length) :
                        if int(s_length) in schema_variability[core_name] :
                            contig_id = values[1]
                            gene_start = values[9]
                            gene_end = values[10]
                            contig_id_start = str(contig_id + '_'+ gene_start)
                            ## skip the allele found in the 100% identity and 100% alignment
                            if not contig_id_start in allele_found :
                                paralog_found[contig_id_start] = values
                    if len(paralog_found) == 0 :
                    # exact match found
                        qseqid = allele_found[contig_id_start][0]
                        sseqid = allele_found[contig_id_start][1]
                        sstart = allele_found[contig_id_start][9]
                        send = allele_found[contig_id_start][10]
                        samples_matrix_dict[sample_value].append(qseqid)
                        if not sseqid in matching_genes_dict[sample_value] :
                            matching_genes_dict[sample_value][sseqid] = []
                        # store the matching genes in forward order
                        if sstart > send :
                            matching_genes_dict[sample_value][sseqid].append([core_name, sstart,send,'-','EXACT'])
                        else:
                            matching_genes_dict[sample_value][sseqid].append([core_name, sstart,send,'+','EXACT'])
                        
                        continue
                    else:
                        # paralog has been found
                        paralog_matrix = {}
                        samples_matrix_dict[sample_value].append('NIPH')
                        if not sample_value in paralog_dict :
                            paralog_dict[sample_value] = {}
                        if not core_name in paralog_dict[sample_value] :
                            paralog_dict[sample_value] [core_name]= []
                        # merging the 2 dictionary
                        paralog_matrix[sample_value] = {**allele_found, **paralog_found}
                        
                        for paralog in paralog_matrix[sample_value] :
                            sstart = paralog_matrix[sample_value][paralog][9]
                            send = paralog_matrix[sample_value][paralog] [10]
                            sseqid = paralog_matrix[sample_value][paralog] [1]
                            qseqid = paralog_matrix[sample_value][paralog] [0]
                            bitscore = paralog_matrix[sample_value][paralog] [8]
                            sseq = paralog_matrix[sample_value][paralog] [13]
                            paralog_dict[sample_value][core_name].append([qseqid,sseqid,bitscore,sstart, send, sseq])
                            if not sseqid in matching_genes_dict[sample_value] :
                                matching_genes_dict[sample_value][sseqid] = []
                            if sstart > send :
                                matching_genes_dict[sample_value][sseqid].append([core_name, sstart,send,'-', 'NIPH'])
                            else:
                                matching_genes_dict[sample_value][sseqid].append([core_name, sstart,send,'+', 'NIPH'])

                        continue

            #print('blast len is 0  or not full length was matched')
            cline = NcbiblastnCommandline(db=blast_db_name, evalue=0.001, perc_identity = 80, outfmt= blast_parameters, max_target_seqs=1, max_hsps=1,num_threads=1, query=reference_query)
            #cline = NcbiblastnCommandline(db=Gene_Blast_DB_name, evalue=0.001, outfmt=5, max_target_seqs=10, max_hsps=10,num_threads=1, query='/srv/project_wgmlst/seqSphere_listeria_cgMLST_test/targets/lmo0001.fasta')
            out, err = cline()
            out_lines = out.splitlines( )
            bigger_bitscore = 0
            if len (out_lines) == 0:
                samples_matrix_dict[sample_value].append('LNF')
                logger.info('Locus not found at sample %s, for gene  %s', sample_value, core_name)
                continue
            for line in out_lines :
                values = line.split('\t')
                if  int(values[8]) > bigger_bitscore:
                    qseqid , sseqid , pident ,  qlen , s_length , mismatch , gapopen , evalue , bitscore , sstart , send , qstart , qend ,sseq , qseq = values
                    #print('q len seq is : ', len(qseq), ' s len seq is : ', len(sseq))
                    bigger_bitscore = int(bitscore)
            #cline = NcbiblastnCommandline(db=Gene_Blast_DB_name, evalue=0.001, outfmt=5, max_target_seqs=10, max_hsps=10,num_threads=1, query='/srv/project_wgmlst/seqSphere_listeria_cgMLST_test/targets/lmo0001.fasta')= values
            #print ( 'number of matches is : ', len(out_lines))
            #print ('qlen is: ',qlen, ' seq_len is : ', length , 'query_reference_length is : ', query_length)
            #print('mistmatch is : ', mismatch, 'gaps is : ', gapopen)
            #print('q start : ', qstart, ' q end : ', qend )
            #print ('s start : ', sstart, ' s end', send)
            
            
            if int(s_length) in schema_variability[core_name] :
            
            #if int(s_length) == int(query_length) : # if equal then a new allele has been detected
                logger.info('Found new allele for core gene %s ', core_name)
                if not sample_value in inf_dict :
                    inf_dict[sample_value] = {}
                if not core_name in inf_dict[sample_value] :
                    inf_dict[sample_value] [core_name]= []
                ### adding new allele to the  inferred allele list if it is not already included                
                if not core_name in inferred_alleles_dict :
                    inferred_alleles_dict[core_name] = []
                if not sseq in inferred_alleles_dict[core_name] :
                    inferred_alleles_dict[core_name].append(sseq)
                ### find the index to include in the sample matrix dict
                index_inferred = inferred_alleles_dict[core_name].index(sseq)
                inferred_allele = 'INF_' + core_name + '_' + str(index_inferred)
                samples_matrix_dict[sample_value].append(inferred_allele)
                inf_dict[sample_value][core_name].append([qseqid,sseqid,bitscore,sstart, send, sseq])
                
                if not sseqid in matching_genes_dict[sample_value] :
                    matching_genes_dict[sample_value][sseqid] = []
                if sstart > send :
                    matching_genes_dict[sample_value][sseqid].append([core_name, sstart,send,'-',inferred_allele])
                else:
                    matching_genes_dict[sample_value][sseqid].append([core_name, sstart,send,'+',inferred_allele])
                continue
            
            #if int(s_length) / int(query_length) < ( 1- percentlength/100)  or int(s_length) / int(query_length) > (1 + percentlength/100) :
            #    samples_matrix_dict[sample_value].append('LNF')
            #    logger.info('Locus found at sample %s, for gene  %s', sample_value, core_name)
                
            #    continue
            
            alleles_in_gene = list (SeqIO.parse(reference_query, "fasta"))
            for allele_item in alleles_in_gene :
                if allele_item.id == qseqid :
                    break
            allele_sequence = allele_item.seq
                    
                    
                    
            #if int(s_length) < min(schema_variability[core_name]) : 
            if int(s_length) < int(query_length) :
                ## check if the blast alignment could be clasified as PLOT
                seq_id_split = sseqid.split('_')
                length_sseqid = seq_id_split[3]
                if sstart == length_sseqid or send == length_sseqid:
                    samples_matrix_dict[sample_value].append('PLOT')
                    logger.info('PLOT found at sample %s, for gene  %s', sample_value, core_name)
                    if sample_value not in plot_dict :
                        plot_dict[sample_value] = {}
                    if not core_name in plot_dict[sample_value] :
                        plot_dict[sample_value][core_name] = []
                    plot_dict[sample_value][core_name].append([qseqid,sseqid,bitscore,sstart, send, sseq])
                    continue
                else:
                    # print ('There is a deletion of ', gapopen,'gaps', 'or shorter mapping')
                    # print ('qlen is: ',qlen, ' seq_len is : ', length,  'query_reference_length is : ', query_length)
                    # print('mistmatch is : ', mismatch, 'gaps is : ', gapopen)
                    # print('q start : ', qstart, ' q end : ', qend )
                    # print ('s start : ', sstart, ' s end', send) 

                    
                    query_direction = check_sequence_order(allele_sequence, logger)
                    contig_file = os.path.join(inputdir,str(sample_value + '.fasta'))
                    records = list (SeqIO.parse(contig_file, "fasta"))
                    #accession_sequence = records[accession]
                    for record in records:
                        if record.id == sseqid :
                            break
                    accession_sequence = record.seq
                    
                    if allele_sequence.endswith ('TGA') or  allele_sequence.startswith ('TCA') :
                        tga_stop_codon = True
                    else:
                        tga_stop_codon = False
                        #it is assume that reference query is in reverse complement. 
                        #tga_stop_codon = allele_sequence.startswith('TCA')
                        #bassed_added_start = int(qlen) -len(qseq) - int(qstart)
                    if query_direction == 'reverse' :
                        if int(send) > int (sstart): ## increasin the number of nucleotides to check if getting  longer protein
                            sample_gene_sequence = accession_sequence[int(sstart) - 51 :  int(send)  ]
                            sample_gene_sequence = sample_gene_sequence.reverse_complement()
                        else:
                                sample_gene_sequence = accession_sequence[int(send) -1 : int(sstart)  + 51]
                    else:
                        if int(sstart) > int (send):
                            sample_gene_sequence = accession_sequence[int(send) - 51 :  int(sstart)  ]
                            sample_gene_sequence = sample_gene_sequence.reverse_complement()
                        else:
                            sample_gene_sequence = accession_sequence[int(sstart) -1 : int(send)  + 51]
                    #else:
                    #    tga_stop_codon = allele_sequence.endswith ('TGA') or  allele_sequence.startswith ('TCA')
                    #    bassed_added_end = int(qlen) -len(qseq) 
                    #    if query_direction == 'forward' :
                            # increasing the number of nucleotides in the sequence to find the stop codon when becomes longer because of the deletion
                    #        sample_gene_sequence = accession_sequence[int(send) - bassed_added_end - 51: int(sstart)  ]
                    #        sample_gene_sequence = sample_gene_sequence.reverse_complement()
                    #    else:
                            # increasing the number of nucleotides in the sequence to find the stop codon when becomes longer because of the deletion
                    #        sample_gene_sequence = accession_sequence[int(send) -1: int(sstart) + bassed_added_end + 51]
                    stop_index = get_stop_codon_index(sample_gene_sequence, tga_stop_codon, int(qlen)- int(qstart))
                    if stop_index != False:
                        new_sequence_lenght = stop_index +3                      
                        new_sseq = str(sample_gene_sequence[0:new_sequence_lenght])
                        
                        ### adding ASM allele to the asm_allele_matrix if it is not already include
                        if not core_name in deletions_dict :
                            deletions_dict[core_name] = []
                        if not new_sseq in deletions_dict[core_name] :
                            deletions_dict[core_name].append(new_sseq)
                        ### find the index of ASM  to include it in the sample matrix dict
                        index_delete = deletions_dict[core_name].index(new_sseq)
                        if new_sequence_lenght < query_length :
                            delete_allele = 'ASM_DELETE_' + core_name + '_' + str(index_delete)
                        elif new_sequence_lenght == query_length:
                            delete_allele = 'AEM_DELETE_' + core_name + '_' + str(index_delete)
                        else:
                            delete_allele = 'ALM_DELETE_' + core_name + '_' + str(index_delete)
                        samples_matrix_dict[sample_value].append(delete_allele)
                        
                        if not sseqid in matching_genes_dict[sample_value] :
                            matching_genes_dict[sample_value][sseqid] = []
                        if sstart > send :
                            matching_genes_dict[sample_value][sseqid].append([core_name, str(int(sstart)-new_sequence_lenght -1), sstart,'-', delete_allele])
                        else:
                            matching_genes_dict[sample_value][sseqid].append([core_name, sstart,str(int(sstart)+ new_sequence_lenght),'+', delete_allele])
                        
                        ### add the deletion into deletion list
                        if not core_name in list_deletions :
                            list_deletions [core_name] = {}
                        if not sample_value in list_deletions[core_name] :
                            list_deletions[core_name][sample_value] = {}
                        list_deletions[core_name][sample_value][delete_allele] = [qseqid, sseqid,  bitscore, str(query_length) , s_length, str(new_sequence_lenght), mismatch , gapopen, sstart, send,  new_sseq ]
                        
                        if check_sequence_order(qseq, logger) == 'reverse' :
                            qseq = str(allele_sequence.reverse_complement())
                        else:
                            qseq = str(allele_sequence)
                        # get the SNP for the  delection
                        if not core_name in snp_dict :
                            snp_dict[core_name] = {}
                        if not sample_value in snp_dict[core_name] :
                            snp_dict[core_name][sample_value] = []
                        snp_dict[core_name][sample_value] = get_snp(new_sseq, qseq)
                        
                        
                        # execute again blast with the reference query the previous query found to get the aligment format to get the SNPs
                        if not core_name in match_alignment_dict :
                            match_alignment_dict[core_name] = {}
                            if not sample_value in match_alignment_dict[core_name] :
                                match_alignment_dict[core_name][sample_value] = get_aligments_for_deletions (new_sseq,  str(qseq))
                                               
                        
                        
                        
                        
                        
                        # convert the sequence to protein
                        if not core_name in protein_dict :
                            protein_dict[core_name] = {}
                        if not sample_value in protein_dict[core_name] :
                            protein_dict[core_name][sample_value] = []
                        protein_dict[core_name][sample_value] = nucleotide_to_protein_aligment(new_sseq, qseq ) 
                    else:
                        logger.error('ERROR : Stop codon was not found for the core %s and the sample %s', core_name, sample_value)
                        
            #if int(s_length) > int(query_length) :
            elif int(s_length) > max(schema_variability[core_name]) :
            #elif int(s_length) > int(query_length) :   
                #print ('there is a insertion of  ', gapopen ,' bases in the sequence')
                #print ('qlen is: ',qlen, ' seq_len is : ', length,  'query_reference_length is : ', query_length)
                #query_seq = Seq.Seq(qseq)
                tga_stop_codon = qseq.endswith('TGA')
                sseq = sseq.replace('-','')
                stop_index = get_stop_codon_index(sseq, tga_stop_codon, qseq.find('-'))
                
                if stop_index != False:
                    new_sequence_lenght = stop_index +3
                    ### adding ASM allele to the asm_allele_matrix if it is not already include
                    new_sseq = sseq[0:new_sequence_lenght]
                    
                    if not core_name in insertions_dict :
                        insertions_dict[core_name] = []
                    if not new_sseq in insertions_dict[core_name] :
                        insertions_dict[core_name].append(new_sseq)
                    ### find the index of ASM  to include it in the sample matrix dict
                    index_insert = insertions_dict[core_name].index(new_sseq)
                    #if new_sequence_lenght < query_length :
                    if new_sequence_lenght < min(schema_variability[core_name]) :
                        insert_allele = 'ASM_INSERT_' + core_name + '_' + str(index_insert)
                    else:
                        insert_allele = 'AEM_INSERT_' + core_name + '_' + str(index_insert)
                    samples_matrix_dict[sample_value].append(insert_allele)
                else:
                    samples_matrix_dict[sample_value].append('ALM_INSERT_')
                if not sseqid in matching_genes_dict[sample_value] :
                    matching_genes_dict[sample_value][sseqid] = []
                if sstart > send :
                    matching_genes_dict[sample_value][sseqid].append([core_name, sstart,send,'-', insert_allele])
                else:
                    matching_genes_dict[sample_value][sseqid].append([core_name, sstart,send,'+', insert_allele])
                ### add the insertion into insertion list
                if not core_name in list_insertions :
                    list_insertions [core_name] = {}
                if not sample_value in list_insertions[core_name] :
                    list_insertions[core_name][sample_value] = {}
                list_insertions[core_name][sample_value][insert_allele] = [qseqid, sseqid,  bitscore, str(query_length) , s_length, str(new_sequence_lenght), mismatch , gapopen, sstart, send,  new_sseq ]
            
                if check_sequence_order(qseq, logger) == 'reverse' :
                    qseq = str(allele_sequence.reverse_complement())
                else:
                    qseq = str(allele_sequence)
                # get the SNP for the  delection
                if not core_name in snp_dict :
                    snp_dict[core_name] = {}
                if not sample_value in snp_dict[core_name] :
                    snp_dict[core_name][sample_value] = []
                snp_dict[core_name][sample_value] = get_snp(new_sseq, qseq)
                    
                '''    
                cline = NcbiblastnCommandline(db=blast_db_name, evalue=0.001, perc_identity = 80, outfmt= 5, max_target_seqs=10, max_hsps=10,num_threads=3)
                out, err = cline(stdin = qseq)
                psiblast_xml = StringIO(out)
                blast_records = NCBIXML.parse(psiblast_xml)   
                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        for match in alignment.hsps:
                            match_alignment = [['sample', match.sbjct],['match', match.match], ['schema',match.query]]
                ''' 
                if not core_name in match_alignment_dict :
                    match_alignment_dict[core_name] = {}
                if not sample_value in match_alignment_dict[core_name] :
                    match_alignment_dict[core_name][sample_value] = get_alignment_for_indels (blast_db_name, qseq) 
                # index_not_match = [m.start() for m in re.finditer(' ', match.match)]
                
                # convert the sequence to protein
                if not core_name in protein_dict :
                    protein_dict[core_name] = {}
                if not sample_value in protein_dict[core_name] :
                    #protein_dict[core_name][sample_value] = []
                    protein_dict[core_name][sample_value] = nucleotide_to_protein_aligment(new_sseq, qseq )
                
                # get the SNP from the alignment
                
            
            else:
                samples_matrix_dict[sample_value].append('ERROR ')
                
                print ('ERROR when looking the allele match for core gene ', core_name, 'at sample ', sample_value )

    logger.debug ('matching genes =  %s', matching_genes_dict)
    logger.debug ('---------------------------------------------------')
    logger.debug ('sample matrix  = %s', samples_matrix_dict)
    logger.debug ('---------------------------------------------------')
    logger.debug ('inferred alleles = %s', inferred_alleles_dict)
    logger.debug ('---------------------------------------------------')
    logger.debug ('inferred in the sequence  = %s' , inf_dict)
    logger.debug ('---------------------------------------------------')
    logger.debug ('insertions = %s',  insertions_dict)
    logger.debug ('---------------------------------------------------')
    logger.debug ('list of insertions = %s ', list_insertions)
    logger.debug ('---------------------------------------------------')
    logger.debug ('deletions = %s', deletions_dict)
    logger.debug ('---------------------------------------------------')
    logger.debug ('list of deletions = %s ', list_deletions)
    logger.debug ('---------------------------------------------------')
    logger.debug ('list of SNPs = %s', snp_dict)
    logger.debug ('---------------------------------------------------')
    logger.debug ('list of proteins = %s' , protein_dict)
    logger.debug ('---------------------------------------------------')

            
    

    # print ( 'valor de retorno ', samples_matrix_dict)
    result_file = os.path.join ( outputdir, 'result.tsv')
    # saving the reult information to file
    logger.info('Saving result information to file..')
    with open (result_file, 'w') as out_fh:
        out_fh.write ('Sample Name\t'+'\t'.join( full_gene_list) + '\n')
        for key in sorted (samples_matrix_dict):
            out_fh.write (key + '\t' + '\t'.join(samples_matrix_dict[key])+ '\n')
    
    # saving paralog sequence to file
    logger.info('Saving paralog information to file..')
    paralog_file =  os.path.join(outputdir, 'paralog.tsv')
    with open (paralog_file , 'w') as paralog_fh :
        paralog_fh.write('\t'.join(header_paralogs) + '\n')
        for sample in sorted (paralog_dict) :
            for core in sorted (paralog_dict[sample]):
                for paralog in paralog_dict[sample][core] :
                    paralog_fh.write(sample + '\t' + core + '\t' + '\t'.join (paralog) + '\n')
    
    # saving inferred alleles to file
    logger.info('Saving inferred alleles information to file..')
    inferred_file =  os.path.join(outputdir, 'inferred_alleles.tsv')
    with open (inferred_file , 'w') as infer_fh :
        infer_fh.write('\t'.join(header_inferred) + '\n')
        for key in sorted (inf_dict) :
            seq_in_inferred_allele = '\t'.join (inf_dict[key])
            infer_fh.write(key + '\t' + seq_in_inferred_allele + '\n')
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
    
    # saving insertions bases in contigs to file
    logger.info('Saving insert bases  information to file..')
    insertions_file =  os.path.join(outputdir, 'insertions.tsv')
    with open (insertions_file , 'w') as insertions_fh :
        insertions_fh.write('\t'.join(header_insertions )+ '\n')
        for c_gene in list_insertions :
            for sample in list_insertions[c_gene] :
                for insert in list_insertions[c_gene][sample]:
                    insertions_fh.write(c_gene + '\t' + sample + '\t' + insert + '\t' + '\t'.join(list_insertions[c_gene][sample][insert]) + '\n')
    
    # saving deletions bases in contigs to file
    logger.info('Saving deleted bases information to file..')
    deletions_file =  os.path.join(outputdir, 'deletions.tsv')
    with open (deletions_file , 'w') as deletions_fh :
        deletions_fh.write('\t'.join(header_deletions) + '\n')
        for c_gene in list_deletions :
            for sample in list_deletions[c_gene] :
                for delete in list_deletions[c_gene][sample]:
                    deletions_fh.write(c_gene + '\t' + sample + '\t' + delete + '\t' + '\t'.join(list_deletions[c_gene][sample][delete]) + '\n')
         
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
                for snp in snp_dict[core][sample] :
                    snp_fh.write(core + '\t' + sample + '\t' + '\t'.join (snp) + '\n')
    
    
    match_alignment_dict
    
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
    
        
            
    
    return True

if __name__ == '__main__' :
    version = ' Taranis  0.0.3'
    if sys.argv[1] == '-v' or sys.argv[1] == '--version':
        print( version, '\n')
        exit (0)
    arguments = check_arg(sys.argv[1:])
    start_time = datetime.now()
    # open log file
    logger = open_log ('taranis_wgMLST.log')
    # check additional programs are installed in your system
    if not check_prerequisites (logger):
        print ('your system does not fulfill the pre-requistes to run the script ')
        exit(0)
    ##############################################
    # Check that given directories contatin fasta files
    ##############################################
    valid_core_gene_files = get_fasta_file_list(arguments.coregenedir, logger)
    if not valid_core_gene_files :
        print ('There are not valid  fasta files in ',  arguments.coregenedir , ' directory. Check log file for more information ')
        exit(0)
    
    valid_sample_files = get_fasta_file_list(arguments.inputdir, logger)
    if not valid_sample_files :
        print ('There are not valid  fasta files in ',  arguments.inputdir , ' directory. Check log file for more information ')
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

    core_gene_dict_files , core_first_alleles_files, schema_variability , schema_statistics = prepare_core_gene (valid_core_gene_files , tmp_core_gene_dir , logger)
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
            logger.info ( 'Temporary folder %s  has been created again', tmp_samples_dir)
        except:
            logger.info('Unable to create again the temporary directory %s', tmp_samples_dir)
            shutil.rmtree(os.path.join(arguments.outputdir, 'tmp'))
            logger.info('Cleaned up temporary directory ', )
            print('Cannot create temporary directory on ', tmp_samples_dir, 'Check the log file to get more information \n')
            exit(0)
    sample_dict_files = prepare_samples (valid_sample_files, tmp_samples_dir, logger)
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
    #reference_query_directory = os.path.join(tmp_core_gene_dir,'first_alleles')
    ##########   Modified to get all alleles instead of the first one  #############
    reference_query_directory = arguments.coregenedir
    blast_db_directory = os.path.join(tmp_samples_dir,'blastdb')
    if not allele_call_nucleotides( core_gene_dict_files, reference_query_directory, sample_dict_files,  blast_db_directory, arguments.inputdir, arguments.outputdir,  int(arguments.cpus), int(arguments.percentlength) , schema_variability, logger):
        print('There is an error while processing the allele calling. Check the log file to get more information \n')
        exit(0)



print ('script ends ')
 