#!/usr/bin/env python3

from datetime import datetime
from progressbar import ProgressBar
from Bio.Blast.Applications import NcbiblastnCommandline
import numpy as np
from utils.taranis_utils import *


# · * · * · * · * · * · * · * · * · * · * · * · * · #  
# Get reference alleles for one locus in the schema #
# · * · * · * · * · * · * · * · * · * · * · * · * · #

def get_reference_allele(locus_quality, fasta_file, store_dir, evalue, perc_identity, reward, penalty, gapopen, gapextend, num_threads, logger):

    logger.info('Searching reference alleles')
 
    ## Create mash directory where to store temporal files
    f_name = os.path.basename(fasta_file).split('.')
    #full_path_reference_allele = os.path.join(store_dir, 'reference_alleles')
    #full_path_mash = os.path.join(full_path_reference_allele, 'mash')
    full_path_mash = os.path.join(store_dir, 'mash')
    full_path_locus_mash = os.path.join(full_path_mash, f_name[0])
    
    if not os.path.exists(full_path_locus_mash):
        try:
            os.makedirs(full_path_locus_mash)
            logger.info('Directory %s has been created', full_path_locus_mash)
        except:
            print ('Cannot create the directory ', full_path_locus_mash)
            logger.info('Directory %s cannot be created', full_path_locus_mash)
            exit (0)


    ## Split locus multifasta into fastas containing one allele sequence each
    alleles_in_locus_number = 0 # Get alleles in locus number to set max_target_seqs value in final BLAST search
    alleles_in_locus = [] # List to store alleles in locus IDs to intersect final BLAST search results
    for record in list(SeqIO.parse(fasta_file, "fasta")):
        alleles_in_locus.append(str(record.id))
        split_fasta_path = os.path.join(full_path_locus_mash, str(record.id) + ".fasta")
        alleles_in_locus_number += 1
        with open (split_fasta_path, 'w') as out_fh:
            out_fh.write ('>' + str(record.id) + '\n' + str(record.seq))


    ## Get mash sketch file to get pairwise sequences distances at a time
    sketch_path = os.path.join(full_path_locus_mash, "reference.msh")
    mash_sketch_command = ["mash", "sketch", "-o", sketch_path]
    
    # Get file paths to include in mash sketch file
    split_multifasta_files_list = get_fasta_file_list(full_path_locus_mash, logger)

    for fasta_path in split_multifasta_files_list:
        mash_sketch_command.append(fasta_path)

    mash_sketch_result = subprocess.run(mash_sketch_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


    ## Get pairwise allele sequences mash distances
    mash_distance_command = ["mash", "dist", sketch_path, sketch_path]
    mash_distance_result = subprocess.Popen(mash_distance_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    out, err = mash_distance_result.communicate()
    out = out.decode('UTF-8').split('\n')

    comp_dist_list = []
    for n in range(len(out)-1):
        comp = out[n].split('\t')
        comp_dist = float(comp[2])
        comp_dist_list.append(comp_dist)


    ## Get distances matrix and mean distances matrix
    comp_dist_list_per_allele = []
    alleles_number = len(split_multifasta_files_list)

    for index_distance in range(0, len(comp_dist_list), alleles_number):  
        dist_per_allele = comp_dist_list[index_distance : index_distance + alleles_number] 
        comp_dist_list_per_allele.append(dist_per_allele)

    comp_dist_arr_per_allele = np.asarray(comp_dist_list_per_allele)
    allele_mean_distance = np.mean(comp_dist_arr_per_allele, 0)
    

    ## Get reference allele (centroid): max average distance allele tagged as 'good_quality'
    min_mean = max(allele_mean_distance)
    ref_allele_id = str()

    for mean_index in range(len(allele_mean_distance)):
        if allele_mean_distance[mean_index] <= min_mean:
            allele_path = split_multifasta_files_list[mean_index]
            allele_id = os.path.basename(split_multifasta_files_list[mean_index]).split('.')[0]
            if locus_quality[allele_id] == 'good_quality':
                min_mean = allele_mean_distance[mean_index]
                ref_allele_id = allele_id


    ## Check that chosen reference allele represents every allele in the locus
    
    # Create local BLAST database for all alleles in the locus
    db_name = os.path.join(store_dir, 'locus_blastdb')
    if not create_blastdb(fasta_file, db_name, 'nucl', logger):
        print('Error when creating the blastdb for locus %s. Check log file for more information. \n ', f_name[0])
        return False

    locus_db_name = os.path.join(db_name, f_name[0], f_name[0])

    # All alleles in locus VS reference allele chosen (centroid) BLAST 
    blast_parameters = '"6 , qseqid , sseqid , pident ,  qlen , length , mismatch , gapopen , evalue , bitscore , sstart , send , qstart , qend , sseq , qseq"'
    cline = NcbiblastnCommandline(db=locus_db_name, evalue=evalue, perc_identity=perc_identity, reward=reward, penalty=penalty, gapopen=gapopen, gapextend=gapextend, outfmt=blast_parameters, max_target_seqs=alleles_in_locus_number, max_hsps=alleles_in_locus_number, num_threads=num_threads, query=allele_path)

    out, err = cline()
    out_lines = out.splitlines()

    # Keep not represented alleles along with the centroid as reference alleles of themselves
    alleles_in_blast = [] 

    for line in out_lines:
        values = line.split('\t') 
        alleles_in_blast.append(values[1])

    ids_intersect = list(set(alleles_in_locus) - set(alleles_in_blast))
    ids_intersect.insert(0, ref_allele_id)

    reference_file_path = os.path.join(store_dir, os.path.basename(fasta_file))
    with open (reference_file_path, 'w') as out_fh:
        for record in list(SeqIO.parse(fasta_file, "fasta")):
            if record.id in ids_intersect:
                out_fh.write ('>' + str(record.id) + '\n' + str(record.seq) + '\n')

    shutil.rmtree(full_path_locus_mash)
    shutil.rmtree(db_name)

    return True, full_path_mash


# · * · * · * · * · * · * · * · * · * · * · * · * · * #  
# Get reference alleles for every locus in the schema #
# · * · * · * · * · * · * · * · * · * · * · * · * · * #

def processing_reference_alleles (arguments) :
    '''
    Description:
        This is the main function for getting reference alleles.
        With the support of additional functions it will obtain reference alleles for each locus in the schema.
    Input:
        arguments   # Input arguments given on command line 
    Functions:
        
    Variables:  ?????

    Return: ?????
    '''

    start_time = datetime.now()
    print('Start the execution at :', start_time )
    # Open log file
    logger = open_log ('taranis_ref_alleles.log')
    print('Checking the pre-requisites.') 

    ############################################################
    ## Check additional programs are installed in your system ##
    ############################################################
    pre_requisite_list = [['blastp', '2.5'], ['makeblastdb' , '2.5'], ['mash', '1.1']]
    if not check_prerequisites (pre_requisite_list, logger):
        print ('your system does not fulfill the pre-requistes to run the script ')
        exit(0)

    ######################################################
    ## Check that given directories contain fasta files ##
    ######################################################
    print('Validating schema fasta files in ' , arguments.coregenedir , '\n')
    core_gene_files_list = get_fasta_file_list(arguments.coregenedir, logger)
    if not core_gene_files_list :
        print ('There are not valid fasta files in ',  arguments.coregenedir , ' directory. Check log file for more information ')
        exit(0)

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
            logger.info ('Output folder %s  has been created again', arguments.outputdir)
        except:
            logger.info('Unable to create again the output directory %s', arguments.outputdir)
            print('Cannot create output directory on ', arguments.outputdir)
            exit(0)

    #######################################################
    ## Obtain reference alleles for each locus in schema ##
    #######################################################
    logger.info('Getting reference alleles for each locus in schema')
    print('Getting reference alleles...')
    
    pbar = ProgressBar ()
    for fasta_file in pbar(core_gene_files_list):

        # Get core gene alleles quality
        locus_quality = check_core_gene_quality(fasta_file, logger)

        # Get core gene reference alleles
        complete_reference_alleles, full_path_mash = get_reference_allele(locus_quality, fasta_file, arguments.outputdir, float(arguments.evalue), float(arguments.perc_identity), int(arguments.reward), int(arguments.penalty), int(arguments.gapopen), int(arguments.gapextend), int(arguments.num_threads), logger)
        if not complete_reference_alleles: 
            print('There is an error while processing reference alleles. Check the log file to get more information \n')
            logger.info('Deleting the directory to clean up the temporary files created')
            shutil.rmtree(os.path.join(arguments.outputdir))
            exit(0)

    shutil.rmtree(full_path_mash)

    end_time = datetime.now()
    print('Completed execution at :', end_time )

    return True