#import analyze_schema
import sys
import argparse
import allele_calling
#import analyze_schema
import create_schema
import reference_alleles
from taranis_configuration import *

def check_arg (args=None) :
    '''
    The function is used for parsing the input parameters form the command line
    using the standard python package argparse. The package itself is handling
    the validation and the return errors messages
    Input:
        args    # Contains the arguments from the command line
    Variables:
        allele_calling_parser   # It is used for the allele calling input parameters
        evaluate_schema_parser  # It is used for the schema evaluation input parameters
        compare_schema_parser   # It is used for schema comparison input parameters
        create_schema_parser    # It is used for create an schema input parameters
    Return:
        parser.parse_args()     # The variable contains the valid parameters
    '''
    parser = argparse.ArgumentParser(prog = 'taranis.py',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description= 'Taranis is a set of utilities related to cgMSLST ')

    parser.add_argument('--version', action='version', version='%(prog)s 0.3.5')

    subparser = parser.add_subparsers(help = 'allele calling /interactive/schema '
                                      + 'are the available actions to execute taranis',
                                      dest = 'chosen_action')

    ### Input parameters for allele calling option
    allele_calling_parser = subparser.add_parser('allele_calling',
                                    help = 'Allele calling way to downloads the  schema locus')
    allele_calling_parser.add_argument('-coregenedir', required= True,
                                    help = 'Directory where the core gene files are located ')
    allele_calling_parser.add_argument('-refalleles', required= True,
                                    help = 'Directory where the core gene references files are located ') ### cambiando/modificando: añadiendo path a alelos de referencia de cada locus del esquema
    allele_calling_parser.add_argument('-inputdir', required= True,
                                    help ='Directory where are located the sample fasta files')
    allele_calling_parser.add_argument('-refgenome', required= True,
                                    help = 'Reference genome file for genes prediction') ### cambiando/modificiando: introduciendo genoma de referencia para predicción de genes con prodigal
    allele_calling_parser.add_argument('-outputdir', required= True,
                                    help = 'Directory where the result files will be stored')
    allele_calling_parser.add_argument('-cpus', required= False,
                                    help = 'Number of CPUS to be used in the program. Default is 1.',
                                    default = 1)
    allele_calling_parser.add_argument('-percentlength', required=False,
                                    help = 'Allowed length percentage to be considered as INF. '
                                    + 'Outside of this limit it is considered as ASM or ALM. Default is SD.',
                                    default = 'SD') ### c/m: percentlength por defecto SD
    allele_calling_parser.add_argument('-coverage', required=False,
                                    help = 'Coverage threshold to exclude found sequences. '
                                    + 'Outside of this limit it is considered LNF. Default is 50 %.',
                                    default = 50) ### c/m: incluyendo -coverage como argumento. 50% por defecto de momento
    allele_calling_parser.add_argument('-evalue', required=False,
                                    help = 'E-value in BLAST searches. Default is 0.001. ',
                                    default = 0.001) ### c/m: introduciendo evalue como argumento
    allele_calling_parser.add_argument('-perc_identity_ref', required=False,
                                    help = 'Identity percent in BLAST searches using reference alleles for each locus detection in samples. Default is 90 %. ',
                                    default = 90) ### c/m: introduciendo perc_ident_ref como argumento
    allele_calling_parser.add_argument('-perc_identity_loc', required=False,
                                    help = 'Identity percent in BLAST searches using all alleles in each locus for allele identification in samples. Default is 90 %. ',
                                    default = 90) ### c/m: introduciendo perc_ident_loc como argumento
    allele_calling_parser.add_argument('-reward', required=False,
                                    help = 'Match reward in BLAST searches. Default is 1. ',
                                    default = 1) ### c/m: introduciendo reward como argumento    
    allele_calling_parser.add_argument('-penalty', required=False,
                                    help = 'Mismatch penalty in BLAST searches. Default is -2. ',
                                    default = -2) ### c/m: introduciendo penalty como argumento
    allele_calling_parser.add_argument('-gapopen', required=False,
                                    help = 'Gap open penalty in BLAST searches. Default is 1. ',
                                    default = 1) ### c/m: introduciendo gapopen como argumento
    allele_calling_parser.add_argument('-gapextend', required=False,
                                    help = 'Gap extension penalty in BLAST searches. Default is 1. ',
                                    default = 1) ### c/m: introduciendo gapextend como argumento       
    allele_calling_parser.add_argument('-max_target_seq', required=False,
                                    help = 'max_target_seq in BLAST searches. Default is 10. ', ########## BUSCAR QUÉ ERA ESTO
                                    default = 10) ### c/m: introduciendo max_target_seq como argumento    
    allele_calling_parser.add_argument('-max_hsps', required=False,
                                    help = 'max_hsps in BLAST searches. Default is 10. ', ########## BUSCAR QUÉ ERA ESTO
                                    default = 10) ### c/m: introduciendo max_hsps como argumento
    allele_calling_parser.add_argument('-num_threads', required=False,
                                    help = 'num_threads in BLAST searches. Default is 1. ',
                                    default = 1) ### c/m: introduciendo num_threads como argumento   
    allele_calling_parser.add_argument('-flankingnts' , required=False,
                                    help = 'Number of flanking nucleotides to add to each BLAST result obtained after locus detection in sample using reference allele for correct allele identification. Default is 100. ',
                                    default = 100)  ### c/m: introduciendo flankingnts como argumento  
    allele_calling_parser.add_argument('-updateschema' , required=False,
                                    help = 'Add INF alleles found for each locus to the analysis schema. Default is True.',
                                    default = True)                                          ### Dejar esto así?
    allele_calling_parser.add_argument('-updatenewschema' , required=False,
                                    help = 'Add INF alleles found for each locus to a new schema. Default is False.',
                                    default = False)                                         ### Dejar esto así?
   


    ### Input parameters for schema evaluation options
    evaluate_schema_parser = subparser.add_parser('evaluate_schema',
                                    help = 'Evaluate the schema.')
    evaluate_schema_parser.add_argument('-inputdir',required= True,
                                    help = 'Directory where are the schema files.')
    evaluate_schema_parser.add_argument('-outputdir', required= True,
                                    help = 'Directory where the result files will be stored.')
    evaluate_schema_parser.add_argument('-alt', required = False, action = "store_true" ,
                                    help = 'Set this parameter if alternative start codon should be considered. '
                                    + 'Do not include to accept only ATG as a start codon.',
                                    default = False)

    ### Input parameters for schema comparison options
    compare_schema_parser = subparser.add_parser('compare_schema', help = 'Compare 2 schema.')
    compare_schema_parser.add_argument('-scheme1',
                                       help = 'Directory where are the schema files for the schema 1.')
    compare_schema_parser.add_argument('-scheme2',
                                       help = 'Directory where are the schema files for the schema 2.')

    ### Input parameters for schema creation options
    create_schema_parser = subparser.add_parser('create_schema', help = 'Create a schema.')
    create_schema_parser.add_argument('-xlsfile',
                                      help = 'xls file name which contains the list of the core genes.')
    create_schema_parser.add_argument('-outputdir', help = 'Directory where the core gene files '
                                      + 'will be stored. If directory exists it will be prompt for '
                                      + 'deletion confirmation.')

    ### Input parameters for reference alleles options   ### introduciendo script para obtener alelos de referencia
    reference_alleles_parser = subparser.add_parser('reference_alleles', help = 'Obtain reference allele(s) for each locus.')
    reference_alleles_parser.add_argument('-coregenedir', required= True,
                                    help = 'Directory where the core gene files are located ')
    reference_alleles_parser.add_argument('-outputdir', required= True,
                                    help = 'Directory where the result files will be stored')
        allele_calling_parser.add_argument('-evalue', required=False,
                                    help = 'E-value in BLAST searches. ',
                                    default = 0.001) ### c/m: introduciendo evalue como argumento
    reference_alleles_parser.add_argument('-perc_identity', required=False,
                                    help = 'Identity percent in BLAST searches. ',
                                    default = 90) ### c/m: introduciendo perc_ident como argumento
    reference_alleles_parser.add_argument('-reward', required=False,
                                    help = 'Match reward in BLAST searches. ',
                                    default = 1) ### c/m: introduciendo reward como argumento    
    reference_alleles_parser.add_argument('-penalty', required=False,
                                    help = 'Mismatch penalty in BLAST searches. ',
                                    default = -2) ### c/m: introduciendo penalty como argumento
    reference_alleles_parser.add_argument('-gapopen', required=False,
                                    help = 'Gap open penalty in BLAST searches. ',
                                    default = 1) ### c/m: introduciendo gapopen como argumento
    reference_alleles_parser.add_argument('-gapextend', required=False,
                                    help = 'Gap extension penalty in BLAST searches. ',
                                    default = 1) ### c/m: introduciendo gapextend como argumento       
    reference_alleles_parser.add_argument('-num_threads', required=False,
                                    help = 'num_threads in BLAST searches. ',
                                    default = 1) ### c/m: introduciendo num_threads como argumento   

    return parser.parse_args()


def processing_evaluate_schema (arguments) :
    print ('evaluate_schema')

    return True

def processing_compare_schema (arguments) :
    print ('compare_schema')
    return True

def processing_create_schema (arguments) :
    print ('create_schema')
    create_schema.processing_create_schema(arguments)
    return True

if __name__ == '__main__' :
    version = 'taranis version 0.2.2'
    if len(sys.argv) == 1 :
        print( 'Mandatory parameters are missing to execute the program. \n ' ,'Usage: "taranis.py  -help " for more information \n')
        exit (0)
    
    arguments = check_arg(sys.argv[1:])

    if arguments.chosen_action == 'allele_calling' :
        result = allele_calling.processing_allele_calling(arguments)
    elif arguments.chosen_action == 'evaluate_schema':
        result = analyze_schema.processing_evaluate_schema(arguments)
    elif arguments.chosen_action == 'compare_schema' :
        result = processing_compare_schema(arguments)
    elif arguments.chosen_action == 'create_schema' :
        result = processing_create_schema(arguments)
    elif arguments.chosen_action == 'reference_alleles' : ### añadiendo script para obtener alelos de referencia
        result = reference_alleles.processing_reference_alleles(arguments)
    else:
        print('not allow')
        result = 'Error'
    '''
    if 'Error' in result :
        print('Exiting the code with errors. ')
        print('Check the log for more information')
    else:
    '''
    print('completed')


