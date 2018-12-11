#!/usr/bin/env python3
#import analyze_schema
import sys
import argparse
#import allele_calling
#import analyze_schema
import create_schema
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
    
    parser = argparse.ArgumentParser(prog = 'tara.py',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description= 'Taranis is a set of utilities related to cgMSLST ') 

    parser.add_argument('--version', action='version', version='%(prog)s 0.3.5')
    
    subparser = parser.add_subparsers(help = 'allele calling /interactive/schema '
                                      + 'are the available actions to execute taranis',
                                      dest = 'chosen_action')
    
    ### Input parameters for allele calling option
    allele_calling_parser = subparser.add_parser('allele_calling',
                                    help = 'Allele calle way to downloads the  schema locus')
    allele_calling_parser.add_argument('-coregenedir', required= True,
                                    help = 'Directory where the core gene files are located ')
    allele_calling_parser.add_argument('-inputdir', required= True,
                                    help ='Directory where are located the sample fasta files')
    allele_calling_parser.add_argument('-outputdir', required= True,
                                    help = 'Directory where the result files will be stored')
    allele_calling_parser.add_argument('-cpus', required= False,
                                    help = 'Number of CPUS to be used in the program. Default is 1.',
                                    default = 1)
    allele_calling_parser.add_argument('-updateschema' , required=False,
                                    help = 'Create a new schema with the new locus found. Default is True.',
                                    default = True)
    allele_calling_parser.add_argument('-percentlength', required=False,
                                    help = 'Allowed length percentage to be considered as ASM or ALM. '
                                    + 'Outside of this limit it is considered as LNF Default is 20.',
                                    default = 20)
    
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
    version = 'tara version 0.2.1'
    if len(sys.argv) == 1 :
        print( 'Mandatory parameters are missing to execute the program. \n ' ,'Usage: "tara.py  -help " for more information \n')
        exit (0)
    arguments = check_arg(sys.argv[1:])

    if arguments.chosen_action == 'allele_calling' :
        allele_calling.processing_allele_calling(arguments)
    elif arguments.chosen_action == 'evaluate_schema':
        analyze_schema.processing_evaluate_schema(arguments)
    elif arguments.chosen_action == 'compare_schema' :
        processing_compare_schema(arguments)
    elif arguments.chosen_action == 'create_schema' :
        processing_create_schema(arguments)
    print('completed') 
   
