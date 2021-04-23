import sys
import argparse
import allele_calling
import analyze_schema
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
        analyze_schema_parser  # It is used for the schema evaluation input parameters
        compare_schema_parser   # It is used for schema comparison input parameters
        create_schema_parser    # It is used for create an schema input parameters
    Return:
        parser.parse_args()     # The variable contains the valid parameters
    '''
    parser = argparse.ArgumentParser(prog = 'taranis.py',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description= 'Taranis is a set of utilities related to cgMSLST ')

    parser.add_argument('--version', action='version', version='%(prog)s 0.3.5')

    subparser = parser.add_subparsers(help = 'analyze_schema, reference_alleles, allele_calling'
                                      + 'are the available actions to execute taranis',
                                      dest = 'chosen_action')


    ### Input parameters for analyze schema option
    analyze_schema_parser = subparser.add_parser('analyze_schema',
                                    help = 'Analyze the schema.')
    analyze_schema_parser.add_argument('-inputdir',required = True,
                                    help = 'Directory where are the schema files.')
    analyze_schema_parser.add_argument('-outputdir', required = True,
                                    help = 'Directory where the result files will be stored.')
    analyze_schema_parser.add_argument('-removesubsets', required = False,
                                    help = 'Remove allele subsequences from the schema.'
                                    + 'True: Remove subsets.'
                                    + 'False: Do not remove subsets.'
                                    + 'Default is False.',
                                    default = False)    
    analyze_schema_parser.add_argument('-removeduplicates', required = False,
                                    help = 'Remove duplicated alleles from the schema.'
                                    + 'True: Remove duplicates.'
                                    + 'False: Do not remove duplicates.'
                                    + 'Default is False.',
                                    default = False) 
    analyze_schema_parser.add_argument('-removenocds', required = False,
                                    help = 'Remove no CDS alleles from the schema.'
                                    + 'True: Remove no CDS alleles.'
                                    + 'False: Do not remove no CDS alleles.'
                                    + 'Default is False.',
                                    default = False) 
    analyze_schema_parser.add_argument('-newschema', required = False,
                                    help = 'Filter a copy of the core genes schema preserving the analysis core genes schema.' 
                                            #Create an analysis core genes schema copy for filtering alleles when this option is selected.'
                                    + 'True: Create a copy of the core genes schema for filtering.'
                                    + 'False: Do not create a copy of the core genes schema for filtering.'
                                    + 'Default is False.',
                                    default = False) 
    analyze_schema_parser.add_argument('-genus' , required = False,
                                    help = 'Genus name for Prokka schema genes annotation. Defualt is Genus. ',
                                    default = 'Genus') 
    analyze_schema_parser.add_argument('-species' , required = False,
                                    help = 'Species name for Prokka schema genes annotation. Defualt is species. ',
                                    default = 'species') 
    analyze_schema_parser.add_argument('-usegenus' , required = False,
                                    help = 'Use genus-specific BLAST databases for Prokka schema genes annotation (needs --genus). Defualt is False. ',
                                    default = 'False')
    analyze_schema_parser.add_argument('-cpus', required = False,
                                    help = 'Number of CPUS to be used in the program. Default is 1.',
                                    default = 1)                                


    ### Input parameters for reference alleles options  
    reference_alleles_parser = subparser.add_parser('reference_alleles', help = 'Obtain reference allele(s) for each locus.')
    reference_alleles_parser.add_argument('-coregenedir', required = True,
                                    help = 'Directory where the core gene files are located. ')
    reference_alleles_parser.add_argument('-outputdir', required = True,
                                    help = 'Directory where the result files will be stored. ')
    reference_alleles_parser.add_argument('-evalue', required = False,
                                    help = 'E-value in BLAST searches. Default is 0.001.',
                                    default = 0.001) 
    reference_alleles_parser.add_argument('-perc_identity', required = False,
                                    help = 'Identity percent in BLAST searches. Default is 90 %. ',
                                    default = 90) 
    reference_alleles_parser.add_argument('-reward', required = False,
                                    help = 'Match reward in BLAST searches. Default is 1. ',
                                    default = 1)     
    reference_alleles_parser.add_argument('-penalty', required = False,
                                    help = 'Mismatch penalty in BLAST searches. Default is -2. ',
                                    default = -2) 
    reference_alleles_parser.add_argument('-gapopen', required = False,
                                    help = 'Gap open penalty in BLAST searches. Default is 1. ',
                                    default = 1) 
    reference_alleles_parser.add_argument('-gapextend', required = False,
                                    help = 'Gap extension penalty in BLAST searches. Default is 1. ',
                                    default = 1)        
    reference_alleles_parser.add_argument('-num_threads', required = False,
                                    help = 'num_threads in BLAST searches. Default is 1. ',
                                    default = 1)    
    reference_alleles_parser.add_argument('-cpus', required = False,
                                    help = 'Number of CPUS to be used in the program. Default is 1.',
                                    default = 1)


    ### Input parameters for allele calling option
    allele_calling_parser = subparser.add_parser('allele_calling',
                                    help = 'Gene by gene allele calling') 
    allele_calling_parser.add_argument('-coregenedir', required = True,
                                    help = 'Directory where the core gene files are located ')
    allele_calling_parser.add_argument('-refalleles', required = True,
                                    help = 'Directory where the core gene references files are located ') 
    allele_calling_parser.add_argument('-inputdir', required = True,
                                    help ='Directory where are located the sample fasta files')
    allele_calling_parser.add_argument('-refgenome', required = True,
                                    help = 'Reference genome file for genes prediction') 
    allele_calling_parser.add_argument('-outputdir', required = True,
                                    help = 'Directory where the result files will be stored')
    allele_calling_parser.add_argument('-percentlength', required = False,
                                    help = 'Allowed length percentage to be considered as INF. '
                                    + 'Outside of this limit it is considered as ASM or ALM. Default is SD.',
                                    default = 'SD') 
    allele_calling_parser.add_argument('-coverage', required = False,
                                    help = 'Coverage threshold to exclude found sequences. '
                                    + 'Outside of this limit it is considered LNF. Default is 50 %.',
                                    default = 50) 
    allele_calling_parser.add_argument('-evalue', required = False,
                                    help = 'E-value in BLAST searches. Default is 0.001. ',
                                    default = 0.001) 
    allele_calling_parser.add_argument('-perc_identity_ref', required = False,
                                    help = 'Identity percentage in BLAST searches using reference alleles for each locus detection in samples. Default is 90 %. ',
                                    default = 90) 
    allele_calling_parser.add_argument('-perc_identity_loc', required = False,
                                    help = 'Identity percentage in BLAST searches using all alleles in each locus for allele identification in samples. Default is 90 %. ',
                                    default = 90) 
    allele_calling_parser.add_argument('-reward', required = False,
                                    help = 'Match reward in BLAST searches. Default is 1. ',
                                    default = 1)     
    allele_calling_parser.add_argument('-penalty', required = False,
                                    help = 'Mismatch penalty in BLAST searches. Default is -2. ',
                                    default = -2) 
    allele_calling_parser.add_argument('-gapopen', required = False,
                                    help = 'Gap open penalty in BLAST searches. Default is 1. ',
                                    default = 1) 
    allele_calling_parser.add_argument('-gapextend', required = False,
                                    help = 'Gap extension penalty in BLAST searches. Default is 1. ',
                                    default = 1)        
    allele_calling_parser.add_argument('-max_target_seqs', required = False,
                                    help = 'max_target_seqs in BLAST searches. Default is 10. ', 
                                    default = 10)     
    allele_calling_parser.add_argument('-max_hsps', required = False,
                                    help = 'max_hsps in BLAST searches. Default is 10. ',
                                    default = 10) 
    allele_calling_parser.add_argument('-num_threads', required = False,
                                    help = 'num_threads in BLAST searches. Default is 1. ',
                                    default = 1)    
    allele_calling_parser.add_argument('-flankingnts' , required = False,
                                    help = 'Number of flanking nucleotides to add to each BLAST result obtained after locus detection in sample using reference allele for correct allele identification. Default is 100. ',
                                    default = 100)
    allele_calling_parser.add_argument('-updateschema' , required = False,
                                    help = 'Add INF alleles found for each locus to the core genes schema. ' 
                                    + 'True: Add INF alleles to the analysis core genes schema. ' 
                                    + 'New: Add INF alleles to a copy of the core genes schema preserving the analysis core genes schema. ' 
                                    + 'False: Do not update the core gene schema adding new INF alleles found. '
                                    + 'Default is True. ',
                                    default = True)
    allele_calling_parser.add_argument('-profile' , required = False,
                                    help = 'ST profile file based on core genes schema file to get ST for each sample. Default is empty and Taranis does not calculate samples ST. ',
                                    default = '') 
    allele_calling_parser.add_argument('-cpus', required = False,
                                    help = 'Number of CPUS to be used in the program. Default is 1.',
                                    default = 1)
    allele_calling_parser.add_argument('-genus' , required = False,
                                    help = 'Genus name for Prokka schema genes annotation. Defualt is Genus. ',
                                    default = 'Genus') 
    allele_calling_parser.add_argument('-species' , required = False,
                                    help = 'Species name for Prokka schema genes annotation. Defualt is species. ',
                                    default = 'species') 
    allele_calling_parser.add_argument('-usegenus' , required = False,
                                    help = 'Use genus-specific BLAST databases for Prokka schema genes annotation (needs --genus). Defualt is False. ',
                                    default = 'False')


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


#def processing_analyze_schema (arguments) :
 #   print ('analyze_schema')

  #  return True

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
    elif arguments.chosen_action == 'analyze_schema':
        result = analyze_schema.processing_analyze_schema(arguments)
    elif arguments.chosen_action == 'compare_schema' :
        result = processing_compare_schema(arguments)
    elif arguments.chosen_action == 'create_schema' :
        result = processing_create_schema(arguments)
    elif arguments.chosen_action == 'reference_alleles' :
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


