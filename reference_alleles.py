#!/usr/bin/env python3

from datetime import datetime
from Bio.Blast.Applications import NcbiblastnCommandline
import numpy as np
from utils.taranis_utils import *


### poner como argumento ID, reward etc de BLAST (lo suyo es que sea igual para este paso y para allele calling)
#### cuando meta id blast, etc, como argumento, cambiar los params de get_reference_allele para introducir esta info como param
### arreglar para que se borre el directoro mash al finalizar, ya que ahora mismo borra el interior pero se mantiene la carpeta mash vacía


def get_reference_allele(locus_quality, fasta_file, store_dir, evalue, perc_identity, reward, penalty, gapopen, gapextend, num_threads, logger): ### MODIFICANDO/CAMBIANDO: Función reference_allele añadida para encontrar el alelo de referencia
                                                                                                ### PROBANDO obtención de alelo de referencia con función distance de paquete Levenshtein
                                                                                                ### logger?
    ### logger.info('search reference allele')
    
    # alleles_seqs_dict --> diccionario generado previamente key:id alelo - value:seq alelo
    # fasta_file --> Es el path entero del locus en la carpeta de core_gene que se le mete en la terminal

    ### Me llevo la generación del directorio de alelos de referencia al inicio del código para poder así generar en reference_alleles el directorio mash donde se va a guardar el split del multifasta para cada locus y el sketch de cada locus
#    reference_allele_directory = 'reference_alleles'
 
#    full_path_reference_allele = os.path.join(store_dir, reference_allele_directory)

#    if not os.path.exists(full_path_reference_allele):
#        try:
#            os.makedirs(full_path_reference_allele)
#            logger.info('Directory %s has been created', full_path_reference_allele)
#        except:
#            print ('Cannot create the directory ', full_path_reference_allele)
#            logger.info('Directory %s cannot be created', full_path_reference_allele)
#            exit (0)
 
    ### Aquí metería una función para generar los archivos fasta individuales a partir del multifasta, generando previamente el directorio mash dentro del directorio reference_alleles
    ### y dentro del directorio mash crearía un directorio para el locus en cuestión donde se guardarían los fasta individuales
    ### En esa función metería un subprocess de mash sketch para generar el sketch con todos los fastas individuales y haría otro subproces de mash dist para obtener las distancias y las 
    ### guardaría en forma de matriz que va a sacar como output dentro de la función get_reference_allele y no sé si también la lista con todos los archivos que hay en la carpeta que
    ### debería haberse creado para 

    ### Creacion del directorio mash
#    mash_directory = 'mash'
#    full_path_mash = os.path.join(full_path_reference_allele, mash_directory)
#    if not os.path.exists(full_path_mash):
#        try:
#            os.makedirs(full_path_mash)
#            logger.info('Directory %s has been created', full_path_mash)
#        except:
#            print ('Cannot create the directory ', full_path_mash)
#            logger.info('Directory %s cannot be created', full_path_mash)
#            exit (0)

    ### Creacion del directory mash para cada locus. Si cada vez que se obtiene la referencia para un locuse se borra el directorio mash igual no haría falta generar un directorio
    ### dentro de mash para cada locus

  #  print("fasta_file (path del locus en cuestión): ", fasta_file, '\n')
    
    f_name = os.path.basename(fasta_file).split('.')
    full_path_reference_allele = os.path.join(store_dir, 'reference_alleles')
    full_path_locus_mash = os.path.join(full_path_reference_allele, 'mash', f_name[0])
    
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

    alleles_in_locus_number = 0 # sacando numero de alelos en locus para establecer params de blast
    alleles_in_locus = []     # añadiendo lista alelos en locus para comparar con resultado de ids de alelos tras BLAST alelos locus VS alelo de referencia
    for record in list(SeqIO.parse(fasta_file, "fasta")):
        alleles_in_locus.append(str(record.id))     # añadiendo lista alelos en locus para comparar con resultado de ids de alelos tras BLAST alelos locus VS alelo de referencia
        split_fasta_path = os.path.join(full_path_locus_mash, str(record.id) + ".fasta")
        alleles_in_locus_number += 1
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

 #   print("comp_dist_list: ", comp_dist_list, '\n')
    ### Obteniendo array de distancias y array de means
    comp_dist_list_per_allele = []
    alleles_number = len(split_multifasta_files_list)

    for index_distance in range(0, len(comp_dist_list), alleles_number):  
    
        dist_per_allele = comp_dist_list[index_distance : index_distance + alleles_number] 
        comp_dist_list_per_allele.append(dist_per_allele)

    comp_dist_arr_per_allele = np.asarray(comp_dist_list_per_allele)
    allele_mean_distance = np.mean(comp_dist_arr_per_allele, 0)
    
    #print("allele_mean_distance: ", allele_mean_distance)

 #   print("allele_mean_distance (matriz de distancias medias de un alelo al resto): ", allele_mean_distance, '\n')
    min_mean = max(allele_mean_distance) # inicializando con la mayor mean
    ref_allele_id = str()

    for mean_index in range(len(allele_mean_distance)):
  #      print("Ha entrado al bucle para obtener alelo de ref", '\n')
        if allele_mean_distance[mean_index] <= min_mean: # añado = para casos en los que las distancias sean todas iguales como en el caso del locus AEJV01_03887.fasta de Shigella (que solo tiene 2 alelos), para que cumplan la condición y se pueda obtener el path y el id del alelo y con ello el alelo de referencia
   #         print("Ha entrado al if para obtener path e id del alelo", '\n')
            allele_path = split_multifasta_files_list[mean_index]
            allele_id = os.path.basename(split_multifasta_files_list[mean_index]).split('.')[0]
            #if locus_alleles_ids[mean_index] in locus_quality['good_quality']:
            if locus_quality[allele_id] == 'good_quality':
                min_mean = allele_mean_distance[mean_index]
                ref_allele_id = allele_id

   # print("ref_allele_id: ", ref_allele_id, '\n')
  #  print("min_mean: ", min_mean, '\n')

    # comprobando que el alelo de referencia obtenido representa a todos los alelos del locus. Los alelos no representados se incluyen en el fasta junto al alelo de referencia general

    # creando base de datos BLAST de todos los alelos del locus
    # IGUAL DEBERÍA CREAR LAS BASES DE DATOS DE TODOS LOS LOCUS EN EL MISMO LUGAR DEL CÓDIGO QUE DONDE SE GENERAN LAS BASES DE DATOS DE LAS MUESTRAS AL INICIO?

    #create_blastdb(core_file, db_name, 'nucl', logger) 

    db_name = os.path.join(full_path_reference_allele, 'locus_blastdb') ### path al directorio donde se van a guardar las bases de datos de las secuencias encontradas tras el primer blast a partir del alelo de referencia para esta muestra
    ##create local blast db for sample fasta file
    if not create_blastdb(fasta_file, db_name, 'nucl', logger):
        print('Error when creating the blastdb for locus %s. Check log file for more information. \n ', f_name[0])
        return False

    #def create_blastdb (file_name, db_name,db_type, logger ):
        #f_name = os.path.basename(file_name).split('.') ----> path al fasta con las secuencias obtenidas tras el blast con el alelo de referencia
        #db_dir = os.path.join(db_name,f_name[0]) ----> path a directorio donde se guardan las bases de datos, pero añado el nombre de la muestra donde se tiene que crear no?
        #output_blast_dir = os.path.join(db_dir, f_name[0]) --->

    # path a la base de datos de BLAST creada para las secuencias obtenidas tras BLAST con alelo de referencia
    locus_db_name = os.path.join(db_name, f_name[0], f_name[0])


    # BLAST de todos los alelos frente alelo de referencia
    blast_parameters = '"6 , qseqid , sseqid , pident ,  qlen , length , mismatch , gapopen , evalue , bitscore , sstart , send , qstart , qend , sseq , qseq"'
    cline = NcbiblastnCommandline(db=locus_db_name, evalue=evalue, perc_identity=perc_identity, reward=reward, penalty=penalty, gapopen=gapopen, gapextend=gapextend, outfmt=blast_parameters, max_target_seqs=alleles_in_locus_number, max_hsps=alleles_in_locus_number, num_threads=num_threads, query=allele_path)

    out, err = cline()
    out_lines = out.splitlines()

    alleles_in_blast = [] ### lista de ids de alelos para los que se ha obtenido match con el alelo de referencia

    for line in out_lines:
        values = line.split('\t') ### tener en cuenta que estos values son del segundo BLAST, es decir, no tengo aquí info sobre posición en contig, ni contig como tal si guardo la secuencia como la he guardado (id númeroalelo_contig), pero tengo el pident, s_length, qlen, sseq, qseq, bitscore, etc
        alleles_in_blast.append(values[1])

  #  print("len(alleles_in_blast): ", len(alleles_in_blast), '\n')

    ids_intersect = list(set(alleles_in_locus) - set(alleles_in_blast)) ### alelos del locus no representados por el alelo de referencia

    #all_reference_alleles = ids_intersect # uniendo en una unica lista el alelo de referencia (primer elemento) y el resto de alelos no representados por el alelo de referencia
    #all_reference_alleles.insert(0, ref_allele_id) 

    ids_intersect.insert(0, ref_allele_id) # uniendo en una unica lista el alelo de referencia (primer elemento) y el resto de alelos no representados por el alelo de referencia

  #  print("ids_intersect: ", ids_intersect, '\n')
  #  print("len(ids_intersect): ", len(ids_intersect), '\n')

    # parseando fasta locus para generar fasta de referencia con los ids y alelos correspondientes
    reference_file_path = os.path.join(full_path_reference_allele, os.path.basename(fasta_file)) #### path donde se va a guardar el fasta con los alelos ref y nombre archivo
    with open (reference_file_path, 'w') as out_fh:
        for record in list(SeqIO.parse(fasta_file, "fasta")):
            if record.id in ids_intersect:
                out_fh.write ('>' + str(record.id) + '\n' + str(record.seq))

    # Codigo anterior a introducir BLAST para comprobacion de alelo de referencia
    
#    ### Moviendo fasta que contiene el alelo de referencia al directorio reference
#    #for file_path in split_multifasta_files_list:
#     #   f_name = os.path.basename(file_path).split('.')[0]
#      #  if f_name == ref_allele_id:
#       #     shutil.move(file_path, full_path_reference_allele)

#    # Probando mover fasta metiendo path obtenido generándolo a partir del path del locus en el directorio mash y el id del alelo de referencia obtenido
#    reference_file_path = os.path.join(full_path_locus_mash, ref_allele_id + ".fasta")
#    new_reference_file_path = os.path.join(full_path_reference_allele, os.path.basename(fasta_file))
#    #shutil.move(reference_file_path, full_path_reference_allele)
#    shutil.move(reference_file_path, new_reference_file_path) ### tengo que renombrar el archivo del alelo de referencia poniéndole el nombre del locus
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

    return reference_file_path


def processing_reference_alleles (arguments) :
    '''
    Description:
        This is the main function for getting reference alleles.
        With the support of additional functions it will obtain reference alleles for each locus in the schema.
    Input:
        arguments   # Input arguments given on command line 
    Functions:
        
    Variables:  ?????
        run_metric_processed # True or False if there are some rows in
                            StatsRunSummary for this run
    Return: ?????
        experiment_name if the run is updated. Empty if not
    '''
    #logger = logging.getLogger(__name__)
    #logger.debug ('Starting function check_run_metrics_processed')
    start_time = datetime.now()
    print('Start the execution at :', start_time )
    # open log file
    logger = open_log ('taranis_wgMLST.log')
    print('Checking the pre-requisites./n') #--------> Hace falta checkear los prerrequisitos para este script
    # check additional programs are installed in your system
    pre_requisite_list = [['blastp', '2.5'], ['makeblastdb' , '2.5'], ['mash', '1.1']]
    if not check_prerequisites (pre_requisite_list, logger):
        print ('your system does not fulfill the pre-requistes to run the script ')
        exit(0)

    ####################################################
    # Check that given directories contain fasta files #
    ####################################################
    print('Validating schema fasta files in ' , arguments.coregenedir , '\n')
    core_gene_files_list = get_fasta_file_list(arguments.coregenedir, logger)
    if not core_gene_files_list :
        print ('There are not valid fasta files in ',  arguments.coregenedir , ' directory. Check log file for more information ')
        exit(0)

    #####################################################
    # Obtain reference alleles for each locus in schema #
    #####################################################
    logger.info('getting reference alleles for each locus in schema')
    for fasta_file in core_gene_files_list:

        # obteniendo calidad del locus
        locus_quality = check_core_gene_quality(fasta_file, logger) ### cambiando/modificando: obteniendo diccionario de calidad de cada locus introduciendo path al fasta del locus en cuestión en lugar del diccionario

        #### cuando meta id blast, etc, como argumento, cambiar los params de get_reference_allele para introducir esta info como param
        # llamada a get_reference_allele y ponerlo rollo como se llama a la función de allele_Call_nucleotides y meter dento de get_reference_alleles logger?
        if not get_reference_allele(locus_quality, fasta_file, arguments.outputdir, float(arguments.evalue), float(arguments.perc_identity), int(arguments.reward), int(arguments.penalty), int(arguments.gapopen), int(arguments.gapextend), int(arguments.num_threads), logger): ## si meto todo esto en una fucnión rollo preprocessing allele calling, arguments.outputdir pasará a llamarse como el param en cuesetión en la función 
            print('There is an error while processing reference alleles. Check the log file to get more information \n')
            logger.info('Deleting the directory to clean up the temporary files created')
            shutil.rmtree(os.path.join(arguments.outputdir))
            exit(0)

    end_time = datetime.now()
    print('completed execution at :', end_time )
    return True