import logging

import click
import concurrent.futures
import glob
import os
import rich.console
import rich.logging
import rich.traceback
import sys
import time

import taranis.prediction
import taranis.utils
import taranis.analyze_schema
import taranis.reference_alleles
import taranis.allele_calling

log = logging.getLogger()

# Set up rich stderr console
stderr = rich.console.Console(
    stderr=True, force_terminal=taranis.utils.rich_force_colors()
)


def run_taranis():
    # Set up the rich traceback
    rich.traceback.install(console=stderr, width=200, word_wrap=True, extra_lines=1)

    # Print taranis header
    # stderr.print("\n[green]{},--.[grey39]/[green],-.".format(" " * 42), highlight=False)
    stderr.print(
        "[blue]                ______           ___                     ___    ",
        highlight=False,
    )
    stderr.print(
        "[blue]   \    |-[grey39]-|  [blue] |__--__|   /\    |   \    /\    |\   | | |  ",
        highlight=False,
    )
    stderr.print(
        "[blue]    \   \  [grey39]/ [blue]     ||     /  \   |__ /   /  \   | \  | | |___   ",
        highlight=False,
    )
    stderr.print(
        "[blue]    /  [grey39] / [blue] \      ||    /____\  |  \   /____\  |  \ | |     |",
        highlight=False,
    )
    stderr.print(
        "[blue]   /   [grey39] |-[blue]-|      ||   /      \ |   \ /      \ |   \| |  ___| ",
        highlight=False,
    )

    # stderr.print("[green]                                          `._,._,'\n", highlight=False)
    __version__ = "3.0.0"
    stderr.print(
        "\n" "[grey39]    Taranis version {}".format(__version__), highlight=False
    )

    # Lanch the click cli
    taranis_cli()


# Customise the order of subcommands for --help
class CustomHelpOrder(click.Group):
    def __init__(self, *args, **kwargs):
        self.help_priorities = {}
        super(CustomHelpOrder, self).__init__(*args, **kwargs)

    def get_help(self, ctx):
        self.list_commands = self.list_commands_for_help
        return super(CustomHelpOrder, self).get_help(ctx)

    def list_commands_for_help(self, ctx):
        """reorder the list of commands when listing the help"""
        commands = super(CustomHelpOrder, self).list_commands(ctx)
        return (
            c[1]
            for c in sorted(
                (self.help_priorities.get(command, 1000), command)
                for command in commands
            )
        )

    def command(self, *args, **kwargs):
        """Behaves the same as `click.Group.command()` except capture
        a priority for listing command names in help.
        """
        help_priority = kwargs.pop("help_priority", 1000)
        help_priorities = self.help_priorities

        def decorator(f):
            cmd = super(CustomHelpOrder, self).command(*args, **kwargs)(f)
            help_priorities[cmd.name] = help_priority
            return cmd

        return decorator


@click.group(cls=CustomHelpOrder)
@click.version_option(taranis.__version__)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    default=False,
    help="Print verbose output to the console.",
)
@click.option(
    "-l", "--log-file", help="Save a verbose log to a file.", metavar="filename"
)
def taranis_cli(verbose, log_file):
    # Set the base logger to output DEBUG
    log.setLevel(logging.DEBUG)

    # Set up logs to a file if we asked for one
    if log_file:
        log_fh = logging.FileHandler(log_file, encoding="utf-8")
        log_fh.setLevel(logging.DEBUG)
        log_fh.setFormatter(
            logging.Formatter(
                "[%(asctime)s] %(name)-20s [%(levelname)-7s]  %(message)s"
            )
        )
        log.addHandler(log_fh)


# Analyze schema
# taranis analyze-schema -i /media/lchapado/Reference_data/proyectos_isciii/taranis/documentos_antiguos/pasteur_schema -o /media/lchapado/Reference_data/proyectos_isciii/taranis/test/analyze_schema
# testing data for analyze schema
# taranis analyze-schema -i /media/lchapado/Reference_data/proyectos_isciii/taranis/taranis_testing_data/listeria_testing_schema -o /media/lchapado/Reference_data/proyectos_isciii/taranis/test/analyze_schema


@taranis_cli.command(help_priority=1)
@click.option(
    "-i",
    "--inputdir",
    required=True,
    multiple=False,
    type=click.Path(),
    help="Directory where the schema with the core gene files are located. ",
)
@click.option(
    "-o",
    "--output",
    required=True,
    multiple=False,
    type=click.Path(),
    help="Output folder to save analyze schema",
)
@click.option(
    "--remove-subset/--no-remove-subset",
    required=False,
    default=False,
    help="Remove allele subsequences from the schema.",
)
@click.option(
    "--remove-duplicated/--no-remove-duplicated",
    required=False,
    default=False,
    help="Remove duplicated subsequences from the schema.",
)
@click.option(
    "--remove-no-cds/--no-remove-no-cds",
    required=False,
    default=False,
    help="Remove no CDS alleles from the schema.",
)
@click.option(
    "--output-allele-annot/--no-output-allele-annot",
    required=False,
    default=True,
    help="get extension annotation for all alleles in locus",
)
@click.option(
    "--genus",
    required=False,
    default="Genus",
    help="Genus name for Prokka schema genes annotation. Default is Genus.",
)
@click.option(
    "--species",
    required=False,
    default="species",
    help="Species name for Prokka schema genes annotation. Default is species",
)
@click.option(
    "--usegenus",
    required=False,
    default="Genus",
    help="Use genus-specific BLAST databases for Prokka schema genes annotation (needs --genus). Default is False.",
)
@click.option(
    "--cpus",
    required=False,
    multiple=False,
    type=int,
    default=1,
    help="Number of cpus used for execution",
)
def analyze_schema(
    inputdir: str,
    output: str,
    remove_subset: bool,
    remove_duplicated: bool,
    remove_no_cds: bool,
    output_allele_annot: bool,
    genus: str,
    species: str,
    usegenus: str,
    cpus: int,
):
    schema_files = taranis.utils.get_files_in_folder(inputdir, "fasta")

    """
    schema_analyze = []
    for schema_file in schema_files:
        schema_obj = taranis.analyze_schema.AnalyzeSchema(schema_file, output, remove_subset, remove_duplicated, remove_no_cds, genus, species, usegenus)
        schema_analyze.append(schema_obj.analyze_allele_in_schema())
    import pdb; pdb.set_trace()
    _ = taranis.analyze_schema.collect_statistics(schema_analyze, output, output_allele_annot)
    sys.exit(0)
    # for schema_file in schema_files:
    """
    results = []
    start = time.perf_counter()
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpus) as executor:
        futures = [
            executor.submit(
                taranis.analyze_schema.parallel_execution,
                schema_file,
                output,
                remove_subset,
                remove_duplicated,
                remove_no_cds,
                genus,
                species,
                usegenus,
            )
            for schema_file in schema_files
        ]
        # Collect results as they complete
        for future in concurrent.futures.as_completed(futures):
            results.append(future.result())
    _ = taranis.analyze_schema.collect_statistics(results, output, output_allele_annot)
    finish = time.perf_counter()
    print(f"Schema analyze finish in {round((finish-start)/60, 2)} minutes")


# Reference alleles
@taranis_cli.command(help_priority=2)
@click.option(
    "-s",
    "--schema",
    required=True,
    multiple=False,
    type=click.Path(),
    help="Directory where the schema with the core gene files are located. ",
)
@click.option(
    "-o",
    "--output",
    required=True,
    multiple=False,
    type=click.Path(),
    help="Output folder to save reference alleles",
)
def reference_alleles(
    schema,
    output,
):
    # taranis reference-alleles -s ../../documentos_antiguos/datos_prueba/schema_1_locus/ -o ../../new_taranis_result_code
    # taranis reference-alleles -s /media/lchapado/Reference_data/proyectos_isciii/taranis/taranis_testing_data/listeria_testing_schema/ -o /media/lchapado/Reference_data/proyectos_isciii/taranis/test/reference_alleles
    schema_files = taranis.utils.get_files_in_folder(schema, "fasta")

    # Check if output folder exists
    if taranis.utils.folder_exists(output):
        q_question = (
            "Folder "
            + output
            + " already exists. Files will be overwritten. Do you want to continue?"
        )
        if "no" in taranis.utils.query_user_yes_no(q_question, "no"):
            log.info("Aborting code by user request")
            stderr.print("[red] Exiting code. ")
            sys.exit(1)
    else:
        try:
            os.makedirs(output)
        except OSError as e:
            log.info("Unable to create folder at %s with error %s", output, e)
            stderr.print("[red] ERROR. Unable to create folder  " + output)
            sys.exit(1)
    """Create the reference alleles from the schema """
    for f_file in schema_files:
        ref_alleles = taranis.reference_alleles.ReferenceAlleles(f_file, output)
    _ = ref_alleles.create_ref_alleles()


# Allele calling
#  taranis -l ../../test/taranis.log  allele-calling -s ../../documentos_antiguos/datos_prueba/schema_test/ -r ../../documentos_antiguos/datos_prueba/reference_alleles/ -g ../../taranis_data/listeria_genoma_referencia/listeria.fasta -a ../../taranis_data/listeria_sampled/RA-L2073_R1.fasta -o ../../test/
# taranis allele-calling -s ../../documentos_antiguos/datos_prueba/schema_test/ -r ../../documentos_antiguos/datos_prueba/reference_alleles/ -g ../../taranis_data/listeria_genoma_referencia/listeria.fasta -a ../../taranis_data/listeria_sampled/RA-L2073_R1.fasta -o ../../test/
# taranis allele-calling -s ../../documentos_antiguos/datos_prueba/schema_test/ -r ../../documentos_antiguos/datos_prueba/reference_alleles/ -g ../../taranis_data/listeria_genoma_referencia/listeria.fasta -a ../../taranis_data/muestras_listeria_servicio_fasta/3789/assembly.fasta -o ../../test/


@taranis_cli.command(help_priority=3)
@click.option(
    "-s",
    "--schema",
    required=True,
    multiple=False,
    type=click.Path(),
    help="Directory where the schema with the core gene files are located. ",
)
@click.option(
    "-r",
    "--reference",
    required=True,
    multiple=False,
    type=click.Path(),
    help="Directory where the schema reference allele files are located. ",
)
@click.option(
    "-g",
    "--genome",
    required=True,
    multiple=False,
    type=click.Path(),
    help="Genome reference file",
)
@click.option(
    "-a",
    "--sample",
    required=True,
    multiple=False,
    type=click.Path(),
    help="Sample location file in fasta format. ",
)
@click.option(
    "-o",
    "--output",
    required=True,
    multiple=False,
    type=click.Path(),
    help="Output folder to save reference alleles",
)
def allele_calling(
    schema,
    reference,
    genome,
    sample,
    output,
):
    folder_to_check = [schema, reference]
    for folder in folder_to_check:
        if not taranis.utils.folder_exists(folder):
            log.error("folder %s does not exists", folder)
            stderr.print("[red] Folder does not exist. " + folder + "!")
            sys.exit(1)
    if not taranis.utils.file_exists(sample):
        log.error("file %s does not exists", sample)
        stderr.print("[red] File does not exist. " + sample + "!")
        sys.exit(1)
    schema_files = taranis.utils.get_files_in_folder(schema, "fasta")
    if len(schema_files) == 0:
        log.error("Schema folder %s does not have any fasta file", schema)
        stderr.print("[red] Schema folder does not have any fasta file")
        sys.exit(1)

    # Check if output folder exists
    if taranis.utils.folder_exists(output):
        q_question = (
            "Folder "
            + output
            + " already exists. Files will be overwritten. Do you want to continue?"
        )
        if "no" in taranis.utils.query_user_yes_no(q_question, "no"):
            log.info("Aborting code by user request")
            stderr.print("[red] Exiting code. ")
            sys.exit(1)
    else:
        try:
            os.makedirs(output)
        except OSError as e:
            log.info("Unable to create folder at %s because %s", output, e)
            stderr.print("[red] ERROR. Unable to create folder  " + output)
            sys.exit(1)
    # Filter fasta files from reference folder
    ref_alleles = glob.glob(os.path.join(reference, "*.fasta"))
    # Create predictions
    pred_out = os.path.join(output, "prediction")
    pred_sample = taranis.prediction.Prediction(genome, sample, pred_out)
    pred_sample.training()
    pred_sample.prediction()

    """Analyze the sample file against schema to identify outbreakers
    """
    sample_allele = taranis.allele_calling.AlleleCalling(
        pred_sample, sample, schema, ref_alleles, output
    )
    sample_allele.analyze_sample()