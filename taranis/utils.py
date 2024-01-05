#!/usr/bin/env python
"""
Common utility function used for relecov_tools package.
"""

import glob

import logging
import numpy as np
import questionary
import os
import plotly.graph_objects as go
import rich.console
import shutil
import subprocess


import sys

from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import pdb
log = logging.getLogger(__name__)

def rich_force_colors():
    """
    Check if any environment variables are set to force Rich to use coloured output
    """
    if (
        os.getenv("GITHUB_ACTIONS")
        or os.getenv("FORCE_COLOR")
        or os.getenv("PY_COLORS")
    ):
        return True
    return None
stderr = rich.console.Console(
    stderr=True,
    style="dim",
    highlight=False,
    force_terminal=rich_force_colors(),
)

START_CODON_FORWARD= ['ATG','ATA','ATT','GTG', 'TTG']
start_codon_reverse= ['CAT', 'TAT','AAT','CAC','CAA']

STOP_CODON_FORWARD = ['TAA', 'TAG','TGA']
stop_codon_reverse = ['TTA', 'CTA','TCA']

POSIBLE_BAD_QUALITY = ["Can not be converted to protein", "Start codon not found", "Stop codon not found", "Multiple stop codons found" ,"Duplicate allele", "Sub set allele"]

def check_sequence_order(allele_sequence):
    # check direction
    if allele_sequence[0:3] in START_CODON_FORWARD or allele_sequence[-3:] in STOP_CODON_FORWARD:
        return 'forward'
    if allele_sequence[-3:] in start_codon_reverse or allele_sequence[0:3] in stop_codon_reverse:
        return 'reverse'
    return "Error"

def create_annotation_files(fasta_file, annotation_dir, prefix, genus="Genus", species="species", usegenus=False):
    try:
        _ = subprocess.run (['prokka', fasta_file, '--force', '--outdir', annotation_dir, '--genus', genus, '--species', species, '--usegenus', str(usegenus), '--gcode', '11', '--prefix', prefix, '--quiet'])
    except Exception as e:
        log.error("Unable to run prokka. Error message: %s ", e )
        stderr.print("[red] Unable to run prokka. Given error; " + e)
        sys.exit(1)
    # Check that prokka store files in the requested folder
    # if prokka results are not found in the requested folder then move from the
    # running directory to the right one
    if not folder_exists(annotation_dir):
        try:
            shutil.move(prefix, annotation_dir)
        except Exception as e:
            log.error("Unable to move prokka result folder to %s ", e )
            stderr.print("[red] Unable to move result prokka folder. Error; " + e)
            sys.exit(1)
    return os.path.join(annotation_dir, prefix)


def create_new_folder(folder_name):
    try:
        os.makedirs(folder_name, exist_ok=True)
    except Exception as e:
        log.error("Folder %s can not be created %s", folder_name, e)
        stderr.print("[red] Folder does not have any file which match your request")
        sys.exit(1)
    return


def  create_graphic(out_folder, f_name, mode, x_data, y_data, labels, title ):
    fig = go.Figure()
    # pdb.set_trace()
    if mode == "lines":
        fig.add_trace(go.Scatter(x=x_data, y=y_data, mode=mode, name=title))
    elif mode == "pie":
        fig.add_trace(go.Pie(labels=labels, values=x_data))
        fig.update_layout(title_text= title)
    fig.write_image(os.path.join(out_folder, f_name))


def get_files_in_folder(folder, extension=None):
    """get the list of files, filtered by extension in the input folder. If
    extension is not set, then all files in folder are returned

    Args:
        folder (string): folder path
        extension (string, optional): extension for filtering. Defaults to None.
    
    Returns:
        list: list of files which match the condition
    """
    if not folder_exists(folder):
        log.error("Folder %s does not exists", folder)
        stderr.print("[red] Schema folder does not exist. " + folder + "!")
        sys.exit(1)
    if extension is None:
        extension = "*"
    folder_files = os.path.join(folder , "*." + extension)
    files_in_folder = glob.glob(folder_files)
    if len(files_in_folder) == 0:
        log.error("Folder %s does not have any file which the extension %s", folder, extension)
        stderr.print("[red] Folder does not have any file which match your request")
        sys.exit(1)
    return files_in_folder
    

def file_exists(file_to_check):
    """Checks if input file exists

    Args:
        file_to_check (string): file name  including path of the file

    Returns:
        boolean: True if exists
    """    
    if os.path.isfile(file_to_check):
        return True
    return False

def find_multiple_stop_codons(seq) :
    stop_codons = ['TAA', 'TAG','TGA']
    c_index = []
    for idx in range (0, len(seq) -2, 3) :
        c_seq = seq[idx : idx + 3]
        if c_seq in stop_codons :
            c_index.append(idx)
    if len(c_index) == 1:
        return False
    return True


def find_nearest_numpy_value(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def folder_exists(folder_to_check):
    """Checks if input folder exists

    Args:
        folder_to_check (string): folder name  including path

    Returns:
        boolean: True if exists
    """    
    if os.path.isdir(folder_to_check):
        return True
    return False

def prompt_text(msg):
    source = questionary.text(msg).unsafe_ask()
    return source

def query_user_yes_no(question, default):
    """Query the user to choose yes or no for the query question

    Args:
        question (string): Text message
        default (string): default option to be used: yes or no

    Returns:
        user select: True continue with code 
    """ 
    valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)
    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == "":
            return valid[default]
        elif choice in valid:
            if "y" in choice:
                return "yes"
            else:
                return "no"
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' (or 'y' or 'n').\n")

def read_annotation_file(ann_file, allele_name, only_first_line=True):

    """ example of annotation file
    locus_tag	ftype	length_bp	gene	EC_number	COG	product
    IEKBEMEO_00001	CDS	1344	yeeO_1		COG0534	putative FMN/FAD exporter YeeO
    IEKBEMEO_00002	CDS	1344	yeeO_2		COG0534	putative FMN/FAD exporter YeeO

    """
    ann_data = {}
    with open (ann_file, "r") as fh:
        lines = fh.readlines()
    heading = lines[0].split("\t")
    locus_tag_idx = heading.index("locus_tag")
    gene_idx = heading.index("gene")
    if only_first_line:
        first_line = lines[1].split("\t")
        ann_data[allele_name] = first_line[gene_idx] if first_line[gene_idx] != "" else "Not found by Prokka"
    else:
        # Return all annotation lines
        for line in lines[1:]:
            s_line = line.strip().split("\t")
            allele_key = allele_name + "_" + s_line[locus_tag_idx].split("_")[1]
            ann_data[allele_key] = s_line[gene_idx] if s_line[gene_idx] != "" else "Not found by Prokka"
    return ann_data



def read_fasta_file(fasta_file):
    return SeqIO.parse(fasta_file, "fasta")
    

def write_fasta_file(out_folder, seq_data, allele_name=None, f_name=None):
    try:
        os.makedirs(out_folder, exist_ok=True)
    except OSError as e:
        sys.exit(1)
    if isinstance(seq_data, dict):
        for key, seq in seq_data.items():
            if f_name is None:
                # use the fasta name as file name
                f_name = key + ".fasta"
            f_path_name = os.path.join(out_folder, f_name)
            with open (f_path_name, "w") as fo:
                fo.write(">" + key + "\n")
                fo.write(seq)
    else:
        if f_name is None:
            f_name = allele_name
        f_path_name = os.path.join(out_folder, f_name)
        with open (f_path_name, "w") as fo:
            fo.write(">" + allele_name + "\n")
            fo.write(seq_data)
    return f_name

def write_data_to_file(out_folder, f_name, data, include_header=True, data_type="pandas", extension="csv"):
    f_path_name = os.path.join(out_folder,f_name)
    if data_type == "pandas":
        data.to_csv(f_path_name, sep=",",header=include_header)
        return


