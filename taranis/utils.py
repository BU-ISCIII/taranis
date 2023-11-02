#!/usr/bin/env python
"""
Common utility function used for relecov_tools package.
"""

import glob
# import hashlib
import logging
import questionary
# import json
# import openpyxl
# import yaml
# from itertools import islice

import os
from rich.console import Console
import sys

from Bio import SeqIO

log = logging.getLogger(__name__)


def get_files_in_folder(folder, extension=None):
    """get the list of files, filtered by extension in the input folder. If
    extension is not set, then all files in folder are returned

    Args:
        folder (string): folder path
        extension (string, optional): extension for filtering. Defaults to None.
    
    Returns:
        list: list of files
    """
    if extension is None:
        extension = "*"
    return glob.glob(folder + "*." + extension)


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
            sys.stdout.write("Please respond with 'yes' or 'no' " "(or 'y' or 'n').\n")

def read_fasta_file(fasta_file):
    return SeqIO.parse(fasta_file, "fasta")
    
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

def write_fasta_file(out_folder, seq_data, multiple_files=False, extension=True):
    try:
        os.makedirs(out_folder, exist_ok=True)
    except OSError as e:
        sys.exit(1)
    if isinstance(seq_data, dict):
        for key, seq in seq_data.items():
            if extension:
                f_name = os.path.join(out_folder, key + ".fasta")
            else:
                f_name = os.path.join(out_folder, key)
            with open (f_name, "w") as fo:
                fo.write(">" + key + "\n")
                fo.write(seq)
        