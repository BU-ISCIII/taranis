import logging
import os
import re
import rich.console
import sys
import subprocess

# from Bio import SeqIO
from Bio.Seq import Seq

import taranis.utils

import pdb

log = logging.getLogger(__name__)
stderr = rich.console.Console(
    stderr=True,
    style="dim",
    highlight=False,
    force_terminal=taranis.utils.rich_force_colors(),
)

class ReferenceAlleles:
    def __init__(self, fasta_file, output):
        """
        self.schema_dir = schema_dir
        self.out_dir = out_dir
        if self.schema_dir is None:
            self.schema_dir = taranis.utils.prompt_text("Write the path of the schema`s files")
        if not taranis.utils.folder_exists(self.schema_dir):
            log.error("schema folder %s does not exists", self.schema_dir)
            stderr.print(
               "[red] Schema folder does not exist. " + self.schema_dir + "!"
            )
            sys.exit(1)
        if out_dir is None:
            self.out_dir = taranis.utils.prompt_text("Define the the directory to save results")
        # Check if folder exists
        if taranis.utils.folder_exists(self.out_dir):
            q_question = "Folder " + self.out_dir + " already exists. Files will be overwritten. Do you want to continue?"
            if "no" in taranis.utils.query_user_yes_no(q_question, "no"):
                log.info("Aborting code by user request")
                stderr.print("[red] Exiting code. ")
                sys.exit(1)
        else:
            try:
                os.makedirs(self.out_dir)
            except OSError as e:
                log.info("Unable to create folder at %s", self.out_dir)
                stderr.print("[red] ERROR. Unable to create folder  " + self.out_dir)
                sys.exit(1)
        """
        self.fasta_file = fasta_file
        self.output = output
        self.records = None
        self.locus_quality = {}
        self.selected_locus = {}

    def check_locus_quality(self):
        # START_CODONS_FORWARD = ['ATG', 'ATA', 'ATT', 'GTG', 'TTG', 'CTG']
        # START_CODONS_REVERSE = ['CAT', 'TAT', 'AAT', 'CAC', 'CAA', 'CAG']

        STOP_CODONS_FORWARD = ["TAA", "TAG", "TGA"]
        STOP_CODONS_REVERSE = ["TTA", "CTA", "TCA"]
        for record in self.records:
            # Check if start condon forward
            seq = str(record.seq)
            s_codon_f = re.match(r"^(ATG|ATA|ATT|GTG|TTG|CTG).+(\w{3})$", seq)
            if s_codon_f:
                #  Check if last 3 characters are codon stop forward
                if s_codon_f.group(2) in STOP_CODONS_FORWARD:
                    # Check if multiple stop codon by translating to protein and
                    # comparing length
                    locus_prot = Seq(record.seq).translate()
                    if len(locus_prot) == int(len(seq)/3):
                        self.locus_quality[record.id] = "good quality"
                        self.selected_locus[record.id] = seq
                    else:
                        self.locus_quality[record.id] = "bad quality: multiple_stop"
                else:
                    self.locus_quality[record.id] = "bad quality: no_stop"
            else:
                # Check if start codon reverse
                s_codon_r = re.match(r"^(\w{3}).+ (CAT|TAT|AAT|CAC|CAA|CAG)$", seq)
                if s_codon_r:
                    # Matched reverse start codon
                    if s_codon_f.group(1) in STOP_CODONS_REVERSE:
                        locus_prot = Seq(record.seq).reverse_complement().translate()
                        if len(locus_prot) == int(len(record.seq)/3):
                            self.locus_quality[record.id] = "good quality"
                            self.selected_locus[record.id] = seq
                        else:
                            self.locus_quality[record.id] = "bad quality: multiple_stop"
                    else:
                        self.locus_quality[record.id] = "bad quality: no_stop"
                else:
                    self.locus_quality[record.id] = "bad_quality: no_start"
        return

    def create_matrix_distance(self):
        f_name = os.path.basename(self.fasta_file).split('.')[0]
        mash_folder = os.path.join(self.output, "mash", f_name )
        _ = taranis.utils.write_fasta_file(mash_folder, self.selected_locus, multiple_files=True, extension=False)
        # save directory to return after mash
        working_dir = os.getcwd()
        os.chdir(mash_folder)
        # run mash sketch command
        sketch_path =  "reference.msh"
        mash_sketch_command = ["mash", "sketch", "-o", sketch_path]
        mash_sketch_command += list(self.selected_locus.keys())
       
        mash_sketch_result = subprocess.run(mash_sketch_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # Get pairwise allele sequences mash distances
        mash_distance_command = ["mash", "dist", sketch_path, sketch_path]
        mash_distance_result = subprocess.Popen(mash_distance_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
        out, err = mash_distance_result.communicate()
        with open ("matrix_distance.tsv", "w") as fo:
            fo.write(out.decode("UTF-8"))
        import pandas as pd
        locus_num = len(self.selected_locus)
        
        matrix_df = pd.read_csv("matrix_distance.tsv", sep="\t", header=None)
        list_alleles = matrix_df.iloc[0:locus_num,0]
        values_np = matrix_df.iloc[:,2].to_numpy()

        matrix_np = values_np.reshape(locus_num, locus_num)
        # out = out.decode('UTF-8').split('\n')
        from sklearn.cluster import AgglomerativeClustering
        clusterer = AgglomerativeClustering(n_clusters=3, metric="precomputed", linkage="average", distance_threshold=None)
        clusters = clusterer.fit_predict(matrix_np)
        # clustering = AgglomerativeClustering(affinity="precomputed").fit(matrix_np)
        pdb.set_trace()

        # from sklearn.cluster import AgglomerativeClustering
        # import numpy as np
        # X = np.array([[0, 2, 3], [2, 0, 3], [3, 3, 0]])
        # clustering = AgglomerativeClustering(affinity="precomputed").fit(X)


    def create_ref_alleles(self):
        self.records = taranis.utils.read_fasta_file(self.fasta_file)
        _ = self.check_locus_quality()
        pdb.set_trace()
        # Prepare data to use mash to create the distance matrix
        _ = self.create_matrix_distance()
        

        pass