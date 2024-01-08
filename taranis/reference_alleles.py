import logging
import numpy as np
import os
import re
import rich.console
# import sys
import subprocess

# from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastnCommandline
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
                    if len(locus_prot) == int(len(seq) / 3):
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
                        if len(locus_prot) == int(len(record.seq) / 3):
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
        # f_name = os.path.basename(self.fasta_file).split('.')[0]
        f_name = os.path.basename(self.fasta_file)
        mash_folder = os.path.join(self.output, "mash")
        # _ = taranis.utils.write_fasta_file(mash_folder, self.selected_locus, multiple_files=True, extension=False)
        # save directory to return after mash
        # working_dir = os.getcwd()
        os.chdir(mash_folder)
        # run mash sketch command
        sketch_file = "reference.msh"
        mash_sketch_command = ["mash", "sketch", "-i", "-o", sketch_file, f_name]
        # mash sketch -i -o prueba.msh lmo0003.fasta
        # mash_sketch_command += list(self.selected_locus.keys())
        _ = subprocess.run(
            mash_sketch_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        # Get pairwise allele sequences mash distances
        # mash_distance_command = ["mash", "dist", sketch_path, sketch_path]
        mash_distance_command = ["mash", "triangle", "-i", "reference.msh"]
        mash_distance_result = subprocess.Popen(
            mash_distance_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        # pdb.set_trace()
        out, err = mash_distance_result.communicate()
        with open("matrix_distance.tsv", "w") as fo:
            # adding alleles to create a heading
            # the value are not required to be in order, just only any name and the right length
            fo.write("alleles\t" + "\t".join(list(self.selected_locus.keys())) + "\n")
            fo.write(out.decode("UTF-8"))
        import pandas as pd

        locus_num = len(self.selected_locus)
        # pdb.set_trace()
        matrix_df = pd.read_csv("matrix_distance.tsv", sep="\t").fillna(value=0)
        # remove the first line of the matrix that contain only the number of alleles
        matrix_df = matrix_df.drop(0)
        locus_list = matrix_df.iloc[0:locus_num, 0]
        matrix_np = matrix_df.iloc[:, 1:].to_numpy()
        # convert the triangular matrix to mirror up triangular part
        t_matrix_np = matrix_np.transpose()
        matrix_np = t_matrix_np + matrix_np
        # values_np = matrix_df.iloc[:,2].to_numpy()

        # matrix_np = values_np.reshape(locus_num, locus_num)
        # out = out.decode('UTF-8').split('\n')
        from sklearn.cluster import AgglomerativeClustering

        clusterer = AgglomerativeClustering(
            n_clusters=7,
            metric="precomputed",
            linkage="average",
            distance_threshold=None,
        )
        clusters = clusterer.fit_predict(matrix_np)
        # clustering = AgglomerativeClustering(affinity="precomputed").fit(matrix_np)
        mean_distance = np.mean(matrix_np, 0)
        # std = np.std(matrix_np)
        min_mean = min(mean_distance)
        mean_all_alleles = np.mean(mean_distance)
        max_mean = max(mean_distance)
        # buscar el indice que tiene el minimo valor de media
        min_mean_idx = np.where(mean_distance == float(min_mean))[0][0]
        # create fasta file with the allele
        min_allele = self.selected_locus[locus_list[min_mean_idx]]

        record_allele_folder = os.path.join(os.getcwd(), f_name.split(".")[0])
        min_allele_file = taranis.utils.write_fasta_file(
            record_allele_folder, min_allele, locus_list[min_mean_idx]
        )
        # pdb.set_trace()
        # busca el indice que tiene el valor de la media
        mean_all_closser_value = taranis.utils.find_nearest_numpy_value(
            mean_distance, mean_all_alleles
        )
        mean_all_alleles_idx = np.where(mean_distance == float(mean_all_closser_value))[
            0
        ][0]
        # create fasta file with the allele
        mean_allele = self.selected_locus[locus_list[mean_all_alleles_idx]]
        # record_allele_folder = os.path.join(mash_folder, f_name)
        mean_allele_file = taranis.utils.write_fasta_file(
            record_allele_folder, mean_allele, locus_list[mean_all_alleles_idx]
        )

        # busca el indice con la mayor distancia
        max_mean_idx = np.where(mean_distance == float(max_mean))[0][0]
        # create fasta file with the allele
        max_allele = self.selected_locus[locus_list[max_mean_idx]]
        max_allele_file = taranis.utils.write_fasta_file(
            record_allele_folder, max_allele, locus_list[max_mean_idx]
        )

        # Elijo un outlier lmo0002_185 para ver la distancia
        outlier_allele = self.selected_locus[locus_list[184]]
        outlier_allele_file = taranis.utils.write_fasta_file(
            record_allele_folder, outlier_allele, locus_list[184]
        )

        # elijo un segundo outlier lmo0002_95 que tiene como cluster =1
        outlier2_allele = self.selected_locus[locus_list[95]]
        outlier2_allele_file = taranis.utils.write_fasta_file(
            record_allele_folder, outlier2_allele, locus_list[95]
        )

        # elijo un tercer outlier lmo0002_185 que tiene como cluster =4
        outlier3_allele = self.selected_locus[locus_list[185]]
        outlier3_allele_file = taranis.utils.write_fasta_file(
            record_allele_folder, outlier3_allele, locus_list[185]
        )

        # saca una lista de cuantas veces se repite un valor
        np.bincount(clusters)
        blast_parameters = '"6 , qseqid , sseqid , pident ,  qlen , length , mismatch , gapopen , evalue , bitscore , sstart , send , qstart , qend , sseq , qseq"'

        # Create local BLAST database for all alleles in the locus
        db_name = "/media/lchapado/Reference_data/proyectos_isciii/taranis/new_taranis_result_code/blast/locus_db"
        # db_name = os.path.join("blast", 'locus_blastdb')
        # fasta_file = "/media/lchapado/Reference_data/proyectos_isciii/taranis/documentos_antiguos/datos_prueba/schema_1_locus/lmo0002.fasta"
        # pdb.set_trace()
        # blast_command = ['makeblastdb' , '-in' , fasta_file , '-parse_seqids', '-dbtype',  "nucl", '-out' , db_name]
        # blast_result = subprocess.run(blast_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # taranis.utils.create_blastdb(fasta_file, db_name, 'nucl', logger):
        # locus_db_name = os.path.join(db_name, f_name[0], f_name[0])
        # query_data= self.selected_locus["lmo0002_1"]
        # All alleles in locus VS reference allele chosen (centroid) BLAST

        # ref_query_file="/media/lchapado/Reference_data/proyectos_isciii/taranis/new_taranis_result_code/mash/lmo0002/query.fasta"
        # cline = NcbiblastnCommandline(db=db_name, evalue=0.001, perc_identity=100, reward=1, penalty=-2, gapopen=1, gapextend=1, outfmt=blast_parameters, max_target_seqs=0, max_hsps=0, num_threads=4, query=ref_query_file)

        # minima distancia .
        # min_dist_file="/media/lchapado/Reference_data/proyectos_isciii/taranis/new_taranis_result_code/mash/lmo0002/lmo0002_610"
        # pdb.set_trace()
        min_dist_file = os.path.join(record_allele_folder, min_allele_file)
        cline = NcbiblastnCommandline(
            db=db_name,
            evalue=0.001,
            perc_identity=90,
            reward=1,
            penalty=-2,
            gapopen=1,
            gapextend=1,
            outfmt=blast_parameters,
            max_target_seqs=1100,
            max_hsps=1000,
            num_threads=4,
            query=min_dist_file,
        )
        out, err = cline()
        min_dist_lines = out.splitlines()
        min_dist_alleles = []
        for min_dist in min_dist_lines:
            min_dist_alleles.append(min_dist.split("\t")[1])
        min_np = np.array(min_dist_alleles)
        # pdb.set_trace()
        print("matches con min distancia: ", len(min_dist_lines))
        print("Not coverage using as reference", np.setdiff1d(locus_list, min_np))
        # distancia media. Sale 133 matches
        # mean_dist_file="/media/lchapado/Reference_data/proyectos_isciii/taranis/new_taranis_result_code/mash/lmo0002/lmo0002_870"
        mean_dist_file = os.path.join(record_allele_folder, mean_allele_file)
        cline = NcbiblastnCommandline(
            db=db_name,
            evalue=0.001,
            perc_identity=90,
            reward=1,
            penalty=-2,
            gapopen=1,
            gapextend=1,
            outfmt=blast_parameters,
            max_target_seqs=1100,
            max_hsps=1000,
            num_threads=4,
            query=mean_dist_file,
        )
        out, err = cline()
        mean_dist_lines = out.splitlines()
        mean_dist_alleles = []
        for mean_dist in mean_dist_lines:
            mean_dist_alleles.append(mean_dist.split("\t")[1])
        mean_np = np.array(mean_dist_alleles)
        print("matches con distancia media: ", len(mean_dist_lines))
        print("Not coverage using as reference", np.setdiff1d(locus_list, mean_np))

        # maxima distancia,
        # ref_query_file="/media/lchapado/Reference_data/proyectos_isciii/taranis/new_taranis_result_code/mash/lmo0002/lmo0002_216"
        max_dist_file = os.path.join(record_allele_folder, max_allele_file)
        cline = NcbiblastnCommandline(
            db=db_name,
            evalue=0.001,
            perc_identity=90,
            reward=1,
            penalty=-2,
            gapopen=1,
            gapextend=1,
            outfmt=blast_parameters,
            max_target_seqs=1100,
            max_hsps=1000,
            num_threads=4,
            query=max_dist_file,
        )
        out, err = cline()
        max_dist_lines = out.splitlines()
        max_dist_alleles = []
        for max_dist in max_dist_lines:
            max_dist_alleles.append(max_dist.split("\t")[1])
        max_np = np.array(max_dist_alleles)
        print("matches con max distancia: ", len(max_dist_lines))
        print("Not coverage using as reference", np.setdiff1d(locus_list, max_np))

        # eligiendo uno de los outliers ,
        # outlier_file="/media/lchapado/Reference_data/proyectos_isciii/taranis/new_taranis_result_code/mash/lmo0002/lmo0002_183"
        outlier_file = os.path.join(record_allele_folder, outlier_allele_file)
        cline = NcbiblastnCommandline(
            db=db_name,
            evalue=0.001,
            perc_identity=90,
            reward=1,
            penalty=-2,
            gapopen=1,
            gapextend=1,
            outfmt=blast_parameters,
            max_target_seqs=1100,
            max_hsps=1000,
            num_threads=4,
            query=outlier_file,
        )
        out, err = cline()
        outlier_lines = out.splitlines()
        outlier_alleles = []
        for outlier_line in outlier_lines:
            outlier_alleles.append(outlier_line.split("\t")[1])
        outlier_np = np.array(outlier_alleles)
        print("matches con outliers distancia: ", len(outlier_lines))

        print("Alleles added using outlier as reference", outlier_np)
        new_ref_np = np.unique(np.concatenate((min_np, outlier_np), axis=0))
        print("\n", "remaining alleles ", np.setdiff1d(locus_list, new_ref_np))

        # eligiendo el segundo de los outliers ,
        # outlier_file="/media/lchapado/Reference_data/proyectos_isciii/taranis/new_taranis_result_code/mash/lmo0002/lmo0002_183"
        outlier2_file = os.path.join(record_allele_folder, outlier2_allele_file)
        cline = NcbiblastnCommandline(
            db=db_name,
            evalue=0.001,
            perc_identity=90,
            reward=1,
            penalty=-2,
            gapopen=1,
            gapextend=1,
            outfmt=blast_parameters,
            max_target_seqs=1100,
            max_hsps=1000,
            num_threads=4,
            query=outlier2_file,
        )
        out, _ = cline()
        outlier2_lines = out.splitlines()
        outlier2_alleles = []
        for outlier2_line in outlier2_lines:
            outlier2_alleles.append(outlier2_line.split("\t")[1])
        outlier2_np = np.array(outlier2_alleles)
        print("matches con second outliers distance: ", len(outlier2_lines))
        # print("Alleles added using second outlier as reference" , outlier2_np)
        upd_new_ref_np = np.unique(np.concatenate((new_ref_np, outlier2_np), axis=0))
        print(
            "\n",
            "remaining alleles after second outlier",
            np.setdiff1d(locus_list, upd_new_ref_np),
        )

        # eligiendo el tercero de los outliers ,
        # outlier_file="/media/lchapado/Reference_data/proyectos_isciii/taranis/new_taranis_result_code/mash/lmo0002/lmo0002_183"
        outlier3_file = os.path.join(record_allele_folder, outlier3_allele_file)
        cline = NcbiblastnCommandline(
            db=db_name,
            evalue=0.001,
            perc_identity=90,
            reward=1,
            penalty=-2,
            gapopen=1,
            gapextend=1,
            outfmt=blast_parameters,
            max_target_seqs=1100,
            max_hsps=1000,
            num_threads=4,
            query=outlier3_file,
        )
        out, _ = cline()
        outlier3_lines = out.splitlines()
        outlier3_alleles = []
        for outlier3_line in outlier3_lines:
            outlier3_alleles.append(outlier3_line.split("\t")[1])
        outlier3_np = np.array(outlier3_alleles)
        print("matches con third outliers distance: ", len(outlier3_lines))
        # print("Alleles added using second outlier as reference" , outlier2_np)
        upd2_new_ref_np = np.unique(
            np.concatenate((upd_new_ref_np, outlier3_np), axis=0)
        )
        print(
            "\n",
            "remaining alleles after second outlier",
            np.setdiff1d(locus_list, upd2_new_ref_np),
        )

        print("\n Still missing ", len(np.setdiff1d(locus_list, upd2_new_ref_np)))

        pdb.set_trace()

        # from sklearn.cluster import AgglomerativeClustering
        # import numpy as np
        # X = np.array([[0, 2, 3], [2, 0, 3], [3, 3, 0]])
        # clustering = AgglomerativeClustering(affinity="precomputed").fit(X)

    def create_ref_alleles(self):
        self.records = taranis.utils.read_fasta_file(self.fasta_file)
        _ = self.check_locus_quality()
        # pdb.set_trace()
        # Prepare data to use mash to create the distance matrix
        _ = self.create_matrix_distance()

        pass
