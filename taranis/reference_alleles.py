import logging

import rich.console
from pathlib import Path
import os

import taranis.utils
import taranis.distance
import taranis.clustering
from Bio import SeqIO

log = logging.getLogger(__name__)
stderr = rich.console.Console(
    stderr=True,
    style="dim",
    highlight=False,
    force_terminal=taranis.utils.rich_force_colors(),
)


class ReferenceAlleles:
    def __init__(self, fasta_file: str, output: str):
        """ReferenceAlleles instance creation

        Args:
            fasta_file (str): file name included path for locus
            output (str): output folder
        """
        self.fasta_file = fasta_file
        self.locus_name = Path(fasta_file).stem
        self.output = output
        self.selected_locus = {}

    def create_cluster_alleles(self) -> list:
        """Alleles in fasta file are clustering by using two additional classes:
            DistanceMatrix which creates a matrix of distance using the allele
            sequences, and ClusterDistance which get the matrix and group the
            alleles in clusters. As per result of ClusterDistance methods, the
            reference alleles are saved to file and statistics information is
            returned

        Returns:
            list: two dictionaires are returned, cluster_data having statistics
                and reference_alleles, where keys are cluster number and value
                the reference allele for the cluster
        """
        log.debug("Processing distance matrix for $s", self.fasta_file)
        distance_obj = taranis.distance.DistanceMatrix(self.fasta_file)
        mash_distance_df = distance_obj.create_matrix()
        log.debug(f"Created distance matrix for {self.fasta_file}")
        # fetch the allele position into array
        postition_to_allele = {
            x: mash_distance_df.columns[x] for x in range(len(mash_distance_df.columns))
        }
        # convert the  triangle matrix into full data matrix
        matrix_np = mash_distance_df.to_numpy()
        t_matrix_np = matrix_np.transpose()
        matrix_np = t_matrix_np + matrix_np
        # At this point minimal distance is 0. For clustering requires to be 1
        # the oposite.
        dist_matrix_np = (matrix_np - 1) * -1

        cluster_obj = taranis.clustering.ClusterDistance(
            dist_matrix_np, self.locus_name
        )
        cluster_ptrs, cluster_data = cluster_obj.create_clusters()
        # convert the center pointer to allele name and create list to get
        # sequences
        reference_alleles = []
        for cluster_id, values in cluster_data.items():
            center_allele = postition_to_allele[values["center_id"]]
            values["center_id"] = center_allele
            reference_alleles.append(center_allele)
        alleles_in_cluster = cluster_obj.convert_to_seq_clusters(
            cluster_ptrs, postition_to_allele
        )
        cluster_folder = os.path.join(self.output, "Clusters")
        _ = taranis.utils.create_new_folder(cluster_folder)
        cluster_file = os.path.join(
            cluster_folder, "cluster_alleles_" + self.locus_name + ".txt"
        )
        with open(cluster_file, "w") as fo:
            for cluster_id, alleles in alleles_in_cluster.items():
                fo.write("Cluster number" + str(cluster_id + 1) + "\n")
                fo.write("\n".join(alleles) + "\n")

        return [cluster_data, reference_alleles]

    def save_reference_alleles(self, reference_alleles: list) -> None:
        """From the input list it fetch the allele squence and save it as fasta

        Args:
            reference_alleles (list): list having the allele ids
        """
        record_seq = {}
        with open(self.fasta_file) as fh:
            for record in SeqIO.parse(fh, "fasta"):
                if record.id in reference_alleles:
                    record_seq[record.id] = str(record.seq)
        ref_allele_file = os.path.join(self.output, self.locus_name + ".fa")
        with open(ref_allele_file, "w") as fo:
            for r_id, r_seq in record_seq.items():
                fo.write(r_id + "\n")
                fo.write(r_seq + "\n")
        return

    def create_ref_alleles(self) -> dict:
        """Main method to create the reference alleles

        Returns:
            dict: statistics information for each cluster
        """
        self.records = taranis.utils.read_fasta_file(self.fasta_file)
        # Prepare data to use mash to create the distance matrix
        cluster_data, reference_alleles = self.create_cluster_alleles()
        _ = self.save_reference_alleles(reference_alleles)
        return cluster_data


def collect_statistics(data_alleles: list, out_folder: str) -> None:
    """Collect the individual statistics for each locus to create graphics

    Args:
        data_alleles (list): list having the indiviual statistics data
        out_folder (str): folder to save graphics
    """

    def stats_graphics(stats_folder: str, cluster_alleles: dict) -> None:
        """Create the statistics graphics. Bar graphic for number of cluster
            per alleles

        Args:
            stats_folder (str): folder path to store graphic
            cluster_alleles (dict): contain number of clusters as key and number
                of alleles having the same cluster number as value
        """
        stderr.print("Creating graphics")
        log.info("Creating graphics")
        graphic_folder = os.path.join(stats_folder, "graphics")
        _ = taranis.utils.create_new_folder(graphic_folder)
        cluster, alleles = zip(*cluster_alleles.items())
        _ = taranis.utils.create_graphic(
            graphic_folder,
            "num_genes_per_allele.png",
            "bar",
            cluster,
            alleles,
            ["Number of clusters", "number of genes"],
            "Number of cluster per gene",
        )

    cluster_alleles = {}
    for d_allele in data_alleles:
        cluster_alleles[len(d_allele)] = cluster_alleles.get(len(d_allele), 0) + 1
    _ = stats_graphics(out_folder, cluster_alleles)
    return
