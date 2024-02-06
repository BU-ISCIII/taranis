import logging

import rich.console
from pathlib import Path
import os

from scipy.sparse import coo_matrix
import pdb
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
    def __init__(self, fasta_file, output):
        self.fasta_file = fasta_file
        self.locus_name = Path(fasta_file).stem
        self.output = output
        self.selected_locus = {}

    def create_cluster_alleles(self):
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

        # create a sparse matrix used for summary
        _ = coo_matrix(matrix_np, shape=matrix_np.shape)

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
        cluster_file = os.path.join(
            self.output, "cluster_alleles_" + self.locus_name + ".txt"
        )
        pdb.set_trace()
        with open(cluster_file, "w") as fo:
            for cluster_id, alleles in alleles_in_cluster.items():
                fo.write("Cluster number" + str(cluster_id + 1) + "\n")
                fo.write("\n".join(alleles) + "\n")
        pdb.set_trace()

        return cluster_data, reference_alleles

    def save_reference_alleles(self, reference_alleles: list) -> None:
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

    def create_ref_alleles(self):
        self.records = taranis.utils.read_fasta_file(self.fasta_file)
        # Prepare data to use mash to create the distance matrix
        cluster_data, reference_alleles = self.create_cluster_alleles()
        _ = self.save_reference_alleles(reference_alleles)

        pass
