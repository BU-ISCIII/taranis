import logging

import rich.console
from pathlib import Path

from scipy.sparse import coo_matrix

import taranis.utils
import taranis.distance
import taranis.clustering

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
        """ in alfaclust TO DELETE
        sparse_edge_weight_mtrx = coo_matrix(global_edge_weight_mtrx, shape=global_edge_weight_mtrx.shape)
        """
        # create a sparse matrix used for summary
        _ = coo_matrix(matrix_np, shape=matrix_np.shape)
        pdb.set_trace()
        cluster_seq = taranis.clustering.ClusterDistance(
            dist_matrix_np, self.locus_name
        )
        cluster_ptrs = cluster_seq.create_clusters()
        cluster_seq.convert_to_seq_clusters(cluster_ptrs, postition_to_allele)

    def create_ref_alleles(self):
        self.records = taranis.utils.read_fasta_file(self.fasta_file)
        # _ = self.check_locus_quality()
        # pdb.set_trace()
        # Prepare data to use mash to create the distance matrix
        _ = self.create_cluster_alleles()

        pass
