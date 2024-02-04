import igraph as ig
import leidenalg
import logging
import numpy as np
import rich.console
import os
import taranis.utils

log = logging.getLogger(__name__)
stderr = rich.console.Console(
    stderr=True,
    style="dim",
    highlight=False,
    force_terminal=taranis.utils.rich_force_colors(),
)


class ClusterDistance:
    def __init__(self, dist_matrix: np.array, ref_seq_name: str):
        self.dist_matrix = dist_matrix
        self.num_seq = dist_matrix.shape[0]
        self.ref_seq_name = ref_seq_name
        self.seed = None
        self.res_param = 0.9

    def calculate_mean_cluster(self, src_cluster_row_idxs):
        src_cluster_col_idxs = src_cluster_row_idxs
        src_cluster_mtrx_idxs = np.ix_(src_cluster_row_idxs, src_cluster_col_idxs)
        row_idx_pos = np.argwhere(src_cluster_row_idxs).flatten()
        col_idx_pos = np.argwhere(src_cluster_col_idxs).flatten()
        num_of_diag_elements = np.intersect1d(row_idx_pos, col_idx_pos).size
        num_of_non_diag_elements = (
            row_idx_pos.size * col_idx_pos.size - num_of_diag_elements
        )
        if num_of_non_diag_elements == 0:
            return 1
        return (
            np.sum(self.dist_matrix[src_cluster_mtrx_idxs]) - num_of_diag_elements
        ) / num_of_non_diag_elements

    def convert_to_seq_clusters(seq_cluster_ptrs, seq_id_to_seq_name_map):
        output_seq_clusters = list()

        for cluster_id in range(np.max(seq_cluster_ptrs) + 1):
            seq_cluster = list()
            for seq_id in np.argwhere(seq_cluster_ptrs == cluster_id).flatten():
                seq_cluster.append(
                    "{}{}".format(seq_id_to_seq_name_map[str(seq_id)], os.linesep)
                )

            output_seq_clusters.append(seq_cluster)
        return output_seq_clusters

    def verify_cluster_quality(self, src_cluster_ptrs):
        log.debug(f"Verifying cluster for {self.ref_seq_name}")
        avg_clusters = {}
        for cluster_id in range(np.max(src_cluster_ptrs) + 1):
            log.debug(f"calculating mean for cluster number {cluster_id}")
            src_cluster_bool_ptrs = src_cluster_ptrs == cluster_id
            avg_clusters[cluster_id] = self.calculate_mean_cluster(
                src_cluster_bool_ptrs
            )

        return avg_clusters

    def create_clusters(self):
        comm_graph = ig.Graph.Weighted_Adjacency(
            self.dist_matrix.tolist(), mode=1, loops=False
        )
        graph_clusters = leidenalg.find_partition(
            comm_graph,
            leidenalg.CPMVertexPartition,
            weights="weight",
            n_iterations=-1,
            resolution_parameter=self.res_param,
            seed=self.seed,
        )
        cluster_ptrs = np.array(graph_clusters.membership)

        avg_clusters = self.verify_cluster_quality(cluster_ptrs)
        # check that cluste average values are upper than 0.9
        if not all(x > 0.9 for x in list(avg_clusters.values())):
            log.warning(
                f"There are some cluster below average of 0.9 in locus {self.ref_seq_name} "
            )
            stderr.print(
                f"[red]There are some cluster below average of 0.9 in locus {self.ref_seq_name}"
            )
        return cluster_ptrs
