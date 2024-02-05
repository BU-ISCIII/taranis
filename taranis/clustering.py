import igraph as ig
import leidenalg
import logging
import numpy as np
import rich.console
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

    def calculate_closest_index(
        self, cluster_mtrx_idxs: tuple, cluster_mean: float
    ) -> list:
        cluster_matrix = self.dist_matrix[cluster_mtrx_idxs]
        cluster_flat = cluster_matrix.flatten()
        closest_index = np.argmin(np.abs(cluster_flat - cluster_mean))
        return [np.unravel_index(closest_index, cluster_matrix.shape)]

    def calculate_mean_cluster(self, cluster_mtrx_idxs: tuple, row_idx_pos: np.ndarray):
        col_idx_pos = row_idx_pos
        # src_cluster_mtrx_idxs = np.ix_(src_cluster_row_idxs, src_cluster_col_idxs)
        # row_idx_pos = np.argwhere(src_cluster_row_idxs).flatten()
        # col_idx_pos = np.argwhere(src_cluster_col_idxs).flatten()
        num_of_diag_elements = np.intersect1d(row_idx_pos, col_idx_pos).size
        num_of_non_diag_elements = (
            row_idx_pos.size * col_idx_pos.size - num_of_diag_elements
        )
        if num_of_non_diag_elements == 0:
            return 1
        return (
            np.sum(self.dist_matrix[cluster_mtrx_idxs]) - num_of_diag_elements
        ) / num_of_non_diag_elements

    def convert_to_seq_clusters(
        self, cluster_ids: np.array, id_to_seq_name: dict
    ) -> dict:
        out_clusters = {}
        for cluster_id in range(np.max(cluster_ids) + 1):
            alleles_in_cluster = []
            out_clusters[cluster_id] = [
                id_to_seq_name[seq_id]
                for seq_id in np.argwhere(cluster_ids == cluster_id).flatten()
            ]

        return out_clusters

    def collect_data_cluster(self, src_cluster_ptrs):
        log.debug(f"Collecting data for cluster {self.ref_seq_name}")
        cluster_data = {}
        for cluster_id in range(np.max(src_cluster_ptrs) + 1):
            cluster_data[cluster_id] = {}
            log.debug(f"calculating mean for cluster number {cluster_id}")
            cluster_bool_ptrs = src_cluster_ptrs == cluster_id
            cluster_mtrx_idxs = np.ix_(cluster_bool_ptrs, cluster_bool_ptrs)
            row_idx_pos = np.argwhere(cluster_bool_ptrs).flatten()
            # col_idx_pos = np.argwhere(cluster_bool_ptrs).flatten()
            cluster_mean = self.calculate_mean_cluster(cluster_mtrx_idxs, row_idx_pos)
            # pdb.set_trace()
            cluster_data[cluster_id]["avg"] = cluster_mean
            cluster_data[cluster_id]["closest_idx"] = self.calculate_closest_index(
                cluster_mtrx_idxs, cluster_mean
            )
        return cluster_data

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
        # Convert the partition to a DataFrame
        # df_clusters = pd.DataFrame({'Node': range(len(graph_clusters.membership)), 'Cluster': graph_clusters.membership})
        # Calculate the centroid of each cluster
        # cluster_centers = df_clusters.groupby('Cluster').apply(lambda x: np.mean(self.dist_matrix[x['Node']], axis=0)).values
        clusters_data = self.collect_data_cluster(cluster_ptrs)
        # check that cluste average values are upper than 0.9
        for value in clusters_data.values():
            if value["avg"] < 0.9 :
                log.warning(
                    f"There are some cluster below average of 0.9 in locus {self.ref_seq_name} "
                )
                stderr.print(
                    f"[red]There are some cluster below average of 0.9 in locus {self.ref_seq_name}"
                )

        return cluster_ptrs
