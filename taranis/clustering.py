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
        """ClusterDistance instance creation

        Args:
            dist_matrix (np.array): distance matrix
            ref_seq_name (str): locus name
        """
        self.dist_matrix = dist_matrix
        self.num_seq = dist_matrix.shape[0]
        self.ref_seq_name = ref_seq_name
        self.seed = None
        self.res_param = 0.9

    def calculate_cluster_center(
        self, cluster_mtrx_idxs: tuple, cluster_mean: float
    ) -> int:
        """Get the center allele for the cluster by selecting the allele closest
            value to cluster mean

        Args:
            cluster_mtrx_idxs (tuple): tuple with the filter indexes to create
                submatrix for each cluster
            cluster_mean (float): cluster mean value to compare

        Returns:
            int: index of the allele which is the center of cluster
        """
        cluster_matrix = self.dist_matrix[cluster_mtrx_idxs]
        row_means = np.mean(cluster_matrix, axis=1)
        return cluster_mtrx_idxs[0][np.argmin(np.abs(row_means - cluster_mean))][0]

    def calculate_mean_cluster(
        self, cluster_mtrx_idxs: tuple, row_idx_pos: np.ndarray
    ) -> float:
        """Calculate the mean of cluster distance values

        Args:
            cluster_mtrx_idxs (tuple): tuple with the filter indexes to create
                submatrix for each cluster
            row_idx_pos (np.ndarray): indexes of matrix belongs to cluster

        Returns:
            float: mean of the cluster
        """
        col_idx_pos = row_idx_pos
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
        """convert the index into the alleale names

        Args:
            cluster_ids (np.array): index of matrix belongs to cluster
            id_to_seq_name (dict): having the index as key and allele name in
                value

        Returns:
            dict: where key is the cluster number and value is the list of
                alleles belongs to the cluster
        """
        out_clusters = {}
        for cluster_id in range(np.max(cluster_ids) + 1):
            out_clusters[cluster_id] = [
                id_to_seq_name[seq_id]
                for seq_id in np.argwhere(cluster_ids == cluster_id).flatten()
            ]

        return out_clusters

    def collect_data_cluster(self, src_cluster_ptrs: np.ndarray) -> dict:
        """Collect the mean, index allele center and number of alleles in
            cluster for each cluster

        Args:
            src_cluster_ptrs (np.ndarray): cluster matrix

        Returns:
            dict: where key is the cluster number and value a list of the
                statistic data
        """
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
            # get the closest distance coordenates to cluster mean value
            cluster_data[cluster_id]["avg"] = cluster_mean
            cluster_data[cluster_id]["center_id"] = self.calculate_cluster_center(
                cluster_mtrx_idxs, cluster_mean
            )
            log.debug(f"Get the cluster center for {cluster_id}")
            # get the number of sequences for the cluster
            cluster_data[cluster_id]["n_seq"] = len(cluster_mtrx_idxs[0])
        return cluster_data

    def create_clusters(self) -> list[dict]:
        """main method to create clustering using the Leiden algorithm

        Returns:
            list: two dictionaries are returned first with the cluster and the
            matrix indexes adn second the statistics data for each cluster
        """
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

        clusters_data = self.collect_data_cluster(cluster_ptrs)
        # check that cluste average values are upper than 0.9
        for value in clusters_data.values():
            if value["avg"] < 0.9:
                log.warning(
                    f"There are some cluster below average of 0.9 in locus {self.ref_seq_name} "
                )
                stderr.print(
                    f"[red]There are some cluster below average of 0.9 in locus {self.ref_seq_name}"
                )

        return [cluster_ptrs, clusters_data]
