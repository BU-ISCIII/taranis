import io
import logging
import numpy as np
import rich.console
import os
import taranis.utils
import taranis.blast
from Bio import SeqIO
import pdb

log = logging.getLogger(__name__)
stderr = rich.console.Console(
    stderr=True,
    style="dim",
    highlight=False,
    force_terminal=taranis.utils.rich_force_colors(),
)


class EvaluateCluster:
    def __init__(self, locus_path: str, locus_name: str, output: str):
        self.locus_path = locus_path
        self.locus_name = locus_name

        self.output = os.path.join(output, "evaluate_cluster")
        taranis.utils.create_new_folder(self.output)
        # locus_blast_dir = os.path.join(self.output, locus_name)
        self.blast_obj = taranis.blast.Blast("nucl")
        _ = self.blast_obj.create_blastdb(locus_path, self.output)
        return

    def find_cluster_from_ref_allele(self, cluster_ref_alleles: dict) -> dict:
        return dict(
            [(value["center_id"], c_id) for c_id, value in cluster_ref_alleles.items()]
        )

    def summary(self, cluster_data: dict) -> list:
        summary_table = [
            "Locus name",
            "result",
            "alleles not found",
            "alleles not in cluster",
        ]
        sorted_cluster = sorted(cluster_data.keys())
        for cluster_id in sorted_cluster:
            row_data = [self.locus_name, str(cluster_id)]
            row_data.append(cluster_data[cluster_id]["result"])
            row_data.append(
                ";".join(cluster_data[cluster_id]["alleles_not_found"])
                if "alleles_not_found" in cluster_data[cluster_id]
                else "-"
            )
            row_data.append(
                ";".join(cluster_data[cluster_id]["alleles_not_in_cluster"])
                if "alleles_not_in_culster" in cluster_data[cluster_id]
                else "-"
            )
            summary_table.append(",".join(row_data))
        return summary_table

    def validate_cluster(self, blast_result: dict, cluster_data: list) -> dict:
        """For cluster validation, the sequence id matched in blast are compared
            with the cluster sequences. Return False validation if there are
            difference between them.

        Args:
            blast_result (dict): _description_
            cluster_data (list): _description_

        Returns:
            dict: _description_
        """
        # index of sequence id
        sseqid = 1
        blast_alleles = []
        alleles_not_in_cluster = []
        for match in blast_result:
            blast_allele = match.split("\t")[sseqid]
            if blast_allele in cluster_data:
                blast_alleles.append(blast_allele)
            else:
                alleles_not_in_cluster.append(blast_allele)

        if len(cluster_data) == len(blast_alleles):
            return {"validation": True}
        result = {"validation": False}
        # convert list to numpy array to find out differences
        c_alleles_np = np.array(list(cluster_data))
        blast_alleles_np = np.array(blast_alleles)
        result["alleles_not_found"] = np.setdiff1d(
            c_alleles_np, blast_alleles_np
        ).tolist()
        result["alleles_not_in_cluster"] = np.setdiff1d(
            blast_alleles_np, c_alleles_np
        ).tolist()
        return result

    def evaluate_clusters(
        self, cluster_alleles: dict, cluster_ref_alleles, ref_alleles_file: str
    ):
        reference_alleles = {}
        evaluation_alleles = {}
        ref_allele_in_cluster = self.find_cluster_from_ref_allele(cluster_ref_alleles)
        with open(ref_alleles_file, "r") as fh:
            for record in SeqIO.parse(fh, "fasta"):
                reference_alleles[record.id] = str(record.seq)

        for r_id, r_seq in reference_alleles.items():
            # create file in memory to increase speed
            query_file = io.StringIO()
            query_file.write(">" + r_id + "\n" + r_seq)
            query_file.seek(0)
            blast_result = self.blast_obj.run_blast(
                query_file.read(), perc_identity=90, query_type="stdin"
            )
            # Close object and discard memory buffer
            query_file.close()

            cluster_id = ref_allele_in_cluster[r_id]
            # pdb.set_trace()
            result_eval = self.validate_cluster(
                blast_result, cluster_alleles[cluster_id]
            )
            evaluation_alleles[cluster_id] = {}
            if result_eval["validation"] is False:
                evaluation_alleles[cluster_id]["result"] = "NOK"
                if len(result_eval["alleles_not_found"]) > 0:
                    evaluation_alleles[cluster_id]["alleles_not_found"] = result_eval[
                        "alleles_not_found"
                    ]
                if len(result_eval["alleles_not_in_cluster"]) > 0:
                    evaluation_alleles[cluster_id][
                        "alleles_not_in_cluster"
                    ] = result_eval["alleles_not_in_cluster"]
            else:
                evaluation_alleles[cluster_id]["result"] = "OK"
        return self.summary(evaluation_alleles)
