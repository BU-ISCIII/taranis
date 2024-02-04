import io
import logging
import pandas as pd
import subprocess
import rich
import sys
from pathlib import Path
import taranis.utils

log = logging.getLogger(__name__)
stderr = rich.console.Console(
    stderr=True,
    style="dim",
    highlight=False,
    force_terminal=taranis.utils.rich_force_colors(),
)


class DistanceMatrix:
    def __init__(
        self,
        file_path: str,
    ) -> "DistanceMatrix":
        """DistanceMatrix instance creation

        Args:
            file_path (str): Locus file path

        Returns:
            DistanceMatrix: created instance
        """
        self.file_path = file_path

    def create_matrix(self) -> pd.DataFrame:
        """Create distance matrix using external program called mash

        Returns:
            pd.DataFrame: Triangular distance matrix as panda DataFrame
        """
        allele_name = Path(self.file_path).stem
        mash_distance_command = [
            "mash",
            "triangle",
            "-i",
            self.file_path,
            "-k",
            "17",
            "-s",
            "2000",
        ]
        try:
            mash_distance_result = subprocess.Popen(
                mash_distance_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            out, _ = mash_distance_result.communicate()
            log.debug(f"calculate mash distance for {allele_name}")
        except Exception as e:
            log.error(f"Unable to create distance matrix for {self.file_path}. {e}")
            stderr.print(
                f"[red] Error when creating distance matrix for {self.file_path}"
            )
            stderr.print(f"{e}")
            sys.exit(1)

        out_data = out.decode("UTF-8").split("\n")
        allele_names = [item.split("\t")[0] for item in out_data[1:-1]]
        # create file in memory to increase speed
        dist_matrix = io.StringIO()
        dist_matrix.write("alleles\t" + "\t".join(allele_names) + "\n")
        dist_matrix.write("\n".join(out_data[1:]))
        dist_matrix.seek(0)
        matrix_pd = pd.read_csv(dist_matrix, sep="\t", index_col="alleles").fillna(0)
        # Close object and discard memory buffer
        dist_matrix.close()
        log.debug(f"create distance for {allele_name}")
        return matrix_pd
