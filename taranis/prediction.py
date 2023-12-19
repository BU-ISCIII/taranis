import logging
import os
import rich
import subprocess
import taranis.utils

from pathlib import Path

log = logging.getLogger(__name__)
stderr = rich.console.Console(
    stderr=True,
    style="dim",
    highlight=False,
    force_terminal=taranis.utils.rich_force_colors(),
)


class Prediction():
    def __init__(self, genome_ref, sample_file, out_dir):
        self.genome_ref = genome_ref
        self.sample_file = sample_file
        self.sample_name = Path(sample_file).stem
        self.out_dir = out_dir
        self.train = os.path.join(out_dir, self.sample_name + ".trn")
        self.pred_protein = os.path.join(out_dir, self.sample_name + "_prot.faa")
        self.pred_gene = os.path.join(out_dir, self.sample_name + "_dna.faa")
        self.pred_coord = os.path.join(out_dir, self.sample_name + "_coord.gff")

        if not os.path.exists(self.out_dir):
            try:
                os.makedirs(self.out_dir, exist_ok=True)
                log.debug("Created directory %s for prediction ", self.out_dir)
            except OSError as e:
                log.error("Cannot create %s directory", self.out_dir)
                log.error(e)
                stderr.print (f"[red] Unable to create {self.out_dir} folder")
                exit(1)

    def training(self):
        prodigal_command = ["prodigal" , "-i", self.genome_ref, "-t", self.train]
        try:
            _ = subprocess.run(prodigal_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        except Exception as e:
            log.error("Unable to execute prodigal command for training")
            log.error(e)
            stderr.print (f"[red] Unable to run prodigal commmand. ERROR {e} ")
            exit(1)
        return



    def prediction(self):
        
        prodigal_command = ["prodigal" , "-i", self.sample_file , "-t", self.train, "-f", "gff", "-o", self.pred_coord, "-a", self.pred_protein, "-d", self.pred_gene]
        try:
            _ = subprocess.run(prodigal_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        except Exception as e:
            log.error("Unable to execute prodigal command for training")
            log.error(e)
            stderr.print (f"[red] Unable to run prodigal commmand. ERROR {e} ")
            exit(1)
        return

    def get_sequence(self):
        return