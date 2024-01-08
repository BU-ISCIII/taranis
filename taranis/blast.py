import logging
import os
import rich
import subprocess
import taranis.utils

from pathlib import Path
from Bio.Blast.Applications import NcbiblastnCommandline
import pdb

log = logging.getLogger(__name__)
stderr = rich.console.Console(
    stderr=True,
    style="dim",
    highlight=False,
    force_terminal=taranis.utils.rich_force_colors(),
)


class Blast:
    def __init__(self, db_type):
        self.db_type = db_type

    def create_blastdb(self, file_name, blast_dir):
        self.f_name = Path(file_name).stem
        db_dir = os.path.join(blast_dir, self.f_name)
        self.out_blast_dir = os.path.join(db_dir, self.f_name)

        blast_command = [
            "makeblastdb",
            "-in",
            file_name,
            "-parse_seqids",
            "-dbtype",
            self.db_type,
            "-out",
            self.out_blast_dir,
        ]
        try:
            _ = subprocess.run(
                blast_command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            )
        except Exception as e:
            log.error("Unable to create blast db for %s ", self.f_name)
            log.error(e)
            stderr.print(
                f"[red] Unable to create blast database for sample %s", self.f_name
            )
            exit(1)
        return

    def run_blast(
        self,
        query,
        evalue=0.001,
        perc_identity=90,
        reward=1,
        penalty=-2,
        gapopen=1,
        gapextend=1,
        max_target_seqs=10,
        max_hsps=10,
        num_threads=1,
    ):
        """_summary_
            blastn -outfmt "6 , qseqid , sseqid , pident ,  qlen , length , mismatch , gapopen , evalue , bitscore , sstart , send , qstart , qend , sseq , qseq" -query /media/lchapado/Reference_data/proyectos_isciii/taranis/documentos_antiguos/pasteur_schema/lmo0002.fasta -db /media/lchapado/Reference_data/proyectos_isciii/taranis/test/blastdb/RA-L2073_R1/RA-L2073_R1 -evalue 0.001 -penalty -2 -reward 1 -gapopen 1 -gapextend 1  -perc_identity 100 > /media/lchapado/Reference_data/proyectos_isciii/taranis/test/blast_sample_locus002.txt

        Args:
            query (_type_): _description_
            evalue (float, optional): _description_. Defaults to 0.001.
            perc_identity (int, optional): _description_. Defaults to 90.
            reward (int, optional): _description_. Defaults to 1.
            penalty (int, optional): _description_. Defaults to -2.
            gapopen (int, optional): _description_. Defaults to 1.
            gapextend (int, optional): _description_. Defaults to 1.
            max_target_seqs (int, optional): _description_. Defaults to 10.
            max_hsps (int, optional): _description_. Defaults to 10.
            num_threads (int, optional): _description_. Defaults to 1.
        """
        blast_parameters = '"6 , qseqid , sseqid , pident ,  qlen , length , mismatch , gapopen , evalue , bitscore , sstart , send , qstart , qend , sseq , qseq"'
        pdb.set_trace()
        # db=self.blast_dir, evalue=evalue, perc_identity=perc_identity_ref, reward=reward, penalty=penalty, gapopen=gapopen, gapextend=gapextend, outfmt=blast_parameters, max_target_seqs=max_target_seqs, max_hsps=max_hsps, num_threads=num_threads, query=core_reference_allele_path)
        cline = NcbiblastnCommandline(
            db=self.out_blast_dir,
            evalue=evalue,
            perc_identity=perc_identity,
            reward=reward,
            penalty=penalty,
            gapopen=gapopen,
            gapextend=gapextend,
            outfmt=blast_parameters,
            max_target_seqs=max_target_seqs,
            max_hsps=max_hsps,
            num_threads=num_threads,
            query=query,
        )
        try:
            out, _ = cline()
        except Exception as e:
            log.error("Unable to run blast for %s ", self.out_blast_dir)
            log.error(e)
            stderr.print(
                f"[red] Unable to run blast for database %s", self.out_blast_dir
            )
            exit(1)
        return out.splitlines()
