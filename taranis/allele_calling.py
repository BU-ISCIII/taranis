import logging
import os
import rich.console


import taranis.utils
import taranis.blast
# import numpy
import pandas as pd
from pathlib import Path


import pdb
log = logging.getLogger(__name__)
stderr = rich.console.Console(
    stderr=True,
    style="dim",
    highlight=False,
    force_terminal=taranis.utils.rich_force_colors(),
)

class AlleleCalling:
    def __init__(self, prediction, sample_file, schema, reference_alleles, out_folder):
        self.prediction = prediction
        self.sample_file = sample_file
        self.schema = schema
        self.ref_alleles = reference_alleles
        self.out_folder = out_folder
        self.s_name = Path(sample_file).stem
        self.blast_dir = os.path.join(out_folder,"blastdb")
        self.blast_sample = os.path.join(self.blast_dir, self.s_name)
        self.blast_heading = ["qseqid", "sseqid", "pident", "qlen", "length", "mismatch", "gapopen", "evalue", "bitscore", "sstart", "send", "qstart", "qend", "sseq", "qseq"]

    def assign_allele_type(self, query_seq, allele_name, sample_contig, schema_gene):
        """_summary_

        Args:
            query_seq (_type_): _description_
            allele_name (_type_): _description_
            sample_contig (_type_): _description_
            schema_gene (_type_): _description_
        """        
        s_alleles_blast = taranis.blast.Blast("nucl")
        ref_allele_blast_dir = os.path.join(self.blast_dir, "ref_alleles")
        query_path = os.path.join(self.out_folder, "tmp", allele_name)
        # Write to file the sequence to find out the loci name that fully match 
        f_name = taranis.utils.write_fasta_file(query_path, query_seq, allele_name)
        query_file = os.path.join(query_path, f_name)
        _ = s_alleles_blast.create_blastdb(schema_gene, ref_allele_blast_dir)
        # Blast with sample sequence to find the allele in the schema 
        seq_blast_match = s_alleles_blast.run_blast(query_file, perc_identity=100)
        pdb.set_trace()
        if len(seq_blast_match) >= 1:
            # allele is named as NIPHEM 
           
            # Hacer un blast con la query esta secuencia y la database del alelo
            # Create  blast db with sample file
            

            pass
        elif len(seq_blast_match) == 1:
            pass
        else:
            pass


    def search_alleles (self, ref_allele):
        allele_name = Path(ref_allele).stem
        schema_gene = os.path.join(self.schema, allele_name + ".fasta")  
        allele_name = Path(ref_allele).stem
        # run blast with sample as db and reference allele as query
        sample_blast_match = self.sample_blast.run_blast(ref_allele)
        if len(sample_blast_match) > 0 :
            pd_lines = pd.DataFrame([item.split("\t") for item in sample_blast_match])
            pd_lines.columns = self.blast_heading
            pd_lines["pident"] = pd_lines["pident"].apply(pd.to_numeric)
            sel_max = pd_lines.loc[pd_lines["pident"].idxmax()]
            query_seq = sel_max["qseq"]
            # np_lines = numpy.array(s_lines)
            # convert to float the perc_identity to find out the max value
            # max_val = numpy.max(np_lines[:,2].astype(float))
            # mask = np_lines[:, 2] ==str(max_val)
            # Select rows that match the percent identity. Index 2 in blast results
            # sel_row = np_lines[mask, :] = np_lines[mask, :]
            # query_seq = sel_row[0,14]
            sample_contig = sel_max["sseqid"]
            abbr = self.assign_allele_type(query_seq, allele_name, sample_contig, schema_gene)
        else:
            # Sample does not have a reference allele to be matched
            # Keep LNF info
            # ver el codigo de espe
            #lnf_tpr_tag()
            pass
        pdb.set_trace()


    def analyze_sample(self):
        # Create  blast db with sample file
        self.sample_blast = taranis.blast.Blast("nucl")
        _ = self.sample_blast.create_blastdb(self.sample_file, self.blast_dir)
        result = {}
        # pdb.set_trace()
        for ref_allele in self.ref_alleles:
            # schema_alleles = os.path.join(self.schema, ref_allele)
            # parallel in all CPUs in cluster node
            result[ref_allele] = self.search_alleles(ref_allele)

        pdb.set_trace()
        return

