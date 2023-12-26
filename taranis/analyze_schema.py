import logging
import pandas as pd
import os
import rich.console
import statistics
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import taranis.utils

import pdb

log = logging.getLogger(__name__)
stderr = rich.console.Console(
    stderr=True,
    style="dim",
    highlight=False,
    force_terminal=taranis.utils.rich_force_colors(),
)


class AnalyzeSchema:
    def __init__(
        self,
        schema_allele,
        output,
        remove_subset,
        remove_duplicated,
        remove_no_cds,
        genus,
        species,
        usegenus,
    ):
        self.schema_allele = schema_allele
        self.allele_name = Path(self.schema_allele).stem
        self.output = output
        self.remove_subset = remove_subset
        self.remove_duplicated = remove_duplicated
        self.remove_no_cds = remove_no_cds
        self.genus = genus
        self.species = species
        self.usegenus = usegenus

    def check_allele_quality(self):
        a_quality = {}
        allele_seq = {}
        bad_quality_record = []
        with open(self.schema_allele) as fh:
            for record in SeqIO.parse(self.schema_allele, "fasta"):
                a_quality[record.id] = {"quality": "Good quality", "reason": "-"}
                allele_seq[record.id] = str(record.seq)
                a_quality[record.id]["length"] = len(str(record.seq))
                if len(record.seq) % 3 != 0:
                    a_quality[record.id]["quality"] = "Bad quality"
                    a_quality[record.id]["reason"] = "Can not be converted to protein"
                    a_quality[record.id]["order"] = "-"
                    continue
                sequence_order = taranis.utils.check_sequence_order(str(record.seq))
                if sequence_order == "Error":
                    a_quality[record.id]["quality"] = "Bad quality"
                    a_quality[record.id]["reason"] = "Start or end codon not found"
                    a_quality[record.id]["order"] = "-"
                    continue
                elif sequence_order == "reverse":
                    record_sequence = str(record.seq.reverse_complement())
                else:
                    record_sequence = str(record.seq)
                a_quality[record.id]["order"] = sequence_order
                if record_sequence[0:3] not in taranis.utils.START_CODON_FORWARD:
                    a_quality[record.id]["quality"] = "Bad quality"
                    a_quality[record.id]["reason"] = "Start codon not found"
                    continue
                if record_sequence[-3:] not in taranis.utils.STOP_CODON_FORWARD:
                    a_quality[record.id]["quality"] = "Bad quality"
                    a_quality[record.id]["reason"] = "Stop codon not found"
                    continue
                if taranis.utils.find_multiple_stop_codons(record_sequence):
                    a_quality[record.id]["quality"] = "Bad quality"
                    a_quality[record.id]["reason"] = "Multiple stop codons found"
                    continue
                if (
                    self.remove_no_cds
                    and a_quality[record.id]["quality"] == "Bad quality"
                ):
                    bad_quality_record.append(record.id)

        if self.remove_duplicated:
            # get the unique sequences and compare the length with all sequences
            unique_seq = list(set(list(allele_seq.values())))
            if len(unique_seq) < len(allele_seq):
                tmp_dict = {}
                for rec_id, seq_value in allele_seq.items():
                    if seq_value not in tmp_dict:
                        tmp_dict[seq_value] = 0
                    else:
                        bad_quality_record.append(rec_id)
        if self.remove_subset:
            unique_seq = list(set(list(allele_seq.values())))
            for rec_id, seq_value in allele_seq.items():
                unique_seq.remove(seq_value)
                if seq_value in unique_seq:
                    bad_quality_record.append(rec_id)
        new_schema_folder = os.path.join(self.output, "new_schema")
        _ = taranis.utils.create_new_folder(new_schema_folder)
        new_schema_file = os.path.join(new_schema_folder, self.allele_name + ".fasta")
        with open(self.schema_allele, "r") as fh:
            with open(new_schema_file, "w") as fo:
                for record in SeqIO.parse(self.schema_allele, "fasta"):
                    if len(bad_quality_record) > 0:
                        if record.id not in bad_quality_record:
                            SeqIO.write(record, fo, "fasta")
                    else:
                        SeqIO.write(record, fo, "fasta")
        # update the schema allele with the new file
        self.schema_allele = new_schema_file
        return a_quality

    def fetch_statistics_from_alleles(self, a_quality):
        record_data = {}
        bad_quality_reason = {}
        a_length = []
        bad_quality_counter = 0
        for record_id in a_quality.keys():
            record_data["allele_name"] = self.allele_name
            a_length.append(a_quality[record_id]["length"])
            if a_quality[record_id]["quality"] == "Bad quality":
                bad_quality_counter += 1
            bad_quality_reason[a_quality[record_id]["reason"]] = (
                bad_quality_reason.get(a_quality[record_id]["reason"], 0) + 1
            )
        total_alleles = len(a_length)
        record_data["min_length"] = min(a_length)
        record_data["max_length"] = max(a_length)
        record_data["num_alleles"] = total_alleles
        record_data["mean_length"] = round(statistics.mean(a_length), 2)
        record_data["good_percent"] = round(
            100 * (total_alleles - bad_quality_counter) / total_alleles, 2
        )
        return record_data

    def analyze_allele_in_schema(self):
        allele_data = {}
        # Perform quality
        a_quality = self.check_allele_quality()
        # run annotations
        prokka_folder = os.path.join(self.output, "prokka", self.allele_name)
        anotation_files = taranis.utils.create_annotation_files(
            self.schema_allele, prokka_folder, self.allele_name
        )
        allele_data["annotation_gene"] = taranis.utils.read_annotation_file(
            anotation_files + ".tsv", self.allele_name
        ).get(self.allele_name)
        allele_data.update(self.fetch_statistics_from_alleles(a_quality))
        return allele_data


def prueba_paralelizacion(
    schema_allele,
    output,
    remove_subset,
    remove_duplicated,
    remove_no_cds,
    genus,
    species,
    usegenus,
):
    schema_obj = AnalyzeSchema(
        schema_allele,
        output,
        remove_subset,
        remove_duplicated,
        remove_no_cds,
        genus,
        species,
        usegenus,
    )
    return schema_obj.analyze_allele_in_schema()


def collect_statistics(stat_data):
    stats_df = pd.DataFrame(stat_data)
    print(stats_df)
