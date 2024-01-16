import logging
import pandas as pd
import os
import rich.console
import statistics
from pathlib import Path
import Bio.Data.CodonTable

from Bio import SeqIO

# from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
from typing import Self

import taranis.utils


log = logging.getLogger(__name__)
stderr = rich.console.Console(
    stderr=True,
    style="dim",
    highlight=False,
    force_terminal=taranis.utils.rich_force_colors(),
)


class AnalyzeSchema:
    def __init__(
        self: Self,
        schema_allele: str,
        output: str,
        remove_subset: bool,
        remove_duplicated: bool,
        remove_no_cds: bool,
        genus: str,
        species: str,
        usegenus: str,
    ) -> "AnalyzeSchema":
        """AnalyzeSchema instance creation

        Args:
            self (Self): Self
            schema_allele (str): Folder path where schema files are located
            output (str): Out folder to save result
            remove_subset (bool): Remove subset sequences if True
            remove_duplicated (bool): Remove duplicated sequences if True
            remove_no_cds (bool): Removing non coding sequences if True
            genus (str): Genus name for Prokka schema genes annotation
            species (str): Species name for Prokka schema genes annotation
            usegenus (str): genus-specific BLAST databases for Prokka

        Returns:
            AnalyzeSchema: Instance of the created class
        """
        self.schema_allele = schema_allele
        self.allele_name = Path(self.schema_allele).stem
        self.output = output
        self.remove_subset = remove_subset
        self.remove_duplicated = remove_duplicated
        self.remove_no_cds = remove_no_cds
        self.genus = genus
        self.species = species
        self.usegenus = usegenus

    def check_allele_quality(self: Self, prokka_annotation) -> OrderedDict:
        a_quality = OrderedDict()
        allele_seq = {}
        bad_quality_record = []
        with open(self.schema_allele) as fh:
            for record in SeqIO.parse(fh, "fasta"):
                try:
                    prokka_ann = prokka_annotation[record.id]
                except Exception:
                    prokka_ann = "Not found in prokka"
                a_quality[record.id] = {
                    "allele_name": self.allele_name,
                    "quality": "Good quality",
                    "reason": "-",
                    "direction": "forward",
                    "start_codon_alt": "standard",
                    "protein_seq": "",
                    "cds_coding": prokka_ann,
                }
                allele_seq[record.id] = str(record.seq)
                a_quality[record.id]["length"] = str(len(str(record.seq)))
                a_quality[record.id]["dna_seq"] = str(record.seq)
                sequence_direction = taranis.utils.get_seq_direction(str(record.seq))

                if sequence_direction == "reverse":
                    record.seq = record.seq.reverse_complement()
                    a_quality[record.id]["direction"] = sequence_direction
                elif sequence_direction == "Error":
                    a_quality[record.id]["direction"] = "-"
                try:
                    a_quality[record.id]["protein_seq"] = str(
                        record.seq.translate(table=1, cds=True)
                    )

                except Bio.Data.CodonTable.TranslationError as e:
                    if "not a start codon" in str(e):
                        try:
                            # Check if sequence has an alternative start codon
                            # for protein coding
                            a_quality[record.id]["protein_seq"] = str(
                                record.seq.translate(table=2, cds=True)
                            )
                            a_quality[record.id]["start_codon_alt"] = "alternative"
                        except Bio.Data.CodonTable.TranslationError as e_2:
                            if "stop" in str(e_2):
                                a_quality[record.id]["reason"] = str(e_2).replace(
                                    "'", ""
                                )
                            else:
                                a_quality[record.id]["reason"] = str(e).replace("'", "")
                            a_quality[record.id]["quality"] = "Bad quality"
                    else:
                        a_quality[record.id]["quality"] = "Bad quality"
                        a_quality[record.id]["reason"] = str(e).replace("'", "")

                if (
                    self.remove_no_cds
                    and a_quality[record.id]["quality"] == "Bad quality"
                ):
                    bad_quality_record.append(record.id)

        # check if there are duplicated alleles
        # get the unique sequences and compare the length with all sequences
        unique_seq = list(set(list(allele_seq.values())))
        if len(unique_seq) < len(allele_seq):
            tmp_dict = {}
            for rec_id, seq_value in allele_seq.items():
                if seq_value not in tmp_dict:
                    tmp_dict[seq_value] = 0
                else:
                    a_quality[rec_id]["quality"] = "Bad quality"
                    a_quality[rec_id]["reason"] = "Duplicate allele"
                    if self.remove_duplicated:
                        bad_quality_record.append(rec_id)

        for rec_id, seq_value in allele_seq.items():
            unique_seq.remove(seq_value)
            if seq_value in unique_seq:
                a_quality[rec_id]["quality"] = "Bad quality"
                a_quality[rec_id]["reason"] = "Sub set allele"
                if self.remove_subset:
                    bad_quality_record.append(rec_id)

        new_schema_folder = os.path.join(self.output, "new_schema")
        _ = taranis.utils.create_new_folder(new_schema_folder)
        new_schema_file = os.path.join(new_schema_folder, self.allele_name + ".fasta")
        with open(self.schema_allele, "r") as _:
            with open(new_schema_file, "w") as fo:
                for record in SeqIO.parse(self.schema_allele, "fasta"):
                    if len(bad_quality_record) > 0:
                        if record.id not in bad_quality_record:
                            SeqIO.write(record, fo, "fasta")
                    else:
                        SeqIO.write(record, fo, "fasta")
        # update the schema allele with the new file
        self.schema_allele = new_schema_file

        """
        if self.output_allele_annot:
            # dump allele annotation to file
            ann_heading = ["gene", "allele", "allele direction","nucleotide sequence", "protein sequence", "nucleotide sequence length", "star codon", "CDS coding", "allele quality", "bad quality reason" ]
            ann_fields = ["direction", "dna_seq", "protein_seq", "length", "start_codon_alt","cds_coding", "quality", "reason"]
            f_name = os.path.join(self.output, self.allele_name +"_allele_annotation.csv")
            with open (f_name, "w") as fo:
                fo.write(",".join(ann_heading) + "\n")
                for allele in a_quality.keys():
                    data_field = [a_quality[allele][field] for field in ann_fields]
                    fo.write(self.allele_name + "," + allele + "," + ",".join(data_field) + "\n")
        """

        return a_quality

    def fetch_statistics_from_alleles(self, a_quality):
        # POSIBLE_BAD_QUALITY = ["not a start codon", "not a stop codon", "Extra in frame stop codon", "is not a multiple of three", "Duplicate allele", "Sub set allele"]
        record_data = {}
        bad_quality_reason = {}
        a_length = []
        bad_quality_counter = 0
        for record_id in a_quality.keys():
            record_data["allele_name"] = self.allele_name
            a_length.append(int(a_quality[record_id]["length"]))
            if a_quality[record_id]["quality"] == "Bad quality":
                bad_quality_counter += 1
                for reason in taranis.utils.POSIBLE_BAD_QUALITY:
                    if reason in a_quality[record_id]["reason"]:
                        bad_quality_reason[reason] = (
                            bad_quality_reason.get(reason, 0) + 1
                        )
        total_alleles = len(a_length)
        record_data["min_length"] = min(a_length)
        record_data["max_length"] = max(a_length)
        record_data["num_alleles"] = total_alleles
        record_data["mean_length"] = round(statistics.mean(a_length), 2)
        record_data["good_percent"] = round(
            100 * (total_alleles - bad_quality_counter) / total_alleles, 2
        )
        for item in taranis.utils.POSIBLE_BAD_QUALITY:
            record_data[item] = (
                bad_quality_reason[item] if item in bad_quality_reason else 0
            )

        return record_data

    def analyze_allele_in_schema(self):
        allele_data = {}
        # run annotations
        prokka_folder = os.path.join(self.output, "prokka", self.allele_name)
        anotation_files = taranis.utils.create_annotation_files(
            self.schema_allele, prokka_folder, self.allele_name
        )
        prokka_annotation = taranis.utils.read_annotation_file(anotation_files + ".gff")

        # Perform quality
        a_quality = self.check_allele_quality(prokka_annotation)
        allele_data = self.fetch_statistics_from_alleles(a_quality)
        return [allele_data, a_quality]


def parallel_execution(
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


def collect_statistics(data, out_folder, output_allele_annot):
    def stats_graphics(stats_folder):
        allele_range = [0, 300, 600, 1000, 1500]

        graphic_folder = os.path.join(stats_folder, "graphics")
        _ = taranis.utils.create_new_folder(graphic_folder)
        # create graphic for alleles/number of genes
        # genes_alleles_df = stats_df["num_alleles"].value_counts().rename_axis("alleles").sort_index().reset_index(name="genes")
        group_alleles_df = stats_df.groupby(
            pd.cut(stats_df["num_alleles"], allele_range)
        ).count()
        _ = taranis.utils.create_graphic(
            graphic_folder,
            "num_genes_per_allele.png",
            "bar",
            allele_range[1:],
            group_alleles_df["num_alleles"].to_list(),
            ["Allele", "number of genes"],
            "Number of alleles per gene",
        )
        # _ = taranis.utils.create_graphic(graphic_folder, "num_genes_per_allele.png", "lines", genes_alleles_df["alleles"].to_list(), genes_alleles_df["genes"].to_list(), ["Allele", "number of genes"],"title")
        # create pie graph for good quality

        """
        good_percent = [round(stats_df["good_percent"].mean(),2)]
        good_percent.append(100 - good_percent[0])
        labels = ["Good quality", "Bad quality"]
        # pdb.set_trace()
        _ = taranis.utils.create_graphic(graphic_folder, "quality_of_locus.png", "pie", good_percent, "", labels, "Quality of locus")
        # create pie graph for bad quality reason. This is optional if there are
        # bad quality alleles
        """
        sum_all_alleles = stats_df["num_alleles"].sum()

        labels = []
        values = []
        for item in taranis.utils.POSIBLE_BAD_QUALITY:
            labels.append(item)
            values.append(stats_df[item].sum())
        labels.append("Good quality")
        values.append(sum_all_alleles - sum(values))
        _ = taranis.utils.create_graphic(
            graphic_folder,
            "quality_percent.png",
            "pie",
            values,
            "",
            labels,
            "Quality percent",
        )
        # create box plot for allele length variability
        _ = taranis.utils.create_graphic(
            graphic_folder,
            "allele_variability.png",
            "box",
            "",
            stats_df["mean_length"].to_list(),
            "",
            "Allele length variability",
        )

    summary_data = []
    a_quality = []
    for idx in range(len(data)):
        # pdb.set_trace()
        summary_data.append(data[idx][0])
        a_quality.append(data[idx][1])

    stats_df = pd.DataFrame(summary_data)
    # a_quality = data[1]
    stats_folder = os.path.join(out_folder, "statistics")
    _ = taranis.utils.create_new_folder(stats_folder)
    _ = taranis.utils.write_data_to_file(stats_folder, "statistics.csv", stats_df)
    # pdb.set_trace()
    stats_graphics(stats_folder)

    if output_allele_annot:
        # dump allele annotation to file
        ann_heading = [
            "gene",
            "allele",
            "allele direction",
            "nucleotide sequence",
            "protein sequence",
            "nucleotide sequence length",
            "star codon",
            "CDS coding",
            "allele quality",
            "bad quality reason",
        ]
        ann_fields = [
            "direction",
            "dna_seq",
            "protein_seq",
            "length",
            "start_codon_alt",
            "cds_coding",
            "quality",
            "reason",
        ]
        # f_name = os.path.join(self.output, self.allele_name +"_allele_annotation.csv")
        ann_data = ",".join(ann_heading) + "\n"
        for gene in a_quality:
            for allele in gene.keys():
                data_field = [gene[allele][field] for field in ann_fields]
                ann_data += (
                    gene[allele]["allele_name"]
                    + ","
                    + allele
                    + ","
                    + ",".join(data_field)
                    + "\n"
                )

        _ = taranis.utils.write_data_to_compress_filed(
            out_folder, "allele_annotation.csv", ann_data
        )
    return
