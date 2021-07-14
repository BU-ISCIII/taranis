# Taranis

- [Introduction](#introduction)
- [Dependencies](#dependencies)
- [Installation](#installation)
  - [Install from source](#install-from-source)
  - [Install using conda](#install-using-conda)
- [Quick usage](#quick-usage)
- [Usage](#usage)
- [Output](#output)
- [Illustrated pipeline](#illustrated-pipeline)



## Introduction

**Taranis** is a computational stand-alone pipeline for **gene-by-gene allele calling analysis** based on BLASTn using  whole genome (wg) and core genome (cg) multilocus sequence typing (MLST) schemas on complete or draft genomes resulting from de novo assemblers, while tracking helpful and informative data among the process.

Taranis includes four main functionalities: MLST **schema analysis**, gene-by-gene **allele calling**, **reference alleles** obtainment for allele calling analysis and the final **distance matrix** construction.



## Dependencies

* Python >=3.6
* NCBI_blast >= v2.9
* prokka >=1.14
* prodigal v2.6.3
* mash >=2
* biopython v1.72
* pandas v1.2.4
* progressbar v2.5
* openpyxl v3.0.7
* plotly v5.0.0
* numpy v1.20.3



## Installation

#### Install from source

Install all dependencies and add them to $PATH.

`git clone https://github.com/BU-ISCIII/taranis.git`

Add taranis and ./bin to $PATH.


#### Install using conda

This option is recomended.

Install Anaconda3.

`conda install -c conda-forge -c bioconda -c defaults taranis`

Wait for the environment to solve. <br>
Ignore warnings/errors.



## Quick usage

- **analyze_schema mode:**

  Schema analysis:

```
taranis analyze_schema \
-inputdir schema_dir \
-outputdir YYYY-MM-DD_taranis_analyze_schema_dir
```

  Schema analysis and duplicated alleles, alleles subsequences and no CDS alleles filtering:

```
taranis analyze_schema \
-inputdir schema_dir \
-outputdir YYYY-MM-DD_taranis_analyze_schema_dir \
-removesubsets True \
-removeduplicates True \
-removenocds True
```


- **reference_alleles mode:**

  Get reference alleles:

```
taranis reference_alleles \
-coregenedir schema_dir \
-outputdir YYYY-MM-DD_taranis_reference_alleles_dir
```


- **allele_calling mode:**

  Run allele calling:

```
taranis allele_calling \
-coregenedir schema_dir \
-refalleles YYYY-MM-DD_taranis_reference_alleles_dir \
-inputdir samples_dir \
-refgenome reference_genome.fasta \
-outputdir YYYY-MM-DD_taranis_allele_calling_dir
```

  Run allele calling getting ST profile:

```
taranis allele_calling \
-coregenedir schema_dir \
-refalleles YYYY-MM-DD_taranis_reference_alleles_dir \
-inputdir samples_dir \
-refgenome reference_genome.fasta \
-profile profile.csv \
-outputdir YYYY-MM-DD_taranis_allele_calling_dir
```

- **distance_matrix mode:**

  Get distance matrix:

```
taranis distance_matrix \
-alleles_matrix YYYY-MM-DD_taranis_allele_calling_dir/result.tsv -outputdir YYYY-MM-DD_taranis_distance_matrix_dir
```

  <p>Get distance matrix filtering loci and samples which missing values percentage is above specified threshold:

```
taranis distance_matrix\
-alleles_matrix YYYY-MM-DD_taranis_allele_calling_dir/result.tsv\
-locus_missing_threshold 20 \
-sample_missing_threshold 50 \
-outputdir YYYY-MM-DD_taranis_distance_matrix_dir
```



## Usage

- **analyze_schema mode:**

```
usage: taranis.py analyze_schema [-h] -inputdir INPUTDIR -outputdir OUTPUTDIR [-removesubsets REMOVESUBSETS] [-removeduplicates REMOVEDUPLICATES] [-removenocds REMOVENOCDS] [-newschema NEWSCHEMA]
                                 [-genus GENUS] [-species SPECIES] [-usegenus USEGENUS] [-cpus CPUS]

optional arguments:
  -h, --help            show this help message and exit
  -inputdir INPUTDIR    Directory where are the schema files.
  -outputdir OUTPUTDIR  Directory where the result files will be stored.
  -removesubsets REMOVESUBSETS
                        Remove allele subsequences from the schema.True: Remove subsets.False: Do not remove subsets.Default is False.
  -removeduplicates REMOVEDUPLICATES
                        Remove duplicated alleles from the schema.True: Remove duplicates.False: Do not remove duplicates.Default is False.
  -removenocds REMOVENOCDS
                        Remove no CDS alleles from the schema.True: Remove no CDS alleles.False: Do not remove no CDS alleles.Default is False.
  -newschema NEWSCHEMA  Filter a copy of the core genes schema preserving the analysis core genes schema.True: Create a copy of the core genes schema for filtering.False: Do not create a copy of the
                        core genes schema for filtering.Default is False.
  -genus GENUS          Genus name for Prokka schema genes annotation. Default is Genus.
  -species SPECIES      Species name for Prokka schema genes annotation. Default is species.
  -usegenus USEGENUS    Use genus-specific BLAST databases for Prokka schema genes annotation (needs --genus). Default is False.
  -cpus CPUS            Number of CPUS to be used in the program. Default is 1.
```


- **reference_alleles mode:**

```
usage: taranis.py reference_alleles [-h] -coregenedir COREGENEDIR -outputdir OUTPUTDIR
                                    [-evalue EVALUE] [-perc_identity PERC_IDENTITY]
                                    [-reward REWARD] [-penalty PENALTY] [-gapopen GAPOPEN]
                                    [-gapextend GAPEXTEND] [-num_threads NUM_THREADS] [-cpus CPUS]

optional arguments:
  -h, --help            show this help message and exit
  -coregenedir COREGENEDIR
                        Directory where the core gene files are located.
  -outputdir OUTPUTDIR  Directory where the result files will be stored.
  -evalue EVALUE        E-value in BLAST searches. Default is 0.001.
  -perc_identity PERC_IDENTITY
                        Identity percent in BLAST searches. Default is 90.
  -reward REWARD        Match reward in BLAST searches. Default is 1.
  -penalty PENALTY      Mismatch penalty in BLAST searches. Default is -2.
  -gapopen GAPOPEN      Gap open penalty in BLAST searches. Default is 1.
  -gapextend GAPEXTEND  Gap extension penalty in BLAST searches. Default is 1.
  -num_threads NUM_THREADS
                        num_threads in BLAST searches. Default is 1.
  -cpus CPUS            Number of CPUS to be used in the program. Default is 1.
```


- **allele_calling mode:**

```
usage: taranis.py allele_calling [-h] -coregenedir COREGENEDIR -refalleles REFALLELES -inputdir
                                 INPUTDIR -refgenome REFGENOME -outputdir OUTPUTDIR
                                 [-percentlength PERCENTLENGTH] [-coverage COVERAGE]
                                 [-evalue EVALUE] [-perc_identity_ref PERC_IDENTITY_REF]
                                 [-perc_identity_loc PERC_IDENTITY_LOC] [-reward REWARD]
                                 [-penalty PENALTY] [-gapopen GAPOPEN] [-gapextend GAPEXTEND]
                                 [-max_target_seqs MAX_TARGET_SEQS] [-max_hsps MAX_HSPS]
                                 [-num_threads NUM_THREADS] [-flankingnts FLANKINGNTS]
                                 [-updateschema UPDATESCHEMA] [-profile PROFILE]
                                 [-updateprofile UPDATEPROFILE] [-cpus CPUS] [-genus GENUS]
                                 [-species SPECIES] [-usegenus USEGENUS]

optional arguments:
  -h, --help            show this help message and exit
  -coregenedir COREGENEDIR
                        Directory where the core gene files are located
  -refalleles REFALLELES
                        Directory where the core gene references files are located
  -inputdir INPUTDIR    Directory where are located the sample fasta files
  -refgenome REFGENOME  Reference genome file for genes prediction
  -outputdir OUTPUTDIR  Directory where the result files will be stored
  -percentlength PERCENTLENGTH
                        Allowed length percentage to be considered as INF. Outside of this limit it
                        is considered as ASM or ALM. Default is SD.
  -coverage COVERAGE    Coverage threshold to exclude found sequences. Outside of this limit it is
                        considered LNF. Default is 50.
  -evalue EVALUE        E-value in BLAST searches. Default is 0.001.
  -perc_identity_ref PERC_IDENTITY_REF
                        Identity percentage in BLAST searches using reference alleles for each
                        locus detection in samples. Default is 90.
  -perc_identity_loc PERC_IDENTITY_LOC
                        Identity percentage in BLAST searches using all alleles in each locus for
                        allele identification in samples. Default is 90.
  -reward REWARD        Match reward in BLAST searches. Default is 1.
  -penalty PENALTY      Mismatch penalty in BLAST searches. Default is -2.
  -gapopen GAPOPEN      Gap open penalty in BLAST searches. Default is 1.
  -gapextend GAPEXTEND  Gap extension penalty in BLAST searches. Default is 1.
  -max_target_seqs MAX_TARGET_SEQS
                        max_target_seqs in BLAST searches. Default is 10.
  -max_hsps MAX_HSPS    max_hsps in BLAST searches. Default is 10.
  -num_threads NUM_THREADS
                        num_threads in BLAST searches. Default is 1.
  -flankingnts FLANKINGNTS
                        Number of flanking nucleotides to add to each BLAST result obtained after
                        locus detection in sample using reference allele for correct allele
                        identification. Default is 100.
  -updateschema UPDATESCHEMA
                        Add INF alleles found for each locus to the core genes schema. True: Add
                        INF alleles to the analysis core genes schema. New: Add INF alleles to a
                        copy of the core genes schema preserving the analysis core genes schema.
                        False: Do not update the core gene schema adding new INF alleles found.
                        Default is True.
  -profile PROFILE      ST profile file based on core genes schema file to get ST for each sample.
                        Default is empty and Taranis does not calculate samples ST.
  -updateprofile UPDATEPROFILE
                        Add new ST profiles found to the ST profile file. True: Add new ST profiles
                        to the analysis ST profile file. New: Add Add new ST profiles to a copy of
                        the ST profile file preserving the analysis ST file. False: Do not update
                        the ST profile file adding new ST profiles found. Default is True.
  -cpus CPUS            Number of CPUS to be used in the program. Default is 1.
  -genus GENUS          Genus name for Prokka schema genes annotation. Default is Genus.
  -species SPECIES      Species name for Prokka schema genes annotation. Default is species.
  -usegenus USEGENUS    Use genus-specific BLAST databases for Prokka schema genes annotation
                        (needs --genus). Default is False.
```


- **distance_matrix mode:**

```
usage: taranis.py distance_matrix [-h] -alleles_matrix ALLELES_MATRIX [-locus_missing_threshold LOCUS_MISSING_THRESHOLD] [-sample_missing_threshold SAMPLE_MISSING_THRESHOLD]
                                  [-paralog_filter PARALOG_FILTER] [-lnf_filter LNF_FILTER] [-plot_filter PLOT_FILTER] -outputdir OUTPUTDIR

optional arguments:
  -h, --help            show this help message and exit
  -alleles_matrix ALLELES_MATRIX
                        Alleles matrix file from which to obtain distances between samples
  -locus_missing_threshold LOCUS_MISSING_THRESHOLD
                        Missing values percentage threshold above which loci are excluded for distance matrix creation. Default is 100.
  -sample_missing_threshold SAMPLE_MISSING_THRESHOLD
                        Missing values percentage threshold above which samples are excluded for distance matrix creation. Default is 100.
  -paralog_filter PARALOG_FILTER
                        Consider paralog tags (NIPH, NIPHEM) as missing values. Default is True
  -lnf_filter LNF_FILTER
                        Consider locus not found tag (LNF) as missing value. Default is True
  -plot_filter PLOT_FILTER
                        Consider incomplete alleles found on the tip of a contig tag (PLOT) as missing value. Default is True
  -outputdir OUTPUTDIR  Directory where the result files will be stored
```



## Output

- **analyze_schema mode:**

  * **FOLDERS:**
  
    * **raw_schema_information:**  General information about each allele of each locus

  * **FILES:**
  
    * **alleles_subsets.tsv:** Report of alleles that are subsequences of other alleles of the same locus
    * **duplicated_alleles.tsv:** Report of duplicate alleles within the same locus
    * **length_statistics.tsv:** Allele length statistics report for each locus
    * **schema_quality.tsv:** Quality report of alleles of each locus


- **reference_alleles mode:**

  * **FILES:**
  
    * **[refalleles_locusX].fasta:** One fasta file for each schema locus containing reference alleles for that locus


- **allele_calling mode:**

  * **FOLDERS:**
    * **alignments:** Nucleotide alignment between sequence found in the sample and allele
    * **proteins:** Protein alignment between sequence found in sample and allele
    * **plots:** Interactive pie charts of allele call results for each sample

  * **FILES:**
    * **alm.tsv:** Sample sequences found x% larger than the locus alleles mean length report
    * **asm.tsv:** Sample sequences found x% shorter than the locus alleles mean length report
    * **exact.tsv:** Exact matches report
    * **inferred_alleles.tsv:** New inferred alleles report
    * **lnf_tpr.tsv:** Locus not found (LNF) and truncated protein (TPR) report
    * **paralog.tsv:** Possible paralogs (NIPHEM (100% ID paralogs) and NIPH (<=100% ID paralogs)) report
    * **plot.tsv:** Possible loci on the tip of the sample contig (PLOT) report
    * **snp.tsv:** SNPs report
    * **matching_contigs.tsv:** Summary report of loci found in samples
    * **result.tsv:** Allele calling main results
    * **summary_result.tsv:** Allele calling results summary. Count of each tag type found for each sample is indicated
    * **stprofile.tsv:** Sequence type report 


- **distance_matrix mode:**

  * **FILES:**
    * **filtered_result.tsv:** Filtered allele calling matrix filtered
    * **matrix_distance.tsv:** Samples matrix distance
    * **matrix_distance_filter_report.tsv:** Allele calling matrix filtering report
    
 
 
## Illustrated pipeline

Under construction
