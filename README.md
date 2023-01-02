<p align="center"><img src="Paratype.png" alt="paratype" width="400"></p>

# Paratype
**Assigns genotypes to _Salmonella_ Paratyphi A isolates using whole-genome sequencing data.**

### Introduction
This script is a beta version of the **paratype** tool that assigns genotypes to *Salmonella* Paratyphi A genomes. It is written in _**python3**_. 
**Paratype** also detects mutations in the quinolone resistance-determining regions (_gyrA_-83, _gyrA_-87, _parC_-80, _parC_-84) responsible for resistance to ciprofloxacin and _acrB_ gene (_acrB_-717) which can cause azithromycin resistance in _Salmonella_ Typhi and Paratyphi ([Hooda Y et al. 2019](https://doi.org/10.1371/journal.pntd.0007868), [Sajib MSI et al. 2021](https://doi.org/10.1128/mBio.03481-20)).

[Tanmoy AM et al. 2022](https://doi.org/10.1038/s41467-022-35587-6) described the design of the **Paratype** genotyping scheme. 
Inspiration to design such a scheme came from [genotyphi](https://github.com/katholt/genotyphi), a tool that has been used for genotyping of a related serovar, _Salmonella_ Typhi.

### Dependencies
Dependencies are listed below *(tested versions are in the parentheses)*
1. [Python3](https://www.python.org/) (_v3.8.10_)
2. [Biopython](https://biopython.org/wiki/Download) (_v1.79_)
3. [Samtools](https://github.com/samtools/samtools) (_v1.10_ & v1.13)
4. [BCFtools](https://github.com/samtools/bcftools) (_v1.10.2_ & _v1.13_)
5. [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (_v2.3.5.1_) *(required for fastq and fqin modes)*
6. [BWA](http://bio-bwa.sourceforge.net/) (_v0.7.17_) *(required for fasta and nano modes)*

**Python modules:** _os_, _argparse_ 
(Both modules should be present by default. If not, install it using _"pip install modulename"_. Use _"sudo pip install modulename"_ if you require administrative access for installation.)

**Note:** Paratype assumes that all dependencies are already installed in the system, at their default location. User may notice a few warning messages from samtools mpileup (_for options - u, g and I_). Please ignore those messages.


### Input files
Currently, **Paratype** has **six** working modes with **four** different file types: FASTQ, BAM, VCF and FASTA.
It can run with _fastq_ files from both **Illumina** and **Nanopore** platforms.Use of different modes are as follows _(you can use any of them)_:
```
--id Sample --mode bam --bam $folder\/Sample.bam
```
**bam is the default running mode.** It requires a mapped *.bam* file, mapped against the reference **AKU_12601** genome.

```
--id Sample --mode fastq --fastq $folder\/Sample_1.fastq.gz $folder\/Sample_2.fastq.gz
```
**fastq** mode is the **slowest** but also the **most accurate**  mode and requires two paired-end raw illumina read files. 

```
--id Sample --mode fqin --fqin $folder\/Sample.fastq.gz
```
**fqin** mode requires one paired-end interleaved raw illumina read file, but runs exactly the same as the _fastq_ mode.
```
--id Sample --mode nano --nano $folder\/Sample.fastq.gz
```
**nano** mode requires one raw nanopore read file, generated using MinION. This mode uses _bwa mem_ to map the reads. Rest of the codes are exactly the same as the fastq mode. 
```
--id Sample --mode fasta --fasta $folder\/Sample.fasta
```
**fasta** mode requires one assembled fasta file. As this mode is highly dependent on quality of the genome assembly, user should be careful while choosing the assembly program. 
```
--id Sample --mode vcf --vcf $folder\/Sample.vcf
```
**VCF** mode is the **fastest** mode, but also the **least accurate**. Thus, it is not recommended unless you are highly confident about your SNP data. Moreover, the **Paratype** script requires a **VCF** file of genome-wide locations for the referene AKU_12601 genome, not only the SNP-occurring genomic location. 

**Reference sequence has to be **AKU_12601** genome ([NC_011147.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_011147.1)).** Fasta file is provided with the script here. Genbank accession [FM200053.1](https://www.ncbi.nlm.nih.gov/nuccore/FM200053.1) can also be used.

**The reference accession for the BAM and VCF alignments _(--mode bam or vcf)_ must be NC_011147.1.** It is also the default accession for this script, so you do not need to use the *--ref_id* option. If you manually change the accession to *NC_011147_1* or, *NC_011147-1* or, *NC_011147*, please use the option *--ref_id* to provide the changed accession. Also, provide your reference fasta file with changed accession using *--ref* option. For example, if you change the accession to *NC_011147_1*, you can use the following option:

```
--id Sample --ref NC_011147_1.fasta --ref_id NC_011147_1
```


**The paratype script needs the allele definition and gene region (codon) files** _(provided with the script)_. You do not need to use the option, *--allele*  or, *--gene* to provide the files, if you download the folder and keep all the files at the directory of **paratype** script. 
However, if you found a new genotype or want to detect a new mutation, you can either edit the designated text file for allele definitions _(SParatyphiA_genotype_specific_alleles.txt)_ and gene_regions _(SParatyphiA_gene_mutation_codons.txt)_, or you can prepare files with the specific format used here. In that case, please use the following options: 
```
--id Sample --allele new_alele_definition.txt --gene new_gene_codons_definition.txt
```
However, if you add a new genotype to the provided allele_definition file or, make a new one, please follow the numbered_nomenclature we followed here. For example, you can use N.N.N format (e.g. 2.4.1), but not N.N.TEXT format (e.g. 2.4.Ab). _(Otherwise, the script will show errors)_.

### Usage
Paratype assumes that **python3** is the default **_python_** in your system. If it is not, you should use _python3_ instead of _python_ in the following commands.

Provide the sample name or ID using _--id_ option _(mandatory)_. Paratype will use this to name all necessary files. 

If the reference accession is exactly **NC_011147.1**, the use of **_--ref_id_** is not required. If all provided files with the script are in the same folder, the use of **_--allele_**, **_--gene_** and **_--ref_** options are not required either. Use of **_--output_** is also optional. If not used, a file named _paratype_results.txt_ will be generated. 

#### BAM mode _(default)_
Use of *--mode* is not required.
```
python paratype.py --id Sample --bam Sample.bam --output Sample_paratype.txt
```

#### FASTQ and FQIN mode _(Recommended for illumina reads)_
Use of *--threads* is recommended (default: 1). 
```
python paratype.py --id Sample --mode fastq --fastq Sample_R1.fastq.gz Sample_R2.fastq.gz --threads 8 --output Sample_paratype.txt
```
```
python paratype.py --id Sample --mode fqin --fqin Sample.fastq.gz --threads 8 --output Sample_paratype.txt
```
#### NANO mode _(Recommended for nanopore reads)_
Use of *--threads* is recommended (default: 1). 
```
python paratype.py --id Sample --mode nano --nano Sample.fastq.gz Sample_R2.fastq.gz --threads 8 --output Sample_paratype.txt
```
#### FASTA mode _(Recommended for quick and moderately accurate results)_
Use of *--threads* is recommended (default: 1). 
```
python paratype.py --id Sample --mode fasta --fasta Sample.fasta --threads 8 --output Sample_paratype.txt
```

#### VCF mode _(Not recommended unless SNP data is highly trusted)_
```
python paratype.py --id Sample --mode vcf --vcf Sample.vcf --output Sample_paratype.txt
```


### Options and details

```
usage: paratype.py [-h] --id ID [--mode MODE] [--fastq FASTQ [FASTQ ...]] [--fqin FQIN] [--bam BAM] [--vcf VCF] [--fasta FASTA] [--nano NANO] [--ref REF] [--ref_id REF_ID] [--phrd_cutoff PHRD_CUTOFF]
                   [--read_cutoff READ_CUTOFF] [--threads THREADS] [--allele ALLELE] [--genes GENES] [--output OUTPUT]

Genotyping of Salmonella Paratyphi A using fastq or fasta or bam or vcf files, against the strain AKU_12601 as reference.

optional arguments:
  -h, --help            show this help message and exit
  --id ID               Sample ID
  --mode MODE           Mode to run in based on input files (fastq, fastq interleaved, bam, vcf, fasta, and nanopore). Default: bam
  --fastq FASTQ [FASTQ ...]
                        Raw fastq read files (paired-end).
  --fqin FQIN           Raw fastq read files (paired-end interleaved).
  --bam BAM             Mapped BAM file against the AKU_12601 reference genome.
  --vcf VCF             Mapped VCF file against the AKU_12601 reference genome.
  --fasta FASTA         Assembled fasta files (not recommended unless the contigs are highly trusted).
  --nano NANO           Raw nanopore fastq read files.
  --ref REF             Fasta Reference sequence of AKU_12601 (default file is provided with the script)
  --ref_id REF_ID       Reference sequence id (default: NC_011147.1).
  --phrd_cutoff PHRD_CUTOFF
                        Minimum phred quality score to consider a variant call as a true allele (default: 20).
  --read_cutoff READ_CUTOFF
                        Minimum proportion of reads required to call a true allele (default: 0.75).
  --threads THREADS     Number of threads to use for Bowtie2 or bwa mapping (only for "fastq" mode). (default: 1)
  --allele ALLELE       Allele definition in tab-delimited format (default file is provided with the script).
  --genes GENES         List of codons to find targeted gene mutation (tab-delimited format; default file is provided with the script).
  --output OUTPUT       output file.

```

#### Required options
```
  --id  ID          Sample ID
  --mode    MODE    Mode to run in based on input files (fastq, fastq interleaved, bam, vcf, fasta, and nanopore). Default: bam
```


#### Mode-specific options

##### --mode bam
Requires [SAMtools](http://samtools.sourceforge.net/) and [BCFtools](https://samtools.github.io/bcftools/)
```
--bam   BAM    Mapped BAM file against the AKU_12601 reference genome.
```


##### --mode fastq and fqin
Requires [Bowtie2]( http://bowtie-bio.sourceforge.net/bowtie2/), [SAMtools](http://samtools.sourceforge.net/) and [BCFtools](https://samtools.github.io/bcftools/)
```
--fastq FASTQ   [FASTQ ...]   Raw fastq read files (paired-end).
```
```
--fqin  FQIN    Raw fastq read file (paired-end interleaved).
```

##### --mode nano
Requires [bwa](http://bio-bwa.sourceforge.net/), [SAMtools](http://samtools.sourceforge.net/) and [BCFtools](https://samtools.github.io/bcftools/)
```
--nano  NANO Raw nanopore fastq read files.
```

##### --mode fasta
Requires [bwa](http://bio-bwa.sourceforge.net/), [SAMtools](http://samtools.sourceforge.net/) and [BCFtools](https://samtools.github.io/bcftools/)
```
--fasta FASTA   Assembled fasta files (not recommended unless the contigs are highly trusted).
```

##### --mode vcf
```
--vcf   VCF Mapped VCF file against the AKU_12601 reference genome. 
```


#### Other options
```
  --ref REF             Fasta Reference sequence of AKU_12601 (default file is provided with the script)
  --ref_id REF_ID       Reference sequence id (default: NC_011147.1).
  --phrd_cutoff PHRD_CUTOFF
                        Minimum phred quality score to consider a variant call as a true allele (default: 20).
  --read_cutoff READ_CUTOFF
                        Minimum proportion of reads required to call a true allele (default: 0.75).
  --threads THREADS     Number of threads to use for Bowtie2 mapping (only for "fastq" mode). (default: 1)
  --allele ALLELE       Allele definition in tab-delimited format (default file is provided with the script).
  --genes GENES         List of codons to find targeted gene mutation (tab-delimited format; default file is provided with the script).
  --output OUTPUT       output file.
```
#
### Batch mode
Paratype does not have a batch mode yet. However, a bash script is added here that can be used to run multiple isolate data. All input files need to be in one folder. 
```
bash Batch_run_paratype.sh <Input file directory> <Results directory> <Paratype directory>
```
The script has commands for all six different modes. Please unmute the mode you want to run. By default, the **_bam_** mode is unmuted in the script. 

### Citation
If you use this tool or the scheme, please cite the **paratype** article on _Nature Commnications_([Tanmoy AM et al.](https://doi.org/10.1038/s41467-022-35587-6)).


### Python 2.7 version
Paratype is no longer maintained on python2.7. Please download the original codes (release: *Original paratype codes, v1_beta*), if you have no other option than working with python2.7. 
