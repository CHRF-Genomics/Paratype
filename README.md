<p align="center"><img src="Paratype.png" alt="paratype" width="400"></p>

# Paratype
**Assigns genotypes to _Salmonella_ Paratyphi A isolates using whole-genome sequencing data.**

### Introduction
This script is a beta version of the **paratype** tool that assigns genotypes to *Salmonella* Paratyphi A genomes. It is written in *python2.7*. In addition, **paratype** detects mutations in the quinolone resistance-determining regions (_gyrA_-83, _gyrA_-87, _parC_-80, _parC_-84) responsible for resistance to ciprofloxacin and _acrB_ gene (_acrB_-717) which can cause azithromycin resistance in _Salmonella_ Typhi and Paratyphi ([Hooda Y et al. 2019](https://doi.org/10.1371/journal.pntd.0007868), [Sajib MSI et al. 2021](https://doi.org/10.1128/mBio.03481-20)).

The article that will describe the design of the genotyping scheme will be published soon. Inspiration to write such a script came from [genotyphi](https://github.com/katholt/genotyphi), used for typing isolates from the related serovar _Salmonella_ Typhi.


### Dependencies
Dependencies are listed below *(tested versions are in the parentheses)*
1. python2.7 (_v2.7.18_)
2. samtools (_v1.10_ & v1.13)
3. bcftools (_v1.10.2_ & _v1.13_)
4. bowtie2 (_v2.3.5.1_) *(required for fastq mode only)*

Python libraries: os, argparse 
*(Both libraries should be present by default. If not, install it using "pip install libraryname". Use "sudo pip install libraryname" if you require administrative access for installation.)*

*Note: You may see some warning messages from samtools mpileup (for options - u, g and I). Please ignore those messages.*

**Paratype has been tested with illumina paired-end reads only and it assumes that all dependencies are already installed in the system, at their default location. Also, it does not have a _batch_ running mode yet (under development).**


### Input files
Currently, **paratype** has three working modes with three different file types: FASTQ, BAM, and VCF. FASTA mode will be added in the next version of the script. 
You can use them like this:
```
--mode bam --bam $folder\/Sample.bam
```
or, 
```
--mode fastq --fastq $folder\/Sample_1.fastq.gz $folder\/Sample_2.fastq.gz
```
or, 
```
--mode fastq --fastq $folder\/Sample_*R1*.fastq.gz $folder\/Sample_*R2*.fastq.gz
```
or, 
```
--mode vcf --vcf $folder\/Sample.vcf
```
**BAM is the default and recommended running mode.** It requires a mapped *.bam* file, sorted and preferably indexed. However, if it is not indexed, **paratype** will index the bam file before using it. 

**FASTQ** mode is the slowest mode and requires two paired-end raw read files. Default file extensions are - *Name_1.fastq.gz & Name_2.fastq.gz*. (Please rename your files to these file extensions, or, you can change the number in line 309 of the *paratpe.py* script. For example, “_1.fastq.gz” has 11 characters, so we used “[:-11]” in line 309 of the script. For “_1.fq.gz” (8 characters), you can use “[:-8]” instead. (**paratype** is tested with illumina short-read fastq files.)

**VCF** mode is faster than the other modes, but it is not recommended unless you are highly confident about your SNP data. Also, the script requires a VCF file of genome-wide locations, not only the SNP-occurring genomic location. 


**Both BAM and VCF files have to be mapped against the _S._ Paratyphi A AKU_12601 reference genome, [FM200053.1](https://www.ncbi.nlm.nih.gov/nuccore/FM200053.1).** The reference sequence is provided with this script. You need to use the *--ref* option to refer to this reference fasta file (the one we used is provided with the script).

```
--ref SParatyphiAKU12601.fasta
```


**The reference accession for the BAM and VCF alignments must be FM200053.1.** It is also the default accession for this script, so you do not need to use the *--ref_id* option. If you manually change the accession to *FM200053_1* or, *FM200053-1* or, *FM200053*, please use the option *--ref_id* to provide the changed accession. Also, provide your reference fasta file with changed accession using *--ref* option. For example – if you change the accession to *FM200053_1*, you can use the following options:

```
--ref FM200053_1.fasta --ref_id FM200053_1
```


**The paratype script needs the allele definition and gene region (codon) files** _(provided with the script)_. You do not need to use the option, *--allele*  or, *--gene* to provide the files, if you download the folder and keep all the files at the directory of paratype script. 
However, if you found a new genotype or want to detect a new mutation, you can either edit the designated text file for allele definitions (SParatyphiA_genotype_specific_alleles_v1_b1.txt) and gene_regions (SParatyphiA_gene_mutation_codons_v1_b2.txt), or you can prepare files with specific format. In that case, please use following options: 
```
--allele new_alele_definition.txt --gene new_gene_codons_definition.txt
```
However, if you add a new genotype to the provided allele_definition file or, make a new one, please follow the numbered_nomenclature we followed here. For example, you can use N.N.N format (e.g. 2.4.1), but not N.N.TEXT format (e.g. 2.4.Ab). _(Otherwise, the script will show errors)_.

### Usage
If the reference accession is exactly **FM200053.1**, the use of **_--ref_id_** is not required. If all provided files with the script are in the same folder, the use of **_--allele_**, **_--gene_** and **_--ref_** options are not required either. Use of **_--output_** is also optional. If not used, a file named _paratype_results.txt_ will be created. 

#### Basic Usage – BAM mode (recommended)
Use of *--mode* is not required.
```
python paratype_v1_b2.py --bam Sample.bam --output Sample_paratype.txt
```

#### Basic Usage – FASTQ mode
Use of *--threads* is recommended (default: 1). 
```
python paratype_v1_b2.py --mode fastq --fastq Sample_R1.fastq.gz Sample_R2.fastq.gz --threads 8 --output Sample_paratype.txt
```

#### Basic Usage – VCF mode
```
python paratype_v1_b2.py --mode vcf --vcf Sample.vcf --output Sample_paratype.txt
```


### Options and details

```
usage: paratype_v1_b2.py [-h] [--mode MODE] [--fastq FASTQ [FASTQ ...]]
                         [--bam BAM] [--vcf VCF] [--ref REF] [--ref_id REF_ID]
                         [--phrd_cutoff PHRD_CUTOFF]
                         [--read_cutoff READ_CUTOFF] [--threads THREADS]
                         [--allele ALLELE] [--genes GENES] [--output OUTPUT]

Genotyping of Salmonella Paratyphi A using fastq or bam or vcf files, against
AKU_12601 as reference.

optional arguments:
  -h, --help            show this help message and exit
  --mode MODE           Mode to run in based on input files (fastq or, bam or,
                        vcf)
  --fastq FASTQ [FASTQ ...]
                        Raw fastq read files (paired-end).
  --bam BAM             Mapped BAM file against the AKU_12601 reference
                        genome.
  --vcf VCF             Mapped VCF file against the AKU_12601 reference
                        genome.
  --ref REF             Fasta Reference sequence of AKU_12601 (default file is
                        provided with the script)
  --ref_id REF_ID       Reference sequence id (default: FM200053.1).
  --phrd_cutoff PHRD_CUTOFF
                        Minimum phred quality score to consider a variant call
                        as a true allele (default: 20).
  --read_cutoff READ_CUTOFF
                        Minimum proportion of reads required to call a true
                        allele (default: 0.75).
  --threads THREADS     Number of threads to use for Bowtie2 mapping (only for
                        "fastq" mode). (default: 1)
  --allele ALLELE       Allele definition in tab-delimited format (default
                        file is provided with the script).
  --genes GENES         File for Gene mutation finding (tab-deleimited format;
                        provided with the script).
  --output OUTPUT       output file.

```

#### Required options
```
--mode    Mode to run the script based on the type of input files (fast or, bam or, vcf). BAM is set as the default option. So, this option is not required for BAM mode. 
```


#### Mode-specific options

##### --mode bam
Requires [SAMtools](http://samtools.sourceforge.net/) and [BCFtools](https://samtools.github.io/bcftools/)
```
--bam   Path to one BAM file, generated by mapping reads to the Salmonella Paratyphi A AKU_12601 reference genome (FM200053.1). Note that the SNP coordinates used here for genotyping are specifically relative to Salmonella Paratyphi A AKU_12601 (FM200053.1). So the input MUST be a BAM obtained via mapping to this reference sequence. Also, the BAM file has to be sorted. 
```


##### --mode fastq
Requires [Bowtie2]( http://bowtie-bio.sourceforge.net/bowtie2/), [SAMtools](http://samtools.sourceforge.net/) and [BCFtools](https://samtools.github.io/bcftools/)
```
--fastq	  Path to two paired-end raw fastq read files. The default file extension is ‘fastq.gz’. Rename your files accordingly if needed. If you use another file extension, your output file name can get a bit messy, but it will not affect the results. 
```


##### --mode vcf
```
--vcf	  Path to one VCF file, generated by mapping reads to the Salmonella Paratyphi A AKU_12601 reference genome (FM200053.1). Note that the paratype use alleles, not SNP to define the genotypes. So, it requires a VCF with both alleles (REF & ALT) for all genome coordinates. An alternative way – you can get the genomic coordinates from the allele file (provided with the script) and only extract alleles for those from BAM file using bcftools. 
```


#### Other options
```
--ref_id      Accession of the Paratyphi A AKU_12601 chromosome reference used for mapping.
--phrd_cutoff Minimum Phred quality score to consider a variant call as a true allele. Default is set to 20.
--read_cutoff Minimum proportion of reads required to call a true allele. Default is set to 75% (0.75).
--output      Specify the location of the output text file. Default is set to *Sample_paratype_results.txt*. The output file format is – (tab-separated) StrainID, Primary_clade, Secondary_clade, Subclade, Genotype, Support.
--allele      Allele definition file in tab-delimited format (provided with the script). If you think you found a new genotype, you can add its unique allele location, nucleotide, locus_tag, and the new genotype designation in this file. You can then run the paratype script with that file to know the accuracy of your findings.  
--gene      Gene mutation file with codon location and nucleotide details in tab-delimited format. 
--ref       Fasta Reference sequence of Paratyphi AKU_12601 used for mapping (provided with the script).
--ref_id    Fasta ID of the reference sequence. Default is FM200053.1.
```
