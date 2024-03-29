#!/usr/bin/env python3
'''
This script assigns genotypes to Salmonella Paratyphi A isolates using whole-genome sequencing data.
It uses FASTQ or BAM (recommended) or VCF (if highly trusted SNP data) or FASTA (assembled contigs) files relative to Paratyphi A AKU_12601 (NC_011147.1). Also, it reads allele definitions of the genotypes from a text file (given with this script).

Authors: Arif Mohammad Tanmoy (arif.tanmoy@chrfbd.org) wrote the script. He and Yogesh Hooda (yhooda@chrfbd.org) defined the genotype-specific alleles.

Last modified - 25 May, 2023
Version: 1.1
'''
import os
from argparse import ArgumentParser
from Bio.Seq import Seq


def parse_args():
    "Parse the input arguments, use '-h' for help"
    commands = ArgumentParser(
        description='Genotyping of Salmonella Paratyphi A using fastq or fasta or bam or vcf files, against the strain AKU_12601 as reference.')
    commands.add_argument('--id', type=str, required=True,
                          help='Sample ID')
    commands.add_argument('--mode', required=False, default='bam',
                          help='Mode to run in based on input files (fastq, fastq interleaved, bam, vcf, fasta, and nanopore). Default: bam')
    commands.add_argument('--fastq', nargs='+', required=False,
                          help='Raw fastq read files (paired-end).')
    commands.add_argument('--fqin', type=str, required=False,
                          help='Raw fastq read files (paired-end interleaved).')
    commands.add_argument('--bam', type=str, required=False,
                          help='Mapped BAM file against the AKU_12601 reference genome.')
    commands.add_argument('--vcf', type=str, required=False,
                          help='Mapped VCF file against the AKU_12601 reference genome.')
    commands.add_argument('--fasta', type=str, required=False,
                          help='Assembled fasta files (not recommended unless the contigs are highly trusted).')
    commands.add_argument('--nano', type=str, required=False,
                          help='Raw nanopore fastq read files.')
    commands.add_argument('--ref', type=str, required=False, default='SParatyphiAKU12601.fasta',
                          help='Fasta Reference sequence of AKU_12601 (default file is provided with the script)')
    commands.add_argument('--ref_id', type=str, required=False,	default='NC_011147.1',
                          help='Reference sequence id (default: NC_011147.1).')
    commands.add_argument('--mapq_cutoff', type=float, required=False, default=20,
                          help='Minimum mapping quality (by phred score) to consider a variant call as a true allele (default: 20).')
    commands.add_argument('--phrd_cutoff', type=float, required=False, default=20,
                          help='Minimum base quality (by phred score) to consider a variant call as a true allele (default: 20).')
    commands.add_argument('--read_cutoff', type=float, required=False, default=0.75,
                          help='Minimum proportion of reads required to call a true allele (default: 0.75).')
    commands.add_argument('--threads', type=int, required=False, default=1,
                          help='Number of threads to use for mapping and variant calling. (default: 1).')
    commands.add_argument('--allele', type=str, required=False,	default='SParatyphiA_genotype_specific_alleles.txt',
                          help='Allele definition in tab-delimited format (default file is provided with the script).')
    commands.add_argument('--genes', type=str, required=False,	default='SParatyphiA_gene_mutation_codons.txt',
                          help='List of codons to find targeted gene mutation (tab-delimited format; default file is provided with the script).')
    commands.add_argument('--output', type=str, required=False,
                          default='paratype_results.txt', help='output file.')
    return commands.parse_args()


args = parse_args()

version = "1.1"

# Define allele, gene region and reference files
# Let's define a simple function to check and decide

def define_files(path, filename, argss):
    if (argss == filename) or ("/" not in argss):
        final_file = '/'.join([path, argss])
    else:
        final_file = argss
    return final_file


# the path from this script is running
parapath = os.path.dirname(os.path.realpath(__file__))

genotype_allele_file = define_files(
    parapath, "SParatyphiA_genotype_specific_alleles.txt", args.allele)
ref_fasta_file = define_files(parapath, "SParatyphiAKU12601.fasta", args.ref)
gene_regions_file = define_files(
    parapath, "SParatyphiA_gene_mutation_codons.txt", args.genes)

# Define genotypes

def define_genotypes(allele_file):
    clades, alleles, loci = ([] for l in range(3))
    for line in open(allele_file, 'r'):
        if "#" not in line:
            rec = line.strip('\n').split('\t')
            clades.append(rec[0])
            loci.append(rec[1])
            alleles.append(rec[2])
    temp_bed = open((args.ref_id+'.bed'), 'w')
    for locus in sorted(loci):
        temp_bed.write(
            args.ref_id + '\t' + str(int(locus)-1) + '\t' + locus + '\n')
    temp_bed.close()
    return clades, loci, alleles

# calculate read proportion

def calculate_read_proportion(INFO, FORMAT):
    x = INFO.split('DP4=')[1].split(';')[0].split(',')
    if x != None:
        alt_read_count = float(int(x[2]) + int(x[3]))
        total_read_count = alt_read_count + float(int(x[0]) + int(x[1]))
        # print total_read_count
        # print alt_read_count
        if total_read_count != 0:
            snp_proportion = float(alt_read_count / total_read_count)
        else:
            snp_proportion = float(-1)
    else:
        try:
            ad = FORMAT.split(':')[1].split(',')  # get the AD ratio
            alt_read_count = float(ad[1])
            total_read_count = float(ad[0]) + alt_read_count
            # print total_read_count
            # print alt_read_count
            snp_proportion = float(alt_read_count / total_read_count)
        except IndexError:
            snp_proportion = float(-1)
    return snp_proportion

# check vcf for loci

def check_allele_from_vcf(vcf_file, clades, loci, alleles):
    type_list, propr = [], []
    for line in open(vcf_file, 'r'):
        if not line.startswith('#'):
            CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, FILE = line.rstrip().split()
            if (float(QUAL) > args.phrd_cutoff) and (POS in loci):
                a = loci.index(str(POS))
                if (ALT != '.') and (str(alleles[a]) == ALT):
                    read_ratio = calculate_read_proportion(INFO, FORMAT)
                    #print ('\t'.join([POS, ID, REF, ALT, QUAL, str(read_ratio)])) + '\n'
                    if read_ratio > args.read_cutoff:
                        type_list.append(str(clades[a]))
                        propr.append(str(read_ratio))
                elif (ALT != '.') and (str(alleles[a]) == REF):
                    read_ratio = calculate_read_proportion(INFO, FORMAT)
                    #print ('\t'.join([POS, ID, REF, ALT, QUAL, str(read_ratio)])) + '\n'
                    if read_ratio < (1 - args.read_cutoff):
                        type_list.append(str(clades[a]))
                        propr.append(str(1-read_ratio))
                elif (ALT == '.') and (str(alleles[a]) == REF):
                    # ALT is empty, means all reads are for REF allele.
                    read_ratio = 1.0
                    # print '\t'.join([POS, ID, REF, ALT, QUAL, str(read_ratio)])
                    type_list.append(str(clades[a]))
                    propr.append(str(read_ratio))
    return type_list, propr

# classify type_list inside the length_condition loop

def classify_each_type(type1):
    if type1.endswith('0') or type1.startswith('0'):
        decision = "Primary_Clade"
    elif len(type1) <= 3:
        decision = "Secondary_Clade"
    else:
        decision = "Subclade"
    return decision, type1

# generate printable results (to use in the designate_genotypes function

def generate_print_results(final_type, clades):
    prim, sec, sub, geno = [], [], [], []
    for i in range(0, len(final_type)):
        if "Prim" in str(final_type[i]):
            prim.append(str(clades[i]))
        if "Sec" in str(final_type[i]):
            sec.append(str(clades[i]))
        if "Sub" in str(final_type[i]):
            sub.append(str(clades[i]))

    if not prim:  # check if subclade clade is missing
        prim.append("missing")
    if not sec:  # check if secondary clade is missing
        sec.append("missing")
    if not sub:  # check if subclade is missing
        sub.append("missing")
        if "missing" in sec:
            geno = prim
        else:
            geno = sec
    else:
        geno = sub

    print_result = '\t'.join([str(','.join(prim)), str(
        ','.join(sec)), str(','.join(sub)), str(','.join(geno))])
    return print_result

# calculate list mean

def Average(lst):
    return sum(lst) / len(lst)

# detect the primary, secondary and subclades (genotypes)

def designate_genotypes(type_list, propr):
    final_type, clades = [], []
    propr = [float(i) for i in propr]		# Convert strings to float numbers

    if len(type_list) == 0:
        print("We found no genotype. Please check if your data is truly from Paratyphi A. If you are certain about the serovar, you may have found an new genotype. Please report it in our github page.")
        final_result = "No result"
        ratio = "NA"

    elif len(type_list) == 1:
        typee, name = classify_each_type(str(type_list[0]))
        print(' '.join(['We only found allele for a', typee, '.']))
        final_type.append(typee)
        clades.append(name)
        final_result = generate_print_results(final_type, clades)
        ratio = str(propr[0])

    elif len(type_list) == 2:
        char = (type_list[0].split('.')[0])
        if (str(type_list[0]).startswith(char)) and (str(type_list[-1]).startswith(char)):
            if (str(type_list[0]).endswith('0')) or str(type_list[-1]).endswith('0'):
                for i in range(0, len(type_list)):
                    typee, name = classify_each_type(str(type_list[i]))
                    final_type.append(typee)
                    clades.append(name)
            else:
                print("We found no allele for primary clade.")
                for i in range(0, len(type_list)):
                    typee, name = classify_each_type(str(type_list[i]))
                    final_type.append(typee)
                    clades.append(name)
        else:
            for i in range(0, len(type_list)):
                typee, name = classify_each_type(str(type_list[i]))
                final_type.append(typee)
                clades.append(name)
                print("We found two different genotypes for this isolate. Please check the data for probable contamination. If you are certain about quality, Please report this result in our github page.")
        final_result = generate_print_results(final_type, clades)
        ratio = str(round(Average(propr), 2))

    elif len(type_list) == 3:
        char = type_list[0].split('.')[0]
        if (str(type_list[0]).startswith(char)) and (str(type_list[1]).startswith(char)) and (str(type_list[-1]).startswith(char)):
            for i in range(0, len(type_list)):
                typee, name = classify_each_type(str(type_list[i]))
                final_type.append(typee)
                clades.append(name)
        else:
            print("We found alleles of more than one genotypes for this isolate. Please check the data for probable contamination. If you are certain about quality, please report this result in our github page.")
            for i in range(0, len(type_list)):
                typee, name = classify_each_type(str(type_list[i]))
                final_type.append(typee)
                clades.append(name)
        final_result = generate_print_results(final_type, clades)
        ratio = str(round(Average(propr), 2))

    else:
        char = type_list[0].split('.')[0]
        for i in range(0, len(type_list)):
            typee, name = classify_each_type(str(type_list[i]))
            final_type.append(typee)
            clades.append(name)
        print("We found alleles of at least two genotypes for this isolate. Please check the data for probable contamination. If you are certain about quality, please report this result in our github page.")
        final_result = generate_print_results(final_type, clades)
        ratio = str(round(Average(propr), 2))
    return final_result, ratio

# Define strain ID
def generate_strainID(idd):
    if idd !="":
        strainID = idd
    else:
        strainID = "Strain ID is missing. Please add."
    return strainID

########################
# DETECT GENE MUTATION
# Define mutation position

def define_genes(gene_file):
    mutation, strand, location, codon, base, codon_pos, region = (
        [] for l in range(7))
    for line in open(gene_file, 'r'):
        if not line.startswith('#'):
            rec = line.strip('\n').split('\t')
            mutation.append(rec[0])
            strand.append(rec[1])
            location.append(int(rec[2]))
            #end.append(int(rec[2]))
            codon.append(rec[3])
            base.append(rec[4])
            codon_pos.append(int(rec[5]))
    temp_bed = open((args.ref_id+'.bed'), 'a')
    for x in range(0, len(location)):
        temp_bed.write(
            args.ref_id + '\t' + str(location[x] - 1) + '\t' + str(location[x]) + '\n')
    temp_bed.close()
    return mutation, strand, location, codon, base, codon_pos, region

# check vcf for high-quality SNP for mutations

def select_snp_from_vcf(vcf_file, location):
    snp_pos, snp_alt = [], []
    for line in open(vcf_file, 'r'):
        if not line.startswith('#'):
            CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, FILE = line.rstrip().split()
            if (int(POS) in location) and (float(QUAL) > args.phrd_cutoff) and (ALT != "") and (len(ALT) == 1):
                read_ratio = calculate_read_proportion(INFO, FORMAT)
                if read_ratio > args.read_cutoff:
                    snp_pos.append(POS)
                    snp_alt.append(ALT)
    return snp_pos, snp_alt

# Detect mutations

def classify_mutations(snp_pos, snp_alt, mutation, strand, location, codon, codon_pos):
    final_mut_list = []
    for i in range(0, len(snp_pos)):
        try:
            j = location.index(int(snp_pos[i]))
            triplet = []
            triplet[:0] = codon[j]
            triplet[(codon_pos[j]-1)] = snp_alt[i]
            alt_triplet = ''.join(triplet)
            if strand[j] == 'for':
                alt_amino = str(Seq(alt_triplet).translate(
                    table='11', stop_symbol='*'))
            else:
                alt_amino = str(Seq(alt_triplet).reverse_complement(
                ).translate(table='11', stop_symbol='*'))
            final_mut_list.append(mutation[j] + alt_amino)
        except ValueError:
            pass
    if len(final_mut_list) == 0:
        final_mut_list.append('No_mutation')
    return final_mut_list


# Run Samtools And BCFtools

def run_samtools(ref_fasta_file, mapq_cutoff, ref_id, threads, sam, bam, sorted_bam, vcf_file):
    os.system(
        ' '.join(['samtools view -ubS --threads', str(args.threads), sam, '>', bam]))
    os.system(' '.join(['samtools sort --threads',
                        str(args.threads), bam, '-o', sorted_bam]))
    os.system(' '.join(['samtools index', sorted_bam]))
    os.system(' '.join(['samtools mpileup -q', str(args.mapq_cutoff), '-uBf', ref_fasta_file,
                        '-l', str(args.ref_id + '.bed'), sorted_bam, ' |', 'bcftools call -c', '-o', vcf_file]))

    os.system('rm ' + sam + ' ' + bam + ' ' + sorted_bam + ' ' + sorted_bam + '.bai ' + args.ref_id + '.bed')
    return print("Samtools and bcftools run completed")


# THE DRIVER
def main():
    print("\nYou are using Paratype version " + version + '\n')
    if args.ref and args.allele and args.id and args.genes:
        # setup strain ID
        strainID = generate_strainID(args.id)

        # setup the reference fasta file
        # index ref.fasta if not indexed
        if os.path.exists(ref_fasta_file + '.fai') == False:
            print('\nreference fasta file is not indexed. Indexing now ...\n')
            os.system(' '.join(['samtools faidx', ref_fasta_file]))

        # set ref_index for bowtie2
        # assuming the extension of reference file is 'fasta', not 'fa'. If it is 'fa', change the number -6 to -3.
        ref_index = '_'.join([(ref_fasta_file[:-6]), 'index'])

        # set up sam, bam, vcf files
        sam = '.'.join([strainID, 'sam'])
        bam = '.'.join([strainID, 'bam'])
        sorted_bam = '.'.join([strainID, 'sorted.bam'])
        vcf_file = '.'.join([strainID, 'vcf'])

        # Setup output file
        if args.output == 'paratype_results.txt':
            outfile = open((strainID + '_' + args.output), 'w')
        else:
            outfile = open(args.output, 'w')

        # setup the lists
        clades, loci, alleles = define_genotypes(genotype_allele_file)

        # setup gene region list
        mutation, strand, location, codon, base, codon_pos, region = define_genes(
            gene_regions_file)

        # If mode is set to fastq:
        if args.mode == "fastq" and args.fastq:
            fastq1 = args.fastq[0]
            fastq2 = args.fastq[1]

            # Build ref_index for bowtie2
            os.system(' '.join(['bowtie2-build -q', ref_fasta_file, ref_index]))

            # Run bowtie2 for mapping
            print("\nBowtie2 mapping is starting....\n")
            os.system(' '.join(['bowtie2 -x', ref_index, '-1', fastq1,
                                '-2', fastq2, '-S', sam, '-p', str(args.threads)]))
            print("\nMapping is complete.\n")

            # Run samtools view, sort, mpileup and bcftools call to generate vcf for a fixed number of loci
            run_samtools(ref_fasta_file, args.mapq_cutoff, args.ref_id, args.threads, sam, bam, sorted_bam, vcf_file)

        # If mode is set to fastq_interleaved:
        if args.mode == "fqin" and args.fqin:
            fastq = args.fqin

            # Build ref_index for bowtie2
            os.system(' '.join(['bowtie2-build -q', ref_fasta_file, ref_index]))

            # Run bowtie2 for mapping
            print("\nBowtie2 mapping is starting....\n")
            os.system(' '.join(['bowtie2 -x', ref_index, '--interleaved', fastq,
                                '-S', sam, '-p', str(args.threads)]))
            print("\nMapping is complete.\n")

            # Run samtools view, sort, mpileup and bcftools call to generate vcf for a fixed number of loci
            run_samtools(ref_fasta_file, args.mapq_cutoff, args.ref_id, args.threads, sam, bam, sorted_bam, vcf_file)

        # If mode is set to nanopore fastq:
        if args.mode == "nano" and args.nano:
            fastq = args.nano

            # Build ref_index for bwa
            os.system(' '.join(['bwa index ', ref_fasta_file]))

            # Run bwa mapping
            print("\nbwa mapping is starting....\n")
            os.system(' '.join(['bwa mem -t ', str(args.threads), ref_fasta_file, fastq, ' > ', sam]))
            print("\nMapping is complete.\n")

            # Run samtools view, sort, mpileup and bcftools call to generate vcf for a fixed number of loci
            run_samtools(ref_fasta_file, args.mapq_cutoff, args.ref_id, args.threads, sam, bam, sorted_bam, vcf_file)

        # If mode is set to fasta:
        if args.mode == "fasta" and args.fasta:
            fasta = args.fasta

            # Build ref_index for bwa
            os.system(' '.join(['bwa index ', ref_fasta_file]))

            # Run bwa mem mapping
            print("\nbwa mem mapping is starting....\n")
            os.system(' '.join(['bwa mem -t ', str(args.threads), ref_fasta_file, fasta, ' > ', sam]))
            print("\nMapping is complete.\n")

            # Run samtools view, sort, mpileup and bcftools call to generate vcf for a fixed number of loci
            run_samtools(ref_fasta_file, args.mapq_cutoff, args.ref_id, args.threads, sam, bam, sorted_bam, vcf_file)

        # If mode is set to bam:
        elif args.mode == "bam" and args.bam:
            # setup bam file
            bam = args.bam

            # Run samtools mpileup and bcftools call to generate vcf for a fixed number of loci
            print("\nSamtools mpilieup is running....\n")
            os.system(' '.join(['samtools sort --threads',
                                str(args.threads), bam, '-o', sorted_bam]))
            os.system(' '.join(['samtools index', sorted_bam]))
            os.system(' '.join(['samtools mpileup -q', str(args.mapq_cutoff), '-ugBf', ref_fasta_file,
                                '-l', str(args.ref_id + '.bed'), sorted_bam, '-I |', 'bcftools call -c', '-o',
                                vcf_file]))

            os.system('rm ' + sorted_bam + ' ' + sorted_bam + '.bai ' + args.ref_id + '.bed')

        # If mode is set to vcf:
        elif args.mode == "vcf" and args.vcf:
            vcf_file = args.vcf
            
            os.system('rm ' + args.ref_id + '.bed')

        # detect type/ clade list from vcf file
        type_list, propor = check_allele_from_vcf(
            vcf_file, clades, loci, alleles)

        # Classify primary, secondary and subclades (genotypes)
        final_result, avg_ratio = designate_genotypes(type_list, propor)

        # Detect position-specific high-quality SNPs for mutations
        snp, nucl = select_snp_from_vcf(vcf_file, location)

        # Classify Mutations
        mutation_list = classify_mutations(
            snp, nucl, mutation, strand, location, codon, codon_pos)
        print_mutations = ','.join(mutation_list)
        print(strainID + '\t' + str(type_list) + '\t' + print_mutations)
        
        if args.mode != "vcf":
            os.system('rm ' + vcf_file)
		
        # print results
        outfile.write('Strain\tPrimary_clade\tSecondary_clade\tSubclade\tGenotype\tSupport\tMutations\n' +
                      strainID + '\t' + final_result + '\t' + avg_ratio + '\t' + print_mutations + '\n')

    else:
        print('Please check if you have provided right allele file and reference fasta.')
    print("Thank you for using Paratype.")


# call main function
if __name__ == '__main__':
    main()
