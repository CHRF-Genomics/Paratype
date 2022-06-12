## Author: Arif Mohammad Tanmoy
## USAGE: bash Batch_run_paratype.sh <input file directory> <results directory> <Paratype directory>
echo "You are using Batch_run_paratype script to run Paratype with multiple isolates."
path=$1
results=$2
paratype=$3

mkdir -p $results

for file in $path\/*
do
	## Considering the 1st part before dot(.) as strain name.
	name=`basename $file | cut -f 1 -d '.'`
	echo $name

	## Considering python3 as the default python in user system.
	## The followings are the commands required for each paratype mode. Unmute the line you need as run accordingly.
	## Mode: bam
	time(python $paratype\/paratype.py --id $name --bam $file --output $results\/$name\.paratype --threads 8)
	## Mode: fastq interleaved
	#time(python $paratype\/paratype.py --id $name --mode fqin --fqin $file --output $results\/$name\.paratype --threads 8)
	## Mode: nano
	#time(python $paratype\/paratype.py --id $name --mode nano --nano $file --output $results\/$name\.paratype --threads 8)
	## Mode: fasta
	#time(python $paratype\/paratype.py --id $name --mode fasta --fasta $file --output $results\/$name\.paratype --threads 8)
	## Mode: vcf
	#time(python $paratype\/paratype.py --id $name --mode vcf --vcf $file --output $results\/$name\.paratype --threads 8)
	
	## Mode: fastq
	## Running fastq mode without knowing the filename extension of paired-end fastq files. Assuming the file name is (Name_1.fastq.gz). Unmute both the lines below.
	#newname=`cut -f 1 -d '_'`
	#time(python $paratype --id $name --mode fastq --fastq $path\/$newname\*1.fastq.gz $path\/$newname\*2.fastq.gz --output $results\/$newname\.paratype --threads 8)
done
