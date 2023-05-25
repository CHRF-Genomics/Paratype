# Batch_run_paratype.sh version-2.0
# Author: Preonath Chandrow Dev, Anisur Rahman, and Arif Mohammad Tanmoy
echo "# Batch_run_paratype.sh"
echo "# version: 2.0"
echo "# You are using Batch_run_paratype.sh script to run Paratype for multiple isolates."$'\n'



### Display Argument ####

display_help() {
    
	echo $'\n'"# Authors: Preonath Chandrow Dev, Anisur Rahman, and Arif Mohammad Tanmoy."$'\n'
	echo "# USAGE: bash Batch_run_paratype.sh -i Input_directory | -o Result_directory |-p Paratype_directory | -m Run_mode" 

	echo $'\n'"Mandatory Arguments:"

	echo "-i, --input           Directory of Input files"
	echo "-o, --output          Directory of Result files"
	echo "-p, --paratype        Paratype directory"
	echo "-m, --mode            Running modes: fasta, fastq, fqin, bam, nano, vcf"

	echo $'\n'"Optional Arguments:"
	echo "-t, --threads         Number of threads to run Paratype"

	echo                        $'\n'"Other Arguments:"
	echo "-h, --help 	      Print the usage message"
 
}



### Argument Handling ####


while getopts ":m:i:o:p:t:h" flag
do  
    

    case "${flag}" in
    m) run_type=${OPTARG};;
    i) input=${OPTARG};;
    o) results=${OPTARG};;
    p) paratype=${OPTARG};;
    t) threads=${OPTARG};;
    h) display_help && exit 0;;
    *) echo "Invalid option -$OPTARG" >&2 && exit 1;;
    esac

done


if [ -z "$run_type" ] || [ -z "$input" ] || [ -z "$results" ] || [ -z "$paratype" ]; then
    echo "Error: Required argument is missing." >&2
    display_help >&2
    exit 1
fi



declare array=()

for file in $input\/*
do  
    name=$(basename "$file")
    extension=""
    extension_2=""

    if [[ $name =~ \.gz$ ]]; then
        name=$(basename "$file" | sed 's/\.gz$//')
        if [[ $name =~ \.fastq$ ]]; then
            extension="fastq"
        elif [[ $name =~ \.fq$ ]]; then
            extension="fq"
        else
            extension="unknown"
        fi
    else
        extension="${name##*.}"
    fi

    case $extension in
        fastq)
            extension_2="fastq"
            ;;
        fq)
            extension_2="fq"
            ;;
        bam)
            extension_2="bam"
            ;;
        fasta)
            extension_2="fasta"
            ;;
        vcf)
            extension_2="vcf"
            ;;
        *)
            extension_2="unknown"
            ;;
    esac




	if [[ "$run_type" == "bam" && "$extension_2" == "bam" ]]
	then
        ## Considering the 1st part before dot(.) as strain name.
        
        echo $name
		
        ## Mode: bam
		echo "Running bam mode........"
		
		mkdir -p -v $results
		time (python $paratype\/paratype.py --id $name --mode bam --bam $file --output $results\/$name\_paratype.tsv --threads $threads)
    
    elif [[ "$run_type" == "fqin" && "$extension_2" == "fastq" ]]
	then
        
        echo $name
    
		echo "Running fastq_interleaved mode......"
		
		mkdir -p -v $results
		## Mode: fastq interleaved
		time (python $paratype\/paratype.py --id $name --mode fqin --fqin $file --output $results\/$name\_paratype.tsv --threads $threads)
    
    elif [[ "$run_type" == "nano" && "$extension_2" == "fastq" || "$extension_2" == "fq" ]]
	then

        
        echo $name
		
		mkdir -p -v $results
		## Mode: nano
		echo "Running nanopore mode........"
		time (python $paratype\/paratype.py --id $name --mode nano --nano $file --output $results\/$name\_paratype.tsv --threads $threads)
    
    elif [[ "$run_type" == "fasta" && "$extension_2" == "fasta" ]]
	then    
        
        echo $name

		## Mode: fasta
		
		mkdir -p -v $results
		echo "Running fasta mode......"
		time (python $paratype\/paratype.py --id $name --mode fasta --fasta $file --output $results\/$name\_paratype.tsv --threads $threads)

	elif [[ "$run_type" == "vcf" && "$extension_2" == "vcf" ]]
	then	
		
        
        echo $name
        
        ## Mode: vcf
        
		mkdir -p -v $results
		echo "Running vcf mode...."
		time (python $paratype\/paratype.py --id $name --mode vcf --vcf $file --output $results\/$name\_paratype.tsv --threads $threads)

    elif [ "$run_type" == "fastq" ]
	then
        
		## Considering the 1st part before dot(_) as strain name.
        name=`basename $file | cut -f 1 -d '_'`
        
        if [[ ${array[*]} != ${name}  ]]; 
        then
            ## Mode: fastq
            echo "Running fastq mode....."
            
			mkdir -p -v $results
            time (python $paratype\/paratype.py --id $name --mode fastq --fastq $input\/$name\*1.fastq.gz $input\/$name\*2.fastq.gz --output $results\/$name\_paratype.tsv --threads $threads)
                # Add new element at the end of the array
                array+=(${name})
        fi
    
    elif [ "$run_type" == "fastq" ]
	then
        
		## Considering the 1st part before dot(_) as strain name.
        name=`basename $file | cut -f 1 -d '_'`
        extension="${name##*.}"
        if [[ ${array[*]} != ${name}  ]]; 
        then
            ## Mode: fastq
            echo "Running fastq mode....."
            
			mkdir -p -v $results
			if ["$extension" == "fastq"]
			then
				time (python $paratype\/paratype.py --id $name --mode fastq --fastq $input\/$name\*1.fastq.gz $input\/$name\*2.fastq.gz --output $results\/$name\_paratype.tsv --threads $threads)
					# Add new element at the end of the array
					array+=(${name})
			fi

			elif [ "$extension" == "fq" ]
			then
				time (python $paratype\/paratype.py --id $name --mode fastq --fastq $input\/$name\*1.fq.gz $input\/$name\*2.fq.gz --output $results\/$name\_paratype.tsv --threads $threads)
					# Add new element at the end of the array
					array+=(${name})
		

        fi
    
	else	
		echo "Error: Invalid run type specified. Please use one of the following: fasta, fastq, fqin, bam, nano, vcf mode for your specific files"
		exit 1
	fi
	
done

