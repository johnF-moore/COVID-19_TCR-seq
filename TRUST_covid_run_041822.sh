#!/bin/bash

# The goal of this script is to run TRUST4 to pull all TCR sequences with V gene allele information from the
# TCR fastq output from Cellranger. Cellranger as of 10/31/21 doesn't support v gene allele information.

N=5
dir=/stor/work/Ehrlich_COVID19/analysis/cellranger/T_cell_activated
run=run_covid-reactive_041822
cd $dir
mkdir_dir=/stor/work/Ehrlich_COVID19/John/VDJ_data/TRUST4_out/
if [ ! -d $mkdir_dir/$run ]; then
    mkdir $mkdir_dir/$run
fi
seq_jobs=$(ls | grep JA)


for job in $seq_jobs; do
    echo "Inside first loop"
    echo $job
    cd $dir/$job/fastq_symlinks/TCR
     # 4 lanes (L001, L002, L003, & L004) with two reads each.
     # Read 1 is the barcode.
     # Read 2 is the sequencing data.
     # TRUST4 requires both reads for 10X Genomics data.
     # TRUST4 also recommends using the barcode whitelist, but that doesn't work as of 11/7/21.

     # There are multiple patients run in each job for these samples.
    patients=$(ls | gawk -F_ '{ print $1 }' | grep TCR | sort | uniq)
    echo $patients    
     ## ls and then select everything prior to the first "_",
     ## sorting is necessary because uniq only removes duplicate neighboring values.
     ## Together, this gives the unique patient jobs
    echo $patients | gawk -F'-TCR ' '{print $2}'
    patients=$(echo $patients | gawk -F'-TCR ' '{print $2}')

    if [ "$job" = JA21414 ]; then 
	echo "skipping JA21414"
	printf "\n"
	continue
	
    for patient in $patients; do
	cd $dir/$job/fastq_symlinks/TCR
	pat_files=$(ls | grep $patient)
	echo "inside second loop"
	echo $patient
	echo $pat_files
	 
	L001_R1=($(find $(ls $pat_files) -regextype sed -regex ".*L001_R1.*" | sed 's|^./||'))
	L001_R2=($(find $(ls $pat_files) -regextype sed -regex ".*L001_R2.*" | sed 's|^./||'))
	L002_R1=($(find $(ls $pat_files) -regextype sed -regex ".*L002_R1.*" | sed 's|^./||'))
	L002_R2=($(find $(ls $pat_files) -regextype sed -regex ".*L002_R2.*" | sed 's|^./||'))
	L003_R1=($(find $(ls $pat_files) -regextype sed -regex ".*L003_R1.*" | sed 's|^./||'))
	L003_R2=($(find $(ls $pat_files) -regextype sed -regex ".*L003_R2.*" | sed 's|^./||'))
	L004_R1=($(find $(ls $pat_files) -regextype sed -regex ".*L004_R1.*" | sed 's|^./||'))
	L004_R2=($(find $(ls $pat_files) -regextype sed -regex ".*L004_R2.*" | sed 's|^./||')) 

	## find gets specific files from the given files in ls $pat_files, which has files that are specific for one patient
	 ## -regextype is what it sounds like
	 ## -regex is saying that as long as it has that string somewhere in the file name, it's a match
	 ## I don't understand what the sed is saying, and I'm not sure if it's doing anything for this script.
	 ## There was a reason I used it in a previous script, but I forget. 

	# Paths for TRUST4 code
	trust_path=/stor/work/Ehrlich_COVID19/John/TCR_code/TRUST4
	reads_path=$dir/$job/fastq_symlinks/TCR
			
	# Combining lanes to simplify TRUST4 call                                                                                         
	# Read 1 and 2 have to have the same order b/c they correspond to each other.                                                     
	cat $reads_path/$L001_R1 \
	    $reads_path/$L002_R1 \
	    $reads_path/$L003_R1 \
	    $reads_path/$L004_R1 > R1_reads

	cat $reads_path/$L001_R2 \
	    $reads_path/$L002_R2 \
	    $reads_path/$L003_R2 \
	    $reads_path/$L004_R2 > R2_reads
       
	#### Running TRUST4 ####
	# See TRUST4 github: https://github.com/liulab-dfci/TRUST4
	trust_out_path=/stor/work/Ehrlich_COVID19/John/VDJ_data/TRUST4_out
	cd $trust_out_path/$run
	echo "Starting TRUST4"
	($trust_path/run-trust4 \
            -f $trust_path/hg38_bcrtcr.fa \
	    --ref $trust_path/human_IMGT+C.fa \
	    -u $reads_path/R2_reads \
	    --barcode $reads_path/R1_reads \
	    --barcodeRange 0 15 + \
	    -t 4 \
	    -o $job\_$patient\_combined) &

	if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
                # now there are $N jobs already running, so wait here for any job
                # to be finished so there is a place to start next one.
                wait -n
        fi
    done

    else
	for patient in $patients; do
	cd $dir/$job/fastq_symlinks/TCR    
        pat_files=$(ls | grep $patient)
        echo "inside second loop"
	echo $patient
        echo $pat_files
	
        L001_R1=($(find $(ls $pat_files) -regextype sed -regex ".*L001_R1.*" | sed 's|^./||'))
	L001_R2=($(find $(ls $pat_files) -regextype sed -regex ".*L001_R2.*" | sed 's|^./||'))
        L002_R1=($(find $(ls $pat_files) -regextype sed -regex ".*L002_R1.*" | sed 's|^./||'))
        L002_R2=($(find $(ls $pat_files) -regextype sed -regex ".*L002_R2.*" | sed 's|^./||'))

	trust_path=/stor/work/Ehrlich_COVID19/John/TCR_code/TRUST4
        reads_path=$dir/$job/fastq_symlinks/TCR
        echo "L001_R1 files: ${L001_R1}"                                                                                                      
        cat $reads_path/$L001_R1 \
            $reads_path/$L002_R1 > R1_reads

        cat $reads_path/$L001_R2 \
            $reads_path/$L002_R2  > R2_reads

	trust_out_path=/stor/work/Ehrlich_COVID19/John/VDJ_data/TRUST4_out
        cd $trust_out_path/$run
        echo "Starting TRUST4"
        ($trust_path/run-trust4 \
            -f $trust_path/hg38_bcrtcr.fa \
            --ref $trust_path/human_IMGT+C.fa \
            -u $reads_path/R2_reads \
            --barcode $reads_path/R1_reads \
            --barcodeRange 0 15 + \
            -t 4 \
            -o $job\_$patient\_combined) &

	
        if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
                # now there are $N jobs already running, so wait here for any job                        
                # to be finished so there is a place to start next one.                                       
                wait -n
        fi
	done
    fi
    printf "\n"
done
