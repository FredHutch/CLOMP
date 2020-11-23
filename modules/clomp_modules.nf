/*
 * Define the parameters used in the processes below
 */

params.SECOND_PASS = false
params.BWT_SECOND_PASS_OPTIONS = false
params.BWT_DB_PREFIX = false
params.SEQUENCER = false
params.TRIMMOMATIC_OPTIONS = ':2:30:10 HEADCROP:10 SLIDINGWINDOW:4:20 CROP:65 MINLEN:65'
params.SNAP_OPTIONS = false
params.BLAST_CHECK = false
params.WRITE_UNIQUES = false
params.BLAST_EVAL = false
params.ADD_HOST_FILTERED_TO_REPORT = true
params.HOST_FILTER_TAXID = 9606
params.H_STRICT = false
params.H_TAXID = 9606
params.FILTER_LIST = "[12908,28384,48479,99802]"
params.LOGIC = "strict"
params.INCLUSION_TAXID = 2759
params.EXCLUSION_TAXID = 9604
params.ENTREZ_EMAIL = "uwvirongs@gmail.com"
params.BASE_DELIMITER = "_"
params.MIN_READ_CUTOFF = 10
params.SAM_NO_BUILD_LIST = "[2759,77133]"
params.EDIT_DISTANCE_OFFSET = 6
params.BUILD_SAMS = false
params.TIEBREAKING_CHUNKS = 16
params.FASTQ_MIN_READ = 10

/*
 * Define the processes used in this workflow
 */

// Validate that the input data is GZIP-compressed FASTQ - paired-end
process validate_paired {
    container "ubuntu:18.04"

    input:
      tuple val(prefix), file(R1), file(R2)

    output:
      tuple file("${prefix}.R1.fastq.gz"), file("${prefix}.R2.fastq.gz")

"""
#!/bin/bash

set -e

# Check for gzip-compressed input
gzip -t ${R1} || (echo "${R1} is not gzip-compressed" && exit 1)
gzip -t ${R2} || (echo "${R2} is not gzip-compressed" && exit 1)

# Check that both read pairs have the same number of lines
[[ \$(gunzip -c ${R1} | wc -l) != \$(gunzip -c ${R2} | wc -l) ]] & echo "${R1} and ${R2} have different numbers of lines" && exit 1

# Rename the input files
mv ${R1} ${prefix}.R1.fastq.gz
mv ${R2} ${prefix}.R2.fastq.gz
"""
}

// Validate that the input data is GZIP-compressed FASTQ - single-end
process validate_single {
    container "ubuntu:18.04"

    input:
      tuple val(prefix), file(R1)

    output:
      file "${prefix}.R1.fastq.gz"

"""
#!/bin/bash

set -e

# Check for gzip-compressed input
gzip -t ${R1} || (echo "${R1} is not gzip-compressed" && exit 1)

# check number of lines in fastq
if [[ `zcat ${R1} | wc -l` -le ${params.FASTQ_MIN_READ} ]] ; then echo "${R1} contains less than ${params.FASTQ_MIN_READ} reads" ; exit 1 ; fi

# Rename the input file
mv ${R1} ${prefix}.R1.fastq.gz


"""
}

process filter_human_paired {
    // Retry at most 3 times
    errorStrategy 'retry'
    maxRetries 3

    // Define the Docker container used for this step
    container "quay.io/fhcrc-microbiome/bowtie2:bowtie2-2.2.9-samtools-1.3.1"

    // Define the input files and set their filename in the execution folder
    input:
      tuple file(R1), file(R2)
      file "*"

    // Define the output files
    output:
      tuple file("${r1}"), file("${r2}")
      file "*.log"

    // Code to be executed inside the task
    script:
      """
#!/bin/bash

set -e

# For logging and debugging, list all of the files in the working directory
ls -lahtr

# Get the sample name from the input FASTQ name
sample_name=\$(echo ${r1} | sed 's/.R1.fastq.gz//')

echo "Starting the alignment of ${r1} and ${r2}"
bowtie2 \
    ${params.BWT_SECOND_PASS_OPTIONS} \
    --threads ${task.cpus} \
    
    -x ${params.BWT_DB_PREFIX} \
    -q \
    -1 <(gunzip -c ${r1}) \
    -2 <(gunzip -c ${r2}) |\
    samtools view -Sb - > \${sample_name}_mappedBam 2>&1 | \
    tee -a \${sample_name}.log

# Delete the input R1 and R2 so that we don't have to worry
# that they are being used (in error) as outputs
rm ${r1} ${r2}

#echo "Extracting the BAM alignments"
#samtools view -Sb -@ 16 \${sample_name}_mappedSam > \${sample_name}_mappedBam

# Extract the R1
echo "Extracting the R1 FASTQ"
samtools view -@ ${task.cpus} -F 2 \${sample_name}_mappedBam | \
samtools view -@ ${task.cpus} -f 64 - | \
    awk \'{if(\$3 == \"*\") print \"@\" \$1 \"\\n\" \$10 \"\\n\" \"+\" \$1 \"\\n\" \$11}\' | \
    gzip -c > \${sample_name}.R1.fastq.gz

echo "Extracting the R2 FASTQ"
# Extract the R2
samtools view -@ ${task.cpus} -F 2 \${sample_name}_mappedBam | \
samtools view -@ ${task.cpus} -f 128 - | \
    awk \'{if(\$3 == \"*\") print \"@\" \$1 \"\\n\" \$10 \"\\n\" \"+\" \$1 \"\\n\" \$11}\' | \
    gzip -c > \${sample_name}.R2.fastq.gz
      """
}


process filter_human_single {

    // Retry at most 3 times
    errorStrategy 'retry'
    maxRetries 3
    
    // Define the Docker container used for this step
    container "quay.io/fhcrc-microbiome/bowtie2:bowtie2-2.2.9-samtools-1.3.1"

    // Define the input files
    input:
      tuple val(base), file (r1)
      file "*"

    // Define the output files
    output:
      file "${r1}"
      file "*.log"
      tuple env(base), env(difference)
    // Code to be executed inside the task
    script:
      """
#!/bin/bash

set -e

# For logging and debugging, list all of the files in the working directory
ls -lahtr




base=`basename ${base} .R1`
echo \$base

# Get the sample name from the input FASTQ name
sample_name=\$(echo ${r1} | sed 's/.R1.fastq.gz//')

readsbefore=`gzip -dc ${r1} | awk 'NR%4==2{c++} END{ print c;}'`

echo "\$readsbefore Reads before host filtering "





echo "Starting the alignment of ${r1}"
bowtie2 \
    ${params.BWT_SECOND_PASS_OPTIONS} \
    --threads ${task.cpus} \
    -x ${params.BWT_DB_PREFIX} \
    -q \
    -U <(gunzip -c ${r1}) | \
    samtools view -@ 4 -Sb - > \${sample_name}_mappedBam 2>&1 | \
    tee -a \${sample_name}.log

# Delete the input R1 so that we don't have to worry
# that it is being used (in error) as output
rm ${r1}

# echo ${workDir} > \${sample_name}.log
#echo "Extracting the BAM alignments"
#samtools view -Sb -@ 16 \${sample_name}_mappedSam > \${sample_name}_mappedBam

# Extract the R1
echo "Extracting the FASTQ"
samtools view -@ ${task.cpus} -f 4 \${sample_name}_mappedBam | \
    awk \'{if(\$3 == \"*\") print \"@\" \$1 \"\\n\" \$10 \"\\n\" \"+\" \$1 \"\\n\" \$11}\' | \
    gzip -c > \${sample_name}.R1.fastq.gz

readsafter=`gzip -dc ${r1} | awk 'NR%4==2{c++} END{ print c;}'`
echo "\$readsafter Reads after host filtering "

# Reads lost in host filtering
difference=\$((\$readsbefore-\$readsafter))

echo "host filtered \$difference reads"

      """
}

// Depracated because of java memory errors with big files. Using bbduk instead. 
process trimmomatic_single {

    // Retry at most 3 times
    errorStrategy 'retry'
    maxRetries 3
    
    // Define the Docker container used for this step
    container "quay.io/fhcrc-microbiome/openjdk:jre-8"

    // Define the input files
    input:
      file r1
      file TRIMMOMATIC_JAR
      file TRIMMOMATIC_ADAPTER 

    // Define the output files
    output:
      file "${r1}"

    // Code to be executed inside the task
    script:
      """
#!/bin/bash

set -e

# For logging and debugging, list all of the files in the working directory
ls -lahtr

echo "Reads before adapter trimming (Trimmomatic) "

gzip -dc ${r1} | awk 'NR%4==2{c++} END{ print c;}'


echo "Starting to trim ${r1} with options ${params.TRIMMOMATIC_OPTIONS}"

# Rename the file to prevent collision
mv ${r1} INPUT.${r1}

java -jar \
    ${TRIMMOMATIC_JAR} \
    SE \
    -threads ${task.cpus} \
    INPUT.${r1} \
    ${r1} \
    ${params.SEQUENCER}${TRIMMOMATIC_ADAPTER}${params.TRIMMOMATIC_OPTIONS}


echo 'Reads after adapter trimming '
gzip -dc ${r1} | awk 'NR%4==2{c++} END{ print c;}'

"""
}



process trimmomatic_paired {

    // Retry at most 3 times
    errorStrategy 'retry'
    maxRetries 3
    
    // Define the Docker container used for this step
    container "quay.io/fhcrc-microbiome/openjdk:jre-8"

    // Define the input files
    input:
      tuple file(r1), file(r2)
      file TRIMMOMATIC_JAR
      file TRIMMOMATIC_ADAPTER 

    // Define the output files
    output:
      tuple file("${r1}"), file("${r2}")

    // Code to be executed inside the task
    script:
      """
#!/bin/bash

set -e

# For logging and debugging, list all of the files in the working directory
ls -lahtr



java -jar \
    ${TRIMMOMATIC_JAR} \
    PE \
    -threads ${task.cpus} \
    ${r1} \
    ${r2} \
    OUTPUT_R1_trimmed.fastq.gz \
    OUTPUT_R1_unpaired.fastq.gz \
    OUTPUT_R2_trimmed.fastq.gz \
    OUTPUT_R2_unpaired.fastq.gz \
    ${params.SEQUENCER}${TRIMMOMATIC_ADAPTER}${params.TRIMMOMATIC_OPTIONS}

mv OUTPUT_R1_trimmed.fastq.gz ${r1}
mv OUTPUT_R2_trimmed.fastq.gz ${r2}

"""
}

process bbMask_Single {

    // Retry at most 3 times
    errorStrategy 'retry'
    maxRetries 3
    
    // Define the Docker container used for this step
    // should build our own docker image for this 
    container "quay.io/biocontainers/bbmap:38.76--h516909a_0"

    // Define the input files
    input:
      file r1
      file TRIMMOMATIC_ADAPTER 

    // Define the output files
    output:
      file("${r1}")

    // Code to be executed inside the task
    script:
      """
#!/bin/bash

set -e

# For logging and debugging, list all of the files in the working directory
ls -lahtr

# Get the sample name from the file name
sample_name=\$(echo ${r1} | sed 's/.R1.fastq.gz//')
echo "Processing \$sample_name"

# Rename the input file to make sure we don't use it as the output
mv ${r1} INPUT.${r1}

echo "Reads before low complexity filtering "
gzip -dc INPUT.${r1} | awk 'NR%4==2{c++} END{ print c;}'

echo "Masking ${r1}"
bbduk.sh \
    in=INPUT.${r1} \
    out=${r1}.lcf.fastq.gz \
    entropy=0.7 \
    entropywindow=50 \
    entropyk=4 
    

    

    mv ${r1}.lcf.fastq.gz ${r1}

echo "Reads after low complexity filtering"
gzip -dc ${r1} | awk 'NR%4==2{c++} END{ print c;}'


"""
}

// Not being used anymore
process deduplicate { 

    // Retry at most 3 times
    errorStrategy 'retry'
    maxRetries 3
    
    // Define the Docker container used for this step
    // should build our own docker image for this 
    container "quay.io/biocontainers/bbmap:38.76--h516909a_0"

    // Define the input files
    input:
      file r1

    // Define the output files
    output:
      file("${r1}")

    // Code to be executed inside the task
    script:
      """
      #!/bin/bash

      set -e

      # For logging and debugging, list all of the files in the working directory
      ls -lahtr

      # Get the sample name from the file name
      sample_name=\$(echo ${r1} | sed 's/.R1.fastq.gz//')
      echo "Processing \$sample_name"

      dedupe2.sh in=${r1} out=${r1}.deduped.fastq.gz -Xmx120g

      mv ${r1}.deduped.fastq.gz ${r1}
      
      """


}
process snap_paired {

    // Retry at most 3 times
    errorStrategy 'retry'
    maxRetries 3
    
    // Define the Docker container used for this step
    container "quay.io/fhcrc-microbiome/snap-no-header:v1.0beta.18--0"

    // Define the input files
    input:
      file paired_fastq_list
      each file(SNAP_DB)

    // Define the output files
    output:
      file("*${SNAP_DB.name}.sam")

    // Code to be executed inside the task
    script:
      """
#!/bin/bash

set -e

# For logging and debugging, list all of the files in the working directory
ls -lahtr

# Iterate over each of the R1 files
for r1 in *.R1.fastq*; do

    # Get the sample name from the file name
    sample_name=\$(echo \${r1} | sed 's/.R1.fastq.gz//')
    echo "Processing \$sample_name"

    r2=\${sample_name}.R1.fastq.gz
    echo "Using \${r2} as the paired-end reads"

    # Check for R2 file existence
    [[ -s \${r2} ]]

    echo "Aligning \${r1} and \${r2}"

    # Decompress the input files
    echo "Decompressing \${r1}"
    gunzip -c \${r1} > R1.fastq && rm \${r1}
    echo "Decompressing \${r2}"
    gunzip -c \${r2} > R2.fastq && rm \${r2}

    echo "Running SNAP"
    snap-aligner paired ${SNAP_DB} R1.fastq R2.fastq -o \${sample_name}__${SNAP_DB.name}.sam -t ${task.cpus} ${params.SNAP_OPTIONS}

    echo "Removing temporary files"
    rm R1.fastq R2.fastq
    
done
"""
}

process snap_single {

   // Retry at most 3 times
    errorStrategy 'retry'
    maxRetries 3
    
    // Define the Docker container used for this step
    container "quay.io/fhcrc-microbiome/snap-no-header:v1.0beta.18--0"

    // Define the input files
    input:
      file r1_list
      each file(SNAP_DB)

    // Define the output files
    output:
      file("*${SNAP_DB.name}.bam")

    // Code to be executed inside the task
    script:
      """
#!/bin/bash
set -e
# For logging and debugging, list all of the files in the working directory
ls -lahtr

# Free space in directory 
df -h 

#Crashes and retries if samples aren't downloaded to worker correctly 
for fp in ${r1_list}; do
  echo Checking to make sure that \$fp was downloaded to the worker
  [[ -s \$fp ]]
done

# Logging ls 
ls -lah ${SNAP_DB}/

# Makes sure snap databases are downloaded correctly
echo Checking to make sure that the full database is available at ${SNAP_DB}
#[[ -f ${SNAP_DB}/GenomeIndexHash ]]
#[[ -f ${SNAP_DB}/OverflowTable ]]
[[ -f ${SNAP_DB}/Genome ]]
[[ -f ${SNAP_DB}/GenomeIndex ]]


echo "Aligning ${r1_list}"
echo "snap-aligner " | tr -d "\n" > snap_cmd.sh
# Iterate over each of the input files
for r1 in ${r1_list}; do
    # Get the sample name from the file name
    sample_name=\$(echo \${r1} | sed 's/.R1.fastq.gz//')
    echo "Processing \$sample_name"
    # Decompress the input files
    echo "Decompressing \${r1}"
    gunzip -c \${r1} > \${sample_name}.fastq
    rm \${r1}
    filename=\${sample_name}__${SNAP_DB.name}.bam
    temp_cmd=\$(echo "single ${SNAP_DB} \${sample_name}.fastq -t ${task.cpus} ${params.SNAP_OPTIONS} -o \$filename")
    temp_cmd="\$temp_cmd"" ,  "
    echo \$temp_cmd | tr "\n" " ">> snap_cmd.sh
    
    #echo \$cmd
    #echo "Running SNAP"
    #snap-aligner single ${SNAP_DB} R1.fastq -t ${task.cpus} ${params.SNAP_OPTIONS} -o -bam - `> \${sample_name}__${SNAP_DB.name}.bam`
    echo "Removing temporary files"
    # rm R1.fastq
done
bash snap_cmd.sh
"""
}

process collect_snap_results {

    // Retry at most 3 times
    //errorStrategy 'retry'
    //maxRetries 3
    
    // Define the Docker container used for this step
    container "quay.io/vpeddu/clomp_containers:latest"
    containerOptions = "--user root"

    // Define the input files
    input:
      tuple val(base), file(bam_list)
      file SAM_SPLIT 

    // Define the output files
    output:
      tuple val(base), file("${base}*.sam")
      file "sortedsam.sam"
    // Code to be executed inside the task
    script:
      """
#!/bin/bash

set -e

# For logging and debugging, list all of the files in the working directory
ls -lahtr


echo "Merging BAM files for ${base}"

for i in ${bam_list}; do /usr/bin/samtools view \$i >> ${base}.sam; done

echo "sorting ${base} pseudosam by readname"

sort -k 1 --buffer-size=200G ${base}.sam > ${base}.sorted.sam

echo "ls after sort" 
ls -lah

echo "splitting ${base} pseudosam"

#linenum=`cat ${base}.sorted.sam | wc -l`
#splitnum=`echo \$(( \$linenum / ${params.TIEBREAKING_CHUNKS} ))`

echo "lines to split: "\$splitnum 

# split files by splitnum number of lines
#split -d -a 3  -l \$splitnum ${base}.sorted.sam ${base}

echo "ls after split"
ls -latr

python3 ${SAM_SPLIT} ${base}.sorted.sam ${params.TIEBREAKING_CHUNKS} ${base}

# MICHELLE SPLIT REALLY BROKEN 
#basename=${base}

#files=(\$basename*)

#echo \$files

# number of files made in split
#total=`ls ${base}0* | wc -l`

# total=\${#files[@]}
#total=\$((total-2))

#echo "basename "\$basename
#echo "total "\$total

# Michelle's split command
# Checks head/tail of each file to check record boundaries 
#fixSplit() {
#        inum=\$(printf %03d \$i)
#        last_record=\$(tail -1 \$basename\$inum | cut -f1)
#        j=\$((i + 1))
#        j=\$(printf %03d \$j)
#        echo "Last record in "\$basename\$inum" is "\$last_record". Grepping in "\$basename\$j"."
#        head -n 30000 \$basename\$j | grep \$last_record >> \$basename\$inum
#        echo "Removing "\$last_record" from "\$basename\$j"."
#        vim -e +:g/\$last_record/d -cwq \$basename\$j
#        echo "Done."
#}

#for i in \$(seq 0 \$total); do fixSplit & done

#wait


# removing to avoid output collisions 
mv ${base}.sorted.sam sortedsam.sam
rm ${base}.sam

# remove files of size zero. This happens when all of the records from the last file should be in the second to last file. 
#find . -name "${base}0*" -size  0 -print -delete

"""
}

process generate_coverage {

    // Retry at most 3 times
    //errorStrategy 'retry'
    //maxRetries 3
    
    // Define the Docker container used for this step
    container "quay.io/vpeddu/clomp_containers:latest"

    // Define the input files
    input:
      tuple val(base), file(bam_list)

    // Define the output files
    output:
      tuple val(base), file("*depth.txt")
      file "${base}.depth.txt"

    // Code to be executed inside the task
    script:
      """
      #!/bin/bash
      

      # logging statement  
      ls -lah 

      echo "sorting"

      for i in *.bam 
      do
      samtools sort -@ ${task.cpus} -o \$i.sorted.bam \$i
      done

      #echo "merging"
      #samtools merge ${base}.merged.bam *.sorted.bam 

        
      echo "creating depth file"
      for i in *.sorted.bam 
      do
      samtools depth \$i > \$i.depth.txt
      done



      echo "creating depth file"
      cat *.depth.txt > ${base}.depth.txt

"""
}


process CLOMP_summary {

    // Retry at most 3 times
    errorStrategy 'retry'
    maxRetries 3
    
    // Define the Docker container used for this step
    container "quay.io/fhcrc-microbiome/clomp:v0.1.3"

    // Define the input files
    input:
      tuple val(base), file(bam_file)
      file BLAST_CHECK_DB
      file "kraken_db/"
      file TAXONOMY_DATABASE
    output:
      tuple val(base), file("${base}.*.temp_kraken.tsv")
      tuple val(base), file("*unassigned.txt")
      tuple val(base), file("*assignments.txt")
   

    // Code to be executed inside the task
    script:
    """
#!/usr/bin/env python3

print("Processing BAM file: ${bam_file}")

import uuid
import ast 
import subprocess
import pysam
import glob
import argparse 
import os
import operator
from collections import Counter
from ete3 import NCBITaxa
import timeit
from collections import defaultdict

subprocess.call('ls -latr', shell = True)



ncbi = NCBITaxa(dbfile = 'taxa.sqlite')

# add print statement for getlineage 9606
# print('here', flush = True)
# print(ncbi.get_lineage('9606'), flush = True)


# Make a function to run a shell command and catch any errors
def subprocess_call(cmd):
    return_code = subprocess.call(cmd, shell = True)
    assert return_code == 0, "Exit code %d for %s" % (
        return_code,
        cmd
    )


# The tie-breaking function takes a list of taxid,edit_distances in the form [[taxid,edit_distance],[taxid,edit_distance],..
# references the global variable LOGIC to control the underlying tiebreaking code, breaks ties and 
# returns a single taxid as a string and a Boolean value (as to whether the read needs to re-BLASTEd against host).
def tie_break(taxid_list):
	score_list = [] 
	actual_taxid_list = []
	for id in taxid_list:
		score_list.append(id[1])
	
	# this can filter out any hits that are sufficiently worse than the edit distance for the 'best'
	# alignment - set to zero to only hold scores that have the best edit distance , and increase to a 
	# number greater than the snap option -d to hold all hits
	best_edit_distance = min(score_list) + ${params.EDIT_DISTANCE_OFFSET}
	
	# Keep taxids that have an edit distance less than the acceptable edit distance defined above 
	#i = 0
	#total = len(taxid_list)
	for id in taxid_list:
		#i += 1
		#percent = (i / total) * 100
		#if(percent % 2 == 0):
			#print(percent)
		if id[1] <= best_edit_distance and str(id[0]) != str('4558') and str(id[0]) != str('99802') and str(id[0]) != str('1945254') and str(id[0]) != str('1945390') and str(id[0]) != str('32630'):
		#if id[1] <= best_edit_distance and str(id[0]) not in ${params.IGNORE_TAXA}:
			actual_taxid_list.append(id[0])
	#No longer holding edit distances		
	taxid_list = actual_taxid_list
	lineage_list = []
	
	for id in taxid_list:
		# Not all taxids have valid lineages 
		try:
			#Not every taxid has a valid NCBI lineage, so this statement has to be encased in a try statement.
			lineage = ncbi.get_lineage(id)
			# filters out lineages that contain taxids in the FILTER_LIST variable
			# commonly, this is 'other sequences', 'artificial sequences' or 'environmental samples' 

		except:
			lineage = []
					
		if any(x in ${params.FILTER_LIST} for x in lineage):
			lineage = []
    
		if lineage:
			lineage_list.append(lineage)
	
	if not lineage_list:
		return ['*',False]
	
	# controls if use any alignment to the human genome as grounds for classification as human source
	if "${params.H_STRICT}" == "true":
		# check if H_TAXID ever shows up in 
		if any(int(${params.H_TAXID}) in sl for sl in lineage_list):
			return [${params.H_TAXID},False]
			

	# count all taxids in all lineages 
	taxid_to_count_map = {}
	for each_lineage in lineage_list:
		for each_taxid in each_lineage:
			if each_taxid in taxid_to_count_map:
				taxid_to_count_map[each_taxid] += 1
			else:
				taxid_to_count_map[each_taxid] = 1
	
	#Set the threshold according to the pre-specified LOGIC in the initialization file
	num_assignments = len(lineage_list)
	if "${params.LOGIC}" == 'strict':
		threshold = num_assignments
	elif "${params.LOGIC}" == '90':
		threshold = num_assignments - ((num_assignments /10) + 1) 
	elif "${params.LOGIC}" == 'oneoff':
		threshold = num_assignments - 1
	else:
		print('invalid logic threshold: defaulting to strict')
		threshold = num_assignments
	
	#Now we will find all the taxids that meet threshold/LOGIC specified above.
	surviving_taxids = []
	for taxid_key in taxid_to_count_map:
		# main filtering - everything that passes this step gets a list intersection and the 
		# most specific taxid left is returned 
		if taxid_to_count_map[taxid_key] >= threshold:
			surviving_taxids.append(taxid_key)
			
	if len(surviving_taxids) == 0:
		return ['*',False]

	d = {}

	for the_value in surviving_taxids:
		d[the_value] = len(ncbi.get_lineage(the_value))
	#Find the remaining taxid with the longest lineage.
	#The longest lineage is defined as the longest list.  #nicetohave: this is not pulling the taxonomic rank at any point.
	assigned_hit = max(d.items(), key=operator.itemgetter(1))[0]
	recheck = False 
	
	#Assign a Boolean value for each read as to whether it needs to be searched against a custom BLAST database
	#Here, we are just assigning whether it needs to be searched or not.  The custom BLAST database would need to be made separately.
	#All reads downstream of INCLUSION_TAXID and but not downstream of EXCLUSION_TAXID will be assigned a True value.
	#if "${params.BLAST_CHECK}" == "true":
	#	assigned_lineage = ncbi.get_lineage(assigned_hit)
	#	if ${params.INCLUSION_TAXID} in assigned_lineage and ${params.EXCLUSION_TAXID} not in assigned_lineage:
	#		recheck = True
	
	return [assigned_hit, recheck]


#Every read has a taxid assignment or is unassigned at this point.

# wrapper for a kraken script that converts tab seperated taxid\tcount file and writes a 
# Pavian output file for it. Requires a copy of ncbi's taxonomy database and some blank files
# This is a map of assigned taxids to number of occurrences as well as the total number of unassigned reads.
def new_write_kraken(basename, final_counts_map, num_unassigned):
	print('Preparing output for ${base}')
	# we write a file in the form taxid\tcount 
  # Name the file using a random string in the filename
  # to prevent using the same file name when combining multiple shards
  # in the step immediately after this
	l = open('${base}.%s.temp_kraken.tsv' % str(uuid.uuid4()), 'w')
	
	# initialize with the number of unassigned, we'll need to add human host filtering in earlier
	# because some reads will get tie broken to human 
	
	l.write('0\t' + str(num_unassigned))
	# write the rest of the taxids to the file
	for key in final_counts_map.keys():
		l.write('\\n' + str(key) + '\t' + str(final_counts_map[key]))
	
	# this close is required so we get an EOF before passing the file to kraken-report 
	l.close()
	
	# kraken-report creates a file that Pavian likes - we name the file base_final_report.tsv
	#kraken_report_cmd = '/usr/local/miniconda/bin/krakenuniq-report --db kraken_db --taxon-counts ${base}_temp_kraken.tsv > ${base}_final_report.tsv'
	#subprocess_call(kraken_report_cmd)
	

# takes a list of finished output files and builds sam files for species level assignments 
# for each sample sam files are wrote directly to disk 
def build_sams(input_list):
	
	# only try to import this if we're trying to build sam files 
	from Bio import Entrez 
	Entrez.email = ${params.ENTREZ_EMAIL}
	
	# go through each file report output files 
	for file_name in glob.glob('*report.tsv'):
		base = file_name.split("${params.BASE_DELIMITER}")[0]
		taxid_to_assemble = []
		
		# grab a list of taxids that are at a species level and also have greater than MIN_READ_CUTOFF assigned to them 
		for line in open(file_name):
			line_list = line.split('\t')
			if line_list[3] == 'S' and int(line_lsit[2]) >= ${params.MIN_READ_CUTOFF}:
				lineage = ncbi.get_lineage(line_list[4])
				# filter out any taxids that have lineages that include anything from the blacklist
				if not any(x in "${params.SAM_NO_BUILD_LIST}" for x in lineage):
					taxid_to_assemble.append(line_list[4])
		
		# go through each taxid that we pulled in the last loop and parse the sam file for the accession numbers of the entries that each read aligned to
		# Each read can align to more than one NT entry across all SNAP databases so we grab every accession number that is this taxid 
		for taxid in taxid_to_assemble:
			taxid_search_list = [str(taxid)]
			taxid_search_list = taxid_search_list + ncbi.get_descendant_taxa(taxid, intermediate_nodes=True)
			list_of_reads_to_pull = []
			# this gets a list of every read that was assigned a given taxid 
			for a_line in open(base + '_assignments.txt'):
				a_line_list = a_line.split('\t')
				if a_line_list[1]  in taxid_search_list:
					list_of_reads_to_pull.append(a_line_list[0])
			acc_num_list = []
			# this gets us all accession numbers that all these reads aligned to 
			for s_file in glob.glob(base + '*' + '.sam'):
				for line in open(s_file):
					sam_line_list = sam_line.split('\t')
          # Hard coded for bad klebsiella reference
          print(sam_line_list[2])
					if sam_line_list[0] in list_of_reads_to_pull and 'complete_genome' in sam_line_list[2] and 'CP025541.1' not in sam_line_list[2]:
            print("appended")
						acc_num_list.append(sam_line_list[2].split('.')[0])
			if len(acc_num_list) == 0:
				print('No complete genome reference found, not assembling taxid: ' + str(taxid))
				break
				
			# now we figure out the most common accession number that was assigned to this taxid
			most_common_acc_num = max(set(acc_num_list), key = acc_num_list.count)
			
			taxid_lineage = ncbi.get_lineage(taxid)
			# we walk up the taxonomy tree ASSEMBLY_NODE_OFFSET nodes and pull all reads that were taxonomically assigned at or below that node
			taxid_to_pull = taxid_lineage[ASSEMBLY_NODE_OFFSET]
			taxid_search_list = taxid_to_pull + ncbi.get_descendant_taxa(taxid_to_pull, intermediate_nodes = True)
			
			header_list = []
			seq_list = []
			g = open(base + '_' + taxid + '.fasta', 'w')
			# then we write all the reads that are at or below the taxid_to_pull variable and write them into a fasta file 
			for line in open(base + '_assignments.txt'):
				line_list = line.split('\t')
				if int(line_list[1]) in taxid_search_list:
					g.write('>' + line_list[0] + '\\n')
					g.write(line_list[2])
			g.close()
			
			# now we download the reference fasta file for the most common acession number 
			print('Searching NCBI for Accession number:' + most_common_acc_num + ' for taxid ' + str(taxid))
			record = Entrez.read(Entrez.esearch(db='nucleotide', term=most_common_acc_num))
			try:
				h2 = Entrez.efetch(db='nucleotide', id=record['IdList'][0], rettype='fasta', retmode='text')
			except:
				print(str(taxid) + ' did not return hits - not assembling')
				break
			
			# build some file names for the bowtie index and our reference fasta 
			ref_fasta = base + '_' + str(taxid) + '_ref.fasta'
			ref_db = base + '_' + str(taxid) + '_bwt_db'
			g = open(ref_fasta, 'w')
			g.write(h2.read())
			g.close()
			print('building bowtie2 index') 
			# build the bowtie2 index 
			subprocess_call('bowtie2-build ' + ref_fasta + ' ' + ref_db + ' > /dev/null 2>&1 ')
			print('Done with index build. Aligning...')
			# aling and output the sam file 
			subprocess_call('bowtie2 -x ' + ref_db + ' -@ ' + THREADS + ' -f -U ' + base + '_' + str(taxid) + '.fasta --no-unal > ' + base + '_' + str(taxid) + '.sam')
			subprocess_call('rm ' + ref_db)
				

base_start_time = timeit.default_timer()
read_to_taxids_map = {}
reads_seq_map = {}
#For every SAM file for a given sample, read in the SAM files.
file_start_time = timeit.default_timer()

bam_file = "${bam_file}"
print("Starting to iterate over every line in ${bam_file}")
#For every line in the BAM file
line_count = 0
for line in  open(bam_file):
    line_count += 1
    if line_count > 0:
        #For each read, pull the SAM information for that read.
        line_list = line.split('\t')
        current_read = line_list[0]
        snap_assignment_of_current_read = line_list[2]
        sequence_of_current_read = line_list[9]
        
        #If read is unassigned, call it unassigned and throw it out.  Unassigned reads do not have an edit distance, assign it 100.
        if snap_assignment_of_current_read == '*':
            # higher edit distance than anything else makes sure this gets parsed out 
            current_read_taxid = [snap_assignment_of_current_read,100]
        else:
            #Pull the taxid and the edit distance from each line.
            #print(line)
            current_read_taxid = [snap_assignment_of_current_read.split('#')[-1],
                int(line_list[17].split(':')[-1])]
        #Create map for each sample.
        #The key in each map is the read ID and the values are lists.
        #For every read, append a list of taxid and edit distance from each SAM file.
        if current_read in read_to_taxids_map:
            read_to_taxids_map[current_read].append(current_read_taxid)
        else: 
            # if this is the first time we've seen this read add the sequence to the list
            # and initialize the read -> taxid_list map  
            read_to_taxids_map[current_read] = [current_read_taxid]
            # also store the read and the sequence, this does need to be in a map 
            reads_seq_map[current_read] = sequence_of_current_read

per_base_runtime = str(timeit.default_timer() - base_start_time)
print("${base}" + ' took ' + per_base_runtime + ' in total to read')
final_assignment_counts = defaultdict(int)

#Now we've read all the reads for all the SAM files for one sample.  We are still within the sample For loop here.
#We have all the reads with all the taxids and edit distances that have been assigned.
print('Breaking ties for ' + "${base}")
tie_break_start = timeit.default_timer()
# now we're done with the loop and we have a map with reads to list of taxids assigned to them
g = open("${base}_" + str(uuid.uuid4()) + '_assignments.txt', 'w')
e = open("${base}_" + str(uuid.uuid4()) + '_unassigned.txt','w')
unass_count = 0
taxid_to_read_set = {}

#if "${params.BLAST_CHECK}" == "true":
#    z = open("${base}_recheck.txt", 'w')


for read_key in read_to_taxids_map.keys():
    # Now we need to assign only one taxid to each read, which means we need to run the tie-breaking function. 
    loaded_read = reads_seq_map[read_key]
    
    #Create a results list which is a taxid and a boolean, which answers whether I should re-BLAST this or not.
    r_list = tie_break(read_to_taxids_map[read_key])
    tax_assignment = r_list[0]
    
    # This will only ever be true if BLAST_CHECK is set to True
    if r_list[1]:
        z.write('>' + read_key + '\\n' + loaded_read + '\\n')
        
        
    
    #If the read is unassigned, write it to the unassigned file.
    if tax_assignment == '*':
        e.write('>' + read_key + '\\n' + loaded_read + '\\n')
        unass_count += 1
    # otherwise write it out to the read-by-read assignment file 
    else:
        g.write(read_key + '\t' + str(tax_assignment) + '\t' + loaded_read + '\\n')
        # create a mapping of taxid -> unique reads.  Unique reads are defined as reads without the exact same sequence.  This can help in debugging.
        if "${params.WRITE_UNIQUES}" == "true":
            if str(tax_assignment) in taxid_to_read_set:
                taxid_to_read_set[str(tax_assignment)].add(loaded_read,)
            else:
                taxid_to_read_set[str(tax_assignment)] = set([loaded_read])
            
        if tax_assignment in final_assignment_counts:
            final_assignment_counts[tax_assignment] += 1
        else:
            final_assignment_counts[tax_assignment] = 1

g.close()
e.close()
#if "${params.BLAST_CHECK}" == "true":
#    z.close()
#    subprocess_call('blastn -db ${BLAST_CHECK_DB} -task blastn -query ${base}_recheck.txt -num_threads 20 -evalue ${params.BLAST_EVAL} -outfmt "6 qseqid" -max_target_seqs 1 -max_hsps 1 > blast_check.txt')
#    redo_taxid_list = []
#    for line in open('blast_check.txt'):
#        redo_taxid_list.append(line.split())
#    n = open('new_assignments.txt', 'w')
#    for line in open("${base}" + '_assignments.txt'):
#        ll = line.split('\t')
#        if ll[0] in redo_taxid_list:
#            n.write(ll[0] + '\t' + DB_TAXID + '\t' + ll[2].strip() + '\\n')
#        else:
#            n.write(line)
#    n.close()
#    for item in redo_taxid_list:
#        final_assignment_counts[read_to_taxids_map[item]] += 1
#        final_assignment_counts[DB_TAXID] += 1
    

#For each sample, we make a folder and for every taxid, we create a FASTA file that are named by their taxid.  We lose the read ID in this file.  #nicetohave would be hold the read ID here.
#Here we will write a FASTA of unique reads
if "${params.WRITE_UNIQUES}" == "true":
    subprocess_call('mkdir ' + "${base}".split('.')[0])
    for id in taxid_to_read_set.keys():
        f = open("${base}".split('.')[0] + '/' + str(id) + '_uniques.txt', 'w')
        count = 0
        for item in taxid_to_read_set[id]:
            f.write('>' + str(count) + '\\n')
            f.write(item + '\\n')
            count += 1
        f.close()
    
tie_break_time = str(timeit.default_timer() - tie_break_start)
#print('Tie breaking ' + "${base}" + ' took ' + tie_break_time)


new_write_kraken("${base}" + '_with_host', final_assignment_counts, unass_count)


if "${params.BUILD_SAMS}" == "true":
    build_sams(sam_list)

    """
}

process generate_report {


    //Retry at most 3 times
    errorStrategy 'retry'
    maxRetries 5
    
    // Define the Docker container used for this step
    container "quay.io/fhcrc-microbiome/clomp:v0.1.3"

    // Define the input files
    input:
      tuple val(base), file(kraken_tsv_list), file(unassigned_txt_list), file(assigned_txt_list), val(HOST_FILTERED_READS)
      file BLAST_CHECK_DB
      file "kraken_db/"
      file r1
      file NODES
      file MERGED
      file TAXONOMY_DATABASE
      // val host_filtered_count
     // file r1
    // Define the output files
    output:
      file "${base}.final_report.tsv"
      file "${base}_unassigned.txt"
      file "${base}_assigned.txt"
      file "*metagenome.fastq.gz"
      file "*.clompviz.tsv"
      file "${base}.with_host_final_report.tsv"
      tuple val(base), file("${base}_unassigned.txt")
     // file "*metagenome.fastq.gz"
    // Code to be executed inside the task
    script:
      """
#!/usr/bin/env python3

import glob
import csv
import uuid
import ast 
import subprocess
import pysam
import argparse 
import os
import operator
from collections import Counter
from ete3 import NCBITaxa
import timeit
import sys
import string

from collections import defaultdict

#literally only exists to make test work 
#taxa.sqlite is too big to have in github actions. Workaround is to have it download a new one during tests but use real one in production
ncbi = NCBITaxa(dbfile = 'taxa.sqlite')


host_lineage = ncbi.get_lineage(${params.H_TAXID})
print(host_lineage)
print(type(host_lineage))

# Combine all of the input TSVs into a single file
print('base')
print("${base}")
print('difference')
print(type(${HOST_FILTERED_READS}))
print(${HOST_FILTERED_READS}[0])
host_filtered_reads = ${HOST_FILTERED_READS}[0]
input_files = "${kraken_tsv_list}".split(" ")

unassigned_generate = 'cat ${unassigned_txt_list} > ${base}_unassigned.txt'
assigned_generate = 'cat ${assigned_txt_list} > ${base}_assigned.txt'

subprocess.call(unassigned_generate, shell = True)
subprocess.call(assigned_generate, shell = True)


def tsv_line_to_lst(line):
    row_lst = []
    for number in row.split("\t"):
        row_lst.append(int(number))
    return row_lst


# Reading kraken_temp file to initialize lists of known taxids and associated counts 

path = input_files[0]
first_lst = open(path)  
main_tax_id_lst = []
main_number_lst = []

for row in first_lst:
    row_lst = tsv_line_to_lst(row)
    main_tax_id_lst.append(row_lst[0]) 
    main_number_lst.append(row_lst[1])


# Reading all other kraken_temp files to se
for i in range(1,len(input_files)):
    path = input_files[i]
    lst = open(path)
    for row in lst:
        other_row_lst = tsv_line_to_lst(row)
        if other_row_lst[0] in main_tax_id_lst:
            location = -1
            for i in main_tax_id_lst:
                location = location + 1
                if i == other_row_lst[0]:
                    sum = other_row_lst[1] + main_number_lst[location]
                    main_number_lst[location] = sum      
        else:
            main_tax_id_lst.append(other_row_lst[0])
            main_number_lst.append(other_row_lst[1])


temp_filename = "${base}" + "_kraken_temp_merged.tsv"
print(temp_filename)
with open(temp_filename,'w') as output_file:
    tsv_writer = csv.writer(output_file, delimiter='\t')
    tsv_writer.writerows(zip(main_tax_id_lst, main_number_lst))
output_file.close() 

final_filename = "${base}" + ".final_report.tsv"
kraken_report_cmd = '/usr/local/miniconda/bin/krakenuniq-report --db kraken_db --taxon-counts ' + temp_filename + ' > ' + final_filename
subprocess.call(kraken_report_cmd, shell = True)
subprocess.call("echo FILES  ;ls -latr",shell = True)
# subprocess.call(" mv ${base}.fastq.gz ${base}.metagenome.fastq.gz",shell = True)
subprocess.call(' x=`basename -s ".fastq.gz" *.fastq.gz` ; mv *.fastq.gz \$x.metagenome.fastq.gz',shell = True)


#Adding in phylum designations from nodes.dmp (replaces s/d/k/- with species/domain/kingdom/subspecies in *.clompviz.tsv report)
nodes = open ('nodes.dmp')

taxa = {} 

taxa[0] = 'root'
taxa[1] = 'cellular_organisms'

for line in nodes: 
	accession =  line.strip().split('\t|\t')
	ranks = accession[2].replace(" ", "_")
	taxa[int(accession[0])] = ranks 

# not all taxa are in nodes.dmp anymore 
# some have been merged with other taxa. These are in merged.dmp
merged_taxa = {}
merged_nodes = open('merged.dmp')
for line in merged_nodes: 
	accession =  line.split('\t|')
	merged_taxa[int(accession[0])] = int(accession[1].strip())


files = glob.glob('*.final_report.tsv')
print('files')
print(${params.H_TAXID})
new_host_filter_count = False
for file in files: 
	currentFile = open(file)
	newFileName = os.path.split(file)[1].split('.tsv')[0] + '.clompviz.tsv'
	newFile = open(newFileName, "w")
	for line in currentFile:
		if int(line.split('\t')[4]) in host_lineage: 
			print('host filtered')
			print(host_filtered_reads)
			print('previous count')
			print(int(line.split('\t')[1]))
			new_host_filter_count = int(host_filtered_reads) + int(line.split('\t')[1])
			print('added count')
			print(new_host_filter_count)
		# checks if in nodes.dmp
		try:
			new_phylum = taxa[int(line.split('\t')[4])]
		#checks if in merged.dmp
		except: 
			try: 
				new_phylum = merged_taxa[int(line.split('\t')[4])]
			except: 
				print("warning: "+  str(int(line.split('\t')[4])) +  ' not found')
				new_phylum = 'not_found'
		if new_host_filter_count:

			newline = line.replace(line.split('\t')[1], str(new_host_filter_count))

			new_host_filter_count = False
		else:
			newline = line.replace(line.split('\t')[3], str(new_phylum))
		newFile.write(newline)


# Writing just host filtered output 
with_host_filename = "${base}" + ".with_host_final_report.tsv"
new_host_filter_count = False
for file in files: 
	with_host_file = open(with_host_filename, "w")
	currentFile = open(file)
	for line in currentFile:
		if int(line.split('\t')[4]) in host_lineage: 
				print('host filtered')
				print(host_filtered_reads)
				print('previous count')
				print(int(line.split('\t')[1]))
				new_host_filter_count = int(host_filtered_reads) + int(line.split('\t')[1])
				print('added count')
				print(new_host_filter_count)
		if new_host_filter_count:
			print('here')
			newline = line.replace(line.split('\t')[1], str(new_host_filter_count))
			new_host_filter_count = False
			line = newline
		with_host_file.write(line)



			

"""
}

process blast_unassigned { 

    //Retry at most 3 times
    errorStrategy 'retry'
    maxRetries 3
    
    // Define the Docker container used for this step
    container "quay.io/fhcrc-microbiome/clomp:v0.1.3"

    // Define the input files
    input:
      tuple val(base), file(unassigned_file)
      path BLAST_DB
      file BLAST_UNASSIGNED_SCRIPT
      file "kraken_db/"
      file TAXONOMY_DATABASE

    // Define the output files
    output:
      file "${base}_unassigned_report.tsv" 
      file 'blast_check.txt'

    // Code to be executed inside the task
    script:
      """
      #!/bin/bash
      
      ls -lah

      # script to BLAST unassigned reads
      python3 ${BLAST_UNASSIGNED_SCRIPT} ${base} ${unassigned_file} ${BLAST_DB} ${task.cpus} 1e-4

      # report generation command 
      /usr/local/miniconda/bin/krakenuniq-report --db kraken_db --taxon-counts unassigned_temp_kraken.tsv > ${base}_unassigned_report.tsv



      """


}


process collect_results { 

    //Retry at most 3 times
    errorStrategy 'retry'
    maxRetries 3
    

    // Publish results to params.OUTDIR
    publishDir "${params.OUTDIR}"


    // Define the Docker container used for this step
    container "quay.io/vpeddu/rgeneratesummary:latest"

    // Define the input files
    input:
      file final_reports
      file unassigneds
      file assigneds
      file metagenomes
      file clompviz_tsvs
      file with_host_final_reports
      file depth_files
      
    // Define the output files
    output:
      file "pavian_input/"
      file "metagenomes/"
      file "clompviz/"
      file "assigned"
      file "depth_files/"
      file "unassigned/"

    // Code to be executed inside the task
    script:
      """
      #!/bin/bash
      
      ls -lah

      mkdir pavian_input
      mkdir pavian_input/with_host
      mkdir pavian_input/without_host
      mkdir metagenomes
      mkdir clompviz
      mkdir assigned
      mkdir unassigned
      mkdir depth_files

      mv *with_host_final_report.tsv pavian_input/with_host

      mv *_report.tsv pavian_input/without_host

      mv *metagenome.fastq.gz metagenomes

      mv *.clompviz.tsv clompviz

      mv *_assigned.txt assigned

      mv *_unassigned.txt unassigned

      mv *.depth.txt depth_files

      """
}

process collect_results_with_unassigned { 

    //Retry at most 3 times
    errorStrategy 'retry'
    maxRetries 3
    
    // Define the Docker container used for this step
    container "quay.io/vpeddu/rgeneratesummary:latest"

    // Publish results to params.OUTDIR
     publishDir "${params.OUTDIR}"


    // Define the input files
    input:
      file final_reports
      file unassigneds
      file assigneds
      file metagenomes
      file clompviz_tsvs
      file with_host_final_reports
      file unassigned_tsv
      file depth_files
      
    // Define the output files
    output:
      file "pavian_input/"
      file "metagenomes/"
      file "clompviz/"
      file "assigned/"
      file "unassigned/"
      file "unassigned_blast/"
      file "depth_files/"
    // Code to be executed inside the task
    script:
      """
      #!/bin/bash
      
      ls -lah

      mkdir pavian_input
      mkdir pavian_input/with_host
      mkdir pavian_input/without_host
      mkdir metagenomes
      mkdir clompviz
      mkdir assigned
      mkdir unassigned
      mkdir unassigned_blast
      mkdir depth_files


      mv *_unassigned_report.tsv unassigned_blast

      mv *with_host_final_report.tsv pavian_input/with_host

      mv *_report.tsv pavian_input/without_host

      mv *metagenome.fastq.gz metagenomes

      mv *.clompviz.tsv clompviz

      mv *_assigned.txt assigned

      mv *_unassigned.txt unassigned

      mv *.depth.txt depth_files

      """
}
