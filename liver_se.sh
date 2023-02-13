#!/bin/bash
##some maintenance stuff
mkdir -p ../adapters mkdir -p ../aligned ../aligned/liver_se ../aligned/liver_se/counts/ 
mkdir -p ../trimmed ../trimmed/liver_se
##specify file locations - not necessary for binaries that are installed already
reads_path='/nfs1/CGRB/Bionaz_Lab/RNAseq/Feb2020_RNASeq/data/liver_se'
aligned_path='/nfs1/CGRB/Bionaz_Lab/RNAseq/Feb2020_RNASeq/aligned/liver_se'
#ensure nothing else is in this folder
for file in $(ls ${reads_path}) 
do
	filename=$(echo $file | egrep -o '^[^.]+')
	echo "working on file $filename"
	echo $filename | tr '-' '\t' | awk '{print ">" $4"_""SE" "\n" "AAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"$4"ATCTCGTATGCCGTCTTCTGCTTG"}' > ../adapters/${filename}-adapters.fa
	echo $filename | tr '-' '\t' | awk '{print ">" $4"_""SE_rc" "\t""AAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"$4"ATCTCGTATGCCGTCTTCTGCTTG"}' | python /nfs1/CGRB/Bionaz_Lab/RNAseq/programs/reverse.py >> ../adapters/${filename}-adapters.fa
	#Trim and remove adapters --works
	java -jar /nfs1/CGRB/Bionaz_Lab/RNAseq/programs/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 -phred33 \
	/nfs1/CGRB/Bionaz_Lab/RNAseq/Feb2020_RNASeq/data/liver_se/${filename}.fastq.gz /nfs1/CGRB/Bionaz_Lab/RNAseq/Feb2020_RNASeq/trimmed/liver_se/${filename}_trimmed.fastq.gz \
	ILLUMINACLIP:../adapters/${filename}-adapters.fa:2:30:10 \
	LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:3
	echo "$filename trimmed"
	#unzip files
	new_path=/nfs1/CGRB/Bionaz_Lab/RNAseq/Feb2020_RNASeq/trimmed/liver_se/${filename}_trimmed.fastq
	echo $new_path
	gunzip $new_path.gz
	#align with hisat2
	hisat2 --threads 4 --dta -x /nfs1/CGRB/Bionaz_Lab/RNAseq/Feb2020_RNASeq/refs/bt_idx_final/bt_idx_final -U $new_path -S $aligned_path/${filename}_trimmed.sam
	#convert to BAM
	samtools sort -@ 4 -o $aligned_path/${filename}_trimmed.bam $aligned_path/${filename}_trimmed.sam
	#count with stringtie
	stringtie -p 4 -e -G /nfs1/CGRB/Bionaz_Lab/RNAseq/Feb2020_RNASeq/refs/Bos_taurus.ARS-UCD1.2.99.gtf -o $aligned_path/counts/${filename}_trimmed.gtf $aligned_path/${filename}_trimmed.bam
	#rezip the file
	gzip $new_path 
	#remove sam and bam files 
	rm -rf $aligned_path/${filename}_trimmed.bam
	rm -rf $aligned_path/${filename}_trimmed.sam
done
