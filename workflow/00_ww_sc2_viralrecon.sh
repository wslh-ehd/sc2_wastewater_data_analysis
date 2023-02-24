#!/bin/bash

### Prerequires
# nextflow
# docker
# docker pull staphb/freyja
# docker pull staphb/ivar

## Help
usage() { echo "Usage: $0 [-s <string> example: ""'Seq078'""]" 1>&2; exit 1; }


## Arguments
while getopts ":s:" o; do
    case "${o}" in
        s)
            seq_folder=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${seq_folder}" ]; then
    usage
fi



## Check seq_folder does not already exist
if [ -d /scratch/projects/SARS-CoV-2/$seq_folder/ ]; then
	echo "The directory $seq_folder exists already! Are you sure you want to overwrite it? If so, rename the existing folder. This script WON'T overwrite it :-)"
	exit 1
fi



## Prepare the working directory
mkdir /scratch/projects/SARS-CoV-2/$seq_folder/
mkdir /scratch/projects/SARS-CoV-2/$seq_folder/fastq/


## Move to Seq### Working directory
cd /scratch/projects/SARS-CoV-2/$seq_folder/


## Import the fastq, config, samplesheet and this script
mv /scratch/projects/SARS-CoV-2/CopyFastqHere/*.gz ./fastq/.
# Move samplesheet.csv to the working directory
mv /scratch/projects/SARS-CoV-2/CopyFastqHere/samplesheet.csv samplesheet_$seq_folder.csv
# Update the path within the file with the output directory 
sed -i 's/seq_folder/'${seq_folder}'/' samplesheet_$seq_folder.csv; sed -i 's/seq_folder/'${seq_folder}'/' samplesheet_$seq_folder.csv
curl -o custom_ww_viralrecon.config https://raw.githubusercontent.com/wslh-ehd/sc2_wastewater_data_analysis/main/resources/custom_ww_viralrecon.config
# Download kraken2 (human + phiX) database  - workflow to obtain the files described here: entire workflow to generate these files: https://hackmd.io/@AstrobioMike/kraken2-read-filtering#Download-database-as-built-on-11-Sept-2020-LATEST
curl -L -o kraken2_human_and_phiX_db.tar.gz https://ndownloader.figshare.com/files/24658262 
tar -xzvf kraken2_human_and_phiX_db.tar.gz
# Save the current script into the output directory
cp "$0" script_prep_viralrecon.sh


## Get ready for the next run
curl -o /scratch/projects/SARS-CoV-2/CopyFastqHere/samplesheet.csv https://raw.githubusercontent.com/wslh-ehd/sc2_wastewater_data_analysis/main/data/samplesheet.csv


## Run viralrecon
printf "Sequences processing using viralrecon (nextflow) workflow\n\n" 
nextflow pull nf-core/viralrecon # Make sure viralrecon is up to date 
nextflow run nf-core/viralrecon --input samplesheet_$seq_folder.csv \
	--outdir /scratch/projects/SARS-CoV-2/$seq_folder/ \
	--platform illumina \
	--protocol amplicon \
	--genome MN908947.3 \
	--primer_set qiaseq \
	--primer_set_version '1' \
	--kraken2_db /scratch/projects/SARS-CoV-2/$seq_folder/kraken2_human_and_phiX_db/ \
	--kraken2_db_name 'human_phiX' \
	--kraken2_variants_host_filter true \
	-profile docker,test \
	-c custom_ww_viralrecon.config


## Remove kraken2 (human + phiX) database
rm -rf kraken2_human_and_phiX_db.tar.gz
rm -rf kraken2_human_and_phiX_db


## Copy/Extract depth and all mutations info and put all results in one place
printf 'Viralrecon workflow is done!\n\n\nNow, it is time to gather all important files into the "output" directory\n\n'
mkdir output
mkdir output/freyja


## Copy variant (filtered info) into the output directory
printf 'Copy variant calls files (after filtering, from viralrecon) into the "output" directory\n\n'
for tsv in variants/ivar/*.tsv; do
	out=${tsv/.tsv/_filtered.tsv}
	out=${out/variants\/ivar\//}
	cp $tsv ./output/$out
done

## Extract depth and copy into the output directory
printf 'Extract the depth (before filtering) into the "output" directory\n\n'
for mpileup in variants/ivar/*.mpileup; do
	out=${mpileup/.mpileup/_depth.tsv}; out=${out/variants\/ivar\//};
	cat $mpileup | cut -f1-4 > ./output/$out
done

## Extract gg and fasta into the output directory
printf 'Copy the files used to do the alignment/variant calling into the "output" directory\n\n'
find ./work/ -type f -name '*GCA_009858895.3_ASM985889v3_genomic.200409.gff' |  while read P; do cp "$P" ./output/; done
find ./work/ -type f -name '*nCoV-2019.reference.fasta' |  while read P; do cp "$P" ./output/; done


printf 'Save ivar and freyja versions (used below) into "pipeline_info" directory\n\n'
docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) staphb/freyja:latest \
    freyja --version >> ./pipeline_info/freyja_version.log 2>&1
docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) staphb/ivar:latest \
    ivar -v >> ./pipeline_info/ivar_version.log 2>&1



## Extract all SNPs informations
printf 'Variant calls files (no filtering) into the "output" directory'
docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) staphb/ivar:latest \
  /bin/bash -c 'for mpileup in variants/ivar/*.mpileup; do  out=${mpileup/.mpileup/_notfiltered.tsv}; out=${out/variants\/ivar\//}; cat $mpileup | ivar variants -t 0.0001 -q 20 -m 0 -g ./output/*.gff -r ./output/*.fasta -p ./output/$out; done'



## Freyja: Identify SNPs
printf 'Create BAM files for Freyja into the "output/freyja" directory\n'
for sample in ./variants/bowtie2/*.ivar_trim.sorted.bam; do
out=${sample/.ivar_trim.sorted.bam/}; out=${out/variants\/bowtie2\//}
docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) staphb/freyja:latest \
    freyja variants $sample --variants ./output/freyja/$out-variant --depths ./output/freyja/$out-depth --ref ./output/nCoV-2019.reference.fasta
done



printf "\n\n\nThe processing of the fastq from $seq_folder is done.\n\n\n"
printf "Please check the QC of the samples in the file located: scratch > projects > SARS-CoV-2 > $seq_folder > multiqc > multiqc_report.html\n"
printf "Copy and paste the content of the csv file (scratch > projects > SARS-CoV-2 > $seq_folder > multiqc > summary_variants_metrics_mqc.csv) into SampleList.xlsx."
printf "\n\n\nGreat job, you rock!\n\n"

