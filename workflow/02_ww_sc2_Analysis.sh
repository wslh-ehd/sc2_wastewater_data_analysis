#!/bin/bash


# Prerequisites
workflow="all"


# Help
usage() { printf "Usage: \n[-w <string> Workflow to use (default: 'all'): ""'freyja' (process all demix files in SARS-CoV-2 directory), 'database' (process all samples in the SARS-CoV-2 directory), 'all' (freyja & database),  'database_freyja_run' (process all samples in the SARS-CoV-2 directory and process demix files only in the specified SequencingID run), 'final_analysis' (perform final analysis for both database and freyja)""]\n[-o <string> Output directory in 'Results' (example: ""'2023-01-03'"")]\n[-s <string> (if 'freyja_run' specified) SequencingID (example: ""'Seq078'"")]\n" 1>&2; exit 1; }



# Arguments
while getopts "w:s:o:" o; do
    case "${o}" in
        w)
            workflow=${OPTARG}
            ;;
        s)
            seq_folder=${OPTARG}
            ;;
        o)
            output=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))



# Check the specified directory exists
if [[ ! -d /scratch/projects/SARS-CoV-2/Results/$output/ ]] ; then
            printf "\n'The directory $output does not exist! Please indicate an output directory (-o)\n\n"
            usage
    fi


# Check a workflow is correctly defined
if [[ $workflow != "all" ]] && [[ $workflow != "database" ]] && [[ $workflow != "freyja" ]] && [[ $workflow != "database_freyja_run" ]] && [[ $workflow != "final_analysis" ]] ; then
    printf "\nWorkflow not (properly) defined\n\n"
    usage
fi



# Check that all arguments are correctly defined if 'freyja_run' is selected
if [[ $workflow == "database_freyja_run" ]]; then
    if [ -z $seq_folder ]; then
        printf "\n'freyja_run' workflow has been selected, but no SequencingID run (-s) have been defined. Sad :-(\n\n"
        usage
    fi
fi


printf "Parameters used:-w $workflow -o $output -s $seq_folder\n"




#########################################################################################
##### Log
#########################################################################################
printf $(date +%Y%m%d_%H%M%S)"\nParameters used:\n-w $workflow\n-o $output\n-s $seq_folder\n" >> /scratch/projects/SARS-CoV-2/Results/$output/archive/02_DataAnalysis.log
cp "$0" /scratch/projects/SARS-CoV-2/Results/$output/archive/.




#########################################################################################
##### Database (to generate mutation plots)
#########################################################################################


if [[ $workflow == "database" ]] || [[ $workflow == "all" ]] ||  [[ $workflow == "database_freyja_run" ]] ; then

    # Update data
    cd /scratch/projects/SARS-CoV-2/

    awk 'BEGIN{OFS="\t"} {print FILENAME,$0}' Seq*/variants/ivar/*_depth.tsv > ./Results/$output/databases/CallDepthCompiled.tsv
    sed -i 's/\/variants\/ivar\//\t/' ./Results/$output/databases/CallDepthCompiled.tsv
    sed -i 's/_depth.tsv//' ./Results/$output/databases/CallDepthCompiled.tsv

    awk 'BEGIN{OFS="\t"} {print FILENAME,$0}' Seq*/variants/ivar/*_notfiltered.tsv > ./Results/$output/databases/CallVariantALLCompiled.tsv
    sed -i 's/\/variants\/ivar\//\t/' ./Results/$output/databases/CallVariantALLCompiled.tsv
    sed -i 's/_notfiltered.tsv//' ./Results/$output/databases/CallVariantALLCompiled.tsv
fi




if [[ $workflow == "database" ]] || [[ $workflow == "all" ]] || [[ $workflow == "final_analysis" ]] ; then
   
    # Analysis
    cd /scratch/projects/SARS-CoV-2/Results/$output/databases/

    cp ../archive/ListSamples.xlsx .

    docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) -w /data r/dashboard:lastest Rscript Database_2_*.R | tee ../archive/R_database_2.log
    docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) -w /data r/dashboard:lastest Rscript Database_3_*.R | tee ../archive/R_database_3.log

fi





#########################################################################################
##### Freyja 
#########################################################################################




##### Freyja - process all samples in Seq* folders

if [[ $workflow == "freyja" ]] || [[ $workflow == "all" ]]; then


    cd /scratch/projects/SARS-CoV-2/
    mkdir ./Results/$output/freyja/bootstraps/
    
    # Identify SNPs/barcodes
    
    # Less1Year_SeqRuns=$(find Seq*/fastq -type d -mtime -365 | cut -f1 -d"/" | uniq | tr " " "\n" )
    # echo $Less1Year_SeqRuns > Results/2022-01-01_FreyjaFiles_DoNotTouch/Less1Year_SeqRuns.txt
    # sed -i 's/\s\+/\n/g' Results/2022-01-01_FreyjaFiles_DoNotTouch/Less1Year_SeqRuns.txt
    
    while read -r line; do
        for variant in $line/freyja/*.tsv; do
            depth=${variant/-variant.tsv/-depth}
            out=${variant/-variant.tsv/}; out=${out/\/freyja\//@};
            echo $out
            docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) staphb/freyja:latest freyja boot $variant $depth --nt 15 --nb 10 --output_base ./Results/$output/freyja/bootstraps/$out --barcodes ./Results/$output/freyja/usher_barcodes_withRecombinantXBBonly.csv
        done
    done < ./Results/2022-01-01_FreyjaFiles_DoNotTouch/Less1Year_SeqRuns.txt

    
    
    
#    ### Run Freyja analysis on all variants incl. all recombinants - 200 last samples
#    mkdir ./Results/$output/freyja/bootstraps_X/
#    
#    # Identify SNPs/barcodes
#    ls Seq*/freyja/*.tsv | tail -n 200 | while read variant; do 
#        depth=${variant/-variant.tsv/-depth}
#        out=${variant/-variant.tsv/}; out=${out/\/freyja\//@};
#        echo $out
#        docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) staphb/freyja:latest freyja boot $variant $depth --nt 15 --nb 10 --output_base ./Results/$output/freyja/bootstraps_X/$out --barcodes ./Results/$output/freyja/usher_barcodes.csv
#    done

fi




##### Freyja - if just want to process data for a specific sequencing run

if [[ $workflow == "database_freyja_run" ]] ; then

    cd /scratch/projects/SARS-CoV-2/
    mkdir ./Results/$output/freyja/bootstraps/

    # Identify SNPs/barcodes
    for variant in $seq_folder/freyja/*.tsv; do
        depth=${variant/-variant.tsv/-depth}
        out=${variant/-variant.tsv/}; out=${out/\/freyja\//@};
        echo $out
        docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) staphb/freyja:latest freyja boot $variant $depth --nt 15 --nb 10 --output_base ./Results/$output/freyja/bootstraps/$out --barcodes ./Results/$output/freyja/usher_barcodes_withRecombinantXBBonly.csv
        #freyja demix $variant $depth --output ./Results/$output/freyja/demix/$out-results.txt
    done
fi




##### Freyja - curation data and visual preparations

if [[ $workflow == "freyja" ]] || [[ $workflow == "all" ]] || [[ $workflow == "final_analysis" ]] ; then

    # Process data for visualization
    cd /scratch/projects/SARS-CoV-2/Results/$output/freyja/

    cp ../archive/ListSamples.xlsx .

    docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) -w /data r/dashboard:lastest Rscript Freyja_1_*.R | tee ../archive/R_freyja_1.log
    docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) -w /data r/dashboard:lastest Rscript Freyja_2_*.R | tee ../archive/R_freyja_2.log
    docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) -w /data r/dashboard:lastest Rscript Freyja_3_*.R | tee ../archive/R_freyja_3.log
    docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) -w /data r/dashboard:lastest Rscript Freyja_4_*.R | tee ../archive/R_freyja_4.log
    docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) -w /data r/dashboard:lastest Rscript Freyja_5_*.R | tee ../archive/R_freyja_5.log

fi


