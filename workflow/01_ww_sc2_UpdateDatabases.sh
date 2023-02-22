#!/bin/bash


# Initialize the analysis
output="$(date +"%Y-%m-%d")"
workflow="all"


# Help
usage() { printf "Usage: $0 \n[-w <string> Workflow to run between: ""'freyja', 'database', 'all' (default: 'all')""][-o <string> (optional) Output directory in 'Results' (example: ""'2023-01-03'"")]\n" 1>&2; exit 1; }



# Arguments
while getopts "w:o:" o; do
    case "${o}" in
        w)
            workflow=${OPTARG}
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



# Check a workflow is correctly defined
if [ $workflow != "all" ] && [ $workflow != "database" ] && [ $workflow != "freyja" ] ; then
    printf "\nWorkflow not (properly) defined\n\n"
    usage
fi

# Check that the directory does not already exist
if [ $workflow != "freyja_run" ]; then
    if [ -d /scratch/projects/SARS-CoV-2/Results/$output/ ]; then
        echo "The directory $output exists already! Are you sure you want to overwrite it? If so, rename the existing folder. This script WON'T overwrite it :-)"
        exit 1
    fi
fi


printf "Parameters used:-w $workflow -o $output\n"


#########################################################################################
##### Set output directories
#########################################################################################
cd /scratch/projects/SARS-CoV-2/Results/
mkdir -f $output
mkdir -f $output/freyja/
mkdir -f $output/databases/
mkdir -f $output/archive/



#########################################################################################
##### Log
#########################################################################################
printf $(date +%Y%m%d_%H%M%S)"\nParameters used:\n-w $workflow\n-o $output\n-s $seq_folder\n" >> /scratch/projects/SARS-CoV-2/Results/$output/archive/01_UpdateDatabase.log
cp "$0" /scratch/projects/SARS-CoV-2/Results/$output/archive/.




#########################################################################################
##### Database (to generate mutation plots)
#########################################################################################


if [ $workflow == "all" ] ||  [ $workflow == "database" ]; then


	## Set working directory
	cd /scratch/projects/SARS-CoV-2/Results/$output/databases/


	## Import resources
	curl -o Database_0_ImportDatabase.R https://raw.githubusercontent.com/wslh-ehd/sc2_wastewater_data_analysis/main/resources/Database_0_ImportDatabase.R
	curl -o Database_1_Usher_to_NextstrainWHO.R https://raw.githubusercontent.com/wslh-ehd/sc2_wastewater_data_analysis/main/resources/Database_1_Usher_to_NextstrainWHO.R
	curl -o Database_2_CovariantMutations.R https://raw.githubusercontent.com/wslh-ehd/sc2_wastewater_data_analysis/main/resources/Database_2_CovariantMutations.R
	curl -o Database_3_Visual_prep.R https://raw.githubusercontent.com/wslh-ehd/sc2_wastewater_data_analysis/main/resources/Database_3_Visual_prep.R
	cp /scratch/projects/SARS-CoV-2/Workflow/Database_Outbreak_VoC_to_explore.txt .

	docker run --rm=True -ti -v $PWD:/data -u $(id -u):$(id -g) r/dashboard:lastest Rscript Database_0_*.R | tee ../archive/R_database_0.log

	git clone https://github.com/hodcroftlab/covariants.git
	awk 'BEGIN{OFS="\t"} {print FILENAME,$0}' ./covariants*/defining_mutations/*.tsv > covariant_mutations.tsv 
	sed -i 's/\.\/covariants\/defining_mutations\///' covariant_mutations.tsv
	sed -i 's/.tsv//' covariant_mutations.tsv
	rm -rf ./covariants*/


	## Generate databases
	docker run --rm=True -ti -v $PWD:/data -u $(id -u):$(id -g) r/dashboard:lastest Rscript Database_1_*.R | tee ../archive/R_database_1.log



fi







#########################################################################################
##### Freyja 
#########################################################################################


if [ $workflow == "all" ] ||  [ $workflow == "freyja" ]; then


	## Set working directory
	cd /scratch/projects/SARS-CoV-2/Results/$output/freyja/
	mkdir ./bootstraps/


	# Log
	docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) staphb/freyja:latest \
	    freyja --version >> ../archive/freyja_version.log 2>&1


	# Import resources
	curl -o Freyja_0_GetLineagesOutbreak.R https://raw.githubusercontent.com/wslh-ehd/sc2_wastewater_data_analysis/main/resources/Freyja_0_GetLineagesOutbreak.R
	curl -o Freyja_1_Aggregate.R https://raw.githubusercontent.com/wslh-ehd/sc2_wastewater_data_analysis/main/resources/Freyja_1_Aggregate.R
	curl -o Freyja_2_Usher_to_NextstrainWHO.R https://raw.githubusercontent.com/wslh-ehd/sc2_wastewater_data_analysis/main/resources/Freyja_2_Usher_to_NextstrainWHO.R
	curl -o Freyja_3_Prep_visual.R https://raw.githubusercontent.com/wslh-ehd/sc2_wastewater_data_analysis/main/resources/Freyja_3_Prep_visual.R
	curl -o Freyja_4_VisualInternalUse.txt https://raw.githubusercontent.com/wslh-ehd/sc2_wastewater_data_analysis/main/resources/Freyja_4_VisualInternalUse.txt
	curl -o Freyja_5_Dashboard.txt https://raw.githubusercontent.com/wslh-ehd/sc2_wastewater_data_analysis/main/resources/Freyja_5_Dashboard.txt
	curl -o outbreak_nomenclature.json -XGET -L https://raw.githubusercontent.com/hodcroftlab/covariants/master/web/data/nameTable.json
	wget https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json
	cp outbreak_nomenclature.json ../databases/.

	docker run --rm=True -ti -v $PWD:/data -u $(id -u):$(id -g) r/dashboard:lastest Rscript Freyja_0_*.R | tee ../archive/R_freyja_0.log


	# Update freyja reference database
	docker run --rm=True -v $PWD:/data -u $(id -u):$(id -g) staphb/freyja:latest \
	    freyja update --outdir . # update lineage database https://github.com/andersen-lab/Freyja
	# Remove recombinants, except XBB
	awk '{ gsub("XBB","_X_B_B",$1); print $1 }' usher_barcodes.csv | grep -v "^X" | awk '{ gsub("_X_B_B","XBB",$1); print $1 }' > usher_barcodes_withRecombinantXBBonly.csv

fi

