#!/bin/bash


# Help
usage() { printf "Usage: $0 \n\n[-o <string> Output directory in 'Results' (example: ""'2023-01-03'"")]\n" 1>&2; exit 1; }



# Arguments
while getopts ":s:o:" o; do
    case "${o}" in
        o)
            output=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [[ (-z "${seq_folder}") || (-z "${output}") ]]; then
    usage
fi




cd /scratch/projects/SARS-CoV-2/Results/$output/freyja/


### Display dashboard's visuals
git clone https://github.com/wslh-ehd/sc2-wastewater-data-freyja-visual


# Test github plots before uploading the data to Github
cd sc2-wastewater-data-freyja-visual

for dir in ./*/; do
cp ../InternalUseData.RData $dir/.
done


# Run docker
docker build -t visual_freyja_local .
docker run -p 8080:3838 visual_freyja_local
#http://localhost:8080/

