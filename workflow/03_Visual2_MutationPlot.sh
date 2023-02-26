#!/bin/bash


# Help
usage() { printf "Usage: $0 \n\n[-o <string> Output directory in 'Results' (example: ""'2023-01-03'"")]\n" 1>&2; exit 1; }



# Arguments
while getopts ":o:" o; do
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

if [[ (-z "${output}") ]]; then
    usage
fi




cd /scratch/projects/SARS-CoV-2/Results/$output/database/


### Display dashboard's visuals
git clone https://github.com/wslh-ehd/sc2-wastewater-data-mut-visual


# Test github plots before uploading the data to Github
cd sc2-wastewater-data-mut-visual

for dir in ./*/; do
cp ../plot.RData $dir/.
done


# Run docker
docker build -t visual_database_local .
docker run -p 8080:3838 visual_database_local
#http://localhost:8080/

