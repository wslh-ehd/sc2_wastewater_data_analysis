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




cd /scratch/projects/SARS-CoV-2/Results/$output/


### Display dashboard's visuals
git clone https://github.com/wslh-data/sc2-wastewater-data-dashboard
rm -rf sc2-wastewater-data-dashboard/data

# Test github plots before uploading the data to Github
cd sc2-wastewater-data-dashboard
for dir in ./*/; do
cp ../freyja/DashboardData.RData $dir/.
done

# Change the path to DashboardData.RData
sed -i 's/url(file_url)/"DashboardData.RData"/' ./ww*/app.R

# Run docker
docker build -t dashboard_local .
docker run -p 8080:3838 dashboard_local
#http://localhost:8080/

