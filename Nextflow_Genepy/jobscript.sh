#!/bin/bash
#SBATCH --mem=12g
#SBATCH --ntasks=5
#SBATCH --job-name="iman_genepy"
#SBATCH --ntasks-per-node=5
#SBATCH --nodes=1
#SBATCH --time=48:00:00

module load conda
module load apptainer
source activate Genepy
# cd /iridisfs/hgig/private/in1f24/singularity_genepy/Nextflow/
SCRIPT=$(readlink -f .)
# Extract the full VCF file path
vcf_path=$(grep -oP '^vcf\s*=\s*"\K[^"]+' nextflow.config)
vcf_filename=$(basename "$vcf_path" .vcf.gz)
echo $vcf_filename

mkdir -p $vcf_filename
chmod +x $SCRIPT/templates/pre_1.sh
chmod +x $SCRIPT/templates/pre_local1.sh
chmod +x $SCRIPT/templates/pre_local.sh
chmod +x $SCRIPT/templates/genepy.py
output1=$vcf_filename/Results
out_=$SCRIPT/$output1
mkdir -p $out_
p=$(grep -oP '^vcf\s*=\s*"\K[^"]+' nextflow.config | xargs zcat | cut -f 1 | sort -u | grep --color=auto "^chr" | sed 's/^chr//' | grep -E '^(1?[0-9]|2[0-2]|X|Y|M)$' | \
   awk '{ printf "%s%s", (NR>1?",":""), ($1 == "X" || $1 == "Y" ? "\"" $1 "\"" : $1) } END { printf "\n" }')
 echo $p

CONFIG_FILE="nextflow.config"

LOCAL_CONTAINER_KEY1="cadd_"
LOCAL_CONTAINER_KEY2="vep_"
LOCAL_CONTAINER_KEY3="GATK4"
LOCAL_CONTAINER_KEY4="pyR"

LOCAL_CONTAINER_KEYS=("$LOCAL_CONTAINER_KEY1" "$LOCAL_CONTAINER_KEY2" "$LOCAL_CONTAINER_KEY3" "$LOCAL_CONTAINER_KEY4")

check_local_container() {
    local key="$1"
    if grep -qE "${key} *= *\"/.+\"" "$CONFIG_FILE"; then
        return 0  # Found and not empty
    else
        return 1  # Not found or empty
    fi
}

all_containers_available=true

for key in "${LOCAL_CONTAINER_KEYS[@]}"; do
    if ! check_local_container "$key"; then
        all_containers_available=false
        break
    fi
done

if $all_containers_available; then
    PROFILE="local"
else
    PROFILE="dockerhub"
fi

echo "Selected PROFILE: $PROFILE"

cd $vcf_filename
nextflow run $SCRIPT/nextflow_G.nf --chr $p --output $out_ -work-dir $SCRIPT/$vcf_filename/work -resume --enable report.overwrite -with-dag $out_/dag.png -profile $PROFILE --basedir $SCRIPT

rm -r work
