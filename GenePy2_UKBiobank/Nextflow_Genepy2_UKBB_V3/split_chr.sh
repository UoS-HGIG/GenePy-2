#!/bin/bash
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --job-name="iman_genepy"
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --time=48:00:00
###################################################### sbatch split_chr.sh /iridisfs/hgig/private/in1f24/singularity_genepy/Nextflow/Ukbb200k/vcf_files/chr22.23.24 chr22
# v2=x.p2
# v5=x.p1
# v4=x.v1_p1.vep.vcf
# v3=x.v1_p12.vcf
# v6=x.v1_p12_f1.vcf
# v7=x.v1_p12_filterNew_onTarget.vcf.gz


module load conda
module load apptainer
source activate Genepy

SCRIPT=$(readlink -f .)

chmod 770 $SCRIPT/templates/pre_1.sh
chmod 770 $SCRIPT/templates/genepy.py 
chmod 770 $SCRIPT/templates/pre_2.sh
chmod 770 $SCRIPT/templates/duplic.sh


# Directory where your VCF files are located
vcf_directory=$1
# Directory for results
output_directory="$SCRIPT/output_results"
mkdir -p $output_directory

chrom=$2
echo "Processing chromosome: $chrom"


#######################
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


###########################
# Get all VCF files related to the current chromosome
vcf_files=$(ls $vcf_directory/*_${chrom}_*.vcf.gz)
echo "$vcf_files"
chr_output_directory="$output_directory/$chrom"
mkdir -p $chr_output_directory
#echo "$(realpath $chr_output_directory)" >> $folder_list_file
cd $chr_output_directory 
     #Run Nextflow pipeline for the current VCF file
nextflow run $SCRIPT/nextflow_G.nf -c $SCRIPT/nextflow.config --chr $chrom --vcf $vcf_directory  --output $chr_output_directory -work-dir $chr_output_directory/work --enable report.overwrite -with-dag $chr_output_directory/dag.png -profile $PROFILE --basedir $SCRIPT -resume
    



