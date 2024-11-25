#!/bin/bash
#SBATCH --mem=48G
#SBATCH --ntasks=1
#SBATCH --job-name="$2"
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1
#SBATCH --time=60:00:00

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
folder_list_file="${chrom}_folder.lst"
> $folder_list_file  # Clear the file if it already exists

counter=0
batch_size=10
for vcf_file in $vcf_files; do
    
    vcf_filename=$(basename "$vcf_file" .vcf.gz)
    #chunk=$(echo "$vcf_filename" | grep -oP '[0-9]+_[0-9]+')

    # Create a directory for the chromosome-specific and chunk-specific results
    chr_chunk_output_directory="$output_directory/$vcf_filename"
    mkdir -p $chr_chunk_output_directory
    echo "$(realpath $chr_chunk_output_directory)" >> $folder_list_file
    cd $chr_chunk_output_directory
     #Run Nextflow pipeline for the current VCF file
    nohup nextflow run $SCRIPT/nextflow_G.nf -c $SCRIPT/nextflow.config --chr $chrom --vcf $vcf_file  --output $chr_chunk_output_directory -work-dir $chr_chunk_output_directory/work --enable report.overwrite -with-dag $chr_chunk_output_directory/dag.png -profile $PROFILE --basedir $SCRIPT -resume &
    cd $SCRIPT
    
    
    # Increment the counter
    counter=$((counter + 1))

    # If the counter reaches the batch size, wait for the current batch to complete
    if [[ $((counter % batch_size)) -eq 0 ]]; then
        echo "Waiting for the current batch of $batch_size jobs to start..."
        wait  # This waits for all background jobs to finish before starting the next batch
    fi
done

# Wait for all parallel processes for the current chromosome to finish
wait
$SCRIPT/templates/duplic.sh $folder_list_file $chrom $output_directory
echo "Completed processing for chromosome: $chrom part1"
output_G=genepy/${chrom}
mkdir -p $output_G/15
mkdir -p $output_G/20
mkdir -p $output_G/ALL

echo "running part2..."

cd $SCRIPT/${output_G}/15
echo $SCRIPT
nextflow run $SCRIPT/Genepy_wf/nextflow_G.nf  -profile $PROFILE --basedir $SCRIPT -c $SCRIPT/Genepy_wf/nextflow.config --output $SCRIPT/${output_G}/15 -work-dir $SCRIPT/${output_G}/15/work --enable report.overwrite -with-dag $SCRIPT/genepy/dag_${chrom}_15.png --input.list ${output_directory}/${chrom}_final15.lst --chr $chrom 


cd $SCRIPT/${output_G}/20
echo $SCRIPT
nextflow run $SCRIPT/Genepy_wf/nextflow_G.nf  -profile $PROFILE --basedir $SCRIPT -c $SCRIPT/Genepy_wf/nextflow.config --output $SCRIPT/${output_G}/20 -work-dir $SCRIPT/${output_G}/20/work --enable report.overwrite -with-dag $SCRIPT/genepy/dag_${chrom}_20.png --input.list ${output_directory}/${chrom}_final20.lst --chr $chrom


cd $SCRIPT/${output_G}/ALL
echo $SCRIPT
nextflow run $SCRIPT/Genepy_wf/nextflow_G.nf  -profile $PROFILE --basedir $SCRIPT -c $SCRIPT/Genepy_wf/nextflow.config --output $SCRIPT/${output_G}/ALL -work-dir $SCRIPT/${output_G}/ALL/work --enable report.overwrite -with-dag $SCRIPT/genepy/dag_${chrom}_ALL.png --input.list ${output_directory}/${chrom}_finalALL.lst --chr $chrom



