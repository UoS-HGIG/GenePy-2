
# Genepy2 Pipeline on phase 2 local UKBB data

This guide assumes you have `apptainer` and `conda` installed on your machine. Before running the Genepy pipeline, we need to prepare the necessary tools and environment.

## Step 1: Downloading the Repository

### Using `curl`
Download the repository as a zip file:
```bash
curl -L -o repo.zip https://github.com/UoS-HGIG/GenePy-2/archive/refs/heads/main.zip

unzip repo.zip
```

### Using `wget`
Alternatively, you can use `wget`:
```bash
wget -O repo.zip https://github.com/UoS-HGIG/GenePy-2/archive/refs/heads/main.zip
unzip repo.zip
```

## Step 2: Extracting the Files
Navigate into the extracted repository:
```bash
cd GenePy-2-main/GenePy2_UKBiobank/Nextflow_Genepy2_UKBB_V3/
```

## Step 3: Create Conda Environment
Create a Conda environment based on the `.yml` file provided in the repository:
```bash
conda env create --file Genepy.yml
source activate Genepy
```

## Step 4: Downloading Annotation Files
Create a data directory and navigate into it:
```bash
mkdir -p data
cd data
```

### Pre-scored Annotation Files
Download pre-scored annotation files in the background:
```bash
wget -c https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz &
tabix -p vcf whole_genome_SNVs.tsv.gz
```

Download additional annotation files:
```bash
wget -c https://kircherlab.bihealth.org/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz &
tabix -p vcf gnomad.genomes.r3.0.indel.tsv.gz
```

### Other Annotation Files
Download and extract the remaining annotation files:
```bash
wget -c https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/annotationsGRCh38_v1.6.tar.gz &
tar -xzf annotationsGRCh38_v1.6.tar.gz
```

### VEP Database for Homo Sapiens
Download the VEP database for Homo Sapiens and extract it:
```bash
nohup curl -O https://ftp.ensembl.org/pub/release-111/variation/indexed_vep_cache/homo_sapiens_vep_111_GRCh38.tar.gz &
tar -xzf homo_sapiens_vep_111_GRCh38.tar.gz
```

## Step 5: Configure the Pipeline

Navigate back to the main directory:
```bash
cd ..
```

Open the `nextflow.config` file with an editor (e.g., nano):
```bash
nano nextflow.config
```

### Configuration Parameters
Edit the following parameters in the `nextflow.config` file as needed:

#### Input VCF File
Specify the path to your VCF file:
```plaintext
vcf = "path/to/your.vcf.gz"
```

#### CADD Annotation
Specify the path to the CADD annotation files:
```plaintext
annotations_cadd = "${basedir}/data/GRCh38_v1.6/"
```

#### Homo Sapiens VEP Database
Specify the path to the Homo Sapiens VEP database:
```plaintext
homos_vep = "${basedir}/data/homo_sapiens/"
```

#### BED Interval File
Specify the path to the BED interval file:
```plaintext
bed_int = "${basedir}/templates/CCDS_hg38_pad25_sorted.bed"
```

#### Containers Location
Keep these fields empty if this is your first time running the pipeline (containers will be downloaded automatically):
```plaintext
cadd_  = ""
vep_   = ""
pyR    = ""
```

#### VEP Plugins
Specify the paths to the VEP plugin files:
```plaintext
plugin1 = "${basedir}/data/whole_genome_SNVs.tsv.gz"
plugin2 = "${basedir}/data/gnomad.genomes.r3.0.indel.tsv.gz"
vep_plugins = "${basedir}/templates/plugins/"
```


#### Other Configuration
Specify the paths to additional configuration files:
```plaintext
header_meta = "${basedir}/header.meta_org"
genepy_py = "${params.basedir}/templates/genepy.py"
```

## Running the Job Script

Finally, run the job script with the Slurm job scheduler:
```bash
chmod +x split_chr.sh 
sbatch ./split_chr.sh /path/input/vcf_file_folder/ chr_number
```

If you use a different job scheduler, update the job script header and the `process.executor` parameter in the `nextflow.config` file accordingly.

## Optimization of Genepy.py Pipeline and Picard Tools for High I/O and RAM Usage

The Genepy.py pipeline and Picard tools in the main workflow can demand significant I/O and RAM, particularly when processing large VCF files. To mitigate computational burden and prevent system overload, please adjust the following parameters:

1. Limiting Parallel Jobs in Nextflow:

In the nextflow_G.nf file located in the /Genepy_wf/ directory, set the maxForks parameter to 5. This limits the number of job scripts that can run in parallel at any given time, reducing the overall load on the system.
// Example snippet in nextflow_G.nf
maxForks 5

2. Adjusting Parallel Chunk Processing in Picard Tools:

In the Nextflow_G.vf script within the main folder, modify the parallel processing of chunks by adjusting the xargs command. Change the -P parameter to 2, which will limit the number of chunks processed in parallel by Picard tools.

find chunks -name '*.vcf.gz' | xargs -n 1 -P 2 -I {} bash -c 'process_chunk "{}"'

Here, -P 2 ensures that only two chunks are processed simultaneously, helping to manage memory and I/O usage more effectively.





By following these steps, you will have set up and run the Genepy pipeline on your system. This pipeline will help you analyze genomic data efficiently using the provided tools and configurations.

### Explanation of the Commands

1. **Downloading the repository**:
   - `curl -L -o repo.zip URL` and `wget URL -O repo.zip`: Download the repository as a zip file.
   - `unzip repo.zip`: Extract the contents of the zip file.

2. **Creating Conda environment**:
   - `conda env create --file conda_lib.yml`: Create a Conda environment using the specified `.yml` file.
   - `source activate Genepy`: Activate the newly created Conda environment.

3. **Downloading necessary files**:
   - `wget -c URL &`: Download files in the background.
   - `tabix -p vcf file`: Index the downloaded VCF files using `tabix`.

4. **Configuring the pipeline**:
   - `nano nextflow.config`: Open the `nextflow.config` file for editing.

5. **Running the job script**:
   - `sbatch jobscript.sh`: Submit the job script to the Slurm job scheduler.
