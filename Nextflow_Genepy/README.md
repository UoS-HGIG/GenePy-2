
# Genepy Pipeline Running

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
cd GenePy-2-main/Nextflow_Genepy/
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

#### If working on Iridis6
Add -f (for fakeroot) to every line starting with containerOptions e.g.
```plaintext
 containerOptions = "-B ${params.annotations_cadd}:/opt/CADD-scripts-CADD1.6/data/annotations/GRCh38_v1.6/ -f"
```

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
GATK4  = ""
pyR    = ""
```

#### VEP Plugins
Specify the paths to the VEP plugin files:
```plaintext
plugin1 = "${basedir}/data/whole_genome_SNVs.tsv.gz"
plugin2 = "${basedir}/data/gnomad.genomes.r3.0.indel.tsv.gz"
vep_plugins = "${basedir}/templates/Plugins"
```

#### CADD Filter
Choose the CADD filter criteria (`15`, `20`, or `ALL`):
```plaintext
cadd_filter = "ALL"
```

#### Other Configuration
Specify the paths to additional configuration files:
```plaintext
header_meta = "${basedir}/header.meta_org"
genepy_py = "${basedir}/templates/genepy.py"
```

Specify the chunk size for gene calculations:
```plaintext
chunk_size = 1000
```

## Running the Job Script

Finally, run the job script with the Slurm job scheduler:
```bash
sbatch jobscript.sh
```

If you use a different job scheduler, update the job script header and the `process.executor` parameter in the `nextflow.config` file accordingly.

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
