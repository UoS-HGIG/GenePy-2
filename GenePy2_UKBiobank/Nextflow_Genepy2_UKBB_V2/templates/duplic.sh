#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 folder_list.txt chromosome_name output_path" 
    exit 1
fi

folder_list="$1"
chromosome_name="$2"
output_path="$3"
# Create a directory for concatenated genes
OUTPUT_FOLDER_15="$output_path/concatenated_genes_${chromosome_name}/metafiles15"
OUTPUT_FOLDER_20="$output_path/concatenated_genes_${chromosome_name}/metafiles20"
OUTPUT_FOLDER_ALL="$output_path/concatenated_genes_${chromosome_name}/metafilesALL"

mkdir -p "$OUTPUT_FOLDER_15"
mkdir -p "$OUTPUT_FOLDER_20"
mkdir -p "$OUTPUT_FOLDER_ALL"

FINAL_LIST_15="${output_path}/${chromosome_name}_final15.lst"
FINAL_LIST_20="${output_path}/${chromosome_name}_final20.lst"
FINAL_LIST_ALL="${output_path}/${chromosome_name}_finalALL.lst"


# Read folder list into an array
mapfile -t FOLDERS < "$folder_list"

# Declare an associative array to track duplicate genes
declare -A gene_files15
declare -A gene_files20
declare -A gene_filesALL
> "$FINAL_LIST_15"
> "$FINAL_LIST_20"
> "$FINAL_LIST_ALL"
# Loop through each folder and find duplicate genes
for folder in "${FOLDERS[@]}"; do
echo "$folder"/metafiles/metafiles15 >> "$FINAL_LIST_15"
    for file in "$folder"/metafiles/metafiles15/*.meta; do
        [ -e "$file" ] || continue  # Skip if no .meta files found
        gene_name=$(basename "$file")
        gene_files15["$gene_name"]+="$file "
    done

echo "$folder"/metafiles/metafiles20 >> "$FINAL_LIST_20"
    for file in "$folder"/metafiles/metafiles20/*.meta; do
        [ -e "$file" ] || continue  # Skip if no .meta files found
        gene_name=$(basename "$file")
        gene_files20["$gene_name"]+="$file "
    done

echo "$folder"/metafiles/metafilesALL >> "$FINAL_LIST_ALL"
    for file in "$folder"/metafiles/metafilesALL/*.meta; do
        [ -e "$file" ] || continue  # Skip if no .meta files found
        gene_name=$(basename "$file")
        gene_filesALL["$gene_name"]+="$file "
    done
done

# Loop through the associative array and concatenate files for duplicate genes
for gene in "${!gene_files15[@]}"; do
    files=(${gene_files15[$gene]})
    if [ ${#files[@]} -gt 1 ]; then
        output_file="$OUTPUT_FOLDER_15/${gene}"
        
        # Extract the header from the first file
        head -n 1 "${files[0]}" > "$output_file"

        # Concatenate the rest of the files excluding their headers
        for file in "${files[@]}"; do
            tail -n +2 "$file" >> "$output_file"
        done
        
        ## Remove original files
        for file in "${files[@]}"; do
           rm -r "$file"
           echo "15"
        done
    fi
done

# Loop through the associative array and concatenate files for duplicate genes
for gene in "${!gene_files20[@]}"; do
    files=(${gene_files20[$gene]})
    if [ ${#files[@]} -gt 1 ]; then
        output_file="$OUTPUT_FOLDER_20/${gene}"
        
        # Extract the header from the first file
        head -n 1 "${files[0]}" > "$output_file"

        # Concatenate the rest of the files excluding their headers
        for file in "${files[@]}"; do
            tail -n +2 "$file" >> "$output_file"
        done
        
        ## Remove original files
        for file in "${files[@]}"; do
          rm -r "$file"
          echo "20"
        done
    fi
done

# Loop through the associative array and concatenate files for duplicate genes
for gene in "${!gene_filesALL[@]}"; do
    files=(${gene_filesALL[$gene]})
    if [ ${#files[@]} -gt 1 ]; then
        output_file="$OUTPUT_FOLDER_ALL/${gene}"
        
        # Extract the header from the first file
        head -n 1 "${files[0]}" > "$output_file"

        # Concatenate the rest of the files excluding their headers
        for file in "${files[@]}"; do
            tail -n +2 "$file" >> "$output_file"
        done
        
        ## Remove original files
        for file in "${files[@]}"; do
           rm -r "$file"
           echo "ALL"
        done
    fi
done

if [ -d "$OUTPUT_FOLDER_15" ]; then
  echo "${output_path}/concatenated_genes_${chromosome_name}/metafiles15" >> "$FINAL_LIST_15"
fi

if [ -d "$OUTPUT_FOLDER_20" ]; then
  echo "${output_path}/concatenated_genes_${chromosome_name}/metafiles20" >> "$FINAL_LIST_20"
fi

if [ -d "$OUTPUT_FOLDER_ALL" ]; then
  echo "${output_path}/concatenated_genes_${chromosome_name}/metafilesALL" >> "$FINAL_LIST_ALL"
fi
