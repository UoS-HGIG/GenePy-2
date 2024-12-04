process ConcatenateDuplicateGenes {
    //publishDir "${params.output}", mode: "copy", overwrite: true

    input:
    tuple val(cadd),val(chromosome_name),val(folder_paths)

    output:
    path("*.lstx"), emit: path_

    shell:
    """
    OUTPUT_FOLDER="metafiles${cadd}"
    mkdir -p "\$OUTPUT_FOLDER"
    FINAL_LIST="${chromosome_name}_${cadd}_final.lst"
    > "\$FINAL_LIST"

    declare -a FOLDERS
    for folder in ${folder_paths}; do
        folder=\$(echo "\$folder" | tr -d '[],') 
        FOLDERS+=("\$folder")
    done
    declare -A gene_files
    for folder in "\${FOLDERS[@]}"; do
        echo "\$folder" >> "\$FINAL_LIST"
        for file in "\${folder}"/*.meta; do
            [ -e "\$file" ] || continue  
            gene_name=\$(basename "\$file")
            gene_files["\$gene_name"]+="\$file "
            echo "\$gene_name"
        done
    done
    echo "half done!"
    for item in "\${gene_files[@]}"; do
        echo "\$item"
    done

    for gene in "\${!gene_files[@]}"; do
        
        echo "this is  \$gene"
        files=(\${gene_files[\$gene]})
        echo "\$files"
        if [ \${#files[@]} -gt 1 ]; then
            echo "this is duplicated  \$gene"
            output_file="\$OUTPUT_FOLDER/\${gene}"

            head -n 1 "\${files[0]}" > "\$output_file"

            for file in "\${files[@]}"; do
                
                tail -n +2 "\$file" >> "\$output_file"
            done

            for file in "\${files[@]}"; do
                rm -f "\$file"
               # echo "scgfx"
            done
        fi
    done

    if [ -d "\$OUTPUT_FOLDER" ]; then
        realpath "\$OUTPUT_FOLDER" >> "\$FINAL_LIST"
    fi
    
    while IFS= read -r line; do
        base_name=\$(basename "\$line")
        echo "\$line" > "\$base_name".lstx
    done < "\$FINAL_LIST"
    
    """
}
