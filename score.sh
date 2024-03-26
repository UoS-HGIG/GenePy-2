#!/bin/bash
#SBATCH --mem=32g
#SBATCH --nodes=1
#SBATCH --job-name="soccer"
#SBATCH --ntasks-per-node=1
#SBATCH --time=47:59:00

##calculate the genepy score of genes per chunk and the chunk id is the 1st vaibale here; x-genes, genes crossing chunks were extracted to x-genes directory for the process

##2nd varaible is the applied CADD phred score filtration, e.g all, phred >= 15, or phred >=20

cd /home/gc1a20/hgig_me/wes_2023010
mkdir -p ${2}
cd ${2}

cp /home/gc1a20/bin/make_scores_mat.py .

##gene.lst is the list of genes in each chunk/bin
cat ../${1} | while read i; do
    cp ../header.meta ${i}.meta
    grep -w "$i" ../meta_${2}.txt |\
        awk -F"\t" '{OFS=FS}{for (i=7;i<=16;i++) if(length($i)<1 || $i==0) $i="3.98e-6"}1' >> ${i}.meta
    python make_scores_mat.py --gene ${i} --cadd $2
    rm ${i}.meta
done




