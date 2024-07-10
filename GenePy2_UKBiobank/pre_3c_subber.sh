cd /home/gc1a20/hgig_me/ukbb/wes_200k/run_202206
ls -1 original/*_v1.vcf.gz |while read i
    do
        v2=$(echo $i | sed 's/v1.vcf\.gz/p2/g')
        v5=$(echo $i | sed 's/v1.vcf\.gz/p1/g')
        v4=$(echo $i | sed 's/v1\.vcf\.gz/v1_p1\.vep\.vcf/g')
        v3=$(echo $v4|sed 's/p1/p12/g' | sed 's/\.vep//g')
        v6=$(echo $v3 |sed 's/p12/p12_f1/g')
        v7=$(echo $v3| sed 's/p12/p12_filterNew_onTarget/g'| sed 's/vcf/vcf\.gz/g')
        sbatch script/pre_3c.sh $i $v2 $v3 $v4 $v5 $v6 $v7
    done
