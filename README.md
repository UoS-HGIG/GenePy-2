GenePy-2.0 is gene- and region-based pathogenecity score for each individual based the carrier status of genomic variants; major updates from GenePy-1.x is 1). handling of multi-allelic loci with number of alternative alleles up to 10; 2). region-based score is encoporated; 3). computational efficiency improved with GPU-based processing option available

Installation pre-requirement:
1). Ensemble VEP 
2). CADD >= 1.6
3). Python3
3). numpy==1.26.4
4). pandas==2.2.1
5). pyarrow==15.0.2
6). numba==0.59.1

Running GenePy is by running the python make_scores_mat.py, and options can be found by -h; the input meta file is the annotated variant file from vcf, an example header of the meta file is:

CHROM   POS     REF     ALT     CONSEQUENCE     GENE    AF1     AF2     AF3     AF4     AF5     AF6     AF7     AF8     AF9     AF10    SCORE1  SCORE2  SCORE3  SCORE4  SCORE5  SCORE6  SCORE7  SCORE8  SCORE9  SCORE10  Sample1-genotype  Sample2-genotype

*AF: allele frequency; AF1: AF of the 1st alternative allele; SCORE: the CADD raw score of the allele

The meta file can be generated using the pre_1.sh which annotate and quality-based filter the vcf file, followed by pre_local.sh which convert the vcf file into the meta file
