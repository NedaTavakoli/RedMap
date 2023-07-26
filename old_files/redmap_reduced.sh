# On hive
ssh ntavakoli6@login-hive-slurm.pace.gatech.edu
salloc -A hive-saluru8 -phive-himem -N1 --ntasks-per-node=1 --mem-per-cpu=1500G -t96:00:00

project_dir=/storage/hive/project/cse-aluru/ntavakoli6/hged
scratch=/storage/hive/scratch/6/ntavakoli6  

cd /storage/hive/project/cse-aluru/ntavakoli6/hged
module load gurobi
module load anaconda3
module load boost
module load cmake

COUNT=10000
LEN=1000
CHRID=chr22
EROR=0.01

samtools=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/samtools-1.12/samtools
bcftools=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/bcftools-1.9/bcftools
bgzip=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/htslib-1.12/bgzip
tabix=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/htslib-1.12/tabix
GraphAligner=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/GraphAligner/bin/GraphAligner
k8=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/k8/k8-linux


# *****************************************************************************
# All variant construction
/usr/bin/time -v ./vg construct -r data/hs37d5.fa -v data/chr22.vcf.gz -R 22 -C > chr22.vg
/usr/bin/time -v ./vg index -x chr22.xg chr22.vg

# *****************************************************************************
#Simulating reads
/usr/bin/time -v ./vg sim -x chr22.xg -n $COUNT -l $LEN -e $EROR -a  > simulated_reads_len_${LEN}_errRate_$EROR.gam
./vg view -X  simulated_reads_len_${LEN}_errRate_$EROR.gam >  simulated_reads_len_${LEN}_errRate_$EROR.fastq
/usr/bin/time -v ./vg surject -s  simulated_reads_len_${LEN}_errRate_$EROR.gam -x chr22.xg > simulated_reads_len_${LEN}_errRate_$EROR.sam #truth on linear genome
# warning:[vg::get_sequence_dictionary] No reference-sense paths available in the graph; falling back to generic paths.

# *****************************************************************************

##### Our reduced graph
# get variant retained vcf file
# do only once
cat data/chr22.vcf | grep '^#' > header_chr22.vcf
cat data/chr22.vcf | grep  -vE '^#' > non_header_chr22.txt

# do for each config
cp header_chr22.vcf retained_variants_ed_ILP_1000_1_chr22.vcf
awk 'NR == FNR {a[$0]; next } $2 in a {print $0} ' retained_variants_ed_ILP_1000_1.txt  non_header_chr22.txt >> retained_variants_ed_ILP_1000_1_chr22.vcf

# Generate vcf.gz file and its index file vcf.gz.tbi
cp retained_variants_ed_ILP_1000_1_chr22.vcf retained_variants_ed_ILP_100_1_chr22_copy.vcf
$bgzip -c retained_variants_ed_ILP_1000_1_chr22.vcf > retained_variants_ed_ILP_1000_1_chr22.vcf.gz
$tabix -p vcf retained_variants_ed_ILP_1000_1_chr22.vcf.gz
# *****************************************************************************

#======================================================================================================================
#***************************  reduced Variants   *************************** 
# on reduced graph -R 22 -C
/usr/bin/time -v ./vg construct -r data/hs37d5.fa -v retained_variants_ed_ILP_1000_1_chr22.vcf.gz -a -f -m 32 > reduced_1000_1.graph.chr22.vg


# # compelte graph (only snp and indels) and mapping [without pruning]
/usr/bin/time -v ./vg index -x reduced_1000_1.graph.chr22.xg -g chr22.vg --gcsa-out $scratch/reduced_1000_1.graph.chr22.gcsa 
/usr/bin/time -v ./vg map -t 24 -x reduced_1000_1.graph.chr22.xg -g $scratch/reduced_1000_1.graph.chr22.gcsa -f simulated_reads_len_${LEN}_errRate_$EROR.fastq > mapped_reads_reduced_1000_1.gam
./vg view -aj mapped_reads_reduced_1000_1.gam | jq '.score' > scores_mapped_reads_reduced_1000_1.txt
awk '{ sum += $1 } END { print "Overall Alignment Score:", sum; accuracy = sum /(1000*NR); print "Alignment Accuracy:", accuracy*100 "%" }' scores_mapped_reads_reduced_1000_1.txt


# ./vg map -s GGTGGGCCTAAGTGCCTCCTCTCGGGACTGGTATGGGGACGGTCATGCAA -x reduced_100_1.graph.chr22.xg -g reduced_100_1.graph.chr22.gcsa > reduced_100_1.read.chr22.gam
# ./vg surject -x reduced_100_1.graph.chr22.xg -b reduced_100_1.read.chr22.gam > reduced_100_1.read.chr22.gam.bam
# $samtools flagstat reduced_100_1.read.chr22.gam.bam

