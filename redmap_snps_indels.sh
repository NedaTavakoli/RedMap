# dobule check the number of threads
# collect time and memory
# research on vg prune


# On hive
ssh ntavakoli6@login-hive-slurm.pace.gatech.edu
salloc -A hive-saluru8 -phive-himem -N1 --ntasks-per-node=1 --mem-per-cpu=1500G -t96:00:00

project_dir=/storage/hive/project/cse-aluru/ntavakoli6/hged
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

#     Repeat only once
$bcftools view -v 'snps,indels' -Oz chr22.vcf.gz > chr22_snps_indels.vcf.gz
$bgzip -c data/chr22_snps_indels.vcf > data/chr22_snps_indels.vcf.gz
$tabix -f -p vcf data/chr22_snps_indels.vcf.gz

# Let's create the complete graph that has onluy snps and indels
#       Vg graph for only snps and indels variants of chr22
/usr/bin/time -v ./vg construct -t 24 -r data/hs37d5.fa -v data/chr22_snps_indels.vcf.gz -R 22 -C > chr22_snps_indels.vg
/usr/bin/time -v ./vg index -x chr22_snps_indels.xg chr22_snps_indels.vg

#Simulating reads
/usr/bin/time -v  ./vg sim -x chr22_snps_indels.xg -n $COUNT -l $LEN -e $EROR -a  > simulated_reads_len_${LEN}_errRate_${EROR}_snps_indels.gam
./vg view -X  simulated_reads_len_${LEN}_errRate_${EROR}_snps_indels.gam >  simulated_reads_len_${LEN}_errRate_${EROR}_snps_indels.fastq
/usr/bin/time -v  ./vg surject -s  simulated_reads_len_${LEN}_errRate_${EROR}_snps_indels.gam -x chr22_snps_indels.xg > s simulated_reads_len_${LEN}_errRate_${EROR}_snps_indels.sam #truth on linear genome
# warning:[vg::get_sequence_dictionary] No reference-sense paths available in the graph; falling back to generic paths.


 
 # 1
# # compelte graph (only snp and indels) and mapping [without pruning]
/usr/bin/time -v ./vg index -x chr22_snps_indels.xg -g chr22_snps_indels.vg --gcsa-out chr22_snps_indels.gcsa
/usr/bin/time -v ./vg map -t 24 -x chr22_snps_indels.xg -g chr22_snps_indels.gcsa -f  simulated_reads_len_${LEN}_errRate_${EROR}_snps_indels.fastq > mapped_reads.chr22_snps_indels.gam
./vg view -aj mapped_reads.chr22_snps_indels.gam | jq '.score' > scores_mapped_reads.chr22_snps_indels.txt
awk '{ sum += $1 } END { print "Overall Alignment Score:", sum; accuracy = sum /(1000*NR); print "Alignment Accuracy:", accuracy*100 "%" }' scores_mapped_reads.chr22_snps_indels.txt
# Overall Alignment Score: 9597622
# Alignment Accuracy: 95.9762%

# Alignment Score = Match Score - Mismatch Penalty - Gap Penalty
# Alignment Accuracy = (Sum of score) / 1000*Total aligned


# ./vg stats -a mapped_reads.chr22_snps_indels.gam 
# # Total alignments: 10000
# # Total primary: 10000
# # Total secondary: 0
# # Total aligned: 10000
# # Total perfect: 1
# # Total gapless (softclips allowed): 9990
# # Total paired: 0
# # Total properly paired: 0
# # Insertions: 81 bp in 14 read events
# # Deletions: 1 bp in 1 read events
# # Substitutions: 99435 bp in 98498 read events
# # Softclips: 453 bp in 37 read events

# Alignment Accuracy = (Total perfect + Total gapless) / Total aligned

# From the given information:

# Total perfect: 1
# Total gapless (softclips allowed): 9990
# Total aligned: 10000

# Using these values in the formula:

# Alignment Accuracy = (1 + 9990) / 10000 = 0.9991 or 99.91%


# 1.2
# pruned graph (ins and del) and mapping and mapping
/usr/bin/time -v ./vg prune -r -p -t 24  chr22_snps_indels.vg > chr22_snps_indels.pruned.vg
# Original graph chr22_snps_indels.vg: 1603268 nodes, 1603267 edges
# Built a temporary XG index
# Removed all paths
# Pruned complex regions: 1603268 nodes, 1603267 edges 
# Removed small subgraphs: 1603268 nodes, 1603267 edges
# Restored graph: 1603268 nodes
# Serialized the graph: 1603268 nodes, 1603267 edges
/usr/bin/time -v ./vg index -t 24 -g chr22_snps_indels_pruned.gcsa chr22_snps_indels.pruned.vg
/usr/bin/time -v ./vg index -t 24 -x chr22_snps_indels_pruned.xg chr22_snps_indels.pruned.vg

/usr/bin/time -v ./vg map -t 24 -x chr22_snps_indels_pruned.xg -g chr22_snps_indels_pruned.gcsa -f   simulated_reads_len_${LEN}_errRate_${EROR}_snps_indels.fastq > mapped_reads.chr22_snps_indels_pruned.gam
/usr/bin/time -v ./vg view -aj mapped_reads.chr22_snps_indels_pruned.gam | jq '.score' > scores_mapped_reads.chr22_snps_indels_pruned.txt
awk '{ sum += $1 } END { print "Overall Alignment Score:", sum; accuracy = sum /(1000*NR); print "Alignment Accuracy:", accuracy*100 "%" }' scores_mapped_reads.chr22_snps_indels_pruned.txt
# Overall Alignment Score: 9602274
# Alignment Accuracy: 96.0227%

./vg stats -a mapped_reads.chr22_snps_indels_pruned.gam 
# Total alignments: 10000
# Total primary: 10000
# Total secondary: 0
# Total aligned: 10000
# Total perfect: 1
# Total gapless (softclips allowed): 9989
# Total paired: 0
# Total properly paired: 0
# Insertions: 46 bp in 13 read events
# Deletions: 1 bp in 1 read events
# Substitutions: 99429 bp in 98494 read events
# Softclips: 268 bp in 30 read events

# Alignment Accuracy = (Total perfect + Total gapless) / Total aligned

# From the given information:

# Total perfect: 1
# Total gapless (softclips allowed): 9989
# Total aligned: 10000

# Using these values in the formula:

# Alignment Accuracy = (1 + 9989) / 10000 = 0.999 or 99.9%




#####################################

# using graphaligner

git clone https://github.com/maickrau/GraphAligner.git
cd GraphAligner
git submodule update --init --recursive
conda env create -f CondaEnvironment_linux.yml or conda env create -f CondaEnvironment_osx.yml
source activate GraphAligner
make bin/GraphAligner


$GraphAligner -g chr22_snps_indels.vg -f  simulated_reads_len_${LEN}_errRate_${EROR}_snps_indels.fastq -a GraphAligner.${CHRID}.aligned.gam -x vg
# GraphAligner Branch master commit 46e03e290280e174d187358e9f13defa804163a9 2023-05-25 14:34:28 +0300
# GraphAligner Branch master commit 46e03e290280e174d187358e9f13defa804163a9 2023-05-25 14:34:28 +0300
# Load graph from chr22_snps_indels.vg
# Build minimizer seeder from the graph
# Minimizer seeds, length 15, window size 20, density 10
# Seed cluster size 1
# Extend up to 5 seed clusters
# Alignment bandwidth 10
# Clip alignment ends with identity < 66%
# X-drop DP score cutoff 14705
# Backtrace from 10 highest scoring local maxima per cluster
# write alignments to GraphAligner..aligned.gam
# Align
# Alignment finished
# Input reads: 10000 (10000000bp)
# Seeds found: 4747696
# Seeds extended: 49998
# Reads with a seed: 10000 (10000000bp)
# Reads with an alignment: 10000 (9999628bp)
# Alignments: 11120 (11118734bp) (4652 additional alignments discarded)
# End-to-end alignments: 10812 (10812000bp)

./vg stats -a GraphAligner.${CHRID}.aligned.gam
# Total alignments: 11120
# Total primary: 11120
# Total secondary: 0
# Total aligned: 11116
# Total perfect: 0
# Total gapless (softclips allowed): 10511
# Total paired: 0
# Total properly paired: 0
# Insertions: 1325 bp in 710 read events
# Deletions: 1228 bp in 713 read events
# Substitutions: 121392 bp in 119687 read events
# Softclips: 0 bp in 0 read events

# Alignment Accuracy = (Total gapless alignments / Total aligned reads) * 100

# Given the information you provided:

# Total gapless alignments: 10,511
# Total aligned reads: 11,116

# Using these values in the formula:

# Alignment Accuracy = (10,511 / 11,116) * 100 ≈ 94.49%



##### Our reduced graph
# get variant retained vcf file
# do only once
cat chr22_snps_indels.vcf | grep '^#' > header_chr22.vcf
cat chr22_snps_indels.vcf | grep  -vE '^#' > non_header_chr22.txt

# do for each config
cp header_chr22.vcf retained_variants_ed_ILP_1000_1_chr22.vcf
awk 'NR == FNR {a[$0]; next } $2 in a {print $0} ' retained_variants_ed_ILP_100_1.txt  non_header_chr22.txt >> retained_variants_ed_ILP_1000_1_chr22.vcf

# Generate vcf.gz file and its index file vcf.gz.tbi
cp retained_variants_ed_ILP_1000_1_chr22.vcf retained_variants_ed_ILP_100_1_chr22_copy.vcf
$bgzip -c retained_variants_ed_ILP_1000_1_chr22.vcf > retained_variants_ed_ILP_1000_1_chr22.vcf.gz
$tabix -p vcf retained_variants_ed_ILP_1000_1_chr22.vcf.gz

# on reduced graph
/usr/bin/time -v ./vg construct -r data/hs37d5.fa -v retained_variants_ed_ILP_1000_1_chr22.vcf.gz -a -f -m 32 > reduced_1000_1.graph.chr22.vg


# # compelte graph (only snp and indels) and mapping [without pruning]
/usr/bin/time -v ./vg index -x reduced_1000_1.graph.chr22.xg -g chr22_snps_indels.vg --gcsa-out reduced_1000_1.graph.chr22.gcsa 
/usr/bin/time -v ./vg map -t 24 -x reduced_1000_1.graph.chr22.xg -g reduced_1000_1.graph.chr22.gcsa -f s simulated_reads_len_${LEN}_errRate_${EROR}_snps_indels.fastq > mapped_reads_reduced_1000_1.gam
./vg view -aj mapped_reads_reduced_1000_1.gam | jq '.score' > scores_mapped_reads_reduced_1000_1.txt
awk '{ sum += $1 } END { print "Overall Alignment Score:", sum; accuracy = sum /(1000*NR); print "Alignment Accuracy:", accuracy*100 "%" }' scores_mapped_reads_reduced_1000_1.txt


# ./vg map -s GGTGGGCCTAAGTGCCTCCTCTCGGGACTGGTATGGGGACGGTCATGCAA -x reduced_100_1.graph.chr22.xg -g reduced_100_1.graph.chr22.gcsa > reduced_100_1.read.chr22.gam
# ./vg surject -x reduced_100_1.graph.chr22.xg -b reduced_100_1.read.chr22.gam > reduced_100_1.read.chr22.gam.bam
# $samtools flagstat reduced_100_1.read.chr22.gam.bam

