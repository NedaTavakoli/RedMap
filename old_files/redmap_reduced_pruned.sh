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
LEN=100
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

# Map the simulated reads to the original graph

LEN=1000
/usr/bin/time -v ./vg map -t 24 -x chr22_pruned.xg -g chr22_pruned.gcsa -f  simulated_reads_len_${LEN}_errRate_$EROR.fastq > mapped_reads.chr22_${LEN}_pruned.gam
/usr/bin/time -v ./vg view -aj mapped_reads.chr22_${LEN}_pruned.gam | jq '.score' > scores_mapped_reads.chr22_${LEN}_pruned.txt
awk '{ sum += $1 } END { print "Overall Alignment Score:", sum; accuracy = sum /(110*NR); print "Alignment Accuracy:", accuracy*100 "%" }' scores_mapped_reads.chr22_${LEN}_pruned.txt


# *****************************************************************************

##### Our reduced graph
# get variant retained vcf file

# ***************************************************************************** DO ONLY ONCE
# Get list of all variant positions using bctfools 
$bcftools query -f '%POS\n' data/chr22.vcf.gz > all_variant_positons_chr22.txt

# Get the list of only snps and indels using bcftools
$bcftools query -i '(TYPE="snp" || TYPE="indel") && GT="alt"' -f '%POS\n'  data/chr22.vcf.gz > snps_indels_variant_positions_chr22.txt

# Get the positions that are not snps and indels (all other positions):
comm -13 <(sort snps_indels_variant_positions_chr22.txt) <(sort all_variant_positions_chr22.txt) > not_snps_indels_positions_chr22.txt

# *****************************************************************************

# For config 100 and 1 
# Get the retained variant positions from the ILP_sol_100_1
# check if it zero shows retained variant positions *******
awk 'NR==FNR {values[NR]=$1; next} values[FNR] == 0 {print $1}' ILP_sol_100_1.txt snps_indels_variant_positions_chr22.txt  > retained_variants_100_1_chr22.txt

 # get reduced vcf file using retained variants

 #substeps:
 # 1 
 #combine no_snps_indels_positions_chr22.txt and retained_variants_100_1_chr22.txt
cat retained_variants_100_1_chr22.txt not_snps_indels_positions_chr22.txt | sort -n > retained_snp_indel_variants_100_1_chr22_plus_other_variants.txt


# 2 
#I want to filter a  vcf file such that it only contains the variants in the file retained_snp_indel_variants_100_1_chr22_plus_other_variants.txt using bcftools 

positions_of_interest='retained_snp_indel_variants_100_1_chr22_plus_other_variants.txt'
$bcftools view -h data/chr22.vcf > retained_snp_indel_variants_100_1_chr22_plus_other_variants.vcf
grep -Fwf $positions_of_interest data/chr22.vcf  >> retained_snp_indel_variants_100_1_chr22_plus_other_variants.vcf

# Generate vcf.gz file and its index file vcf.gz.tbi
cp retained_snp_indel_variants_100_1_chr22_plus_other_variants.vcf retained_snp_indel_variants_100_1_chr22_plus_other_variants_copy.vcf
$bgzip -c retained_snp_indel_variants_100_1_chr22_plus_other_variants.vcf > retained_snp_indel_variants_100_1_chr22_plus_other_variants.vcf.gz
$tabix -p vcf retained_snp_indel_variants_100_1_chr22_plus_other_variants.vcf.gz



# *****************************************************************************

my_new_vcf='retained_snp_indel_variants_100_1_chr22_plus_other_variants.vcf.gz'
my_new_vg='retained_snp_indel_variants_100_1_chr22_plus_other_variants'
/usr/bin/time -v ./vg construct -r data/hs37d5.fa -v ${my_new_vcf} -R 22 -C > ${my_new_vg}.vg
/usr/bin/time -v ./vg index -x ${my_new_vg}.xg ${my_new_vg}.vg

#=======================================================================================================================
#***************************  Pruned all variants   *************************** 
/usr/bin/time -v ./vg prune -r -p -t 24  ${my_new_vg}.vg > ${my_new_vg}.pruned.vg
# Original graph retained_snp_indel_variants_100_1_chr22_plus_other_variants.vg: 1943888 nodes, 2101643 edges
# Built a temporary XG index
# Removed all paths
# Pruned complex regions: 1943888 nodes, 2089709 edges
# Removed small subgraphs: 1933498 nodes, 2082784 edges
# Restored graph: 1940403 nodes
# Serialized the graph: 1940403 nodes, 2090876 edges


/usr/bin/time -v ./vg index -x  ${my_new_vg}_pruned.xg ${my_new_vg}.pruned.vg
/usr/bin/time -v  $project_dir/./vg index -g   ${my_new_vg}_pruned.gcsa  ${my_new_vg}.pruned.vg


/usr/bin/time -v ./vg map -t 24 -x ${my_new_vg}_pruned.xg  -g ${my_new_vg}_pruned.gcsa -f simulated_reads_len_${LEN}_errRate_$EROR.fastq > mapped_reads.${my_new_vg}_pruned.gam
/usr/bin/time -v ./vg view -aj mapped_reads.${my_new_vg}_pruned.gam | jq '.score' > scores_mapped_reads.${my_new_vg}_pruned.txt
awk '{ sum += $1 } END { print "Overall Alignment Score:", sum; accuracy = sum /(110*NR); print "Alignment Accuracy:", accuracy*100 "%" }' scores_mapped_reads.${my_new_vg}_pruned.txt





#####################################

# using graphaligner

/usr/bin/time -v  $GraphAligner  -t 24 -g chr22.vg -f  simulated_reads_len_${LEN}_errRate_$EROR.fastq -a GraphAligner.${CHRID}.aligned.gam -x vg

# GraphAligner Branch master commit 46e03e290280e174d187358e9f13defa804163a9 2023-05-25 14:34:28 +0300
# GraphAligner Branch master commit 46e03e290280e174d187358e9f13defa804163a9 2023-05-25 14:34:28 +0300
# Load graph from chr22.vg
# Build minimizer seeder from the graph
# Minimizer seeds, length 15, window size 20, density 10
# Seed cluster size 1
# Extend up to 5 seed clusters
# Alignment bandwidth 10
# Clip alignment ends with identity < 66%
# X-drop DP score cutoff 14705
# Backtrace from 10 highest scoring local maxima per cluster
# write alignments to GraphAligner.chr22.aligned.gam
# Align
# Alignment finished
# Input reads: 10000 (1000000bp)
# Seeds found: 337137
# Seeds extended: 35779
# Reads with a seed: 9975 (997500bp)
# Reads with an alignment: 9975 (990387bp)
# Alignments: 12386 (1213283bp) (14454 additional alignments discarded)


# (99750/110000 )*100
# We consider a read correctly mapped if \
# its longest alignment overlaps at least 10% with\
#  the genomic position from where it was simulated \
#  and evaluate the number of reads correctly aligned.


# No need
./vg stats -a GraphAligner.${CHRID}.aligned.gam
# Total alignments: 12386
# Total primary: 12386
# Total secondary: 0
# Total aligned: 8540
# Total perfect: 0
# Total gapless (softclips allowed): 7608
# Total paired: 0
# Total properly paired: 0
# Insertions: 1046 bp in 886 read events
# Deletions: 766 bp in 685 read events
# Substitutions: 24503 bp in 22504 read events
# Softclips: 0 bp in 0 read events




