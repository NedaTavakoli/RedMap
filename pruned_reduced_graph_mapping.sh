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
DELTA=1
CHRID=chr22
EROR=0.01
LEN_ADDED_0.1=110

samtools=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/samtools-1.12/samtools
bcftools=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/bcftools-1.9/bcftools
bgzip=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/htslib-1.12/bgzip
tabix=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/htslib-1.12/tabix
GraphAligner=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/GraphAligner/bin/GraphAligner
k8=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/k8/k8-linux


# **************************************************** On Complete graph ***********************
# ==============================================================================================


# *****************************************  DO ONLY ONCE  ************************************
# Constructing all variant graph construction
/usr/bin/time -v ./vg construct -r data/hs37d5.fa -v data/${CHRID}.vcf.gz -R 22 -C > ${CHRID}.vg

#  Prunning all variants graph pruning
/usr/bin/time -v ./vg prune -r -p -t 24  ${CHRID}.vg > ${CHRID}.pruned.vg

# Indexing all variants graph 
/usr/bin/time -v ./vg index -x ${CHRID}_pruned.xg ${CHRID}.pruned.vg
/usr/bin/time -v  $project_dir/./vg index -g  ${CHRID}_pruned.gcsa ${CHRID}.pruned.vg
# *****************************************  DO ONLY ONCE  ************************************




# *****************************************  Run per config  ************************************
#Simulating reads from unpruned graph
/usr/bin/time -v ./vg sim -x ${CHRID}.xg -n $COUNT -l $LEN -e $EROR -a  > simulated_reads_len_${LEN}_errRate_$EROR.gam
./vg view -X  simulated_reads_len_${LEN}_errRate_$EROR.gam >  simulated_reads_len_${LEN}_errRate_$EROR.fastq
/usr/bin/time -v ./vg surject -s  simulated_reads_len_${LEN}_errRate_$EROR.gam -x ${CHRID}.xg > simulated_reads_len_${LEN}_errRate_$EROR.sam #truth on linear genome

# # Mapping the simulated reads to the original graph uisng vg
# /usr/bin/time -v ./vg map -t 24 -x ${CHRID}_pruned.xg -g ${CHRID}_pruned.gcsa -f  simulated_reads_len_${LEN}_errRate_$EROR.fastq > mapped_reads.${CHRID}_pruned.gam
# /usr/bin/time -v ./vg view -aj mapped_reads.${CHRID}_pruned.gam | jq '.score' > scores_mapped_reads.${CHRID}_pruned.txt
# awk '{ sum += $1 } END { print "Overall Alignment Score:", sum; accuracy = sum /(${LEN_ADDED_0.1}*NR); print "Alignment Accuracy:", accuracy*100 "%" }' scores_mapped_reads.${CHRID}_pruned.txt

# # Map the simulated reads to the original graph uisng GraphAligner
# /usr/bin/time -v  $GraphAligner  -t 24 -g ${CHRID}.vg -f  simulated_reads_len_${LEN}_errRate_$EROR.fastq -a GraphAligner.${CHRID}.aligned.gam -x vg
#  # Acccuracy  Reads with an alignment: 990387bp /(LEN+ .1*LEN)*100
# # ************************************************************************************************



##### On reduced graph


# **************************************************** On Reduced graph ***********************
# ==============================================================================================


# *****************************************  DO ONLY ONCE  ************************************
# Get list of all variant positions using bctfools 
$bcftools query -f '%POS\n' data/${CHRID}.vcf.gz > all_variant_positons_${CHRID}.txt

# Get the list of only snps and indels using bcftools
$bcftools query -i '(TYPE="snp" || TYPE="indel") && GT="alt"' -f '%POS\n'  data/${CHRID}.vcf.gz > snps_indels_variant_positions_${CHRID}.txt

# Get the positions that are not snps and indels (all other positions):
comm -13 <(sort snps_indels_variant_positions_${CHRID}.txt) <(sort all_variant_positions_${CHRID}.txt) > not_snps_indels_positions_${CHRID}.txt

# *****************************************************************************



# *****************************************  Run per config  ************************************
# For config LEN and DELTA
# Get the retained variant positions from the ILP_sol_${LEN}_${DELTA}
# Note that ILP_sol_${LEN}_{DELTA} is the output of running our algoritm to reduce the number of variants


mv ILP_sol_${LEN}_{DELTA} ILP_sol_${LEN}_${DELTA}.txt
awk 'NR==FNR {values[NR]=$1; next} values[FNR] == 0 {print $1}' ILP_sol_${LEN}_${DELTA}.txt snps_indels_variant_positions_${CHRID}.txt  > retained_variants_${LEN}_${DELTA}_${CHRID}.txt

 # get reduced vcf file using retained variants

 # 1 
 #combine no_snps_indels_positions_${CHRID}.txt and retained_variants_${LEN}_${DELTA}_${CHRID}.txt
cat retained_variants_${LEN}_${DELTA}_${CHRID}.txt not_snps_indels_positions_${CHRID}.txt | sort -n > retained_snp_indel_variants_${LEN}_${DELTA}_${CHRID}_plus_other_variants.txt

# 2 
#I want to filter a  vcf file such that it only contains the variants in the file retained_snp_indel_variants_${LEN}_${DELTA}_${CHRID}_plus_other_variants.txt using bcftools 
positions_of_interest='retained_snp_indel_variants_${LEN}_${DELTA}_${CHRID}_plus_other_variants.txt'
$bcftools view -h data/${CHRID}.vcf > retained_snp_indel_variants_${LEN}_${DELTA}_${CHRID}_plus_other_variants.vcf
grep -Fwf $positions_of_interest data/${CHRID}.vcf  >> retained_snp_indel_variants_${LEN}_${DELTA}_${CHRID}_plus_other_variants.vcf

# Generate vcf.gz file and its index file vcf.gz.tbi
cp retained_snp_indel_variants_${LEN}_${DELTA}_${CHRID}_plus_other_variants.vcf retained_snp_indel_variants_${LEN}_${DELTA}_${CHRID}_plus_other_variants_copy.vcf
$bgzip -c retained_snp_indel_variants_${LEN}_${DELTA}_${CHRID}_plus_other_variants.vcf > retained_snp_indel_variants_${LEN}_${DELTA}_${CHRID}_plus_other_variants.vcf.gz
$tabix -p vcf retained_snp_indel_variants_${LEN}_${DELTA}_${CHRID}_plus_other_variants.vcf.gz
# ****************************************************************************



my_new_vcf='retained_snp_indel_variants_${LEN}_${DELTA}_${CHRID}_plus_other_variants.vcf.gz'
my_new_vg='retained_snp_indel_variants_${LEN}_${DELTA}_${CHRID}_plus_other_variants'
/usr/bin/time -v ./vg construct -r data/hs37d5.fa -v ${my_new_vcf} -R 22 -C > ${my_new_vg}.vg
/usr/bin/time -v ./vg index -x ${my_new_vg}.xg ${my_new_vg}.vg

#=======================================================================================================================
#***************************  Pruned all variants   *************************** 
# Prune reduced graph
/usr/bin/time -v ./vg prune -r -p -t 24  ${my_new_vg}.vg > ${my_new_vg}.pruned.vg

# Indexing reduced graph
/usr/bin/time -v ./vg index -x  ${my_new_vg}_pruned.xg ${my_new_vg}.pruned.vg
/usr/bin/time -v  $project_dir/./vg index -g   ${my_new_vg}_pruned.gcsa  ${my_new_vg}.pruned.vg

# Mapping the simulated reads on the reduced graph using vg toolkit
/usr/bin/time -v ./vg map -t 24 -x ${my_new_vg}_pruned.xg  -g ${my_new_vg}_pruned.gcsa -f simulated_reads_len_${LEN}_errRate_$EROR.fastq > mapped_reads.${my_new_vg}_pruned.gam
/usr/bin/time -v ./vg view -aj mapped_reads.${my_new_vg}_pruned.gam | jq '.score' > scores_mapped_reads.${my_new_vg}_pruned.txt
awk '{ sum += $1 } END { print "Overall Alignment Score:", sum; accuracy = sum /(110*NR); print "Alignment Accuracy:", accuracy*100 "%" }' scores_mapped_reads.${my_new_vg}_pruned.txt


#  Mapping the simulated reads on the reduced graph using graphaligner
/usr/bin/time -v  $GraphAligner  -t 24 -g ${CHRID}.vg -f  simulated_reads_len_${LEN}_errRate_$EROR.fastq -a GraphAligner.${CHRID}.aligned.gam -x vg
# Acccuracy  Reads with an alignment: 990387bp /(LEN+ .1*LEN)*100





