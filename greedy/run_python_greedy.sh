# On hive
ssh ntavakoli6@login-hive-slurm.pace.gatech.edu
salloc -A hive-saluru8 -phive-himem -N1 --ntasks-per-node=1 --mem-per-cpu=1500G -t120:00:00
salloc -A hive-saluru8 -phive-himem -N1 --ntasks-per-node=1 --mem-per-cpu=3000G -t120:00:00


project_dir=/storage/hive/project/cse-aluru/ntavakoli6/hged
scratch=/storage/hive/scratch/6/ntavakoli6  

cd /storage/hive/project/cse-aluru/ntavakoli6/hged
module load gurobi
module load anaconda3
module load boost
module load cmake


samtools=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/samtools-1.12/samtools
bcftools=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/bcftools-1.9/bcftools
bgzip=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/htslib-1.12/bgzip
tabix=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/htslib-1.12/tabix
GraphAligner=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/GraphAligner/bin/GraphAligner
k8=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/k8/k8-linux


###################################################################################################


# to get start_pos
start_pos=$($bcftools query -f '%POS\n' -i '(TYPE="snp" || TYPE="indels") && GT="alt"' data/chr22.vcf | head -1)
# $bcftools query -f '%POS\n' -i '(TYPE="snp" || TYPE="indels") && GT="alt"' data/chr22.vcf | head -1
# # 16050075

# to get end_pos, very very important and took time to get it
end_pos=$($bcftools query -f '%POS\n' -i '(TYPE="snp" || TYPE="indels") && GT="alt"' data/chr22.vcf | tail -1)
# $bcftools query -f '%POS\n' -i '(TYPE="snp" || TYPE="indels") && GT="alt"' data/chr22.vcf | tail -1
# # 51244237

# to get total number of variants (snps and indels)
total_variants=$($bcftools query -f '%POS\n' -i '(TYPE="snp" || TYPE="indels") && GT="alt"' data/chr22.vcf | wc -l)
# $bcftools query -f '%POS\n' -i '(TYPE="snp" || TYPE="indels") && GT="alt"' data/chr22.vcf | wc -l
# # 1098701

##################################################################


chr_id=22                 #* change this numbers according to your needs
alpha=50 # 50, 75, 100, 150, 500, 10000                #* change this numbers according to your needs
delta=5             
start_pos=16050075        #* change this numbers according to your needs; first variant position (here is for chr22)
end_pos=51244237         #* change this numbers according to your needs; last variant position (here is for chr22)
total_variants=1098701

############ Do this only once
conda activate myenv
python src/data_wrangler_cost_only.py data/hs37d5.fa \
data/chr${chr_id}.vcf.gz ${chr_id} ${start_pos} ${end_pos} ${total_variants} \
${alpha} graph_chr22_all_updated_cost_only.txt greedy_cost_chr${chr_id}_allVariants_${alpha}.txt   # cost is for each position, how many alternates so we have
##################################################################

#---------------------------------------------  
alpha=150   # 50, 75, 100, 150, 500, 10000            #* change this numbers according to your needs
delta=1
conda activate myenv
# pip3 install networkx

#  to run greedy
/usr/bin/time -v python src/greedy.py graph_chr22_all_updated_cost_only.txt  \
greedy_cost_chr22_allVariants_50.txt \
pos_substrings_chr22_all_${alpha}_updated.txt ${alpha} ${delta}  #> greedy_chr22_allVariants_${alpha}_${delta}.txt




