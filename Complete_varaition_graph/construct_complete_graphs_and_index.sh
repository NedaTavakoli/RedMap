
# On hive
# ssh ntavakoli6@login-hive-slurm.pace.gatech.edu
# salloc -A hive-saluru8 -phive-himem -N1 --ntasks-per-node=1 --mem-per-cpu=1500G -t120:00:00
# #  salloc -A hive-saluru8 -phive-himem -N1 --ntasks-per-node=1 --mem-per-cpu=3000G -t120:00:00


# project_dir=/storage/hive/project/cse-aluru/ntavakoli6/hged
# scratch=/storage/hive/scratch/6/ntavakoli6  

# cd /storage/hive/project/cse-aluru/ntavakoli6/hged
# module load gurobi
# module load anaconda3
# module load boost
# module load cmake


# samtools=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/samtools-1.12/samtools
# bcftools=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/bcftools-1.9/bcftools
# bgzip=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/htslib-1.12/bgzip
# tabix=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/htslib-1.12/tabix
# GraphAligner=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/GraphAligner/bin/GraphAligner
# k8=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/k8/k8-linux


###################################################################################################
CHRID=chr22
VCF=data/chr22.vcf.gz
FA=Homo_sapiens.GRCh37.dna.chromosome.22.fa


# construct the graph
/usr/bin/time -v  ./vg construct -r $FA -v $VCF -a > ${CHRID}.vg   # The alt paths were created by the -a option 

# Indexing all variants graph, xg and gbwt indexes
/usr/bin/time -v ./vg index -x ${CHRID}.xg -G ${CHRID}.gbwt -v $VCF  ${CHRID}.vg

#  Prunning all variants graph pruning, Pruning with Haplotypes
/usr/bin/time -v ./vg prune -u -g ${CHRID}.gbwt -m node_mapping ${CHRID}.vg > ${CHRID}.pruned.vg
/usr/bin/time -v ./vg index -g ${CHRID}.pruned.gcsa -f node_mapping ${CHRID}.pruned.vg   # here gcsa index only can obtain if the graph is pruned

# Unfolding the paths with vg prune -u creates duplicated nodes, and if you
#  don't pass the "node mapping" to vg index -g, you get an index that maps to non-existent nodes.


