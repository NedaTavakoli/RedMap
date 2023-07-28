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


#=======================================================================================================================
#***************************  Pruned all variants   *************************** 
/usr/bin/time -v ./vg prune -r -p -t 24  chr22.vg > chr22.pruned.vg
# Original graph chr22.vg: 4328949 nodes, 5514085 edges
# Built a temporary XG index
# Removed all paths
# Pruned complex regions: 4328949 nodes, 5163857 edges
# Removed small subgraphs: 3974904 nodes, 4892512 edges
# Restored graph: 4212309 nodes
# Serialized the graph: 4212309 nodes, 5167904 edges

/usr/bin/time -v ./vg index -x chr22_pruned.xg chr22.pruned.vg

/usr/bin/time -v  $project_dir/./vg index -g  chr22_pruned.gcsa chr22.pruned.vg


/usr/bin/time -v ./vg map -t 24 -x chr22_pruned.xg -g chr22_pruned.gcsa -f  simulated_reads_len_${LEN}_errRate_$EROR.fastq > mapped_reads.chr22_pruned.gam
/usr/bin/time -v ./vg view -aj mapped_reads.chr22_pruned.gam | jq '.score' > scores_mapped_reads.chr22_pruned.txt
awk '{ sum += $1 } END { print "Overall Alignment Score:", sum; accuracy = sum /(1000*NR); print "Alignment Accuracy:", accuracy*100 "%" }' scores_mapped_reads.chr22_pruned.txt

# Overall Alignment Score: 9466581
# Alignment Accuracy: 94.6658%
./vg stats -a mapped_reads.chr22_pruned.gam 
Overall Alignment Score: 9466581
Alignment Accuracy: 94.6658%

(base) [ntavakoli6@atl1-1-01-018-36-0 hged]$ ./vg stats -a mapped_reads.chr22_pruned.gam 
# Total alignments: 10000
# Total primary: 10000
# Total secondary: 0
# Total aligned: 9999
# Total perfect: 0
# Total gapless (softclips allowed): 9524
# Total paired: 0
# Total properly paired: 0
# Insertions: 1769 bp in 512 read events
# Deletions: 19 bp in 10 read events
# Substitutions: 124291 bp in 123429 read events
# Softclips: 4027 bp in 152 read events




#####################################

# using graphaligner

/usr/bin/time -v  $GraphAligner -g chr22.vg -f  simulated_reads_len_${LEN}_errRate_$EROR.fastq -a GraphAligner.${CHRID}.aligned.gam -x vg

#For 1000
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
# Input reads: 10000 (10000000bp)
# Seeds found: 3855004
# Seeds extended: 49998
# Reads with a seed: 10000 (10000000bp)
# Reads with an alignment: 10000 (9998022bp)
# Alignments: 11056 (11029817bp) (8437 additional alignments discarded)
# End-to-end alignments: 10737 (10737000bp)







# Reads with an alignment: 10000 (9998022bp)
# Alignment accuracy: (9998022/11000000)*100
# We consider a read correctly mapped if \
# its longest alignment overlaps at least 10% with\
#  the genomic position from where it was simulated \
#  and evaluate the number of reads correctly aligned.

./vg stats -a GraphAligner.${CHRID}.aligned.gam