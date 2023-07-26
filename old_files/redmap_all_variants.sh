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


#======================================================================================================================
#***************************  All Variants   *************************** 
 # 1
# # compelte graph (all variants) and mapping [without pruning]
cd $scratch
/usr/bin/time -v  $project_dir/./vg index -x $project_dir/chr22.xg -g $project_dir/chr22.vg --gcsa-out chr22.gcsa 
/usr/bin/time -v ./vg map -t 24 -x chr22.xg $scratch/chr22.gcsa -f  simulated_reads_len_${LEN}_errRate_$EROR.fastq --gam-output > mapped_reads.chr22.gam
# /usr/bin/time -v ./vg map -t 24 -x chr22.xg -g $scratch/chr22.gcsa -f  simulated_reads_len_${LEN}_errRate_$EROR.fastq > mapped_reads.chr22.gam
./vg view -aj mapped_reads.chr22.gam | jq '.score' > scores_mapped_reads.chr22.txt
awk '{ sum += $1 } END { print "Overall Alignment Score:", sum; accuracy = sum /(1000*NR); print "Alignment Accuracy:", accuracy*100 "%" }' scores_mapped_reads.chr22.txt


#####################################

# using graphaligner
# install: [remove later]
# git clone https://github.com/maickrau/GraphAligner.git
# cd GraphAligner
# git submodule update --init --recursive
# conda env create -f CondaEnvironment_linux.yml or conda env create -f CondaEnvironment_osx.yml
# source activate GraphAligner
# make bin/GraphAligner


$GraphAligner -g chr22.vg -f  simulated_reads_len_${LEN}_errRate_$EROR.fastq -a GraphAligner.${CHRID}.aligned.gam -x vg
./vg stats -a GraphAligner.${CHRID}.aligned.gam

