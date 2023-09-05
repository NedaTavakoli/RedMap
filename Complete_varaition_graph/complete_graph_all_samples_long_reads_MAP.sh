
###### long reads
##################################################################################################
##################################################################################################

CHRID=chr22
LEN=1000
EROR=0.025
INDEL_error=0.1
COUNT=40 # Total number of reads you want to simulate

simulated_dir=/storage/hive/project/cse-aluru/ntavakoli6/hged/simulated_reads
output_gam_file=${simulated_dir}/simulated_reads_len_${LEN}_errRate_${EROR}_indelRate_${INDEL_error}_all_samples.gam
output_fastq_file=${simulated_dir}/simulated_reads_len_${LEN}_errRate_${EROR}_indelRate_${INDEL_error}_all_samples.fasta
output_sam_file=${simulated_dir}/simulated_reads_len_${LEN}_errRate_${EROR}_indelRate_${INDEL_error}_all_samples.sam
simulated_read_paf_file=${simulated_dir}/simulated_reads_len_${LEN}_errRate_${EROR}_indelRate_${INDEL_error}_all_samples.paf
mapped_reads_gam=${simulated_dir}/mapped_reads.${CHRID}_len_${LEN}_errRate_${EROR}_indelRate_${INDEL_error}_all_samples.gam
scores_mapped_reads=${simulated_dir}/scores_mapped_reads.${CHRID}_len_${LEN}_errRate_${EROR}_indelRate_${INDEL_error}_all_samples.txt
mapped_reads_with_compare_option=${simulated_dir}/mapped_reads.${CHRID}_len_${LEN}_errRate_${EROR}_indelRate_${INDEL_error}_all_samples_with_compare_option.txt


##################################################################################################
####                  Simulating reads from haplotypes
##################################################################################################
cd /storage/hive/project/cse-aluru/ntavakoli6/hged


# Simulate reads using vg sim and append to the output GAM file
sample_list_file="/storage/hive/project/cse-aluru/ntavakoli6/hged/samples.txt"  # Path to your sample names file

sample_names=()
while IFS= read -r line; do
    sample_names+=("$line")
done < "$sample_list_file"



# Randomly shuffle the sample names array
shuf_sample_names=($(shuf -e "${sample_names[@]}"))


while [[ ${#sample_names[@]} -gt 0 ]]; do
    random_index=$((RANDOM % ${#sample_names[@]}))
    random_sample_name="${sample_names[$random_index]}"
    
    ./vg sim -x ${CHRID}.xg -g ${CHRID}.gbwt -m ${random_sample_name} -e ${EROR} -i ${INDEL_error} -n ${COUNT} -l ${LEN} -a  >> ${output_gam_file}

    # Process the randomly selected sample name here
    
    # Remove the selected sample name from the array
    unset "sample_names[$random_index]"
    sample_names=("${sample_names[@]}")
done


/usr/bin/time -v ./vg view -X ${output_gam_file} > ${output_fastq_file}
/usr/bin/time -v ./vg surject -s ${output_gam_file} -x ${CHRID}.xg > ${output_sam_file} #truth on linear genome

##################################################################################################
####                 Mapping reads to the complete graph
##################################################################################################

# Mapping the simulated reads to the original graph uisng vg
LEN=1000
EROR=0.025
INDEL_error=0.1
CHRID=chr22

output_gam_file=${simulated_dir}/simulated_reads_len_${LEN}_errRate_${EROR}_indelRate_${INDEL_error}_all_samples.gam
mapped_reads_with_compare_option=${simulated_dir}/mapped_reads.${CHRID}_len_${LEN}_errRate_${EROR}_indelRate_${INDEL_error}_all_samples_with_compare_option.txt


# Mapping the simulated reads to the original graph uisng vg
/usr/bin/time -v ./vg map -t 24 -x ${CHRID}.xg -g ${CHRID}.pruned.gcsa --gbwt-name ${CHRID}.gbwt -G ${output_gam_file} --compare -j > ${mapped_reads_with_compare_option}  # map with compare option 
cat ${mapped_reads_with_compare_option} | jq .correct | sed s/null/0/ | awk '{i+=$1; n+=1} END {print i/n}'
# 0.808862
/usr/bin/time -v ./vg view -aj ${mapped_reads_gam} | jq '.score' > ${scores_mapped_reads}
awk '{ sum += $1 } END { print "Total Alignment Score:", sum; avg = sum /(NR); print "Average Alignment Score:", avg }' ${scores_mapped_reads}
awk '{ sum += $1 } END { print "Overall Alignment Score:", sum; accuracy = sum /(1010*NR); print "Alignment Accuracy:", accuracy*100 "%" }' ${scores_mapped_reads}


#####################################

# using graphaligner
# install: [remove later]
# git clone https://github.com/maickrau/GraphAligner.git
# cd GraphAligner
# git submodule update --init --recursive
# conda env create -f CondaEnvironment_linux.yml or conda env create -f CondaEnvironment_osx.yml
# source activate GraphAligner
# make bin/GraphAligner

cd /storage/hive/project/cse-aluru/ntavakoli6/hged
CHRID=chr22
LEN=1000
EROR=0.025
INDEL_error=0.1

simulated_dir=/storage/hive/project/cse-aluru/ntavakoli6/hged/simulated_reads
output_fastq_file=${simulated_dir}/simulated_reads_len_${LEN}_errRate_${EROR}_indelRate_${INDEL_error}_all_samples.fasta
GraphAligner=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/GraphAligner/bin/GraphAligner

/usr/bin/time -v $GraphAligner -g ${CHRID}.vg -f ${output_fastq_file} -a GraphAligner.${CHRID}.aligned_complete_graph.gam -x vg
./vg stats -a GraphAligner.${CHRID}.aligned_complete_graph.gam
