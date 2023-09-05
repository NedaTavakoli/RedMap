#!/bin/bash

# for the complete graph
CHRID=chr22
LEN=1000
EROR=0.025
INDEL_error=0.1

mapped_reads_gam=${simulated_dir}/mapped_reads.${CHRID}_len_${LEN}_errRate_${EROR}_indelRate_${INDEL_error}_all_samples.gam
simulated_dir=/storage/hive/project/cse-aluru/ntavakoli6/hged/simulated_reads
node_offset_file="node_offset_complete_graph.txt"

# Extract node_id and offset for the first entry of each path per read
/usr/bin/time -v ./vg view -aj "$mapped_reads_gam" |
jq -r '.path.mapping | select(length > 0) | .[0] | "\(.position.node_id) \(.position.offset)"' > "$node_offset_file"

echo "Node ID and Offset for first entry of each path per read extracted and saved to $node_offset_file"


# for the reduced graph


##################################################################################################

CHRID=chr22
LEN=1000
EROR=0.025
INDEL_error=0.1
alpha=1000
delta=1

mapped_reads_gam=${simulated_dir}/mapped_reads.${CHRID}_len_${LEN}_errRate_${EROR}_indelRate_${INDEL_error}_all_samples_redduced_graph_${alpha}_${delta}.gam
simulated_dir=/storage/hive/project/cse-aluru/ntavakoli6/hged/simulated_reads
node_offset_file="node_offset_reduced_graph.txt"

# Extract node_id and offset for the first entry of each path per read
/usr/bin/time -v ./vg view -aj "$mapped_reads_gam" |
jq -r '.path.mapping | select(length > 0) | .[0] | "\(.position.node_id) \(.position.offset)"' > "$node_offset_file"