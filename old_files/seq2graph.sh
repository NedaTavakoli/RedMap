REF=hs37d5.fa
CHRID=chr22
VARFILE_SNP=ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf
LEN=10000
COUNT=100
k8=$HOME/k8-0.2.5/k8-linux
hisat2=$HOME/hisat-genotype-top/hisat2
vg=$HOME/vg
FORGe=$HOME/FORGe
tabix=$HOME/vg/tabix
bgzip=$HOME/vg/bgzip
GraphAligner=$HOME/GraphAligner

#script assumes vg, tabix, bgzip, GraphAligner symlinked in curent directory

#Construct graph
./vg/vg construct -r $REF -v $VARFILE_SNP -a -f -m 32 > graph.${CHRID}.vg

#./vg/vg construct -r $REF -v $VARFILE_SNP -a -f -m 32 -R chr22  > graph.${CHRID}.vg

#index the graph
./vg/vg index -x graph.${CHRID}.xg -g graph.${CHRID}.gcsa graph.${CHRID}.vg

#Simulating reads
./vg/vg sim -x graph.${CHRID}.xg -n $COUNT -l $LEN -a > reads.${CHRID}.gam
./vg/vg view -X reads.${CHRID}.gam > reads.${CHRID}.fastq
./vg/vg surject -sx graph.${CHRID} reads.${CHRID}.gam -.xg > reads.${CHRID}.sam #truth on linear genome

#Run FORGe
mkdir -p  ${CHRID}_hisat_index 

$FORGe/src/vcf_to_1ksnp.py  --reference $REF --vcf $VARFILE_SNP --out ${CHRID}_snp.1ksnp 
$FORGe/src/rank.py --method popcov --reference $REF  --vars ${CHRID}_snp.1ksnp   --window-size 100  --output ordered_${CHRID}_snp.txt
$FORGe/src/build.py --reference  $REF --vars ${CHRID}_snp.1ksnp --window-size 100 --hisat hisat_input_${CHRID}_snp.snp --sorted ordered_${CHRID}_snp.txt --pct 10
$hisat/hisat2-build --snp hisat_input_${CHRID}_snp.snp $REF ${CHRID}_hisat_index /index_${CHRID}

#map using HISAT2: simulated reads to graph indexes
$hisat/hisat2 -f -x $FORGe/${CHRID}_hisat_index/index_${CHRID} -U reads.${CHRID}.fastq -S hisat_alignment_${CHRID}_snp.sam #reads.${CHRID}.fastq are simulated reads

#map by graph aligner
./GraphAligner/GraphAligner -g graph.${CHRID}.vg -f reads.${CHRID}.fastq -a GraphAligner.${CHRID}.aligned.gam -x vg
./vg/vg surject -s GraphAligner.${CHRID}.aligned.gam -x graph.${CHRID}.xg > GraphAligner.${CHRID}.aligned.sam #linearize

#convert sam to paf
#sam2paf available here https://github.com/lh3/miniasm/tree/master/misc
#had to comment out line #61 to get sam2paf running
$k8 sam2paf.js reads.${CHRID}.sam > reads.${CHRID}.paf
$k8 sam2paf.js GraphAligner.${CHRID}.aligned.sam > GraphAligner.${CHRID}.aligned.paf

