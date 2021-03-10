#!/usr/bin/sh
ucsc_dir="/home/jxyue/Programs/UCSC_Utilities"

query_tag="S288C"
target_tag="ref"

query_assembly="$query_tag.genome.raw.relabel.fa"
target_assembly="$target_tag.genome.raw.relabel.fa"

chr_list="chrI chrII chrIII chrIV chrV chrVI chrVII chrVIII chrIX chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI"

mkdir ${query_tag}_chr
mkdir ${target_tag}_chr

$ucsc_dir/faSplit byname $query_assembly ./${query_tag}_chr/ 
$ucsc_dir/faSplit byname $target_assembly ./${target_tag}_chr/ 

# Split the G2 chromosomes/scaffolds into 3K chunks and make lift files
mkdir lift
mkdir split
for chr in $chr_list
do
  $ucsc_dir/faSplit size ./${target_tag}_chr/${target_tag}_${chr}.fa  3000 ./split/${target_tag}_${chr}.split  -lift=./lift/${target_tag}_${chr}.lft -oneFile
done
 
# run blat
mkdir psl
for chr in $chr_list
do
  $ucsc_dir/blat $query_assembly ./split/${target_tag}_${chr}.split.fa -t=dna -q=dna -tileSize=12 -fastMap -minIdentity=99 -noHead -minScore=100  ./psl/${target_tag}_${chr}.psl
done
 
# Change coordinates of .psl files to parent coordinate system
mkdir liftup
for chr in $chr_list
do
    $ucsc_dir/liftUp -pslQ ./liftup/${target_tag}_${chr}.liftup.psl ./lift/${target_tag}_${chr}.lft warn ./psl/${target_tag}_${chr}.psl
done
 
# Make chain files
mkdir chain_raw
for chr in $chr_list
do
    $ucsc_dir/axtChain -linearGap=medium  -faQ -faT -psl ./liftup/${target_tag}_${chr}.liftup.psl $query_assembly $target_assembly ./chain_raw/$chr.chain
done
 
# Merge and sort chain files
$ucsc_dir/chainMergeSort ./chain_raw/*.chain | $ucsc_dir/chainSplit chain_split stdin
 
$ucsc_dir/faSize $query_assembly  -detailed >$query_tag.chr_length.txt
$ucsc_dir/faSize $target_assembly  -detailed >$target_tag.chr_length.txt
 
# Make alignment nets from chain files
mkdir net
for i in  ./chain_split/*.chain
do
 tag=${i/\.\/chain_split\//}
 $ucsc_dir/chainNet $i ./$query_tag.chr_length.txt ./$target_tag.chr_length.txt ./net/$tag.net /dev/null
done
 
# Create liftOver chain file
mkdir over
for i in ./chain_split/*.chain
do
  tag=${i/\.\/chain_split\//}
  $ucsc_dir/netChainSubset  ./net/$tag.net $i ./over/$tag.chain
done
 
cat ./over/*.chain >${query_tag}_to_${target_tag}.over.chain

rm -r ${query_tag}_chr
rm -r ${target_tag}_chr
rm -r chain_raw
rm -r chain_split
rm -r net
rm -r over
rm -r split
rm -r lift
rm -r liftup
rm -r psl
rm $query_tag.chr_length.txt
rm $target_tag.chr_length.txt
