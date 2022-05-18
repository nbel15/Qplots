#!/bin/bash
set -e
echo "$(tput setaf 2)## Add variables:$(tput sgr0)"
echo "1. Normalization value:Number of reads per samples (numeric) (default: 50000):"
read sample_size
sample_size="${sample_size:=50000}"

echo "2. Minimum size for an OTU as fraction of all OTUs (numeric) (default: 0.001):"
read min_otu_freq
min_otu_freq="${min_otu_freq:=0.001}"

echo "3. Percent identity for OTUs classifier (numeric) (default: 0.97)"
read p_perc_identity
p_perc_identity="${p_perc_identity:=0.97}"



echo "Create Output directory : out/"
mkdir -p out/

echo "$(tput setaf 3)#####################################################################"
echo "00: Assembly"
echo "#####################################################################$(tput sgr0)"
usearch11 -fastq_mergepairs fq/*_R1*.fastq -fastq_maxdiffs 10 \
          -fastq_pctid 20 -fastq_minmergelen 400 -fastq_maxmergelen 520 \
          -fastqout out/merged.fq -relabel @ -log out/merge.log
echo $'\n'
echo "$(tput setaf 3)#####################################################################"
echo "1: Filter sequences"
echo "#####################################################################$(tput sgr0)"
usearch11 -fastq_filter out/merged.fq -fastq_maxee 1.0 -fastaout out/filtered.fa

echo $'\n'
echo "$(tput setaf 3)#######################################################################"
echo '2: Find uniques'
echo "######################################################################$(tput sgr0)"
usearch10 -fastx_uniques out/filtered.fa -fastaout out/uniques.fa \
          -relabel Uniq -sizeout

echo $'\n'
echo "$(tput setaf 3)######################################################################"
echo '3: Cluster OTUs'
echo "######################################################################$(tput sgr0)"
usearch10 -cluster_otus out/uniques.fa -minsize 2 -otus out/otus.fa \
          -relabel Otu

echo $'\n'
echo "$(tput setaf 3)######################################################################"
echo '4: OTU table'
echo "######################################################################$(tput sgr0)"
usearch10 -otutab out/merged.fq -otus out/otus.fa \
          -notmatched out/notmatched.fa -otutabout out/otutab.txt

echo $'\n'
echo "$(tput setaf 3)######################################################################"
echo '04.0: Normalize OTU table'
echo "######################################################################$(tput sgr0)"
usearch10 -otutab_norm out/otutab.txt -sample_size $sample_size \
          -output out/otutab_norm.txt

echo $'\n'
echo "$(tput setaf 3)######################################################################"
echo '04.1: Uncross'
echo "######################################################################$(tput sgr0)"
usearch10 -uncross out/otutab_norm.txt -tabbedout out/out.txt \
          -report out/rep.txt -otutabout out/result.txt

echo $'\n'
echo "$(tput setaf 3)#######################################################################"
echo '04.2: Trim OTU table'
echo "######################################################################$(tput sgr0)"
usearch11 -otutab_trim out/result.txt -min_otu_freq $min_otu_freq \
          -output out/trimmed.txt

col_otu=$(awk '{print NF; exit}' out/otutab_norm.txt)
col_otu_trim=$(awk '{print NF; exit}' out/trimmed.txt)
if [ "$col_otu" != "$col_otu_trim" ]; then
	echo "$(tput bold)$(tput setaf 1)#####################################################################################"
	echo "                         WARNING: Different number of samples !!!!!                          "
	echo "#####################################################################################$(tput sgr0)"
	exit 0
fi
 
echo $'\n'
echo "$(tput setaf 3)#######################################################################"
echo '08: Grep'
echo "#######################################################################$(tput sgr0)"
awk '{print $1 }' out/trimmed.txt > out/otus2trim.txt
usearch10 -fastx_getseqs out/otus.fa -labels out/otus2trim.txt \
          -fastaout out/otus_trimmed.fa

echo $'\n'
echo "$(tput setaf 3)#######################################################################"
echo '09: Import OTUs (fna) to Qiime2'
echo "#######################################################################$(tput sgr0)"
source /Users/gtsiamis/miniconda3/bin/activate qiime2-2019.1
qiime tools import --input-path out/otus_trimmed.fa \
                   --output-path out/otus_trimmed.qza \
                   --type 'FeatureData[Sequence]'

echo $'\n'
echo "$(tput setaf 3)#######################################################################"
echo '10: Assign taxa'
echo "#######################################################################$(tput sgr0)"
time qiime feature-classifier classify-consensus-blast \
           --i-query out/otus_trimmed.qza \
           --i-reference-reads /Users/gtsiamis/SILVA/qiime2/99_otus_16S.qza \
           --i-reference-taxonomy /Users/gtsiamis/SILVA/qiime2/majority_taxonomy_all_levels.qza \
           --o-classification out/taxonomy_blast_silva.qza \
           --p-perc-identity $p_perc_identity \
           --p-min-consensus 0.8 \
           --p-maxaccepts 1

echo $'\n'
echo "$(tput setaf 3)######################################################################"
echo '13: Align Sequences'
echo "######################################################################$(tput sgr0)"
qiime alignment mafft --i-sequences out/otus_trimmed.qza \
                --p-n-threads 0 \
                --o-alignment out/otus_trimmed_aligned

echo $'\n'
echo "$(tput setaf 3)######################################################################"
echo '14: Phylogeny tree'
echo "######################################################################$(tput sgr0)"
qiime phylogeny fasttree --i-alignment out/otus_trimmed_aligned.qza \
                --p-n-threads 0 \
                --o-tree out/Tree

echo $'\n'
echo "$(tput setaf 3)######################################################################"
echo '15: midpoint tree'
echo "######################################################################$(tput sgr0)"
qiime phylogeny midpoint-root --i-tree out/Tree.qza \
      --o-rooted-tree out/Tree_midpoint

qiime tools export --input-path out/Tree_midpoint.qza --output-path out/Tree_midpoint
qiime tools export --input-path out/taxonomy_blast_silva.qza --output-path out/taxonomy

mkdir -p out/qplot_input/
cp out/Tree_midpoint/tree.nwk out/qplot_input/
cp out/taxonomy/taxonomy.tsv out/qplot_input/
cp out/trimmed.txt out/qplot_input/

echo $'\n'
echo "$(tput setaf 5)*********************************************************************"
echo "******* Find Qplots input in $(tput bold)qplot_input folder *********$(tput sgr0)" 
echo $'\n'
awk '{for(i=2;i<=NF;i++) sum[i]+=$i} END{for(i in sum) if(sum[i]<1000) print "Column " i-1 " is " sum[i];}' out/trimmed.txt > out/samples_warning.txt
echo "$(tput setaf 5)*********************************************************************"
echo "********* Find Critical samples in samples_warning.txt file *********$(tput sgr0)"
echo $'\n'
echo "###################################################################################"
echo "###################################################################################"
echo "                                      END SCRIPT                                   "
echo "###################################################################################"
echo "###################################################################################"
