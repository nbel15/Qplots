#!/bin/bash
set -e
echo "$(tput setaf 2)## Add variables:$(tput sgr0)"
echo "3. Percent identity for OTUs classifier (numeric) (default: 0.97)"
read p_perc_identity
p_perc_identity="${p_perc_identity:=0.97}"



echo "Create Output directory : outqiime/"
mkdir -p outqiime/

echo $'\n'
echo "$(tput setaf 3)#######################################################################"
echo '09: Import OTUs (fna) to Qiime2'
echo "#######################################################################$(tput sgr0)"
source /Users/gtsiamis/miniconda3/bin/activate qiime2-2018.11
qiime tools import --input-path otus_trimmed.fa \
                   --output-path outqiime/otus_trimmed.qza \
                   --type 'FeatureData[Sequence]'

echo $'\n'
echo "$(tput setaf 3)#######################################################################"
echo '10: Assign taxa'
echo "#######################################################################$(tput sgr0)"
time qiime feature-classifier classify-consensus-blast \
           --i-query outqiime/otus_trimmed.qza \
           --i-reference-reads /Users/naima/UNITE/unite-ver7-99-seqs-01.12.2017.qza \
           --i-reference-taxonomy /Users/naima/UNITE/unite-ver7-99-tax-01.12.2017.qza \
           --o-classification outqiime/taxonomy_blast_unit.qza \
           --p-perc-identity $p_perc_identity \
           --p-min-consensus 0.8 \
           --p-maxaccepts 1

echo $'\n'
echo "$(tput setaf 3)######################################################################"
echo '13: Align Sequences'
echo "######################################################################$(tput sgr0)"
qiime alignment mafft --i-sequences outqiime/otus_trimmed.qza \
                --p-n-threads 0 \
                --o-alignment outqiime/otus_trimmed_aligned

echo $'\n'
echo "$(tput setaf 3)######################################################################"
echo '14: Phylogeny tree'
echo "######################################################################$(tput sgr0)"
qiime phylogeny fasttree --i-alignment outqiime/otus_trimmed_aligned.qza \
                --p-n-threads 0 \
                --o-tree outqiime/Tree

echo $'\n'
echo "$(tput setaf 3)######################################################################"
echo '15: midpoint tree'
echo "######################################################################$(tput sgr0)"
qiime phylogeny midpoint-root --i-tree outqiime/Tree.qza \
      --o-rooted-tree outqiime/Tree_midpoint

qiime tools export --input-path outqiime/Tree_midpoint.qza --output-path outqiime/Tree_midpoint
qiime tools export --input-path outqiime/taxonomy_blast_unit.qza --output-path outqiime/taxonomy

mkdir -p outqiime/qplot_input/
cp outqiime/Tree_midpoint/tree.nwk outqiime/qplot_input/
cp outqiime/taxonomy/taxonomy.tsv outqiime/qplot_input/

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
