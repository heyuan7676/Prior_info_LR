#!/bin/bash -l
#SBATCH --time 100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue

# =============================

cd /home/yhe23/Predict_target_gene/sCGG



chr="$1"
chr_out="$2"

export chr
export chr_out
python simulate.py



matlab -r "simulate('$chr_out'); exit;"


cd /scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/

rm -f temp.txt
head -n1 GTEx_chr22_pairs_activeSNP_genes_corrected.txt > temp.txt
cat GTEx_chr${chr_out}_pairs_activeSNP_genes_corrected.txt >> temp.txt
mv temp.txt GTEx_chr${chr_out}_pairs_activeSNP_genes_corrected.txt

