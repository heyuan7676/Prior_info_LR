#!/bin/bash -l
#SBATCH --time 100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue

# =============================

cd /home/yhe23/Predict_target_gene/sCGG


for chr in {3..22}
do
    echo chr${chr}
    matlab -r "callGLasso('$chr', '0'); exit;" >> callGLasso.output.txt
    matlab -r "callGLasso('$chr', '1'); exit;" >> callGLasso.output.txt
done

