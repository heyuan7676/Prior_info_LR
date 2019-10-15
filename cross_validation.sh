#!/bin/bash -l
#SBATCH --time 100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --no-requeue

# =============================

cd /home/yhe23/Predict_target_gene/sCGG

chr="$1"
matlab -r "cross_validation('$chr'); exit;"
