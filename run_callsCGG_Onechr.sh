#!/bin/bash -l
#SBATCH --time 100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue

# =============================

cd /home/yhe23/Predict_target_gene/sCGG

chr=24
lambda2=(0.1 0.5)
C0=(5 10)

for l in "${lambda2[@]}"
do
    for c0 in "${C0[@]}"
    do
        c1=$(($c0 / 5))
        echo $l, $c1, $c0
        sbatch callsCGG_Onechr.sh $chr $l $c1 $c0
    done
done

