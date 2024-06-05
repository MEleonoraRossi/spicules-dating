#!/bin/bash

#SBATCH --job-name=JTTG4f
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --tasks-per-node=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=10G

module load apps/paml/4.9i

cd /user/work/yp19290/Time_calibration_July/Cauchy_ht/Partition/JTTG4F


codeml tmp0001.ctl

