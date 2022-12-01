#!/bin/bash
#SBATCH --job-name=med_job_null
#SBATCH --output=med_job_null
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=40G
#SBATCH --time=135:00:00

module load languages/r/4.1.0



cd /mnt/storage/scratch/bu19525/BristolPhD/MVMRproject/adjustedGWAS/Simulation_Analyses/Null_sims/Null_models/Null_med_exp/
Rscript "sims_code_med_nullexp_1.R" med_res_1.txt
Rscript "sims_code_med_nullexp_2.R" med_res_2.txt
Rscript "sims_code_med_nullexp_3.R" med_res_3.txt
Rscript "sims_code_med_nullexp_4.R" med_res_4.txt
Rscript "sims_code_med_nullexp_5.R" med_res_5.txt