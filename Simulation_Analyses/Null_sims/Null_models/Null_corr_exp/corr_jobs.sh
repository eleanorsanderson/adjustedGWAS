#!/bin/bash
#SBATCH --job-name=corr_job_null
#SBATCH --output=corr_job_null
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=40G
#SBATCH --time=135:00:00

module load languages/r/4.1.0



cd /mnt/storage/scratch/bu19525/BristolPhD/MVMRproject/adjustedGWAS/Simulation_Analyses/Null_sims/Null_models/Null_corr_exp/
Rscript "sims_code_corr_nullexp_1.R" corr_res_1.txt
Rscript "sims_code_corr_nullexp_2.R" corr_res_2.txt
Rscript "sims_code_corr_nullexp_3.R" corr_res_3.txt
Rscript "sims_code_corr_nullexp_4.R" corr_res_4.txt
Rscript "sims_code_corr_nullexp_5.R" corr_res_5.txt