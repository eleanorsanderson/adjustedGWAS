#!/bin/bash
#SBATCH --job-name=corr_job_null
#SBATCH --output=corr_job_null
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=40G
#SBATCH --time=35:00:00

module load languages/r/4.1.0



cd /mnt/storage/scratch/bu19525/BristolPhD/MVMRproject/adjustedGWAS/Simulation_Analyses/Null_sims/Null_models/Null_corr_out/
Rscript "sims_code_corr_nullout_1.R" corr_res_1.txt
Rscript "sims_code_corr_nullout_2.R" corr_res_2.txt
Rscript "sims_code_corr_nullout_3.R" corr_res_3.txt
Rscript "sims_code_corr_nullout_4.R" corr_res_4.txt
Rscript "sims_code_corr_nullout_5.R" corr_res_5.txt