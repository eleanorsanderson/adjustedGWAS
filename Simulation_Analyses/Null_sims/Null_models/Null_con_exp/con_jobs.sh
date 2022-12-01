#!/bin/bash
#SBATCH --job-name=con_job_null
#SBATCH --output=con_job_null
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=40G
#SBATCH --time=135:00:00

module load languages/r/4.1.0

cd /mnt/storage/scratch/bu19525/BristolPhD/MVMRproject/adjustedGWAS/Simulation_Analyses/Null_sims/Null_models/Null_con_exp/
Rscript "sims_code_con_nullexp_1.R" con_res_1.txt
Rscript "sims_code_con_nullexp_2.R" con_res_2.txt
Rscript "sims_code_con_nullexp_3.R" con_res_3.txt
Rscript "sims_code_con_nullexp_4.R" con_res_4.txt
Rscript "sims_code_con_nullexp_5.R" con_res_5.txt