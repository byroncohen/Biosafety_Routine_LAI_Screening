#!/bin/bash
#SBATCH -J PaperOne_4_2_2024_freq # A single job name for the array
#SBATCH -c 1 # Number of cores
#SBATCH -p serial_requeue # Partition
#SBATCH --mem 10000 # Memory request (10Gb)
#SBATCH -t 3-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o PaperOne_4_2_2024_freq%A_%a.out # Standard output
#SBATCH -e PaperOne_4_2_2024_freq%A_%a.err # Standard error
module load R/4.2.2-fasrc01 
R CMD BATCH /n/home01/byroncohen/PaperOne_4_2_2024/PaperOne_4_2_2024_freq${SLURM_ARRAY_TASK_ID}.R
