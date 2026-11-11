#!/bin/bash
#SBATCH --job-name=SC_counts
#SBATCH --output=SC_counts-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=400GB
#SBATCH --time=10-23:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user= (email)  

module load cellranger/7.0.0

cellranger count --localcores=48 --localmem=400 --id=Tx430_D29 --transcriptome=Sbicolor_v3.1.1 --fastqs=. --sample=Tx430_D29
cellranger count --localcores=48 --localmem=400 --id=Tx430_D31 --transcriptome=Sbicolor_v3.1.1 --fastqs=. --sample=Tx430_D31
cellranger count --localcores=48 --localmem=400 --id=N609_D29 --transcriptome=Sbicolor_v3.1.1 --fastqs=. --sample=N609_D29
cellranger count --localcores=48 --localmem=400 --id=N609_D31 --transcriptome=Sbicolor_v3.1.1 --fastqs=. --sample=N609_D31
