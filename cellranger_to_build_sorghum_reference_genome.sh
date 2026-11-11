#!/bin/bash
#SBATCH --job-name=Sb_sc__build_ref
#SBATCH --output=Sb_sc__build_ref-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=150GB
#SBATCH --time=96:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user= (email) 

module load cellranger/7.0.0

module load gffread

gffread -E -T -o Sbicolor_454_v3.1.1.gene_exons.gtf Sbicolor_454_v3.1.1.gene_exons.gff3

cellranger mkref --nthreads=48 --memgb=150 --genome=Sbicolor_v3.1.1 --fasta=Sbicolor_454_v3.0.1.softmasked.fa --genes=Sbicolor_454_v3.1.1.gene_exons.gtf
