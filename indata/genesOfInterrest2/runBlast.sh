#!/bin/sh
#SBATCH --ntasks=1           
#SBATCH --nodes=1               
#SBATCH --job-name=blastx

module load blast+

blastx -query COR_genes_S10.fasta \
       -db /mnt/NOBACKUP/mariansc/share/orthos/orthoMCL/allProteomes.db \
	   -out COR_genes_VS_allProteomes.blastx \
	   -outfmt 6 \
	   -evalue 0.01