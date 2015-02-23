#!/bin/sh
#SBATCH -n ${ntasks}
#SBATCH -N ${nodes}
#SBATCH -J ${jobName}
#SBATCH -D ${workdir}
${SBATCHoptions}

echo "[$(date +"%Y-%m-%d %H:%M:%S")] Job started" > ${jobName}.$SLURM_JOB_ID.started

${script}

echo "[$(date +"%Y-%m-%d %H:%M:%S")] Job finished" >> ${jobName}.$SLURM_JOB_ID.started
mv ${jobName}.$SLURM_JOB_ID.started ${jobName}.$SLURM_JOB_ID.finished
