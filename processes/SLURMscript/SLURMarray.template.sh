#!/bin/bash -l
#SBATCH -n ${ntasks}
#SBATCH -N ${nodes}
#SBATCH -J ${jobName}
#SBATCH -D ${workdir}
#SBATCH -a 1-${arraySize}
${SBATCHoptions}

logFile=${jobName}.$SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID.running
failedCommands=failedCommands.$SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID

echo "[$(date +"%Y-%m-%d %H:%M:%S")] Job started" > $logFile

${preScript}

commandList=${commandList}
arraySize=${arraySize}

nCommands=$(cat $commandList | wc -l)
nCommandsThis=$((($nCommands+$arraySize-$SLURM_ARRAY_TASK_ID)/$arraySize))

echo "Distributed command execution started."
echo ""
echo "This is executer task $SLURM_ARRAY_TASK_ID of $arraySize"
echo "Commands to execute: $nCommandsThis of $nCommands"
echo ""

for i in $(seq $SLURM_ARRAY_TASK_ID $arraySize $nCommands);
do
  CMD=$(cat $commandList|awk ' NR=='$i' { print $0 ; }')
  
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] CMD $i: $CMD"
  
  echo $CMD | bash
  
  exitCode=$?
  if [ "$exitCode" = "0" ]; then
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] Command $i finished." >> $logFile
  else
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] Command $i failed! Exit code=$exitCode" >> $logFile
    echo "$CMD" >> $failedCommands
  fi
  
done



if [ -f $failedCommands ];
then
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] Job finished. $(cat $failedCommands | wc -l)/$nCommandsThis commands exited with ERROR!" >> $logFile
  mv $logFile ${jobName}.$SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID.failed
  exit 1
else
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] Job finished" >> $logFile
  mv $logFile ${jobName}.$SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID.finished
fi

