module load trimmomatic
module load fastqc


samplelist=${samplelist}
arraysize=${arraysize}

numberOfSamples=$(cat $samplelist | wc -l)


for i in $(seq 0 $arraysize $numberOfSamples);
do
  echo "i = $i"
  ID=$(echo $SLURM_ARRAY_TASK_ID-1+$i | bc)
  echo "ID = $ID"
  if [ $ID -gt $numberOfSamples ]
  then
    echo "stopped because $ID > $numberOfSamples"
    break
  fi
  
  params=$(cat $samplelist|awk ' NR=='$ID' { print $0 ; }')
  echo "parameters = $params"

  
  bash RUNtrimmo.sh $params
  
done
