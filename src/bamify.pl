#!/bin/bash 
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

for SAM in $DIR/../input/*sam;
do
  echo $SAM
  samtools view -bS $SAM | samtools sort -o ${SAM/.sam/.bam} -
  samtools index ${SAM/.sam/.bam}
done