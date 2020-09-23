#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )/../"
set -euo pipefail

echo " * Building test BAM"
for i in $DIR/input/*.sam
do
  samtools view -bS $i | samtools sort -o ${i/sam/bam} -
  samtools index  ${i/sam/bam}
done

mkdir -p $DIR/{bin,release}/

for SOURCE in covtocounts.nim covtotarget.nim;
do
	echo "[$SOURCE]"
	OUTBIN=$(basename $SOURCE | cut -f1 -d.)
	echo " * Compile (dyn. link)"
	nim c -o:$DIR/bin/${OUTBIN}_debug $DIR/src/$SOURCE  > $DIR/.compile.log 2>&1 

	echo " * Compile (static)"
	$DIR/bin/hts_nim_static_builder -s src/$SOURCE -n src/${SOURCE}ble > $DIR/.docker.log 2>&1 && \
	 sudo chown ubuntu $DIR/${OUTBIN} && \
	 mv $DIR/${OUTBIN} $DIR/release/
done

sudo chown ubuntu:ubuntu $DIR/bin/*
