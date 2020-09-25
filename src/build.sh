#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )/../"
set -euo pipefail

[ -e "$DIR/src/covtotarget" ] && rm "$DIR/src/covtotarget"
[ -e "$DIR/src/covtocounts" ] && rm "$DIR/src/covtocounts"

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
        cp $DIR/src/${SOURCE}ble $DIR/data.nimble
	$DIR/bin/hts_nim_static_builder -s src/$SOURCE -n $DIR/data.nimble > $DIR/.docker.log 2>&1 && \
	 sudo chown ubuntu $DIR/${OUTBIN} && \
	 mv -v $DIR/${OUTBIN} $DIR/release/
	rm $DIR/data.nimble
        if [ ! -e $DIR/release/${OUTBIN} ]; then
           echo "Static build failed."; exit;
        fi

done

sudo chown ubuntu:ubuntu $DIR/bin/*
