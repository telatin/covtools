#!/bin/bash
DEBUG=" -x "
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )/../"
set -euxo pipefail

PLATFORM=""
if [ $(uname) == "Darwin " ]; then
 PLATFORM="_mac"
fi

for BINARY in covtotarget covtocounts covtobed2;
do
	[ -e "$DIR/src/$BINARY" ] && rm "$DIR/src/$BINARY"
done

echo " * Building test BAM (from input/*.sam)"
for i in $DIR/input/*.sam
do
  if [ ! -e "${i/sam/bam}.bai" ]; then
	  samtools view -bS $i | samtools sort -o ${i/sam/bam} -
    samtools index  ${i/sam/bam}
	fi
done


mkdir -p $DIR/{bin,release}/
for SOURCE in covtocounts.nim covtotarget.nim covtobed2.nim covtocounts2.nim;
do
	echo "   $SOURCE"
	OUTBIN=$(basename $SOURCE | cut -f1 -d.)
	echo "   - Compile dyn. link"
	nim c -o:$DIR/bin/${OUTBIN}${PLATFORM}_debug $DIR/src/$SOURCE  > $DIR/.compile.log 2>&1 

	
  if [ $(uname) == "Linux" ];
	then
		 
     if [ -e "/local/miniconda3/lib/libhts.a" ];
		 then
		 	echo "   - Compile static (only ubuntu)"
		   nim c -d:static -d:release -o:$DIR/release/${OUTBIN}${PLATFORM}  $DIR/src/$SOURCE  > $DIR/.compile.log 2>&1 
		 else
		 echo "   - Compile static (docker)"
		   cp $DIR/src/${SOURCE}ble $DIR/data.nimble
			 $DIR/bin/hts_nim_static_builder $DEBUG -s src/$SOURCE -n $DIR/data.nimble -- --opt:speed > $DIR/.docker.log 2>&1 && \
			   sudo chown ubuntu $DIR/${OUTBIN} && \
				 mv -v $DIR/${OUTBIN} $DIR/release/
			 rm $DIR/data.nimble
		 fi
	else
	  echo " Skipping macOS"
	fi

  # if [ $USER == 'ubuntu' ]; then
	# 	cp $DIR/src/${SOURCE}ble $DIR/data.nimble
	# 	$DIR/bin/hts_nim_static_builder $DEBUG -s src/$SOURCE -n $DIR/data.nimble -- --opt:speed > $DIR/.docker.log 2>&1 && \
	# 	sudo chown ubuntu $DIR/${OUTBIN} && \
	# 	mv -v $DIR/${OUTBIN} $DIR/release/
	# 	rm $DIR/data.nimble
	# 				if [ ! -e $DIR/release/${OUTBIN} ]; then
	# 					echo "Static build failed."; exit;
	# 				fi
  # fi
done

if [ $USER == 'ubuntu' ]; then
  sudo chown ubuntu:ubuntu $DIR/bin/*
fi