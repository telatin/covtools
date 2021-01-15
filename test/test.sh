#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )/../"
set -euo pipefail
echo "Dir: $DIR"

echo " - Preparing BAM files (from input/*.sam:"
$DIR/src/bamify

echo -n " - Compiling: "

nim c --verbosity:0 --hints:off --threads:on --warnings:off -o:$DIR/test/covtobed_test $DIR/src/covtobed2.nim && echo " OK " || echo " FAILED ****"
ERRORS=0

CNT=0
for INPUT in $DIR/input/*.bam "$@";
do
  CNT=$((CNT+1))
  echo " [$CNT] Testing $(basename $INPUT)"
  echo -n "       - Standard: "
  TEST=$($DIR/test/covtobed_test -F 0 -Q 0 $INPUT  | md5sum)
  PROD=$(covtobed $INPUT | md5sum)
  if [[ $TEST == $PROD ]]; then
    echo "OK"
  else
    echo "ERROR        <<<<<< "
    ERRORS=$((ERRORS+1))
  fi

  echo -n "       - Physical output: "
  TEST=$($DIR/test/covtobed_test -F 0 -Q 0 --physical $INPUT  | md5sum)
  PROD=$(covtobed --physical $INPUT | md5sum)
  if [[ $TEST == $PROD ]]; then
    echo "OK"
  else
    echo "ERROR        <<<<<< "
    ERRORS=$((ERRORS+1))
  fi

  echo -n "       - Flags: "
  TEST=$($DIR/test/covtobed_test  -Q 1 $INPUT  | md5sum)
  PROD=$(covtobed -d $INPUT | md5sum)
  if [[ $TEST == $PROD ]]; then
    echo "OK"
  else
    echo "ERROR        <<<<<< "
    ERRORS=$((ERRORS+1))
  fi

       
done

echo " Done. $ERRORS errors."