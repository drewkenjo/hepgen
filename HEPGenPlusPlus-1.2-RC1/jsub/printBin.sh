#!/bin/bash

qlist=`seq $1 $2 $3`
wlist=`seq $4 $5 $6`

for q in $qlist
do
  for w in $wlist
  do
    echo $q $w
  done
done

