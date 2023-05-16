#!/bin/bash

MUPI=2.0
MUETA=2.0
MIXANGLE=1.158

DIR=$(realpath $(dirname $0)/preparation)

ARGS="$@"

if [[ -f $1 ]]
then
  ARGS=$(cut -f1 -d\  $1)
fi

if echo $ARGS | grep -qi "newmix"
then
  MIXANGLE=1.01873
else
  MIXANGLE=1.158
fi

if echo $ARGS | grep -q "2.64"
then
  MUPI=2.64
  ln -sf $DIR/mupi.eq.2.64/pi0P $DIR/
  ln -sf $DIR/mupi.eq.2.64/pi0N $DIR/
else
  MUPI=2.0
  ln -sf $DIR/mu.eq.2/pi0P $DIR/
  ln -sf $DIR/mu.eq.2/pi0N $DIR/
fi


if echo $ARGS | grep -q "2.32"
then
  MUETA=2.32
  ln -sf $DIR/mueta.eq.2.32/etaP $DIR/
elif echo $ARGS | grep -q "1.76"
then
  MUETA=1.76
  ln -sf $DIR/mueta.eq.1.76/etaP $DIR/
else
  MUETA=2.0
  ln -sf $DIR/mu.eq.2/etaP $DIR/
fi


sed 's/MUPI/'"$MUPI"'/;s/MUETA/'"$MUETA"'/;s/MIXANGLE/'"$MIXANGLE"'/' libGKPi0.template.cpp > ../libGKPi0/libGKPi0.cpp


if [[ -f $1 ]]
then
  pars=$(awk '{print("  double etNu = "$4", etNd = "$7", etbu = "$10", etbd = "$13", deltaU = "$16", alphastrU = "$19", deltaD = "$22", alphastrD = "$25";")}' $1)
  sed -i 's/^  double etNu.*/'"$pars"'/' ../libGKPi0/libGKPi0.cpp
fi
