#!/bin/sh

st=$1
while [ $st -le $2 ]; do
  if [ ! -e $(printf "result%d.eps" $st) ]; then
    gnuplot -e "i=${st}" loop.plt
  fi
  st=`expr $st + 1`
done
