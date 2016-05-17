#!/bin/sh
N=8192
i=0
while [ $i -lt 14 ]; do
    a=`expr $i \* $N + 1`
    c=`expr $i + 1`
    b=`expr $c \* $N`
    sh middle.sh data.txt $a $b > result$i.dat
    i=`expr $i + 1`
done
