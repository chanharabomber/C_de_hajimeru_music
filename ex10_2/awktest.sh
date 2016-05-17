st=1
while [ $st -le 20 ]; do
   num=`expr $st \* 2047`
   numtwo=`expr \( $st + 1 \) \* 2047`
   echo $num
   echo $numtwo
   awk "NR==$num,NR==$numtwo" data.txt > result${st}.dat
   st=`expr $st + 1`
done
