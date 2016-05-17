st=1
while [ $st -le 20 ]; do
   num=`expr $st \* 2049`
   numtwo=`expr \( $st + 1 \) \* 2049`
   echo $num
   echo $numtwo
   awk "NR==$num,NR==$numtwo" data.txt > result${st}.dat
   st=`expr $st + 1`
   gnuplot << EOF
   set terminal postscript eps color enhanced
   set output "result${st}.eps"
   plot "result${st}.dat" using 1:2 with lines
   set output
EOF
done
