set xrange [0:5]
set yrange [0:3]
plot "-" using 1:2:3 w yerrorbars
#  X     Y     Z 
   1.0   1.2   0.2
   2.0   1.8   0.3
   3.0   1.6   0.2
   4.0   1.2   0.2
end
pause -1
