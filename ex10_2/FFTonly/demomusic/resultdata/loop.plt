reset
input = sprintf("result%d.dat", i)
set term postscript enhanced eps color
#set terminal pngcairo
set output sprintf("result%d.eps", i)
plot input using 1:2 with lines
set output
reset
