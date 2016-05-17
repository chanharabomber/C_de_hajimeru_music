reset
input = sprintf("result%d.dat", i)
option = sprintf("using 0:1 with lines")
set term postscript enhanced eps color
set output sprintf("result%d.eps", i)
plot input using 0:1 with lines 
plot "result20.dat" using 0:1 with lines
set output
reset
