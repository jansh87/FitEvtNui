
set xrange[-0.5:29.5]
set yrange[-3:2]

set xlabel "#Nuisance-Parameter"
set ylabel "Fit value"

set key off
 set xtics 1

f(x) = 0

plot f(x) notitle ls 0
replot "mm2_sideband_nuis.txt" u :2:3 w yerrorbars ls 1 lc 1 lw 2 

set term postscript enhanced color
set output "mm2_sideband_nuis.eps"
replot
