
set xrange[-0.5:4.5]
set yrange[0.88:1.15]

set xlabel "#Stream"
set ylabel "Fit value"

f(x) = 1

plot f(x) notitle ls 0
replot "streamtest_plot.txt" index 0 u 1:2:3 w yerrorbars title "Signal" ls 1 lc 1 lw 2 
replot "streamtest_plot.txt" index 1 u 1:2:3 w yerrorbars title "Xlnu" ls 1 lc 3 lw 2

set term postscript enhanced color
set output "streamtest.eps"
replot
