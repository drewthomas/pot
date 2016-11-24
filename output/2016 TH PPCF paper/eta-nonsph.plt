## eta-nonsph.plt
## gnuplot script to plot nonspherical dust charging data

set terminal postscript eps enhanced font ',25' size 4,4
set output "eta-nonsph.eps"
set size square
set xtics nomirror out
set ytics nomirror out
set encoding iso_8859_15
set ytics -2.9,0.1,-2.5
set grid xtics ytics
set key at 0.6,-2.83 box
set key height +0.5
set xlabel "A" font ',35'
set bars 2
set style line 1 lc rgb 'black' pt 7 lt 1 lw 2
set logscale x
set label "{/Symbol h}_a" font ',35' at graph -0.18, graph 0.67

set xrange [0.008: 126]
set yrange [-2.95:-2.48]
plot 'nonsph-OML-eta.dat' u 1:($2)*(-1) w l lt 1 lw 4 lc rgb "black" t ' OML' smooth bezier,\
     'eta-nonsph-1ev.dat' u 1:($2)*(-1):3 w errorbars ls 1 t 'pot'

