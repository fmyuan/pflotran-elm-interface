set style line 1 lt 1 lw 3 lc rgb "blue"
set style line 2 lt 2 lw 3 lc rgb "red"
set style line 3 lt 3 lw 3 lc rgb "green"
set style line 4 lt 1 lw 1 lc rgb "blue"
set style line 5 lt 2 lw 1 lc rgb "red"
set style line 6 lt 3 lw 1 lc rgb "green"

set xrange[0:100]
set yrange[0:1]
set xlabel 'Time [d]' 
set ylabel 'Concentration [mol/m^3]' 
plot [0:100] exp(-0.07*x) title 'SOM1-Analytical' with points ls 4, \
     0.72*(1-exp(-0.07*x)) title 'SOM2-Analytical' with points ls 5, \
     0.28*(1-exp(-0.07*x)) title 'CO2-Analytical' with points ls 6, \
     'som2-obs-0.tec' using 1:($5) title 'Lit1 (CLM-DECOMP)' with lines ls 1, \
     'som2-obs-0.tec' using 1:($6) title 'SOM1 (CLM-DECOMP)' with lines ls 2, \
     'som2-obs-0.tec' using 1:($4*250) title 'CO2 (CLM-DECOMP)' with lines ls 3

set terminal postscript eps size 4, 3 enhanced color font 'Helvetica, 20' lw 1 
set output 'som2.eps'

replot
set term pop
pause -1

