set style line 1 lt 1 lw 3 lc rgb "blue"
set style line 2 lt 2 lw 3 lc rgb "red"
set style line 3 lt 3 lw 3 lc rgb "green"
set style line 4 lt 4 lw 3 lc rgb "cyan"
set style line 5 lt 5 lw 3 lc rgb "violet"
set style line 6 lt 6 lw 3 lc rgb "black"
set style line 7 lt 7 lw 3 lc rgb "magenta"
set style line 8 lt 8 lw 3 lc rgb "yellow"


set xrange[0:100]
set yrange[0:1]
set xlabel 'Time [d]' 
set ylabel 'Concentration [mol/m^3]' 
plot 'dom1-obs-0.tec' using 1:($6) title 'SOM1' with lines ls 1, \
     'dom1-obs-0.tec' using 1:($7) title 'SOM2/SOM3' with lines ls 2, \
     'dom1-obs-0.tec' using 1:($9) title 'SOM4' with lines ls 3, \
     'dom1-obs-0.tec' using 1:($10) title 'Bacteria/Fungi' with lines ls 4, \
     'dom1-obs-0.tec' using 1:($5*250) title 'DOM' with lines ls 5, \
     'dom1-obs-0.tec' using 1:($4*250) title 'CO2' with lines ls 6, \
     'dom1-obs-0.tec' using 1:($2*250) title 'NH4+' with lines ls 7, \
     'dom1-obs-0.tec' using 1:($3*250) title 'NO3-' with lines ls 8

set terminal postscript eps size 4, 3 enhanced color font 'Helvetica, 20' lw 1 
set output 'dom1.eps'

replot
set term pop
pause -1

