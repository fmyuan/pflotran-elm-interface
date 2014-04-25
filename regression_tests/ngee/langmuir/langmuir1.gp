set style line 1 lt 1 lw 3 lc rgb "blue"
set style line 2 lt 2 lw 3 lc rgb "red"
set style line 3 lt 3 lw 3 lc rgb "green"
set style line 4 lt 4 lw 3 lc rgb "cyan"

set xrange[0:10]
#set yrange[0:1]
set xlabel 'Time [d]' 
set ylabel 'Concentration [mM]' 
plot 'langmuir1-obs-0.tec' using 1:($4*1000) title 'NH4+' with lines ls 1, \
     'langmuir1-obs-0.tec' using 1:($5/250*1000) title 'NH4+(s)' with lines ls 2, \
     'langmuir1-obs-0.tec' using 1:($4*1000+$5/250*1000) title 'total' with lines ls 3

set terminal postscript eps size 4, 3 enhanced color font 'Helvetica, 20' lw 1 
set output 'langmuir1.eps'

replot
set term pop
pause -1

