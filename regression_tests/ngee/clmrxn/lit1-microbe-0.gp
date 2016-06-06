set terminal postscript eps size 8, 3 enhanced color font 'Helvetica, 20' lw 1 
set output 'lit1-microbe-0.eps'
set style line 1 lt 1 lw 3 lc rgb "blue"
set style line 2 lt 2 lw 3 lc rgb "red"
set style line 3 lt 3 lw 3 lc rgb "cyan"
set style line 4 lt 4 lw 3 lc rgb "magenta"

set multiplot layout 1, 2 rowsfirst

set xrange[0:10.0]
set yrange[0:0.3]
set xlabel 'Time [d]' 
set ylabel 'Concentration (mol/m3)' 

plot \
     'lit1-microbe-0-obs-0.tec' using 1:($5) title 'Lit1C' with lines ls 1 axes x1y1, \
     'lit1-microbe-0-obs-0.tec' using 1:($7) title 'Bacteria' with lines ls 2 axes x1y1, \
     'lit1-microbe-0-obs-0.tec' using 1:($8) title 'Fungi' with lines ls 3 axes x1y1, \
     'lit1-microbe-0-obs-0.tec' using 1:($4*250) title 'CO2' with lines ls 4 axes x1y1

set yrange[0:10]
set ylabel 'Concentration (mmol/m3)' 
plot \
     'lit1-microbe-0-obs-0.tec' using 1:($2*1e3*250) title 'NH4+ ' with lines ls 1 axes x1y1, \
     'lit1-microbe-0-obs-0.tec' using 1:($3*1e3*250) title 'NO3- ' with lines ls 2 axes x1y1, \
     'lit1-microbe-0-obs-0.tec' using 1:($6*1e3) title 'Lit1N' with lines ls 3 axes x1y1
unset multiplot

