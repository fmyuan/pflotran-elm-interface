set terminal postscript eps size 8, 9 enhanced color font 'Helvetica, 20' lw 1 
set output 'clmcn3.eps'
set style line 1 lt 1 lw 3 lc rgb "blue"
set style line 2 lt 1 lw 3 lc rgb "red"
set style line 3 lt 1 lw 3 lc rgb "green"
set style line 4 lt 1 lw 3 lc rgb "violetn"

set multiplot layout 3, 2 rowsfirst

set label 1 '(a)' at graph 0.05, graph 0.9 
set xrange[0:400]
set yrange[0:0.5]
set xlabel 'Time [d]' 
set ylabel 'Concentration [mol/m^3]' 
plot \
     'clmcn3-obs-0.tec' using 1:($10) title 'Lit1' with lines ls 1, \
     'clmcn3-obs-0.tec' using 1:($11) title 'Lit2' with lines ls 2, \
     'clmcn3-obs-0.tec' using 1:($12) title 'Li3' with lines ls 3

set label 1 '(b)' at graph 0.05, graph 0.9 
set yrange[0:1.5]
plot \
     'clmcn3-obs-0.tec' using 1:($6) title 'SOM1' with lines ls 1, \
     'clmcn3-obs-0.tec' using 1:($7) title 'SOM2' with lines ls 2, \
     'clmcn3-obs-0.tec' using 1:($8) title 'SOM3' with lines ls 3, \
     'clmcn3-obs-0.tec' using 1:($9) title 'SOM4' with lines ls 4

set label 1 '(c)' at graph 0.05, graph 0.9 
set yrange[0:2]
plot \
     'clmcn3-obs-0.tec' using 1:($2*250+$6+$7+$8+$9+$10+$11+$12) title 'C total' with lines ls 1, \
     'clmcn3-obs-0.tec' using 1:($6+$7+$8+$9+$10+$11+$12) title 'C left' with lines ls 2, \
     'clmcn3-obs-0.tec' using 1:($2*250) title 'CO2' with lines ls 3

set label 1 '(d)' at graph 0.05, graph 0.9 
set yrange[0:0.05]
set ylabel 'Concentration [mol/m^3]' 
plot \
     'clmcn3-obs-0.tec' using 1:($13+$14+$15+$6*0.071429+$7*0.071429+$8*0.085714+$9*0.085714+$3*250+$4*250+$5*250) title 'N total' with lines ls 1, \
     'clmcn3-obs-0.tec' using 1:($13+$14+$15+$6*0.071429+$7*0.071429+$8*0.085714+$9*0.085714) title 'N left' with lines ls 2, \
     'clmcn3-obs-0.tec' using 1:($3*250) title 'NH4+' with lines ls 3, \
     'clmcn3-obs-0.tec' using 1:($4*250) title 'NO3-' with lines ls 4, \
     'clmcn3-obs-0.tec' using 1:($5*250) title 'N2O' with lines ls 5

set label 1 '(e)' at graph 0.05, graph 0.9 
#set yrange[0:0.1]
set ylabel 'Concentration [mol/m^3]' 
plot \
     'clmcn3-obs-0.tec' using 1:($3*250) title 'NH4+' with lines ls 3, \
     'clmcn3-obs-0.tec' using 1:($4*250) title 'NO3-' with lines ls 4

set label 1 '(f)' at graph 0.05, graph 0.9 
#set yrange[0:0.1]
set ylabel 'Concentration [mmol/m^3]' 
plot \
     'clmcn3-obs-0.tec' using 1:($5*1e3*250) title 'N2O' with lines ls 5

unset multiplot

