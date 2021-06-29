set terminal postscript color
set output "MINI.ps"
set palette defined (0 0 0 0.9, 1 0 0 1, 2 0 0.5 1, 3 0 1 1, 4 0.5 1 0.5, 5 1 1 0, 6 1 0.5 0, 7 1 0 0, 8 0.9 0 0)
set lmargin at screen 0.1
set rmargin at screen 0.88
# BUILD FROM BOTTOM TO TOP #
TOP=0.92  #
DY = 0.24 # DIMENSION OF Y
set multiplot
set tics font ' ,20'
set colorbox user origin screen 0.88, screen 0.2 size screen 0.03, screen 0.72
set grid front ls 8
set offset 0, 0, graph 0.05, graph 0.05 
set cbtics ('1' 1, '0' 0, '-1' -1)
set cblabel ('') font ' ,30' offset 3,0,0
set yrange [ -0.01 : 0.01 ]
set xrange [ -1.5 : 0]
set ylabel "{/: E/t }" font 'Linux Biolinum O,25'  offset 0,3,0
set xlabel "{/Symbol m}"  font 'Linux Biolinum O,25' offset 0,-2,0
set xtics (-1.30, -1.20,-0.81,-0.70, -0.58,-0.48, -0.08,-0.02, ' ' 0.05) rotate by 45 right
set ytics (-0.01, 0 , 0.01) 
set tmargin at screen TOP-1.5*DY
set bmargin at screen TOP-3*DY
plot'/home/bento/FINITO200.dat' u 1:2:3  w points pointtype 7 ps 1 lw 2 palette t " "
unset colorbox
unset xtics
set xtics format " "
unset ytics
unset x2tics
unset ylabel
unset xlabel
set yrange [ -0.25 : 0.25]
set xrange [ -1.5 : 0 ]
set grid front ls 8
set tics font ' ,20'
set ytics (-0.25, 0, 0.25 ) 
set xtics ( 0 ) 
set x2tics (-1.5,0,1.5)
set tmargin at screen TOP
set bmargin at screen TOP-1.5*DY
plot'/home/bento/FINITO200.dat' u 1:2:3  w points pointtype 7 ps 0.5 lw 2 palette t " "
unset multiplot
