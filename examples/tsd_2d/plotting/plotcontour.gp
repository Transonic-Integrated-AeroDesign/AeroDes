wpx = 3.8
hpx = 3.3

# PNG
set terminal png medium size 640, 480
set out 'contour.png'

set style line 1 lc 'black'
set multiplot

set pm3d map
set pm3d interpolate 0,0
set contour
set cntrparam levels 25

unset clabel
unset border
unset key
unset colorbox
unset tics
set lmargin screen 0.05
set rmargin screen 0.95
set tmargin screen 0.95
set bmargin screen 0.05

#set xrange [-0.1:0.1]
set xrange [-1:2]
set yrange [-1:2]
#set palette gray
#set palette magma
load 'inferno.pal'

splot '../tsd.cpcon' u 1:2:3

set xrange [-1:2]
set yrange [-1:2]
plot '../tsd.xzmses' u 1:2 w filledcu lc 'black' ,\
     '../tsd.xzmses' u 1:3 w filledcu lc 'black'
