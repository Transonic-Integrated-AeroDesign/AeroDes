wpx = 1280
hpx = 640

# PNG
set terminal png size wpx, hpx
set out 'contourgrid.png'

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
set lmargin screen 0.01
set rmargin screen 0.99
set tmargin screen 0.99
set bmargin screen 0.01

#set xrange [-0.1:0.1]
set xrange [-1:2]
set yrange [-1:2]
#load 'inferno.pal'
load 'jet.pal'

splot '../tsd.cpcon' u 1:2:3 ,\
      '../tsd.cpcon' u 1:2:3 w lines lc 'black'

set xrange [-1:2]
set yrange [-1:2]
plot '../tsd.xzmses' u 1:2 w filledcu lc 'black' ,\
     '../tsd.xzmses' u 1:3 w filledcu lc 'black'
