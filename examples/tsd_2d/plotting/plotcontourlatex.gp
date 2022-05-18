wpx = 4.5
hpx = 2.5

# PNG
#set terminal epslatex size wpx, hpx standalone color colortext 10 header "\\newcommand{\\ft}[0]{\\huge}"
set terminal epslatex size wpx, hpx standalone
set out 'contour.tex'

set multiplot

set pm3d map
set pm3d interpolate 0,0
set contour
set cntrparam levels 29

unset clabel
unset border
unset key
unset colorbox
unset tics
set lmargin screen 0.05
set rmargin screen 0.97
set tmargin screen 0.95
set bmargin screen 0.05

set xrange [-1:2]
set yrange [-1.6:1.9]
load 'jet.pal'

splot '../tsd.cpcon' u 1:2:3 ,\
#      '../tsd.cpcon' u 1:2:3 w lines lc 1

set xrange [-1:2]
set yrange [-1.6:1.9]
plot '../tsd.xzmses' u 1:2 w filledcu lc 'black'

plot '../tsd.xzmses' u 1:3 w filledcu lc 'black'
