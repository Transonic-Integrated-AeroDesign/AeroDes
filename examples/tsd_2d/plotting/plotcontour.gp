wpx = 3.8
hpx = 3.3

# LATEX
#set terminal epslatex size wpx, hpx standalone color colortext 10 header "\\newcommand{\\ft}[0]{\\huge}"
#set out 'contour.tex'

# PNG
set terminal png medium size 640, 480
set out 'contour.png'

set multiplot

set pm3d map
set pm3d interpolate 0,0
set contour
set cntrparam levels 25

set border linewidth 0
unset key
unset colorbox
#unset tics
set lmargin screen 0.05
set rmargin screen 0.95
set tmargin screen 0.95
set bmargin screen 0.05

splot '../tsd.cpcon' u 1:2:3

#set size 0.5, 0.5
#set origin 50, 0.5
#set yrange [-1:1.5]
#set xrange [-2:2]
#splot '../tsd.cpcon' u 1:2:3
