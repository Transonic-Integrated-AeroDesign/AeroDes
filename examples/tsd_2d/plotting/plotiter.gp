wpx = 3.8
hpx = 3.3

# LATEX
#set terminal epslatex size wpx, hpx standalone color colortext 10 header "\\newcommand{\\ft}[0]{\\huge}"
#set out 'cp.tex'

# PNG
set terminal png medium size 640, 480
set out 'iter.png'

set multiplot layout 3,1

#set border linewidth 0
#unset key
#unset colorbox
#unset tics
set lmargin screen 0.15
#set rmargin screen 0.95
#set tmargin screen 0.95
#set bmargin screen 0.05
set grid
set key left bottom

set ylabel 'res' offset -5
set format y '10^{%T}'
set log y
plot '../tsd.iter' u 1:abs(2) w lines title 'residual' ,\

set ylabel 'cl'
plot '../tsd.iter' u 1:6 w lines title 'cl' ,\

set xlabel 'iteration'
set ylabel 'cdw'
plot '../tsd.iter' u 1:7 w lines title 'cdw' ,\
