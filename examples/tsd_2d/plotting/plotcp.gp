wpx = 3.8
hpx = 3.3

# LATEX
#set terminal epslatex size wpx, hpx standalone color colortext 10 header "\\newcommand{\\ft}[0]{\\huge}"
#set out 'cp.tex'

# PNG
set terminal png medium size 640, 480
set out 'cp.png'

#set multiplot

#set border linewidth 0
#unset key
#unset colorbox
#unset tics
#set lmargin screen 0.05
#set rmargin screen 0.95
#set tmargin screen 0.95
#set bmargin screen 0.05

set xlabel 'x'
set ylabel 'C_p'

plot '../tsd.cp' u 2:3 w lines title 'lower' ,\
     '../tsd.cp' u 2:4 w lines title 'upper'
