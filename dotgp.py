#!/usr/bin/env python

from subprocess import Popen, PIPE
from StringIO import StringIO

GNUPLOT = """
############################################################
set terminal png font "Times New Roman,42" size 800,600
set output "kaka.png"
unset key
set xlabel "Energy [eV]"
set ylabel "Oscillator Strength [arb.units]"
set format y ""
set border 15 lw 2.0
############################################################
plot "DOT.dat" using 1:2 with line linewidth 2.0 \
               linecolor rgb 'red'
quit
############################################################
"""

PLOT = Popen(['gnuplot', ], stdin=PIPE)
PLOT.communicate(input=GNUPLOT)

