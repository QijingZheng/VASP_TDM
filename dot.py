#!/usr/bin/env python
"""
This is an auxiliary script, which takes as input the calculated transition
dipole moment (TDM) from 'vasptdm' and outputs the folded oscillator strength
(FOS) in a file called 'DOT.dat'. In addition, figure of FOS is plotted if
either python-matplotlib module or gnuplot is available in the system.
"""

import numpy as np

inp = np.loadtxt('TDM.dat', usecols=(8, 10, 11, 12, 13))

############################################################
hartree = 27.211395655517308
ext = 0.2                   # extra energy range to plot
Np  = 1000                  # Number of points for plotting
sigma = 0.02                # smearing parameter
############################################################
x = np.abs(inp[:,0])
# oscillator strength
y = inp[:,1:]
for i in range(4):
    y[:,i] = y[:,i] * x * 2./3 * hartree

xmin = x.min()
xmax = x.max()
xdif = xmax - xmin

xmin = xmin - xdif * ext
xmax = xmax + xdif * ext

x0 = np.linspace(xmin, xmax, Np)
y0 = np.zeros((Np, 4))

for i in range(x.size):
    en = x[i]
    for j in range(4):
        mag = y[i,j]
        # Lorentzian Function
        # y0[:,j] += mag * sigma / np.pi / ((x0 - en)**2 + sigma**2)
        # Gaussians
        y0[:,j] += mag / (np.sqrt(2*np.pi) * sigma) * np.exp(-(x0-en)**2/2/sigma**2)

np.savetxt('DOT.dat', np.vstack((x0, y0.T)).T, fmt="%12.4f")

############################################################
# plotting
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    mpl.rcParams['axes.linewidth'] = 1.5 
    mpl.rcParams['xtick.labelsize'] = 12
    colors = ['red', 'olive', 'blue', 'black']
    legends =  ['X', 'Y', 'Z', 'TOT']
    fig, ax = plt.subplots()
    # plotting part
    for i in range(4):
        ax.plot(x0, y0[:,i], lw=1.5, color=colors[i], linestyle='-')
    # add legend
    ax.legend(ax.get_lines(), legends)

    if xmin < 0: ax.set_xlim((0, xmax))
    ax.set_xlabel("Photon Energy [eV]",
                   fontsize=15, 
                   labelpad=10)
    ax.set_ylabel("FOS [1/eV]",
                   fontsize=15,
                   labelpad=10)
    ax.set_yticklabels([])
    [tick.label.set_fontweight('bold') for tick in ax.xaxis.get_major_ticks()]
    fig.set_size_inches((8,6))
    plt.savefig('kaka.png', dpi=360)
    plt.show()
except: 
    from subprocess import Popen, PIPE
    from StringIO import StringIO

    GNUPLOT = """
    #
    set terminal png font ",32" size 800,600
    set output "kaka.png"
    # unset key
    set xlabel "Photon Energy [eV]"
    set ylabel "FOS [1/eV]"
    set format y ""
    set border 15 lw 2.0
    #
    plot "DOT.dat" u 1:2 title "X" w line lw 2.0 linecolor rgb 'red',\
                "" u 1:3 title "Y" w line lw 2.0 linecolor rgb 'cyan', \
                "" u 1:4 title "Z" w line lw 2.0 linecolor rgb 'blue', \
                "" u 1:5 title "TOT" w line lw 2.0 linecolor rgb 'black'
    quit
    #
    """

    PLOT = Popen(['gnuplot', ], stdin=PIPE)
    PLOT.communicate(input=GNUPLOT)
