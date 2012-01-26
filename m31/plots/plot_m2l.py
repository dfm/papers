#!/usr/bin/env python

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as pl

import h5py

from hist import reshist2d, colors

# do the plotting

pl.figure(figsize=(6,6))
extent = [[-0.5, 0.9], [-0.7, 0.7]]

# plot the theory (from Larry)
dmin, dmax = -1.0, 2*np.log10(2.6/1.38)-0.001
diskml  = np.linspace(dmin,dmax,10000)      # col1
xdisk   = 1.38/2.6*np.sqrt(10.**diskml)     # col4
bulgev2 = 10.4/2.2/4.8
xbulge2 = bulgev2/2.6**2
bulgeml = np.log10((1-xdisk**2)/xbulge2)    # col3

# Bell
bellc = 'c'
bell, eb = -0.11, 0.1
pl.fill(bell+2*eb*np.array([1,1,-1,-1]), np.append(extent[1], extent[1][::-1]),
        color=bellc, alpha=0.3, zorder=-100)
pl.gca().axvline(bell,                 color=bellc, zorder=-200)
pl.gca().axvline((bell+eb),   ls='--', color=bellc, zorder=-200)
pl.gca().axvline((bell-eb),   ls='--', color=bellc, zorder=-200)
pl.gca().axvline((bell+2*eb), ls=':',  color=bellc, zorder=-200)
pl.gca().axvline((bell-2*eb), ls=':',  color=bellc, zorder=-200)

# excluded region
excol = "m"
xfill = np.append(diskml,  [extent[0][1], extent[0][1], extent[0][0]])
yfill = np.append(bulgeml, [extent[1][0], extent[1][1], extent[1][1]])
pl.fill(xfill, yfill, color=excol, alpha=0.3)
pl.plot(diskml, bulgeml, excol)

# one-to-one
onec = 'k'
pl.plot(extent[0], extent[0], '-', color=onec)
pl.plot(extent[0], np.array(extent[0])+np.log10(2),    ':',  color=onec)
pl.plot(extent[0], np.array(extent[0])-np.log10(2),    ':',  color=onec)
pl.plot(extent[0], np.array(extent[0])+.5*np.log10(2), '--', color=onec)
pl.plot(extent[0], np.array(extent[0])-.5*np.log10(2), '--', color=onec)

# Courteau & Rix
cr = 2*np.log10((0.61 + 0.07*np.random.randn(10000)) * 260./138)
mu, sig = np.mean(cr), np.std(cr)
crc = '#f89406'
pl.fill(mu+sig*np.array([1,1,-1,-1]), np.append(extent[1], extent[1][::-1]),
        color=crc, alpha=0.3, zorder=-300)
# pl.gca().axvline(mu,     color=crc,          zorder=-400)
# pl.gca().axvline(mu+sig, color=crc, ls="--", zorder=-400)
# pl.gca().axvline(mu-sig, color=crc, ls="--", zorder=-400)

# plot the samples
f = h5py.File("m2l.hdf5")

resamp = 1

# naive
burnin = 1750
x,y = np.log10(f["naive"]["m2ld"])[:,burnin::resamp].flatten(),\
        np.log10(f["naive"]["m2lb"])[:,burnin::resamp].flatten()
x = x[y < 0.5]
y = y[y < 0.5]
reshist2d(x, y, extent=extent, color=colors["naive"])

# m2l
burnin = 750
x,y = np.log10(f["m2l"]["m2ld"])[:,burnin::resamp].flatten(),\
        np.log10(f["m2l"]["m2lb"])[:,burnin::resamp].flatten()
reshist2d(x, y, extent=extent, color=colors["m2l"])

f.close()

# Literature
# pl.plot(np.log10(0.93), np.log10(3),   'or', ms=6)
# pl.plot(np.log10(1.2),  np.log10(3.8), 'or', ms=6)

pl.xlabel(r"$\log_{10} \, \Upsilon_\mathrm{d}$", fontsize=16)
pl.ylabel(r"$\log_{10} \, \Upsilon_\mathrm{b}$", fontsize=16)

pl.yticks(rotation=90, family='Times New Roman', fontsize=16.)
pl.xticks(family='Times New Roman', fontsize=16.)

pl.xlim(extent[0])
pl.ylim(extent[1])

pl.savefig("m2l.pdf")

