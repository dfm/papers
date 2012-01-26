#!/usr/bin/env python

import numpy as np

import matplotlib
matplotlib.use("Agg")
# matplotlib.rc('text', usetex=True)
# matplotlib.rc('font',**{'family':'serif', 'serif':['TimesNewRoman']})
import matplotlib.pyplot as pl

import h5py

from hist import reshist2d

# do the plotting

pl.figure(figsize=(6,6))
extent = [[-500, 500], [-500, 500]]

# plot the samples
f = h5py.File("vt.hdf5")

resamp = 1

# m2l
burnin = 750
x,y = 100*f["m2l"]["vw"][:,burnin::resamp].flatten(),\
        100*f["m2l"]["vn"][:,burnin::resamp].flatten()
f.close()
reshist2d(x, y, extent=extent, color='k', bins=25)

# plot loeb
p0 = np.array([-86, -111])
v1 = np.array([ 128., 130.])
v2 = np.array([-91.20855289, 89.80534438])
loeb = np.array([p0, p0+v1, p0+v2+v1, p0+v2, p0])

pl.plot(loeb[:,0], loeb[:,1], 'c')

# plot vdm
vdm = np.array([line.split()[1:] for line in open("vdm.dat")], dtype=float)
pl.errorbar(vdm[:-1,0], vdm[:-1,2], xerr=vdm[:-1,1], yerr=vdm[:-1,3], zorder=20000,
        fmt="+m", mec="m", lw=1.5, capsize=2, barsabove=True)
pl.errorbar(vdm[-1,0], vdm[-1,2], xerr=vdm[-1,1], yerr=vdm[-1,3], zorder=19000,
        fmt="or", mec="r", lw=2., capsize=0, barsabove=True)

# grid
pl.gca().axhline(0, color='k', alpha=0.5)
pl.gca().axvline(0, color='k', alpha=0.5)

# cov eigenvalues
mux, muy = np.mean(x), np.mean(y)
print "mean_vw: ", mux, " mean_vn: ", muy
c = np.cov(x,y)
print "cov: "
print c

eig = np.linalg.eig(c)
print "eig-vals of cov: ", np.sqrt(eig[0])
ex = np.linspace(extent[0][0], extent[0][1], 100)
pl.plot(ex, eig[1][1,1]/eig[1][0,1]*(ex-mux)+muy, '--g')
pl.plot(ex, eig[1][1,0]/eig[1][0,0]*(ex-mux)+muy, '--g')

pl.xlabel(r"$v_W \, (\mathrm{km/s})$", fontsize=16)
pl.ylabel(r"$v_N \, (\mathrm{km/s})$", fontsize=16)

pl.yticks(rotation=90, family='Times New Roman', fontsize=16.)
pl.xticks(family='Times New Roman', fontsize=16.)

pl.xlim(extent[0])
pl.ylim(extent[1])

pl.savefig("vt.pdf")

