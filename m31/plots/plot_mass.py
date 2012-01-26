#!/usr/bin/env python

import cPickle as pickle

import numpy as np
import numpy.ma as ma

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as pl

import h5py

from hist import reshist2d, colors

LOGPLOT = True

# mass profile
pl.figure(figsize=(7, 4.))
pl.axes([0.1, 0.15, 0.85, 0.8])

mass = {}
for bp in ["m2l", "naive"]:
    massrc = h5py.File("%s/massrc.hdf5"%bp)
    rad = massrc['massrad'][...]
    if LOGPLOT:
        rad = np.log10(rad)
    massrc.close()
    pkl_file = open("%s/cache.pkl"%bp, 'rb')
    (mass[bp], rc, sbp) = pickle.load(pkl_file)
    pkl_file.close()

    m = mass[bp][:,:,:]
    if LOGPLOT:
        m = np.log10(ma.masked_array(m, mask=(m==0)))

    mu  = np.array([np.mean(m[:,:,i]) for i in range(m.shape[-1])])
    std = np.array([np.std(m[:,:,i])  for i in range(m.shape[-1])])

    pl.plot(rad, mu, color=colors[bp], zorder=-1)
    # pl.fill(np.append(rad,rad[::-1]),
    #         np.append(mu+2*std,(mu-2*std)[::-1]),
    #         color=colors[bp], alpha=0.1, zorder=-5)
    pl.fill(np.append(rad,rad[::-1]),
            np.append(mu+std,(mu-std)[::-1]),
            color=colors[bp], alpha=0.3, zorder=-5)

# EW
def ewmass(r,M = 12.5,a = 93.):
    return r*M/np.sqrt(r**2.+a**2.)

rs2 = 10.**np.linspace(0., 3.0, 1000)
xfill = np.append(rs2, rs2[::-1])
yfill = np.append(ewmass(rs2,M=12.3+18), ewmass(rs2[::-1],M=12.3-6))
if LOGPLOT:
    xfill, yfill = np.log10(xfill), np.log10(yfill)
pl.fill(xfill, yfill, color='k', alpha=0.1, zorder=-100)

# literature
lit = np.array([line.split() for line in open('mass.dat')], dtype=float)
c = "#f89406"
D = 785.
M = (1+(D - lit[:,0])/D) * lit[:,2]
R = lit[:,1] * D/lit[:,0]
if LOGPLOT:
    R = np.log10(R)
    M = np.log10(M)
    dM = np.zeros((2, len(M)))
    dM[1,:] = np.log10((1+(D - lit[:,0])/D) * (lit[:,2] + lit[:,3])) - M
    dM[0,:] = M - np.log10((1+(D - lit[:,0])/D) * (lit[:,2] - lit[:,3]))
else:
    dM = (1+(D - lit[:,0])/D) * lit[:,3]
pl.errorbar(R, M, yerr=dM, fmt="o", ms=4.5, color=c, mec=c,
        barsabove=True, capsize=0, elinewidth=1.5)

# Klypin
c_k = "#7a43b6"
kx = [3.5, 100, 300, 290]
ky = [0.5*(0.28+0.25+0.36+0.18), 0.5*(8.5+8.1), 16, 14.3]
if LOGPLOT:
    kx, ky = np.log10(kx), np.log10(ky)
pl.plot(kx, ky, 'o', color=c_k)

if LOGPLOT:
    pl.xlim([0.5,np.log10(550)])
    pl.ylim([-0.49, 1.5])
    pl.ylabel(r'$\log_{10}\,M (<R)/\mathrm{10^{11} \, M_\odot}$',fontsize=16.)
    pl.xlabel(r'$\log_{10} \, R/\mathrm{kpc}$',fontsize=16.)
else:
    pl.xlim([0, 310])
    pl.ylim([0, 25])
    pl.ylabel(r'$M (<r) \, (\mathrm{10^{11} \, M_\odot})$',fontsize=16.)
    pl.xlabel(r'$r \, (\mathrm{kpc})$',fontsize=16.)

pl.yticks(family='Times New Roman', fontsize=16.)
pl.xticks(family='Times New Roman', fontsize=16.)

pl.savefig("mass-profile.pdf")

# M300 vs Md

pl.figure(figsize=(6,6))
if LOGPLOT:
    i = np.arange(len(rad))[rad < np.log10(300.)][-1]
else:
    i = np.arange(len(rad))[rad < 300.][-1]

extent = [[1.1, 4.9], [2.1, 12]]
burnin = {"m2l": 750, "naive": 1750}
for bp in ["m2l", "naive"]:
    f = h5py.File("%s/chain.hdf5"%bp)
    x = mass[bp][:,:,i].flatten()/10.
    y = 2.325e-1*10**f['chain'][...][:, 6, burnin[bp]::10].flatten()
    f.close()
    reshist2d(x, y, extent=extent, color=colors[bp], bins=35, alpha=0.3)

# klypin
pl.plot([1.6, 1.43], [7, 9], 'o', color=c_k)

pl.xlabel(r'$M (<300\,\mathrm{kpc}) \, (\mathrm{10^{12} \, M_\odot})$',fontsize=16.)
pl.ylabel(r'$M_d \, (\mathrm{10^{10} \, M_\odot})$',fontsize=16.)

pl.yticks(family='Times New Roman', fontsize=16.)
pl.xticks(family='Times New Roman', fontsize=16.)

pl.xlim(extent[0])
pl.ylim(extent[1])

pl.savefig("m300.pdf")

# histogram

pl.figure(figsize=(6,6))
n, b = {}, {}
for bp in ["m2l", "naive"]:
    x = mass[bp][:,:,i].flatten()/10.
    n[bp], b[bp] = np.histogram(x, 50, range=extent[0])
norm = np.sum(n["m2l"])/np.sum(n["naive"])
pl.plot(0.5*(b["m2l"][1:]+b["m2l"][:-1]), n["m2l"], 'o', color=colors["m2l"])
pl.plot(0.5*(b["naive"][1:]+b["naive"][:-1]), norm*n["naive"], 'o', color=colors["naive"])

pl.xlabel(r'$M (<300\,\mathrm{kpc}) \, (\mathrm{10^{12} \, M_\odot})$',fontsize=16.)
pl.gca().set_yticklabels([])

pl.savefig("m300-hist.pdf")

