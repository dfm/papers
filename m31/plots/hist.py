import numpy as np
import scipy.special as sp

import matplotlib.pyplot as pl
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.cm as cm

colors = {"m2l": "#0064CD", "naive": "#46a546"}

def reshist2d(x, y, *args, **kwargs):
    ax = kwargs.pop('ax', pl.gca())

    extent = kwargs.pop('extent', [[x.min(), x.max()], [y.min(), y.max()]])
    bins   = kwargs.pop('bins', 100)
    color  = kwargs.pop('color', 'k')
    alpha  = kwargs.pop('alpha', 0.1)

    cmap = cm.get_cmap('gray')
    cmap._init()
    cmap._lut[:-3, -1] = np.linspace(1,0,cmap.N) # alpha
    cmap._lut[:-3,:-1] = 0.

    X = np.linspace(extent[0][0], extent[0][1], bins+1)
    Y = np.linspace(extent[1][0], extent[1][1], bins+1)
    H, X, Y = np.histogram2d(x.flatten(), y.flatten(), bins=(X,Y))

    V = sp.erf(np.arange(0.5, 2.1, 0.5)/np.sqrt(2))
    Hflat = H.flatten()
    inds = np.argsort(Hflat)[::-1]
    Hflat = Hflat[inds]
    sm = np.cumsum(Hflat)
    sm /= sm[-1]

    for i, v0 in enumerate(V):
        try:
            V[i] = Hflat[sm <= v0][-1]
        except:
            V[i] = Hflat[0]

    X, Y = 0.5*(X[1:]+X[:-1]), 0.5*(Y[1:]+Y[:-1])

    ax.plot(x, y, 'o', color=color, ms=1, zorder=-1, alpha=alpha,
            rasterized=True)
    ax.contourf(X, Y, H.T, [V[-1], 0.],
        cmap=LinearSegmentedColormap.from_list('cmap',([1]*3,[1]*3),N=2))
    ax.pcolor(X, Y, H.max()-H.T, cmap=cmap)
    ax.contour(X, Y, H.T, V, colors=color)

