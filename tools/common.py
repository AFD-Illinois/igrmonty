from math import pi, sin, cos, floor
from glob import glob

import re
import parse

import numpy  as np
import pandas as pd

from scipy.stats import poisson, norm
import h5py

def filter(self, **kwargs):
    mask = [True] * len(self)
    for k, v in kwargs.items():
        mask &= self[k] == v
    return self[mask]

pd.DataFrame.filter = filter # monkey patch pandas dataframe

def Parabase(fmt, *args, **kwargs):
    pattern = fmt
    for i in range(10):
        #print(i, pattern, args, kwargs)
        try:
            pattern = pattern.format(*args, **kwargs)
            break
        except KeyError as e:
            match   = r'\{' + e.args[0] + ':?.*?\}'
            pattern = re.sub(match, '{}', pattern, 1)
            args    = *args, '*'
    #print(pattern)

    files = sorted(glob(pattern))
    #print(len(files), files[0])

    parser = parse.compile(fmt)
    return pd.DataFrame({'path':f, **parser.parse(f).named} for f in files)

def load_one(file, i=70, di=10):
    ME   = 9.1093897e-28
    CL   = 2.99792458e10
    HPL  = 6.6260755e-27
    LSUN = 3.827e33
    with h5py.File(file, "r") as fp:
        nu    = np.power(10.0, fp["output"]["lnu"])   * ME * CL * CL / HPL
        nuLnu = np.array(      fp["output"]["nuLnu"]) * LSUN
        # Note that nuLnu dimension is [N_TYPEBINS, N_EBINS, N_THBINS]
        # The theta bins range from 0 to 90 deg.

    nth = nuLnu.shape[-1]
    dth = pi / 2 / nth
    mth = list(range(nth))
    mth = mth + mth[::-1]

    ji = (i - di/2) * nth / 90
    jf = (i + di/2) * nth / 90

    W, T = 0, 0
    for j in range(floor(ji), floor(jf)+1):
        w  = abs(+ cos(max(j,  ji ) * dth)
                 - cos(min(jf, j+1) * dth))
        W += w
        T += w * nuLnu[:,:,mth[j]]
    return nu, T.T / W

def load(files, **kwargs):
    nuLnu = []
    for f in files:
        data   = load_one(f, **kwargs)
        nuLnu += [np.column_stack((np.sum(data[1], axis=-1), data[1]))]
        try:
            if all(nu != data[0]):
                raise ValueError('Frequency bins `nu` do not match')
        except NameError:
            nu = data[0]

    nuLnu = np.array(nuLnu)
    return (nu,
            np.mean(nuLnu, axis=0),
            np.std (nuLnu, axis=0))

def interval(avg, std, sigma=1):
    """We model the different realizations of grmonty results follow a
    Poisson distribution, which has mean and variance both equal to
    `mu`.  Therefore,

        avg = dnuLnu * mu
        std = dnuLnu * sqrt(mu)

        mu = (avg/std)**2

    One we obtain `mu`, we can estimate the lower and upper intervals
    according to `sigma`.

    """
    with np.errstate(invalid='ignore', divide='ignore'):
        mu = (avg/std)**2

    lower = poisson.ppf(norm.cdf(-sigma), mu)
    upper = poisson.ppf(norm.cdf(+sigma), mu)
    units = avg / mu

    return lower * units, upper * units

def step_one(ax, nu, avg, std=None, sigma=1,
             shade=True, ylog=True, **kwargs):
    p = ax.step(nu, avg, where='mid', **kwargs)

    if std is not None and shade:
        l, u = interval(avg, std, sigma=sigma)
        ax.fill_between(nu, l, u, step='mid',
                        color=p[0].get_color(), alpha=1/3, linewidth=0)

    # x-axis must be in log scale; otherwise the bin boundaries are wrong
    ax.set_xscale('log')
    # optionally we may set y-axis to log sacle
    if ylog:
        ax.set_yscale('log')

def step(ax, nu, avg, std=None, color=None, shade=None, **kwargs):
    n   = len(nu)
    avg = avg.reshape(n,-1)
    for i in range(avg.shape[-1]):
        stdi   = std.reshape(n,-1)[:,i] if std is not None else None
        colori = color if color is not None else ('k'  if i == 0 else None)
        widthi = 2 if i == 0 else 1
        shadei = shade if shade is not None else (True if i == 0 else False)
        step_one(ax, nu, avg[:,i], stdi,
                 color=colori, linewidth=widthi,
                 shade=shadei, **kwargs)
