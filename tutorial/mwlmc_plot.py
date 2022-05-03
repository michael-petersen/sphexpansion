import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def plot_density(Ym, Zm, dens, **kwargs):
    fig, ax = plt.subplots()
    pc = ax.pcolormesh(Ym, Zm, 
                       dens * 1e9, 
                       shading = 'nearest', 
                       rasterized = True, 
                       **kwargs)


    cbar = fig.colorbar(pc, ax=ax)

    cbar.set_label(r'$\rho\ [\mathrm{M}_\odot \mathrm{kpc}^{-3}]$')
    ax.set_xlabel('y [kpc]')
    ax.set_ylabel('z [kpc]')
    ax.axhline(0, c = 'grey', ls = '--', lw = 1, zorder = 100)
    ax.axvline(0, c = 'grey', ls = '--', lw = 1, zorder = 100)
    ax.set_aspect('equal')
    return fig, ax
  
def plot_potential(Ym, Zm, pot, **kwargs):
    fig, ax = plt.subplots()
    pc = ax.pcolormesh(Ym, Zm, 
                       pot, 
                       shading = 'nearest',
                       rasterized = True, 
                       **kwargs)


    cbar = fig.colorbar(pc, ax=ax)

    cbar.set_label(r'$\Phi\ [\mathrm{km}^2 \mathrm{s}^{-2}]$')
    ax.set_xlabel('y [kpc]')
    ax.set_ylabel('z [kpc]')
    ax.axhline(0, c = 'grey', ls = '--', lw = 1, zorder = 100)
    ax.axvline(0, c = 'grey', ls = '--', lw = 1, zorder = 100)
    ax.set_aspect('equal')
    return fig, ax