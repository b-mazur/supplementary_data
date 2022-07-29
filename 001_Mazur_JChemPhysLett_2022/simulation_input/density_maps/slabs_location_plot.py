#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 12:03:16 2022

@author: bartoszmazur
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, to_rgb, LinearSegmentedColormap


all_slabs = np.load('all_vtk_data.npy')


colors = ['#364B9A', '#4A7BB7', '#6EA6CD', '#98CAE1', '#C2E4EF', '#EAECCC', '#EAECCC', '#FEDA8B', '#FDB366', '#F67E4B', '#DD3D2D', '#A50026']

for i, color in enumerate(colors):
    colors[i] = to_rgb(color)
    
cmap = LinearSegmentedColormap.from_list('custom', colors)

all_slabs[all_slabs == 0] = 1

T = 92
    
zmax = np.max(all_slabs[:, :, :, 3:6])
log_norm = LogNorm(vmin=1, vmax=zmax)

fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(projection='3d')

X, Y = np.meshgrid(np.linspace(0,1,101), np.linspace(0,1,101))


data = np.random.random((100, 100))

zmax = np.max(all_slabs[:, :, :, 3:6])
log_norm = LogNorm(vmin=1, vmax=zmax)
m_r = plt.cm.ScalarMappable(norm=log_norm, cmap=cmap)
m_r.set_array([])

for j in range(5):
    Z = (0.05 * j + 0.025) * np.ones(X.shape)
    zdata = all_slabs[:, :, j, 5]
    fcolors = m_r.to_rgba(zdata)
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=fcolors, shade=False)
    
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_zlim(0, 0.25)
ax.set_xticks([0, 0.5, 1], ['0', '12.916', '25.832'], ha='left', va='center', rotation=0)
ax.set_yticks([0, 0.5, 1], ['0', '12.916', '25.832'], ha='right', va='center', rotation=0)
ax.set_zticks([0, 0.05, 0.1, 0.15, 0.2, 0.25], ['0', '1.292', '2.583', '3.875', '5.166', '6.458'])
ax.tick_params(axis='x', width=10, labelsize=10, pad=-3)
ax.tick_params(axis='y', width=10, labelsize=10, pad=-3)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.azim = 60
ax.dist = 10
ax.elev = 15
ax.set_box_aspect([1, 1, 1])


plt.savefig('fig_slabs_{}K.png'.format(T), dpi=300, bbox_inches='tight')