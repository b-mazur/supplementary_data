#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 21:03:38 2022

@author: bartoszmazur
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import re
from matplotlib.colors import LogNorm, to_rgb, LinearSegmentedColormap



def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

def vtk1Dto3D(filename):
    data1D = np.genfromtxt(filename, skip_header=11)

    # with open(filename) as f:
    #     header = [next(f) for i in range(11)]

    # size_x = int(header[4].split()[1])
    # size_y = int(header[4].split()[2])
    # size_z = int(header[4].split()[3])
    size_x = size_y = size_z = 100
    data3D = np.zeros((size_x, size_y, size_z))

    for x in range(size_x):
        for y in range(size_y):
            for z in range(size_z):
                index = x + y * size_y + z * size_x * size_y
                data3D[x, y, z] = data1D[index]

    data3D = np.concatenate((data3D, data3D[:25, :, :]), axis=0)
    data3D = np.concatenate((data3D, data3D[:, :25, :]), axis=1)
    data3D = np.concatenate((data3D, data3D[:, :, :25]), axis=2)

    data3D = data3D[25:, 25:, 25:]
    
    slabs = np.zeros((size_x, size_y, 5))

    for i in range(5):
        j = i * 5
        k = j + 5
        slab = np.sum(data3D[:, :, j:k], axis=2)
        slabs[:, :, i] = slab
    
    return slabs
    

npy_file_exists = False
list_of_files = []
for i in os.listdir('.'):
    if os.path.isfile(os.path.join('.', i)) and i.endswith('.vtk'):
        list_of_files.append(i)
    elif os.path.isfile(os.path.join('.', i)) and i == 'all_vtk_data.npy':
        npy_file_exists = True
list_of_files = natural_sort(list_of_files)
        
if npy_file_exists == False:
    size_x = size_y = size_z = 100
    all_slabs = np.zeros((size_x, size_y, 5, len(list_of_files)))
    
    for i, file in enumerate(list_of_files):
        all_slabs[:, :, :, i] = vtk1Dto3D(file)
        
    np.save('all_vtk_data.npy', all_slabs)
else:
    all_slabs = np.load('all_vtk_data.npy')
    
# ticktext = np.around(np.linspace(0, 25.832, num=5), 3)
# tickvals = [0, 25, 50, 75, 99]

colors = ['#364B9A', '#4A7BB7', '#6EA6CD', '#98CAE1', '#C2E4EF', '#EAECCC', '#EAECCC', '#FEDA8B', '#FDB366', '#F67E4B', '#DD3D2D', '#A50026']

for i, color in enumerate(colors):
    colors[i] = to_rgb(color)
    
cmap = LinearSegmentedColormap.from_list('custom', colors)

all_slabs[all_slabs == 0] = 1
    
for T in [82, 92, 102, 110]:
    if T == 82:
        zmax = np.max(all_slabs[:, :, :, 0:3])
        log_norm = LogNorm(vmin=1, vmax=zmax)
        fig = plt.figure(figsize=(15, 9))
        ax = fig.add_gridspec(3, 5, wspace=0.02, hspace=0.05).subplots()
        for i, k in enumerate([0, 2, 1]):
            for j in range(5):
                zdata = all_slabs[:, :, j, k]
                ax[i, j].imshow(zdata, cmap=cmap, interpolation='quadric',
                                norm=log_norm)
                ax[i, j].axis('off')
        plt.savefig('fig_log_{}K.png'.format(T), dpi=300, bbox_inches='tight')
        plt.clf()
        
    elif T == 92:
        zmax = np.max(all_slabs[:, :, :, 3:6])
        log_norm = LogNorm(vmin=1, vmax=zmax)
        fig = plt.figure(figsize=(15, 9))
        ax = fig.add_gridspec(3, 5, wspace=0.02, hspace=0.05).subplots()
        for i, k in enumerate([0, 2, 1]):
            for j in range(5):
                zdata = all_slabs[:, :, j, k+3]
                ax[i, j].imshow(zdata, cmap=cmap, interpolation='quadric',
                                norm=log_norm)
                ax[i, j].axis('off')
        plt.savefig('fig_log_{}K.png'.format(T), dpi=300, bbox_inches='tight')
        plt.clf()
        
    elif T == 102:
        zmax = np.max(all_slabs[:, :, :, 6:9])
        log_norm = LogNorm(vmin=1, vmax=zmax)
        fig = plt.figure(figsize=(15, 9))
        ax = fig.add_gridspec(3, 5, wspace=0.02, hspace=0.05).subplots()
        for i, k in enumerate([0, 2, 1]):
            for j in range(5):
                zdata = all_slabs[:, :, j, k+6]
                ax[i, j].imshow(zdata, cmap=cmap, interpolation='quadric',
                                norm=log_norm)
                ax[i, j].axis('off') 
        plt.savefig('fig_log_{}K.png'.format(T), dpi=300, bbox_inches='tight')
        plt.clf()
        
    elif T == 110:
        zmax = np.max(all_slabs[:, :, :, 9:12])
        log_norm = LogNorm(vmin=1, vmax=zmax)
        fig = plt.figure(figsize=(15, 9))
        ax = fig.add_gridspec(3, 5, wspace=0.02, hspace=0.05).subplots()
        for i, k in enumerate([0, 2, 1]):
            for j in range(5):
                zdata = all_slabs[:, :, j, k+9]
                ax[i, j].imshow(zdata, cmap=cmap, interpolation='quadric',
                                norm=log_norm)
                ax[i, j].axis('off')
        plt.savefig('fig_log_{}K.png'.format(T), dpi=300, bbox_inches='tight')
        plt.clf()

# for T in [82, 92, 102, 110]:
#     if T == 82:
#         zmax = np.max(all_slabs[:, :, :, 0:3])
#         fig = plt.figure(figsize=(15, 9))
#         ax = fig.add_gridspec(3, 5, wspace=0.02, hspace=0.05).subplots()
#         for i in range(3):
#             for j in range(5):
#                 zdata = all_slabs[:, :, j, i]
#                 ax[i, j].imshow(zdata, cmap='RdBu_r', interpolation='quadric',
#                                 vmin=1, vmax=zmax)
#                 ax[i, j].axis('off')
#         plt.savefig('fig_{}K.png'.format(T), dpi=300, bbox_inches='tight')
#         plt.clf()
        
#     elif T == 92:
#         zmax = np.max(all_slabs[:, :, :, 3:7])
#         fig = plt.figure(figsize=(15, 12))
#         ax = fig.add_gridspec(4, 5, wspace=0.02, hspace=0.05).subplots()
#         for i in range(4):
#             for j in range(5):
#                 zdata = all_slabs[:, :, j, i+3]
#                 ax[i, j].imshow(zdata, cmap='RdBu_r', interpolation='quadric',
#                                 vmin=1, vmax=zmax)
#                 ax[i, j].axis('off')
#         plt.savefig('fig_{}K.png'.format(T), dpi=300, bbox_inches='tight')
#         plt.clf()
        
#     elif T == 102:
#         zmax = np.max(all_slabs[:, :, :, 7:11])
#         fig = plt.figure(figsize=(15, 12))
#         ax = fig.add_gridspec(4, 5, wspace=0.02, hspace=0.05).subplots()
#         for i in range(4):
#             for j in range(5):
#                 zdata = all_slabs[:, :, j, i+7]
#                 ax[i, j].imshow(zdata, cmap='RdBu_r', interpolation='quadric',
#                                 vmin=1, vmax=zmax)
#                 ax[i, j].axis('off') 
#         plt.savefig('fig_{}K.png'.format(T), dpi=300, bbox_inches='tight')
#         plt.clf()
        
#     elif T == 110:
#         zmax = np.max(all_slabs[:, :, :, 11:15])
#         fig = plt.figure(figsize=(15, 12))
#         ax = fig.add_gridspec(4, 5, wspace=0.02, hspace=0.05).subplots()
#         for i in range(4):
#             for j in range(5):
#                 zdata = all_slabs[:, :, j, i+11]
#                 ax[i, j].imshow(zdata, cmap='RdBu_r', interpolation='quadric',
#                                 vmin=1, vmax=zmax)
#                 ax[i, j].axis('off')
#         plt.savefig('fig_{}K.png'.format(T), dpi=300, bbox_inches='tight')
#         plt.clf()
