#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 11:50:30 2022

@author: bartoszmazur
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.signal import argrelextrema
import sys


def check_for_files(path='.'):
    """Checks for files with transiton probabilites and with pressure values 
    for isotherm calculation"""
    probability_file_exisits = False
    pressure_file_exists = False
    prob_files = []
    
    files = (file for file in os.listdir(path) 
             if os.path.isfile(os.path.join(path, file)))
    for file in files:
        if file.startswith('tmmc-res_'):
            probability_file_exisits = True
            prob_files.append(file)
        elif file == 'pressures.dat':
            pressure_file_exists = True
    
    if probability_file_exisits is False:
        print("File with transition probabilities not found.")
        print("The file should be named " 
              "\"tmmc-res_STRUCTURE_UNITCELLS_TEMPERATURE_INITPRESSURE.dat\".")
        print("Exiting...")
        raise SystemExit 
    
    # checks if file with pressure values for isotherm is in the directiory
    if pressure_file_exists is False: 
        print("File \"pressures.dat\" not found.\n")
        pressures = None
    else:
        print("File \"pressures.dat\" found.\n")
        pressures = np.genfromtxt('pressures.dat').flatten()
    
    return prob_files, pressures 

def construct_Pmatrix(data_file, file_suffix):
    """Constructs a P-matrix (matrix of probabilites of changing state) from 
    a file containing transition acceptance rates"""
    # read file with transition probabilites data
    # the file layout should be:
    # N   C(N->N+1)   C(N->N-1)   C(N->N)
    data = np.genfromtxt(data_file)
    molecules = np.array(data[:, 0], copy=True, dtype=int).flatten()
    
    # maximal simulated number of molecules
    N_max = molecules[-1]
    C_matrix = np.zeros((N_max + 1, N_max + 2))
    
    # update C-matrix values
    for j, i in enumerate(molecules):
        if j == 0:
            C_matrix[i, i+1] = data[j, 1]
            C_matrix[i, i] = data[j, 3]
        else:
            C_matrix[i, i+1] = data[j, 1]
            C_matrix[i, i-1] = data[j, 2]
            C_matrix[i, i] = data[j, 3]
    
    P_matrix_sim = np.hstack((molecules[:, None], 
                              np.zeros((molecules.size, 3))))
    
    for row, N in enumerate(molecules):
        # P(N->N+1) insertion move
        P_ins = C_matrix[N, N+1] / np.sum((C_matrix[N, N-1], C_matrix[N, N], 
                                           C_matrix[N, N+1]))
        P_matrix_sim[row, 1] = P_ins
        # P(N->N-1) deletion move
        P_del = C_matrix[N, N-1] / np.sum((C_matrix[N, N-1], C_matrix[N, N], 
                                           C_matrix[N, N+1]))
        P_matrix_sim[row, 2] = P_del
        # P(N->N) move rejection
        P_rej = C_matrix[N, N] / np.sum((C_matrix[N, N-1], C_matrix[N, N], 
                                         C_matrix[N, N+1]))
        P_matrix_sim[row, 3] = P_rej
    
    # plots P(N) for visual inspection
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.5), layout='constrained')
    axs[0].plot(P_matrix_sim[:, 0], P_matrix_sim[:, 1], '.', color='tab:blue')
    axs[1].plot(P_matrix_sim[:, 0], P_matrix_sim[:, 2], '.', color='tab:blue')
    axs[2].plot(P_matrix_sim[:, 0], P_matrix_sim[:, 3], '.', color='tab:blue')

    axs[0].set_ylabel('$\mathregular{P_{N, N+1}}$', fontsize=13)
    axs[1].set_ylabel('$\mathregular{P_{N, N-1}}$', fontsize=13)
    axs[2].set_ylabel('$\mathregular{P_{N, N}}$', fontsize=13)
    
    fig.savefig('P_simulated{}.png'.format(file_suffix), dpi=300)
    axs[0].set_yscale('log')
    axs[1].set_yscale('log')
    fig.savefig('P_simulated_log_{}.png'.format(file_suffix), dpi=300)

    for ax in axs:
        ax.tick_params(which='both', direction='in', labelsize=13, 
                       top="true", right="true")
        ax.grid(color='grey', linestyle=':', linewidth=0.5)
        ax.set_xlabel('Number of molecules', fontsize=13)
    
    if np.average(np.diff(molecules)) != 1:
        print("Interpolating transition probability values.")
        P_matrix = interpolate_Pmatrix(P_matrix_sim)
        
        axs[0].set_yscale('linear')
        axs[1].set_yscale('linear')

        axs[0].plot(P_matrix[:, 0], P_matrix[:, 1], '-', 
                    color='tab:red', linewidth=1)
        axs[1].plot(P_matrix[:, 0], P_matrix[:, 2], '-', 
                    color='tab:red', linewidth=1)

        fig.savefig('P_interpolated_{}.png'.format(file_suffix), dpi=300)
    else:
        P_matrix = np.array(P_matrix_sim[:, :-1], copy=True)
        
    return P_matrix

def interpolate_Pmatrix(P_matrix_sim):
    """If simulation was conducted with ∆N > 1 interpolates values of 
    transiton probabilites for full macrostate landscape"""
    N_aranged = np.arange(0, P_matrix_sim[-1, 0]+1, dtype=int)
    Pins_interp = np.interp(N_aranged, P_matrix_sim[:, 0], P_matrix_sim[:, 1])
    Pdel_interp = np.interp(N_aranged, P_matrix_sim[:, 0], P_matrix_sim[:, 2])
    P_matrix_interp = np.stack((N_aranged, Pins_interp, Pdel_interp), axis=-1)
    
    return P_matrix_interp

def calculate_logPI(P_matrix):
    """Calculates logPI from matrix of transition probabilites"""
    # calculate logPI basing on P values using equation
    # logPI(N) = logPI(N-1) + log(P(N-1->N)/P(N->N-1))
    # initially logPI(N=0) is set to some arbitrary positive value
    logPI_array = np.zeros(P_matrix.shape[0])
    logPI_array[0] = (1)
    logPI_old = 1

    for N in range(1, P_matrix.shape[0]):
        logPI = logPI_old + np.log(P_matrix[N-1, 1] / P_matrix[N, 2])
        logPI_array[N] = logPI
        logPI_old = logPI

    return logPI_array

def recalculate_logPI(logPI, pressures, init_pressure):
    """Recalculates logPI for subsequent pressure values"""
    
    logPI_recalc = np.zeros((logPI.size, pressures.size))
    
    for i in range(pressures.size):                
        j = (i + 1) / pressures.size
        sys.stdout.write('\r')
        sys.stdout.write("Recalculating logPI...   [%-50s] %d%%" 
                         % ('='*int(50*j), 100*j))
        sys.stdout.flush()
        for N in range(logPI.size):
            if pressures[i] == 0:
                pressures[i] = 1e-31
            
            # recalculates logPI with equation
            # logPI(N, µ) = logPI(N, µ0) + [beta(µ-µ0)N]
            logPI_recalc[N, i] = logPI[N] + N * np.log(pressures[i] 
                                                       / init_pressure)
    print()
                 
    return logPI_recalc

def calculate_isotherm(data, pressures, system_size):
    """
    Calculates isotherm with metastable phases for given pressure values
    
    Arguments:
        data: array with logPI values for next (N, P) pairs
        pressures: array with pressures for isotherm
        system_size: number of unit cells in the system
        
    Returns:
        isotherm: array with pressure values and loading for each phase
        binodal: location of binodal point on the isotherm
        net_isotherm: array with pressure values and loading calculated over 
                      whole logPI landscape
    """
    # NoP == number of phases
    NoP_array = np.zeros(data.shape[1], dtype=int)
    
    # checks for number of phases, if more than two script exits
    # the phase boundary is at the local minimum logPI
    # TODO: add the possibility of more than two phases
    for i in range(data.shape[1]):
        logPI = data[:, i]
        min_location = argrelextrema(logPI, np.less, order=int(logPI.size/5))
        NoP = min_location[0].size + 1
        if NoP > 2:
            fig, ax = plt.subplots(figsize=(5, 3), layout='constrained')
            ax.plot(logPI, '-')
            for j in min_location[0]:
                ax.plot(j, logPI[j], 'o', color='tab:red', markersize=3)
            plt.savefig('phases_warning.png', dpi=300)
            print('!!! WARNING !!!')
            print('THERE ARE MORE THAN TWO PHASES')
            print('RESULTS CALCULATED WITH THIS SCRIPT WON\'T BE CORRECT')
            print('A plot of the location of the minima was created')
            print('Exiting...')
            raise(SystemExit)
        NoP_array[i] = NoP
    
    NoP_max = np.max(NoP_array)
    isotherm = np.zeros((data.shape[1], NoP_max + 1))
    net_isotherm = np.zeros((data.shape[1], 2))
    next_phase = False
    sum_PI_phase = np.zeros((data.shape[1], NoP_max))
    
    # for each pressure value identify number of phases and calculate loading 
    for i in range(data.shape[1]):
        j = (i + 1) / data.shape[1]
        sys.stdout.write('\r')
        sys.stdout.write("Calculating isotherm...  [%-50s] %d%%" 
                         % ('='*int(50*j), 100*j))
        sys.stdout.flush()
        
        logPI = data[:, i]
        min_location = argrelextrema(logPI, np.less, order=int(logPI.size/10))
            
        # number of phases at given pressure
        NoP = min_location[0].size + 1
            
        molecules = np.arange(logPI.size)
        PI = np.exp(logPI, dtype='longdouble')
        PI_N = np.multiply(PI, molecules, dtype='longdouble')
        
        # if there is only one phase at given pressure
        if min_location[0].size == 0:
            sum_PI_N = np.sum(PI_N)
            sum_PI = np.sum(PI)
            uptake = sum_PI_N / sum_PI
            
            # checks if we already left gas phase
            if next_phase is True:
                isotherm[i, NoP_max] = uptake
                sum_PI_phase[i, NoP_max - 1] = sum_PI
            else:
                isotherm[i, 1] = uptake
                sum_PI_phase[i, 0] = sum_PI
                
        else:
            for phase in range(NoP_max):
                if phase == 0:
                    sum_PI_N = np.sum(PI_N[:min_location[0][0]])
                    sum_PI = np.sum(PI[:min_location[0][0]])
                    uptake = sum_PI_N / sum_PI
                    isotherm[i, phase + 1] = uptake
                    sum_PI_phase[i, phase] = sum_PI
                elif phase == 1:
                    sum_PI_N = np.sum(PI_N[min_location[0][phase-1]:])
                    sum_PI = np.sum(PI[min_location[0][phase-1]:])
                    if sum_PI == 0.0:
                        uptake = 0.0
                    else:
                        uptake = sum_PI_N / sum_PI
                    isotherm[i, phase + 1] = uptake
                    sum_PI_phase[i, phase] = sum_PI
                    next_phase = True
        
        net_sum_PI_N = np.sum(PI_N)
        net_sum_PI = np.sum(PI)
        net_uptake = net_sum_PI_N / net_sum_PI
        net_isotherm[i, 1] = net_uptake

    print('\n')
    isotherm /= system_size
    isotherm = np.ma.masked_equal(isotherm, 0)
    isotherm[:, 0] = pressures
    net_isotherm /= system_size
    net_isotherm[:, 0] = pressures
    
    if NoP_max > 1:
        sum_PI_phase = np.ma.masked_equal(sum_PI_phase, 0)
        bi_loc = np.argwhere(np.diff(np.sign(
            sum_PI_phase[:, 0] - sum_PI_phase[:, 1]))).flatten()
        
        if bi_loc.size > 0:
            bino = np.array([isotherm[bi_loc[0], 0], isotherm[bi_loc[0], 1], 
                             isotherm[bi_loc[0], 0], isotherm[bi_loc[0], 2]])
        else:
            bino = None
    else:
        bino = None
    
    return isotherm, bino, net_isotherm

def check_Nmax_probability(logPI):
    """Checks if the simulation sufficiently samples both vapor and liquid 
    phases by veryfing probability of observing state Nmax"""
    min_location = argrelextrema(logPI, np.less, order=int(logPI.size/10))
    if min_location[0].size == 0:
        check = np.max(logPI) - logPI[-1]
    else:
        phase_beginning = min_location[0][-1]
        check = np.max(logPI[phase_beginning:]) - logPI[-1]
    if check < 50.0:
        print("!!! WARNING !!!")
        print("The probability of observing the maximum particle "
              "number is too high.")
        print("Additional calculations with higher N are needed.")
        print("\n")
        
    return None

def calculate_free_energy(temperature, logPI_recalc, pressures):
    """Calculates free energy in kJ/mol for given N-P pair"""
    
    MASS_UNIT = 1.6605402e-27 # kg (atomic mass unit)
    LENGTH_UNIT = 1e-10 # m (anstrom)
    TIME_UNIT = 1e-12 # s (picosecond)
    BOLTZMANN_CONSTANT = 1.380650324e-23 # J K^-1
    AVOGADRO_CONSTANT = 6.0221419947e23 # mol^-1
    MOLAR_GAS_CONSTANT = 8.314464919 # J mol^-1 K^-1

    ENERGY_CONVERSION_FACTOR = MASS_UNIT * LENGTH_UNIT**2 / TIME_UNIT**2 # J
    ENERGY_TO_KELVIN = (ENERGY_CONVERSION_FACTOR * AVOGADRO_CONSTANT 
                        / MOLAR_GAS_CONSTANT) # K
    ENERGY_TO_KJ_PER_MOL=((ENERGY_CONVERSION_FACTOR 
                           * AVOGADRO_CONSTANT) / 1000.0) # kJ mol^-1
    K_B = BOLTZMANN_CONSTANT / ENERGY_CONVERSION_FACTOR
    
    beta_kJmol = 1 / ((K_B * float(temperature)) * ENERGY_TO_KJ_PER_MOL)
    logPI_norm = np.copy(logPI_recalc)
    
    for P in range(logPI_norm.shape[1]):
        logPI_norm[:, P] = logPI_recalc[:, P] - np.max(logPI_recalc[:, P])
    
    free_en = np.copy(logPI_norm) / -beta_kJmol
    # molecules = np.arange(free_en.shape[0])[:, np.newaxis]
    # pressures = np.append([0], pressures)
    # free_en = np.hstack((molecules, free_en))
    # free_en = np.vstack((pressures, free_en))
    
    return free_en


prob_files, pressures = check_for_files()

# for file in prob_files:
file = 'tmmc-res_IRMOF-1_1.1.1_92_500.dat'
f = file.removesuffix('.dat')
f = f.split('_')
structure = f[1]
unitcells = f[2]
temperature = f[3]
initial_pressure = int(f[4])
file_suffix = '{}_{}_{}_{}'.format(
    structure, unitcells, temperature, initial_pressure)
system_size = np.prod(np.array(unitcells.split('.'), dtype=int))
print('Calculating properties for {}...'.format(structure))
P_matrix = construct_Pmatrix(file, file_suffix)
logPI = calculate_logPI(P_matrix)
check_Nmax_probability(logPI)
logPI_recalc = recalculate_logPI(logPI, pressures, initial_pressure)
isotherm, binodal, net_isotherm = calculate_isotherm(logPI_recalc, 
                                                      pressures, 
                                                      system_size)
if binodal is None:
    footer = ''
else:
    footer = str('{}\t{}\n{}\t{}'.format(binodal[0], binodal[1], 
                                         binodal[2], binodal[3]))
np.savetxt('isotherm_{}.dat'.format(file_suffix), isotherm, 
            delimiter='\t', footer=footer, fmt='%.7f')
np.savetxt('net_isotherm_{}.dat'.format(file_suffix), net_isotherm, 
            delimiter='\t', fmt='%.7f')
free_energy = calculate_free_energy(temperature, logPI_recalc, pressures)
free_energy = np.flip(free_energy, 0)
np.savetxt('free_energy_{}.dat'.format(file_suffix), free_energy,
            delimiter='\t', fmt='%.13f')
fig, ax = plt.subplots(figsize=(3.3, 2.3), layout='tight')
cmap = plt.get_cmap('RdBu_r', 19)

plt.imshow(free_energy, cmap=cmap, extent=[min(pressures)/1000, max(pressures)/1000,
                                           0, free_energy.shape[0]/2], 
           aspect='auto', vmin=0.0, vmax=10.0)

cbar = plt.colorbar(extend='max', pad=0.05)
cbar.set_label("Free energy (kJ/mol)", fontsize=10)
cbar.ax.tick_params(labelsize=10)
ax.set_xlabel('Pressure (kPa)', fontsize=10)
ax.set_ylabel('Uptake (molecules/u.c.)', fontsize=10)
ax.tick_params(which='both', direction='in', labelsize=10, top="true", right="true")

plt.savefig('free_en_map_{}.png'.format(file_suffix), dpi=300)
