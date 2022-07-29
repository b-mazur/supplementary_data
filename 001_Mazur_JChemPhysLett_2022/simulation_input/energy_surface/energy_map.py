import numpy as np
import sys
import re
import os
import fileinput


def polar2cart(r, theta, phi):
    """Converts polar coordinates to Cartesian ones"""
    x = r * np.cos(theta) * np.cos(phi)
    y = r * np.cos(theta) * np.sin(phi)
    z = r * np.sin(theta)
    return x, y, z


r_max = 6.458 * np.sqrt(3)
energy = np.zeros((1, 7), float)

for phi in np.radians(np.arange(0, 10, 1)):
    for theta in np.radians(np.arange(-90, 90, 1)):
        minimal_energy = np.zeros((1, 7), float)
        for r in np.linspace(3.45, 7.35, 79):
            x, y, z = polar2cart(r, theta, phi)

            # checks if we are out of pore 
            if max(abs(x), abs(y), abs(z)) > 6.458:
                continue

            x = x + 12.916
            y = y + 12.916
            z = z + 12.916

            # changes position of molecule in restart file
            for line in fileinput.input('RestartInitial/System_0/restart_IRMOF-1_1.1.1_80.000000_0', inplace=True):
                if line.strip().startswith('Adsorbate-atom-position:'):
                    line = 'Adsorbate-atom-position: 0 0 ' + '{: >18.12f} {: >18.12f} {: >18.12f}'.format(x, y, z) + '\n'
                sys.stdout.write(line)

            # runs simulation 
            os.system('simulate simulation.input')

            # looks for minimal energy for given r
            for line in open('Output/System_0/output_IRMOF-1_1.1.1_80.000000_0.data', 'r'):
                if 'Current Host-Adsorbate energy:' in line:
                    current_energy = re.findall('[+-]?\d+\.\d+', line)[0]
                    if float(current_energy) < float(minimal_energy[0, 6]):
                        minimal_energy = np.array([[r, theta, phi, x, y, z, current_energy]])
                    break

        energy = np.append(energy, minimal_energy, axis=0)
        energy = energy.astype(float)
        np.savetxt('map.csv', energy, fmt='%.5f', delimiter=',', header='r, theta, phi, x, y, z, minimal energy', comments='')

energy = np.delete(energy, (0), axis=0)
np.savetxt('map.csv', energy, fmt='%.5f', delimiter=',', header='r, theta, phi, x, y, z, minimal energy')
