
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 58 15:19:04 2014

@author: gabriel
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

'''
Generate plots of contamination index (CI) for each MASSCLEAN cluster versus
their mass, age, distance (extinction) and metallicity values.

Generate plots of CI versus center and radius deltas.

Generate plots of metallicity, age, distance and extinction deltas.
'''


def skip_comments(f):
    '''
    Read lines that DO NOT start with a # symbol.
    '''
    for line in f:
        if not line.strip().startswith('#'):
            yield line


# Read ocaat_output.dat file to obtain the MASSCLEAN clusters actual
# data and parameters the ones estimated by OCAAT.
out_file = 'ocaat_output.dat'
f = open(out_file)
names, params = [], []
for line in skip_comments(f):
    names.append(line.split()[0])
    params.append(line.split()[1:])

# Separate into masses, distances, metallicities and ages of MASSCLEAN
# clusters.
mass, dist, extinc, metal, age = [], [], [], [], []
for clust_str in names:
    first, second = clust_str.split('/')
    first_a, first_b, first_c = first.split('_')
    mass.append(float(first_a))
    dist.append(float(first_b) / 1000.)
    extinc.append(float(first_c) / 3.1)
    metal.append(float('0.' + second[5:8]))
    age.append(float(second[9:]) / 100.)

# Read clusters parameters obtained by OCAAT.
cenx, ceny, rad, ci, prob, metal_ocaat, e_met, age_ocaat, e_age, \
dist_ocaat, e_dist, ext_ocaat, e_ext = [], [], [], [], [], [], [], [], \
[], [], [], [], []
for par_str in params:
    cenx.append(float(par_str[0]))
    ceny.append(float(par_str[1]))
    rad.append(float(par_str[2]))
    ci.append(float(par_str[7]))
    prob.append(float(par_str[10]))
    metal_ocaat.append(float(par_str[13]))
    e_met.append(float(par_str[14]))
    age_ocaat.append(float(par_str[15]))
    e_age.append(float(par_str[16]))
    ext_ocaat.append(float(par_str[17]))
    e_ext.append(float(par_str[18]))
    dist_ocaat.append(float(par_str[19]))
    e_dist.append(float(par_str[20]))


fig = plt.figure(figsize=(20, 35))  # create the top-level container
gs = gridspec.GridSpec(14, 8)  # create a GridSpec object

ax0 = plt.subplot(gs[0:3, 0:3])
plt.xlabel('dist (kpc)', fontsize=12)
plt.ylabel('CI', fontsize=12)
ax0.minorticks_on()
ax0.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
#cm = plt.cm.get_cmap('jet')  # ('RdYlBu')
#import matplotlib.colors as mcolors
#cm = mcolors.ListedColormap([(0, 0, 1),
                               #(0, 1, 0),
                               #(1, 0, 0)])
#categ = []
#for clust in age:
    #if clust == 6.0:
        #categ.append(0)
    #elif clust == 8.0:
        #categ.append(1)
    #elif clust == 9.0:
        #categ.append(2)
    #elif clust == 10.0:
        #categ.append(3)

#colormap = np.array(['b', 'g', 'y', 'r'])
#categories = np.array(categ)

plt.scatter(dist, ci, s=((np.array(mass) / 40.) ** 2), facecolors='none')


ax1 = plt.subplot(gs[4:7, 0:3])
plt.xlim(-0.06, 0.06)
plt.xlabel('metal', fontsize=12)
plt.ylabel('CI', fontsize=12)
delta_met = np.array(metal) - np.array(metal_ocaat)
plt.scatter(delta_met, ci)

ax2 = plt.subplot(gs[4:7, 4:7])
plt.xlim(-4, 4)
plt.xlabel('log(age)', fontsize=12)
plt.ylabel('CI', fontsize=12)
delta_age = np.array(age) - np.array(age_ocaat)
plt.scatter(delta_age, ci, s=mass)


ax3 = plt.subplot(gs[8:11, 0:3])
#plt.xlim(-1.5, 1.5)
plt.xlabel('dist', fontsize=12)
plt.ylabel('CI', fontsize=12)
dist_ocaat_kpc = 10 ** (0.2 * (np.array(dist_ocaat) + 5)) / 1000.
delta_dist = np.array(dist) - dist_ocaat_kpc
cm = plt.cm.get_cmap('RdYlBu_r')
plt.scatter(delta_dist, ci, c=dist_ocaat_kpc, cmap=cm,
    s=((np.array(mass) / 40.) ** 2))
plt.colorbar()

ax4 = plt.subplot(gs[8:11, 4:7])
plt.xlim(-1.5, 1.5)
plt.xlabel('E(B-V)', fontsize=12)
plt.ylabel('CI', fontsize=12)
delta_ext = np.array(extinc) - np.array(ext_ocaat)
plt.scatter(delta_ext, ci)

# Output png file.
plt.savefig('ci_out.png', dpi=150)

print 'End.'