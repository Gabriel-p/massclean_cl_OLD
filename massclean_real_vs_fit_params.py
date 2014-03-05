# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 07:05:46 2013

@author: gabriel
"""

from os.path import join
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#import matplotlib.collections as collections


'''
Plots the real parameter values (z, E(B-V), age, dist) versus the ones obtained
with the automatic isochrone fitting algorithm.

* Files read by this code:
    
** File with real values of cluster's parameters:
   dir_memb_files + algor_analisys_+algor_sr+.data

** File with fitted values of those parameters:
   dir_memb_files + algor_track_fit_+algor_sr+.data


* Files created by this code:
    
** Output PNG image:
   dir_memb_files + massclean_real_vs_fit_params.png
    
'''

# Path to algorithm files.
dir_memb_files = '/home/gabriel/clusters/massclean_KDE-Scott/'
algor_sr = 'KDE-Scott'


# Name of files storing the actual data and the fitted data.
real_vals_file = 'algor_analisys_'+algor_sr+'.data'
fit_vals_file = 'algor_track_fit_'+algor_sr+'.data'
# Join to obtain valid paths.
real_file_path = join(dir_memb_files, real_vals_file)
fit_file_path = join(dir_memb_files, fit_vals_file)

# Load data from files.
real_data = np.genfromtxt(real_file_path, dtype=None)
real_data = zip(*real_data)
fit_data = np.genfromtxt(fit_file_path, dtype=None)
fit_data = zip(*fit_data)

# Store "real" data in variables.
r_name, r_IM, r_zmet, r_age, r_ebv, r_dis_mod, r_dist, r_ci = real_data[0], \
real_data[1], real_data[2], real_data[3], real_data[4], real_data[5], \
real_data[6], real_data[7]
# Store Initial mass as array times 10 for plotting.
r_IM = np.array(r_IM)*10.
marker = []
for dist_clust in r_dist:
    if dist_clust == 0.5:
        marker.append('o')
    elif dist_clust == 1.0:
        marker.append('s')
    elif dist_clust == 1.5:
        marker.append('^')
    elif dist_clust == 2.0:
        marker.append('p')
    elif dist_clust == 2.5:
        marker.append('*')
    elif dist_clust == 3.0:
        marker.append('D')
        

# Store fitted data into vars.
f_name, f_zmet, f_age, f_ebv, f_dis_mod, f_dist= fit_data[0], \
fit_data[1], fit_data[2], fit_data[3], fit_data[4], fit_data[5]


# Store indexes pointing to the position of clusters as stored in 'r_name'
# in 'f_name'
f_indexes = map(f_name.index, r_name)
    
    
# Get index quantifying the difference between real values and fitted values.
zmet_diff, age_diff, ebv_diff, dist_diff = [], [], [], []
for idx, cluster in enumerate(r_name):
    
#    # Store values for metallicity.
#    if r_zmet[idx] != 0:
#        zmet_perc.append((r_zmet[idx]-f_zmet[f_indexes[idx]])/r_zmet[idx])
#    else:
#        zmet_perc.append(f_zmet[f_indexes[idx]])
#
#    # Store values for age.
#    if r_age[idx] != 0:
#        age_perc.append((r_age[idx]-f_age[f_indexes[idx]])/r_age[idx])
#    else:
#        age_perc.append(f_age[f_indexes[idx]])
#        
#    # Store values for extinction.
#    if r_ebv[idx] != 0:
#        ebv_perc.append((r_ebv[idx]-f_ebv[f_indexes[idx]])/r_ebv[idx])
#    else:
#        ebv_perc.append(f_ebv[f_indexes[idx]])
#
#    # Store values for distance.
#    if r_dist[idx] != 0:
#        dist_perc.append((r_dist[idx]-f_dist[f_indexes[idx]])/r_dist[idx])
#    else:
#        dist_perc.append(f_dist[f_indexes[idx]])

    # Calculate differences in the sense real - fitted param value divided by
    # the full range of the parameter.
#    zmet_perc.append((r_zmet[idx]-f_zmet[f_indexes[idx]])/0.02)
#    age_perc.append((r_age[idx]-f_age[f_indexes[idx]])/3.16)
#    ebv_perc.append((r_ebv[idx]-f_ebv[f_indexes[idx]])/1.)
#    dist_perc.append((r_dist[idx]-f_dist[f_indexes[idx]])/3.)

    # Calculate differences in the sense real - fitted param value.
    zmet_diff.append((r_zmet[idx]-f_zmet[f_indexes[idx]]))
    age_diff.append((r_age[idx]-f_age[f_indexes[idx]]))
    ebv_diff.append((r_ebv[idx]-f_ebv[f_indexes[idx]]))
    dist_diff.append((r_dist[idx]-f_dist[f_indexes[idx]]))


# Assign fixed values for clusters where no members where selected, for the
# values obtained dividing the diff by the full param range.
#zmet_perc = [1. if i < -10. else i for i in zmet_perc]
#age_perc = [1. if i < -10. else i for i in age_perc]
#ebv_perc = [1. if i < -10. else i for i in ebv_perc]
#dist_perc = [1. if i < -10. else i for i in dist_perc]

# Assign fixed values for clusters where no members where selected.
zmet_diff = [0.03 if i < -10. else i for i in zmet_diff]
age_diff = [3.16 if i < -10. else i for i in age_diff]
ebv_diff = [1. if i < -10. else i for i in ebv_diff]
dist_diff = [3. if i < -10. else i for i in dist_diff]





# Make plots.

# Set filename of output file.
out_png = dir_memb_files+'massclean_real_vs_fit_params.png'

# Make plots.
fig = plt.figure(figsize=(16, 16)) # create the top-level container
gs = gridspec.GridSpec(11, 11)  # create a GridSpec object   


ax0 = plt.subplot(gs[0:5, 0:5])
plt.ylim(-0.031,0.031)
#plt.ylim(-1.05,1.05)
plt.xlim(0.,1.)
plt.xlabel('CI')
#plt.ylabel(r'$[z_r - z_f] / \Delta$', fontsize=18)
plt.ylabel(r'$[z_r - z_f]$', fontsize=18)
# Plot grid
plt.grid(b=True, which='major', color='gray', linestyle='--', zorder=1)
#x_range, y1_range = [(0.,1.)], (-0.25, 0.5)
#c1 = collections.BrokenBarHCollection(x_range, y1_range, facecolor='red', alpha=0.2)
#ax0.add_collection(c1)
cm = plt.cm.get_cmap('RdYlBu_r')
for indx, clust in enumerate(r_name):
#    print clust, r_ci[indx], zmet_perc[indx], r_IM[indx], r_age[indx], marker[indx]
    plt.scatter(r_ci[indx], zmet_diff[indx], s=r_IM[indx], c=r_age[indx],
                cmap=cm, lw=0.2, vmin=0., vmax=3.16, marker=marker[indx])


ax1 = plt.subplot(gs[0:5, 6:11])
plt.ylim(-3.2, 3.2)
#plt.ylim(-1.05,1.05)
plt.xlim(0.,1.)
plt.xlabel('CI')
#plt.ylabel(r'$[Age_r - Age_f] / \Delta\, (Gyr)$', fontsize=18)
plt.ylabel(r'$[Age_r - Age_f]\, (Gyr)$', fontsize=18)
# Plot grid
plt.grid(b=True, which='major', color='gray', linestyle='--', zorder=1)
#x_range, y1_range = [(0.,1.)], (-0.25, 0.5)
#c1 = collections.BrokenBarHCollection(x_range, y1_range, facecolor='red', alpha=0.2)
#ax1.add_collection(c1)
cm = plt.cm.get_cmap('RdYlBu_r')
for indx, clust in enumerate(r_name):
    plt.scatter(r_ci[indx], age_diff[indx], s=r_IM[indx], c=r_age[indx],
                cmap=cm, lw=0.2, vmin=0., vmax=3.16, marker=marker[indx])


ax2 = plt.subplot(gs[6:11, 0:5])
plt.ylim(-1.1, 1.1)
#plt.ylim(-1.05,1.05)
plt.xlim(0.,1.)
plt.xlabel('CI')
#plt.ylabel(r'$[E_{(B-V)}^r - E_{(B-V)}^f] / \Delta$', fontsize=18)
plt.ylabel(r'$[E_{(B-V)}^r - E_{(B-V)}^f]$', fontsize=18)
# Plot grid
plt.grid(b=True, which='major', color='gray', linestyle='--', zorder=1)
#x_range, y1_range = [(0.,1.)], (-0.25, 0.5)
#c1 = collections.BrokenBarHCollection(x_range, y1_range, facecolor='red', alpha=0.2)
#ax2.add_collection(c1)
cm = plt.cm.get_cmap('RdYlBu_r')
for indx, clust in enumerate(r_name):
    plt.scatter(r_ci[indx], ebv_diff[indx], s=r_IM[indx], c=r_age[indx],
                cmap=cm, lw=0.2, vmin=0., vmax=3.16, marker=marker[indx])


ax3 = plt.subplot(gs[6:11, 6:11])
plt.ylim(-3.1, 3.1)
#plt.ylim(-1.05,1.05)
plt.xlim(0.,1.)
plt.xlabel('CI')
#plt.ylabel(r'$[dist_r - dist_f] / \Delta\, (kpc)$', fontsize=18)
plt.ylabel(r'$[dist_r - dist_f]\, (kpc)$', fontsize=18)
# Plot grid
plt.grid(b=True, which='major', color='gray', linestyle='--', zorder=1)
#x_range, y1_range = [(0.,1.)], (-0.25, 0.5)
#c1 = collections.BrokenBarHCollection(x_range, y1_range, facecolor='red', alpha=0.2)
#ax3.add_collection(c1)
cm = plt.cm.get_cmap('RdYlBu_r')
for indx, clust in enumerate(r_name):
    plt.scatter(r_ci[indx], dist_diff[indx], s=r_IM[indx], c=r_age[indx],
                cmap=cm, lw=0.2, vmin=0., vmax=3.16, marker=marker[indx])


# Generate output plot.
plt.savefig(out_png, dpi=150)