# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 15:11:04 2013

@author: gabriel
"""

from os import listdir, walk
from os.path import join
import numpy as np
import warnings
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator


def members_index(mi_num, N_T, memb_w, not_memb_w):
    '''
    This is the index that indicates how well the algorithm behaved for
    that particular cluster. Values below 0.5 mean the algorithm did a very
    poor job at identifying true cluster members. A value > 1 means something
    went wrong and too many stars were identified as true members.

    N_T is the real number of cluster stars in the MASSCLEAN data file
    fed to the code. memb_list contains the weights of all stars identified
    as members by the algorithm that are members. not_memb_list is equivalent
    but for field stars that were incorrectly identified as members.
    '''

    if mi_num == 0:
        # Equivalent to the TPR_90 index in UPMASK. It's the ratio of
        # true cluster members recovered and the total number of true cluster
        # members.
        n_w = sum(i > 0.9 for i in memb_w)
        memb_index = float(n_w) / float(N_T)
    #elif mi_num == 1:
        #memb_index = sum(memb_w) / N_T
    elif mi_num == 1:
       # Punishes field stars assigned as cluster members according to the
       # weights assigned to them.
        memb_index = (sum(memb_w) - sum(not_memb_w)) / N_T

    return memb_index


def make_plots(clust_CI, clust_MI, clust_MI_r, clust_params):
    '''
    Plot the CI vs MI diagram for each MI.
    '''

    # Color is associated with the dist; size with the initial mass and
    # the marker with the age.
    mrk = {7.: ('o', '$\log(age)=7.$'), 8.: ('s', '$\log(age)=8.$'),
        9.: ('D', '$\log(age)=9.$')}

    # Make plot.
    plt.figure(figsize=(14, 25))  # create the top-level container
    gs = gridspec.GridSpec(4, 3, width_ratios=[1, 1, 0.05])
    xy_font_s = 21

    ax0 = plt.subplot(gs[0])
    ax0.set_title('Decontamination algorithm', fontsize=xy_font_s)
    plt.ylabel('$MI_1$', fontsize=xy_font_s)
    plt.xlim(0., 0.97)
    plt.ylim(-0.01, 0.99)
    # make these tick labels invisible
    plt.setp(ax0.get_xticklabels(), visible=False)
    # Set steps in axis.
    ax0.yaxis.set_major_locator(MultipleLocator(0.2))
    # Plot grid
    plt.grid(b=True, which='major', color='gray', linestyle='--', zorder=1)
    # Add text box with MI equation.
    text = r'$MI_1 = n_m/N_{cl}$' '\n' r'  $(MP >\,0.9)$'
    x_align, y_align = 0.57, 0.85
    plt.text(x_align, y_align, text, transform=ax0.transAxes,
             bbox=dict(facecolor='white', alpha=0.6), fontsize=(xy_font_s + 2))
    # Define color map.
    cm = plt.cm.get_cmap('RdYlBu_r')
    # Order.
    mass, age, dist = clust_params
    order = np.argsort(-np.array(mass))
    z1 = np.take((np.array(mass) / 5.), order)
    z2 = np.take(age, order)
    z3 = np.take(dist, order)
    # Order before plotting.
    x = np.take(clust_CI, order)
    y = np.take(clust_MI[0], order)
    for key, value in sorted(mrk.items()):
        s1 = (z2 == key)
        plt.scatter(x[s1], y[s1],
            marker=value[0], label=value[1],
            s=z1[s1],
            c=z3[s1], cmap=cm, lw=0.2)
    # Plot regression line.
    m, b = np.polyfit(clust_CI, clust_MI[0], 1)
    range_CI = np.linspace(0., 1., 10)
    plt.plot(range_CI, m * range_CI + b, c='k', ls='--')

    #
    # Random MI.
    ax1 = plt.subplot(gs[1])
    ax1.set_title('Random probability', fontsize=xy_font_s)
    plt.xlim(0., 0.97)
    plt.ylim(-0.01, 0.99)
    # make these tick labels invisible
    plt.setp(ax1.get_yticklabels(), visible=False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.yaxis.set_major_locator(MultipleLocator(0.2))
    # Plot grid
    plt.grid(b=True, which='major', color='gray', linestyle='--', zorder=1)
    # Define color map.
    cm = plt.cm.get_cmap('RdYlBu_r')
    # Order.
    mass, age, dist = clust_params
    order = np.argsort(-np.array(mass))
    z1 = np.take((np.array(mass) / 5.), order)
    z2 = np.take(age, order)
    z3 = np.take(dist, order)
    # Order before plotting.
    x = np.take(clust_CI, order)
    y = np.take(clust_MI_r[0], order)
    for key, value in sorted(mrk.items()):
        s1 = (z2 == key)
        SC = plt.scatter(x[s1], y[s1],
            marker=value[0], label=value[1],
            s=z1[s1],
            c=z3[s1], cmap=cm, lw=0.2)
    # Plot regression line.
    m, b = np.polyfit(clust_CI, clust_MI_r[0], 1)
    range_CI = np.linspace(0., 1., 10)
    plt.plot(range_CI, m * range_CI + b, c='k', ls='--')
    # Plot legend.
    legend = plt.legend(loc="upper right", markerscale=0.7, scatterpoints=1,
        fontsize=17)
    for i in range(len(mrk)):
        legend.legendHandles[i].set_color('k')
    # Colorbar
    axp2 = plt.subplot(gs[2])
    cbar = plt.colorbar(SC, cax=axp2)
    cbar.set_ticks([0.5, 1., 3., 5.])
    cbar.set_ticklabels([0.5, 1., 3., 5.])
    cbar.set_label('$dist\,(kpc)$', fontsize=xy_font_s, labelpad=-15, y=0.35)

    #
    # Second MI.
    ax3 = plt.subplot(gs[3])
    plt.xlabel('$CI$', fontsize=xy_font_s)
    plt.ylabel('$MI_2$', fontsize=xy_font_s)
    plt.xlim(0., 0.97)
    plt.ylim(max(min(clust_MI[1]) - 0.1, -2.5), 0.99)
    ax3.yaxis.set_major_locator(MultipleLocator(0.4))
    # Plot grid
    plt.grid(b=True, which='major', color='gray', linestyle='--', zorder=1)
    # Add text box with MI equation.
    text = (r'$MI_2 = \frac{\left(\sum^{n_m}{p_m} - ' +
           r' \sum^{n_f}{p_f}\right)}{N_{cl}}$')
    x_align, y_align = 0.52, 0.86
    plt.text(x_align, y_align, text, transform=ax3.transAxes,
             bbox=dict(facecolor='white', alpha=0.6), fontsize=(xy_font_s + 2))
    plt.axhline(y=0., linestyle='--', color='r', zorder=3)
    # Define color map.
    cm = plt.cm.get_cmap('RdYlBu_r')
    # Order.
    mass, age, dist = clust_params
    order = np.argsort(-np.array(mass))
    z1 = np.take((np.array(mass) / 5.), order)
    z2 = np.take(age, order)
    z3 = np.take(dist, order)
    # Order before plotting.
    x = np.take(clust_CI, order)
    y = np.take(clust_MI[1], order)
    for key, value in sorted(mrk.items()):
        s1 = (z2 == key)
        plt.scatter(x[s1], y[s1],
            marker=value[0], label=value[1],
            s=z1[s1],
            c=z3[s1], cmap=cm, lw=0.2)
    # Plot regression line.
    m, b = np.polyfit(clust_CI, clust_MI[1], 1)
    range_CI = np.linspace(0., 1., 10)
    plt.plot(range_CI, m * range_CI + b, c='k', ls='--')
    plt.axhline(y=0., linestyle='--', color='r', zorder=3)

    #
    # Second random MI.
    ax4 = plt.subplot(gs[4])
    plt.xlabel('$CI$', fontsize=xy_font_s)
    plt.xlim(0., 0.97)
    plt.ylim(max(min(clust_MI[1]) - 0.1, -2.5), 0.99)
    # make these tick labels invisible
    plt.setp(ax4.get_yticklabels(), visible=False)
    ax4.yaxis.set_major_locator(MultipleLocator(0.4))
    # Plot grid
    plt.grid(b=True, which='major', color='gray', linestyle='--', zorder=1)
    plt.axhline(y=0., linestyle='--', color='r', zorder=3)
    # Define color map.
    cm = plt.cm.get_cmap('RdYlBu_r')
    # Order.
    mass, age, dist = clust_params
    order = np.argsort(-np.array(mass))
    z1 = np.take((np.array(mass) / 5.), order)
    z2 = np.take(age, order)
    z3 = np.take(dist, order)
    # Order before plotting.
    x = np.take(clust_CI, order)
    y = np.take(clust_MI_r[1], order)
    for key, value in sorted(mrk.items()):
        s1 = (z2 == key)
        plt.scatter(x[s1], y[s1],
            marker=value[0], label=value[1],
            s=z1[s1],
            c=z3[s1], cmap=cm, lw=0.2)
    # Plot regression line.
    m, b = np.polyfit(clust_CI, clust_MI_r[1], 1)
    range_CI = np.linspace(0., 1., 10)
    plt.plot(range_CI, m * range_CI + b, c='k', ls='--')
    # Colorbar
    #axp4 = plt.subplot(gs[5])
    #cbar = plt.colorbar(SC2, cax=axp4)
    #cbar.set_ticks([0.5, 1., 3., 5.])
    #cbar.set_ticklabels([0.5, 1., 3., 5.])
    #cbar.set_label('$dist\,(kpc)$', fontsize=xy_font_s, labelpad=-15, y=0.35)

    # Save to output png file.
    plt.tight_layout()
    out_png = dir_memb_files + 'MI_analisys.png'
    plt.savefig(out_png, dpi=150)
    print 'Plot done.'


'''
Calculates all member_index defined, for a given decontamination algorithm.
Plots the output as a CI vs MI diagram.

First step is to read the *.DAT data file for the synthetic cluster and recover
the true number of members: N_T.

Then we read the *_memb.dat file output by the cluster analysis code and
get the weights of each star assigned as a member into a given list:
memb_w if it is a true member and not_memb_w if it is not. The ID's of each
star identifies it: 1 means member, 2 means field star.

Third step: get 'memb_index' value for each cluster using the values obtained
above.

Fourth: get the contamination index 'cont_ind' from the 'data_output'
file generated by the cluster analysis code.

Finally store the cluster's parameters in a file.

The results are plotted.


* Files read by this code:

** All files with the cluster's data (.DAT files) to count the number of true
   cluster members, as produced by MASSCLEAN:
   dir_dat_files + subdirs + *.DAT

** All files with membership probabilities assigned to store the number and
   weights of all the stars saved as most probable members:
   dir_memb_files + sub_dir + *_memb.dat

** File that stores all the output from running the 'cluster analysis' code
   over the synthetic clusters using a given decontamination algorithm:
   dir_memb_files + data_output


* Files created by this code:

** Output PNG images for each member index calculated.
   dir_memb_files + MI_+mi_num+_analisys.png

'''


################### CHANGE THESE NAMES ACCORDINGLY ############################
# Location of the *_memb.dat output files. This is also the location of the
# 'data_output' file from where to get the cont_ind value.
dir_memb_files = '/media/rest/github/ocaat/output/massclean/'

# Location of the MASSCLEAN data files.
dir_dat_files = '/media/rest/github/massclean_cl/synth_clusters'

# File that stores all the output from running the 'cluster analysis' code
# over the synthetic clusters using a given decontamination algorithm.
data_out_file = dir_memb_files + 'ocaat_output.dat'
###############################################################################


# Store subdir names [0] and file names [1] inside each subdir.
dir_files = [[], []]
for root, dirs, files in walk(dir_dat_files):
    if dirs:
        for subdir in dirs:
            for name in listdir(join(dir_dat_files, subdir)):
                # Check to see if it's a valid data file.
                if name.endswith(('.DAT')):
                    dir_files[0].append(subdir)
                    dir_files[1].append(name)


# Lists that will hold all the member and contamination index values.
clust_MI, clust_MI_r, clust_CI, clust_params = [[], []], [[], []], \
[], [[], [], []]

# Initialize files remaining counter.
i = 0
# Iterate through all files inside all sub-dirs.
for f_indx, sub_dir in enumerate(dir_files[0]):

    # dir_files[1][f_indx] is the name of the file being processed.
    clust_name = dir_files[1][f_indx][:-4]

    # Get 'cont_ind' index for this cluster.
    data3 = np.genfromtxt(data_out_file, dtype=None)
    # Transpose list since unpack=True apparently doesn't work above.
    data3 = zip(*data3)
    # Search for cluster in output file.
    for ind, cluster in enumerate(data3[0]):
        first, second = cluster.split('/')
        if first == sub_dir and second == clust_name:
            cont_ind = data3[10][ind]
            cenx, ceny = data3[1][ind], data3[2][ind]

    # Only plot clusters closer than this limit to the true center.
    cent_diff = np.sqrt((cenx - 1024.) ** 2 + (ceny - 1024.) ** 2)
    if cent_diff < 90.:

        # Store CI of cluster.
        clust_CI.append(cont_ind)

        # Get N_T value for the cluster.
        # Loads the data in file as a list of N lists where N is the number
        # of columns. Each of the N lists contains all the data for the column.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            data = np.genfromtxt(join(dir_dat_files, sub_dir,
                dir_files[1][f_indx]), dtype=str, unpack=True)
        # Initiate true members counter.
        N_T = 0
        for star_id in data[0]:
            # Increase counter for true members.
            N_T += 1 if star_id[0] == '1' else 0

        # Get weights for the stars within the cluster region and assigned
        # membership probabilities by the code.
        data2 = np.loadtxt(join(dir_memb_files, sub_dir,
                                str(clust_name) + '_memb.dat'), unpack=True)
        # Store weights.
        memb_w, not_memb_w = [], []
        for ind, star_id in enumerate(data2[0]):
            if str(star_id)[0] == '1':
                memb_w.append(data2[7][ind])
            else:
                not_memb_w.append(data2[7][ind])

        # Calculate 'memb_index' for this cluster.
        for mi_num in range(2):
            memb_ind = members_index(mi_num, N_T, memb_w, not_memb_w)
            # Append value to list.
            clust_MI[mi_num].append(memb_ind)
            #if mi_num == 0:
                #print cont_ind, memb_ind

            # Random membership assignation.
            memb_w_r = np.random.uniform(low=0., high=1., size=(len(memb_w),))
            not_memb_w_r = np.random.uniform(low=0., high=1.,
                size=(len(not_memb_w),))
            memb_ind_r = members_index(mi_num, N_T, memb_w_r, not_memb_w_r)
            # Append value to list.
            clust_MI_r[mi_num].append(memb_ind_r)

        # Put cluster's parameters in a list.
        first_a, first_b, first_c = sub_dir.split('_')
        mass = float(first_a)                 # In solar masses.
        dist = float(first_b) / 1000.         # In Kpc
        #extinc = float(first_c) / 3.1
        age = float(clust_name[9:]) / 100.    # In log[age/yr].

        # Store in list for plotting purposes.
        clust_params[0].append(mass)
        clust_params[1].append(age)
        clust_params[2].append(dist)

        print len(dir_files[0]) - i, sub_dir, clust_name, cont_ind
        i += 1

# Generate plots.
make_plots(clust_CI, clust_MI, clust_MI_r, clust_params)

print 'End.'