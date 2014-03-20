# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 13:35:31 2013

@author: gabriel
"""

from os import listdir, getcwd, walk
from os.path import join, realpath, dirname, exists
from os import makedirs
import numpy as np
import random as rd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


'''
Merge cluster and field output files from MASSCLEAN to feed OCAAT.
Calculates errors and applies gaussian noise to data.
'''


# Define exponential function to assign errors.
def func(x):
    '''
    Exponential function.
    '''
    a = 2.e-05
    #a, d = 0., 0.
    b = 0.9
    c = -10.
    d = 0.015  # Minimum error value.
    return a * np.exp(b * x + c) + d


# Path where the code is running
mypath = realpath(join(getcwd(), dirname(__file__)))
# Path where the synth clusters are stored.
sd_path = join(mypath + '/massclean2.013/clusters/')
# Path where the .DAT output files will be stored.
out_path = join(mypath + '/synth_clusters/')

# Store subdir names [0] and file names [1] inside each subdir.
dir_files = [[], []]
for root, dirs, files in walk(sd_path):
    if dirs:
        for subdir in dirs:
            for name in listdir(join(sd_path, subdir)):
                # Check to see if it's a valid data file.
                if name.endswith(('.trek')):
                    dir_files[0].append(subdir)
                    dir_files[1].append(name)

# Iterate through all synthetic clusters inside sub-dirs.
for f_indx, sub_dir in enumerate(dir_files[0]):

    # Store name of file in 'myfile'.
    myfile = dir_files[1][f_indx]

    clust_name = myfile[:-5]
    print sub_dir, clust_name

    # Output subdir.
    #output_subdir = join(out_path, sub_dir)

    # Open cluster data file and process the data.

    # Loads the data in 'myfile' as a list of N lists where N is the number
    # of columns. Each of the N lists contains all the data for the column.
    data = np.loadtxt(join(sd_path, sub_dir, myfile), unpack=True)

    # Store values.
    x_temp, y_temp, mag_temp, col1_temp = data[13], data[14], data[7], \
    data[6] - data[7]

    x_data, y_data, mag_data, e_mag, col1_data, e_col = [], [], [], [], [], []
    # Go through all stars to add errors, re-scale the
    # cluster and discard stars with too large V values.
    for idx in range(len(x_temp)):

        # Upper and lower limit for max V.
        if mag_temp[idx] < 21.0 and mag_temp[idx] > 12.:

            mag_data.append(mag_temp[idx])
            col1_data.append(col1_temp[idx])

            # Re-scale cluster from 2048x2048 px to 250x250 px positioned in
            # the middle of the frame.
            x_data.append(x_temp[idx] * (125. / 512.) + 774.)
            y_data.append(y_temp[idx] * (125. / 512.) + 774.)

            # Add errors.
            e_mag.append(func(mag_temp[idx]))
            e_col.append(func(mag_temp[idx]))

    # An ID that stars with a 1 identifies the star as a cluster member.
    id_clust = [int('1' + str(i)) for i in range(len(x_data))]

    for idx in range(len(mag_data)):
        # Randomly move mag and color through a Gaussian function.
        mag_data[idx] = rd.gauss(mag_data[idx], e_mag[idx])
        col1_data[idx] = rd.gauss(col1_data[idx], e_col[idx])

    for idx in range(len(mag_data)):
        # Randomly increase errors.
        e_mag[idx] = e_mag[idx] + abs(rd.gauss(e_mag[idx] / 7., e_mag[idx]))
        e_col[idx] = e_col[idx] + abs(rd.gauss(e_col[idx] / 7., e_col[idx]))

    # Read data from field stars file.
    data2 = np.loadtxt(join(sd_path, sub_dir, 'field.plot'), unpack=True)

    # Store values.
    x_f_temp, y_f_temp, mag_f_temp, col1_f_temp = data2[8], data2[9], \
    data2[2], data2[1] - data2[2]

    x_field, y_field, mag_field, e_mag_f, col1_field, e_col_f = [], [], [],\
    [], [], []
    # Go through all stars to add errors, reject too bright stars
    # and add errors.
    for idx in range(len(x_f_temp)):

        # Upper and lower limit for V.
        if mag_f_temp[idx] < 21.0 and mag_f_temp[idx] > 12.:

            mag_field.append(mag_f_temp[idx])
            col1_field.append(col1_f_temp[idx])

            x_field.append(x_f_temp[idx])
            y_field.append(y_f_temp[idx])

            # Add errors.
            e_mag_f.append(func(mag_f_temp[idx]))
            e_col_f.append(func(mag_f_temp[idx]))

    # An ID that stars with a 2 identifies the star as a field member.
    id_field = [int('2' + str(i)) for i in range(len(x_field))]

    for idx in range(len(mag_field)):
        # Randomly move mag and color through a Gaussian function.
        mag_field[idx] = rd.gauss(mag_field[idx], e_mag_f[idx])
        col1_field[idx] = rd.gauss(col1_field[idx], e_col_f[idx])

    for idx in range(len(id_field)):
        # Randomly increase errors.
        e_mag_f[idx] = e_mag_f[idx] + \
        abs(rd.gauss(e_mag_f[idx] / 7., e_mag_f[idx]))
        e_col_f[idx] = e_col_f[idx] + \
        abs(rd.gauss(e_col_f[idx] / 7., e_col_f[idx]))

    # Merge cluster and field.
    id_cl_fl = id_clust + id_field
    x_cl_fl = x_data + x_field
    y_cl_fl = y_data + y_field
    mag_cl_fl = mag_data + mag_field
    e_mag_cl_fl = e_mag + e_mag_f
    col1_cl_fl = col1_data + col1_field
    e_col_cl_fl = e_col + e_col_f

    # Check if subdir already exists, if not create it
    out_path_sub = join(out_path, sub_dir)
    if not exists(out_path_sub):
        makedirs(out_path_sub)

    with open(join(out_path_sub, str(clust_name) + '.DAT'), "w") as f_out:
        f_out.write('# id  x  y  V  eV  BV  eBV\n')
        for idx, item in enumerate(id_cl_fl):
            f_out.write("%d  %f  %f  %f  %f  %f  %f\n" % (item, x_cl_fl[idx],
            y_cl_fl[idx], mag_cl_fl[idx], e_mag_cl_fl[idx], col1_cl_fl[idx],
            e_col_cl_fl[idx]))

    # Plot field + cluster.
    # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 =
    # y1/y2 = 2.5
    fig = plt.figure(figsize=(10, 5))  # create the top-level container
    gs = gridspec.GridSpec(2, 4)  # create a GridSpec object

    # x,y finding chart of full frame
    ax0 = plt.subplot(gs[0:2, 0:2])
    # Get max and min values in x,y
    x_min, x_max = 0., 2048.
    y_min, y_max = 0., 2048.
    #Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    #Set axis labels
    plt.xlabel('x (px)', fontsize=12)
    plt.ylabel('y (px)', fontsize=12)
    # Set minor ticks
    ax0.minorticks_on()
    circle = plt.Circle((1024., 1024.), 250., color='r', fill=False)
    fig.gca().add_artist(circle)
    # Count the number of very bright stars in the field.
    range_10_perc = (max(mag_cl_fl) - min(mag_cl_fl)) / 10. + min(mag_cl_fl)
    bright_stars = len([i for i in mag_cl_fl if i < range_10_perc])
    # Set exponential factor for high and low density fields.
    exp_factor = -0.0037
    # Set a lower amplitude for fields with very bright stars.
    amplitude = 500 if bright_stars < 10 else 200
    plt.scatter(x_cl_fl, y_cl_fl, marker='o', c='black',
                s=amplitude * np.exp(exp_factor * np.asarray(mag_cl_fl) ** 2.5))

    # Cluster's stars CMD (stars inside cluster's radius)
    ax1 = plt.subplot(gs[0:2, 2:4])
    #Set plot limits
    #col1_min, col1_max = min(col1_cl_fl), max(col1_cl_fl)
    #mag_min, mag_max = max(mag_cl_fl), min(mag_cl_fl)
    col1_min, col1_max = -0.9, 3.9
    mag_min, mag_max = 21.9, 12.1
    plt.xlim(col1_min, col1_max)
    plt.ylim(mag_min, mag_max)
    #Set axis labels
    plt.xlabel('$(B-V)$', fontsize=18)
    plt.ylabel('$V$', fontsize=18)
    # Set minor ticks
    ax1.minorticks_on()
    # Only draw units on axis (ie: 1, 2, 3)
    #ax1.xaxis.set_major_locator(MultipleLocator(1.0))
    # Set grid
    ax1.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    # Calculate total number of stars whitin cluster's radius.
    tot_stars = len(id_clust)
    plt.text(0.55, 0.93, '$r \leq R_{cl}\,|\,N=%d$' % tot_stars,
             transform=ax1.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=16)
    # Plot stars.
    plt.scatter(col1_field, mag_field, marker='o', c='k', s=2.5)
    plt.scatter(col1_data, mag_data, marker='o', c='r', s=8., lw=0.5, zorder=2)

    fig.tight_layout()
    # Generate output file for each data file.
    plt.savefig(join(out_path_sub, str(clust_name) + '.png'), dpi=100)

    # Close to release memory.
    plt.clf()
    plt.close()

    #raw_input()

print 'End.'