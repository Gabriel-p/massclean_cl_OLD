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


'''
Merge cluster and field output files from MASSCLEAN to feed the cluster analysis
code. Calculates errors and applies gaussian noise to data.
'''


# Define exponential function to assign errors.
def func(x):
    '''
    Exponential function.
    '''
    a = 1.72e-06
    b = 8.58e-01
    c = -7.137
    d = 9.346e-03
    return a * np.exp(b * x + c) + d


# Path where the .DAT files will be stored.
out_path = '/media/rest/Dropbox/GABRIEL/CARRERA/3-POS-DOC/trabajo/codigo/\
massclean/MASSCLEAN_clusters/'


# Path where the code is running
mypath = realpath(join(getcwd(), dirname(__file__)))


# Store subdir names [0] and file names [1] inside each subdir.
dir_files = [[], []]
for root, dirs, files in walk(mypath):
    if dirs:
        for subdir in dirs:
            for name in listdir(join(mypath, subdir)):
                # Check to see if it's a valid data file.
                if name.endswith(('.trek')):
                    dir_files[0].append(subdir)
                    dir_files[1].append(name)


# Iterate through all files inside sub-dirs.
for f_indx, sub_dir in enumerate(dir_files[0]):
    
    # Store name of file in 'myfile'.
    myfile = dir_files[1][f_indx]
    
    clust_name = myfile[:-5]
    print sub_dir, clust_name

    # Output subdir.   
    output_subdir = join(mypath, sub_dir)


    # Open cluster data file and process the data.

    # Loads the data in 'myfile' as a list of N lists where N is the number
    # of columns. Each of the N lists contains all the data for the column.
    data = np.loadtxt(join(output_subdir, myfile), unpack=True)
    
    # Store values.
    x_temp, y_temp, mag_temp = data[13], data[14], data[7]
    col1_temp = []
    for idx in range(len(mag_temp)):
        col1_temp.append(data[6][idx]-data[7][idx])
    
    x_data, y_data, mag_data, e_mag, col1_data, e_col = [], [], [], [], [], []
    # Go through all stars to add errors, re-escalate and move the
    # cluster and discard stars with too large V values.
    for idx in range(len(x_temp)):
        
        # Upper and lower limit for max V.
        if mag_temp[idx] < 21.0 and mag_temp[idx] > 8.5:
            
            mag_data.append(mag_temp[idx])
            col1_data.append(col1_temp[idx])
            
            # Re-escalate cluster from 2000x2000 px to 250x250 px and then
            # move it back to the middle of the frame.
            x_data.append(x_temp[idx]/8. + 875.)
            y_data.append(y_temp[idx]/8. + 875.)
            
            # Add errors.
            e_mag.append(func(mag_temp[idx]))
            e_col.append(func(mag_temp[idx]))
            
    # The number 1.0 in ID identifies the star as a cluster member.
    id_clust = [1.0 for _ in x_data]
            
    for idx in range(len(mag_data)):
        # Randomly move mag and color through a Gaussian function.
        mag_data[idx] = rd.gauss(mag_data[idx], e_mag[idx])
        col1_data[idx] = rd.gauss(col1_data[idx], e_col[idx])

    for idx in range(len(mag_data)):
        # Randomly increase errors.
        e_mag[idx] = e_mag[idx] + abs(rd.gauss(e_mag[idx]/7., e_mag[idx]))
        e_col[idx] = e_col[idx] + abs(rd.gauss(e_col[idx]/7., e_col[idx]))
            
            
            
    # Read data from field stars file.
    distance = sub_dir[3:]
    field_file = str(distance)+'_field.plot'
    
    data2 = np.loadtxt(join(mypath, field_file), unpack=True)
    
    # Store values.
    x_f_temp, y_f_temp, mag_f_temp = data2[8], data2[9], data2[2]
    col1_f_temp = []
    for idx in range(len(mag_f_temp)):
        col1_f_temp.append(data2[1][idx]-data2[2][idx])
    
    x_field, y_field, mag_field, e_mag_f, col1_field, e_col_f = [], [], [], [], [], []
    # Go through all stars to add errors, reject too bright stars
    # and add errors.
    for idx in range(len(x_f_temp)):
        
        # Upper and lower limit for V.
        if mag_f_temp[idx] < 21.0 and mag_f_temp[idx] > 8.5:
            
            mag_field.append(mag_f_temp[idx])
            col1_field.append(col1_f_temp[idx])
            
            x_field.append(x_f_temp[idx])
            y_field.append(y_f_temp[idx])
            
            # Add errors.
            e_mag_f.append(func(mag_f_temp[idx]))
            e_col_f.append(func(mag_f_temp[idx]))        
    
    # The number 2.0 in ID identifies the star as a field member.
    id_field = [2.0 for _ in x_field]
        
    for idx in range(len(mag_field)):
        # Randomly move mag and color through a Gaussian function.
        mag_field[idx] = rd.gauss(mag_field[idx], e_mag_f[idx])
        col1_field[idx] = rd.gauss(col1_field[idx], e_col_f[idx])            
        
    for idx in range(len(id_field)):
        # Randomly increase errors.
        e_mag_f[idx] = e_mag_f[idx] + abs(rd.gauss(e_mag_f[idx]/7., e_mag_f[idx]))
        e_col_f[idx] = e_col_f[idx] + abs(rd.gauss(e_col_f[idx]/7., e_col_f[idx]))


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

    with open(join(out_path_sub, str(clust_name)+'.DAT'), "w") as f_out:
        f_out.write('# id  x  y  V  eV  BV  eBV\n')
        for idx, item in enumerate(id_cl_fl):
          f_out.write("%d  %f  %f  %f  %f  %f  %f\n" % (item, x_cl_fl[idx],\
          y_cl_fl[idx], mag_cl_fl[idx], e_mag_cl_fl[idx], col1_cl_fl[idx], \
          e_col_cl_fl[idx]))

                
print 'End.'