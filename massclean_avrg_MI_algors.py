# -*- coding: utf-8 -*-
"""
Created on Sat Sep  7 11:41:24 2013

@author: gabriel
"""

from os.path import join
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
from itertools import cycle

   
'''
Plots the average member_index values for all the decontamination algorithms
used, for each of the member indexes defined.
'''



###############################################################################
# Names of the algorithms used.

algor_names = ['KDE-Scott', 'KDE-1', 'KDE-2', 'VB-1', 'VB-2', 'Dias-1', 'Dias-2',
               'random']


# TODO: TOMAR ESTOS VALORES DEL ARCHIVO 'algor_analisys_*.MI-CI.data'

# Store average values in list. The position needs to match the position of
# the curves in the list algor_names.
avrg_list = [[0.53, 0.48, 0.54, 0.42, 0.43, 0.34, 0.34, 0.36],
             [0.12, 0.07, 0.14, 0.08, 0.07, -0.04, -0.04, -0.05],\
             [-0.12, -0.17, -0.10, -0.16, -0.17, -0.28, -0.28, -0.30]]


# Location of the data files for ech algorithm.
alg_dir = '/home/gabriel/clusters/massclean_'
alg_dirs = []
for alg_name in algor_names:
    alg_dirs.append(alg_dir+alg_name+'/')
    
               
# Name of the 'algor_analisys_*.data' files where the * identifies the 
# decontamination algorithm.
algor_fil = 'algor_analisys_' # + *.data
###############################################################################



# Store data of MI ids, MI_index and CI_index in lists.
# Each list (3) corresponds to one MI form. Each sublist corresponds to one
# decontamination algorithm. Each sub-sublist (2) corresponds to the x,y
# coordinates, ie: CI, MI.
ci_vs_mi_1 = [[[], []] for _ in range(len(algor_names))]
ci_vs_mi_2 = [[[], []] for _ in range(len(algor_names))]
ci_vs_mi_3 = [[[], []] for _ in range(len(algor_names))]


# Iterate through all dirs.
for indx, dir_name in enumerate(alg_dirs):
    
    # Join path to file.
    path = join(dir_name, algor_fil+algor_names[indx]+'.MI-CI.data')
    # Store MI ids in data[0], CI in data[1] and MI in data[2]
    data = np.loadtxt(path, unpack=True)

    # Store all values in corresponding list.
    for indx2, id_pt in enumerate(data[0]):
        
        if id_pt == 1:
            ci_vs_mi_1[indx][0].append(round(data[1][indx2],2))
            ci_vs_mi_1[indx][1].append(data[2][indx2])
        elif id_pt == 2:
            ci_vs_mi_2[indx][0].append(data[1][indx2])
            ci_vs_mi_2[indx][1].append(data[2][indx2])
        elif id_pt == 3:
            ci_vs_mi_3[indx][0].append(data[1][indx2])
            ci_vs_mi_3[indx][1].append(data[2][indx2])
        

# Store all lists in single list for plotting.
ci_vs_mi = [ci_vs_mi_1, ci_vs_mi_2, ci_vs_mi_3]


# Make plots.
for mi_num, ci_mi_list in enumerate(ci_vs_mi):
    
    # Set filename of output file.
    out_png = '/home/gabriel/clusters/massclean-algors-'+str(mi_num+1)

    # Make plots.
    fig = plt.figure(figsize=(16, 16)) # create the top-level container
    gs = gridspec.GridSpec(11, 11)  # create a GridSpec object    
    
    
    # Plot KDEs.
    ax0 = plt.subplot(gs[0:5, 0:5])
    
#    plt.title('Member index: %d' % (mi_num+1))
    plt.title('Algorithm: KDE')
    plt.xlabel('CI')
    plt.ylabel(r'MI$_{%d}$' % (mi_num+1))
    plt.xlim(0., 1.0)
    
    ax0.yaxis.set_major_locator(MultipleLocator(0.2))
    
    # Plot grid
    plt.grid(b=True, which='major', color='gray', linestyle='--', zorder=1)
    
    # Line styles for all plots except 'random'.
    lines = ["-","--"]
    linecycler = cycle(lines)
   
    # Plot one line for each KDE algorithm used.
    for i in range(3):
        plt.plot(ci_mi_list[i][0], ci_mi_list[i][1], ls=next(linecycler),
             label='%s (%0.2f)' % (algor_names[i], avrg_list[mi_num][i]))

    # Set max upper y limit to 1.        
    plt.ylim(plt.ylim()[0], 1.0)
              
    # Add text box
    if mi_num == 0:
        text = r'$MI_{%d}$ = $n_m/N_T$' % (mi_num+1)
        x_align, y_align = 0.41, 0.92
    elif mi_num == 1:
        text = r'$MI_{%d}$ = $\frac{\left(\sum^{n_m}{w_m} - \sum^{n_f}{w_f}\right)}{N_T}$' % (mi_num+1)
        x_align, y_align = 0.38, 0.9
    elif mi_num == 2:
        text = r'$MI_{%d}$ = $\frac{\left(\sum^{n_m}{w_m} - \sum^{n_f}{w_f}\right) - |N_T - (n_m+n_f)|}{N_T}$' % (mi_num+1)
        x_align, y_align = 0.2, 0.91
    plt.text(x_align, y_align, text, transform=ax0.transAxes,
             bbox=dict(facecolor='white', alpha=0.6), fontsize=12)
             
    # Plot legend.
    leg = plt.legend(loc="upper right", markerscale=30., scatterpoints=40,
                     fontsize=10)
    # set the alpha value of the legend: it will be translucent
    leg.get_frame().set_alpha(0.5)
    # Make gray background.
    ax0.patch.set_facecolor('gray')
    ax0.patch.set_alpha(0.1)
                
                
                
    # Plot Variable Box algorithm.
    ax1 = plt.subplot(gs[0:5, 6:11])
    
    plt.title('Algorithm: Var Box')
    plt.xlabel('CI')
    plt.ylabel(r'MI$_{%d}$' % (mi_num+1))
    plt.xlim(0., 1.0)
    
    ax1.yaxis.set_major_locator(MultipleLocator(0.2))
    
    # Plot grid
    plt.grid(b=True, which='major', color='gray', linestyle='--', zorder=1)
    
    # Line styles for all plots except 'random'.
    lines = ["-","--"]
    linecycler = cycle(lines)
   
    # Plot one line for each algorithm used.
    for i in range(2):
        j = i+3
        plt.plot(ci_mi_list[j][0], ci_mi_list[j][1], ls=next(linecycler),
             label='%s (%0.2f)' % (algor_names[j], avrg_list[mi_num][j]))

    # Set max upper y limit to 1.        
    plt.ylim(plt.ylim()[0], 1.0)
    # Plot legend.
    leg = plt.legend(loc="upper right", markerscale=30., scatterpoints=40,
                     fontsize=10)
    # set the alpha value of the legend: it will be translucent
    leg.get_frame().set_alpha(0.5)
    # Make gray background.
    ax1.patch.set_facecolor('gray')
    ax1.patch.set_alpha(0.1)


    # Plot Dias algorithm.
    ax2 = plt.subplot(gs[6:11, 0:5])
    
    plt.title('Algorithm: Dias')
    plt.xlabel('CI')
    plt.ylabel(r'MI$_{%d}$' % (mi_num+1))
    plt.xlim(0., 1.0)
    
    ax2.yaxis.set_major_locator(MultipleLocator(0.2))
    
    # Plot grid
    plt.grid(b=True, which='major', color='gray', linestyle='--', zorder=1)
    
    # Line styles for all plots except 'random'.
    lines = ["-","--"]
    linecycler = cycle(lines)
   
    # Plot one line for each algorithm used.
    for i in range(2):
        j = i+5
        plt.plot(ci_mi_list[j][0], ci_mi_list[j][1], ls=next(linecycler),
             label='%s (%0.2f)' % (algor_names[j], avrg_list[mi_num][j]))

    # Set max upper y limit to 1.        
    plt.ylim(plt.ylim()[0], 1.0)
    # Plot legend.
    leg = plt.legend(loc="upper right", markerscale=30., scatterpoints=40,
                     fontsize=10)
    # set the alpha value of the legend: it will be translucent
    leg.get_frame().set_alpha(0.5)
    # Make gray background.
    ax2.patch.set_facecolor('gray')
    ax2.patch.set_alpha(0.1)


    # Plot curves with best values for each algorithm.
    ax3 = plt.subplot(gs[6:11, 6:11])
    
    plt.title('Best')
    plt.xlabel('CI')
    plt.ylabel(r'MI$_{%d}$' % (mi_num+1))
    plt.xlim(0., 1.0)
    
    ax3.yaxis.set_major_locator(MultipleLocator(0.2))
    
    # Plot grid
    plt.grid(b=True, which='major', color='gray', linestyle='--', zorder=1)
    
    # Line styles for all plots except 'random'.
    lines = ["-","--"]
    linecycler = cycle(lines)
   
    # Plot one line for each algorithm used. The index 7 plots the random curve,
    # the others plot the best values for each algorithm.
    for i in (2,4,5,7):
        plt.plot(ci_mi_list[i][0], ci_mi_list[i][1], ls=next(linecycler),
             label='%s (%0.2f)' % (algor_names[i], avrg_list[mi_num][i]))

    # Plot random curve.             
#    plt.plot(ci_mi_list[7][0], ci_mi_list[7][1], ls=':',
#             label='%s (%0.2f)' % (algor_names[7], avrg_list[mi_num][7]))

    # Set max upper y limit to 1.        
    plt.ylim(plt.ylim()[0], 1.0)
        
    # Plot legend.
    leg = plt.legend(loc="upper right", markerscale=30., scatterpoints=40,
                     fontsize=10)
    # set the alpha value of the legend: it will be translucent
    leg.get_frame().set_alpha(0.5)
    # Make gray background.
    ax3.patch.set_facecolor('gray')
    ax3.patch.set_alpha(0.1)

                
    # Generate output plot.
    plt.savefig(out_png, dpi=150)


print 'End.'