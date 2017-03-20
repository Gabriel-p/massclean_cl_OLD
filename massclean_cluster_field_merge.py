
from os import listdir, getcwd, walk
from os.path import join, realpath, dirname, exists
from os import makedirs
import numpy as np
import random as rd
import itertools
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def compl_func(x, a):
    '''
    Function that mimics photometric incompleteness.
    '''
    return 1 / (1 + np.exp(x - a))


def compl_removal(region, c_mags):
    '''
    Remove random stars beyond a given magnitude limit according to a
    completeness decreasing function.
    '''

    mag_data = region[3]
    max_mag = max(mag_data)
    # Number of bins.
    bins1 = int((max(mag_data) - min(mag_data)) / 0.2)

    # Histogram of magnitude values.
    mag_hist, bin_edg = np.histogram(mag_data, bins1)
    # Index of maximum magnitude bin, c_mags mags below the max mag value.
    max_indx = min(range(len(bin_edg)),
                   key=lambda i: abs(bin_edg[i] - (max_mag - c_mags)))
    n1, p1 = mag_hist[max_indx], 100.
    # Get completeness percentages.
    a = rd.uniform(2., 4.)
    comp_perc = [compl_func(i, a) * 100.
                 for i in range(len(mag_hist[max_indx:]))]
    # Number of stars that should be removed from each bin.
    di = np.around((abs(mag_hist[max_indx:] - (n1 / p1) *
                    np.asarray(comp_perc))), 0)

    # Store indexes of *all* elements in mag_data whose magnitude
    # value falls between the ranges given.
    c_indx = np.searchsorted(bin_edg[max_indx:], mag_data, side='left')
    N = len(bin_edg[max_indx:])
    mask = (c_indx > 0) & (c_indx < N)
    elements = c_indx[mask]
    indices = np.arange(c_indx.size)[mask]
    sorting_idx = np.argsort(elements, kind='mergesort')
    ind_sorted = indices[sorting_idx]
    x = np.searchsorted(elements, range(N), side='right',
                        sorter=sorting_idx)
    # Indexes.
    rang_indx = [ind_sorted[x[i]:x[i + 1]] for i in range(N - 1)]

    # Pick a number (given by the list 'di') of random elements in
    # each range. Those are the indexes of the elements that
    # should be removed from the three sub-lists.
    rem_indx = []
    for indx, num in enumerate(di):
        if rang_indx[indx].any() and len(rang_indx[indx]) >= num:
            rem_indx.append(np.random.choice(rang_indx[indx],
                            int(num), replace=False))
        else:
            rem_indx.append(rang_indx[indx])

    # Remove items from list.
    # itertools.chain() flattens the list of indexes and sorted()
    # with reverse=True inverts them so we don't change the
    # indexes of the elements in the lists after removing them.
    d_i = sorted(list(itertools.chain(*rem_indx)), reverse=True)
    # Remove those selected indexes from the three sub-lists.
    clust_compl = np.delete(np.asarray(region), d_i, axis=1)

    bins2 = int((max(clust_compl[3]) - min(clust_compl[3])) / 0.2)

    histos = [bins1, bins2, a]

    return clust_compl, histos


def get_clust_names(ver, dir_out):
    '''
    Return names and paths of synthetic clusters.
    '''
    # Path where the code is running
    mypath = realpath(join(getcwd(), dirname(__file__)))
    # Path where the synth clusters are stored.
    sc_path = join(mypath + '/massclean' + ver + '/clusters/')
    # Path where the .DAT output files will be stored.
    out_path = join(mypath + dir_out)

    # Store subdir names [0] and file names [1] inside each subdir.
    dir_files = [[], []]
    for root, dirs, files in walk(sc_path):
        if dirs:
            for subdir in dirs:
                for name in listdir(join(sc_path, subdir)):
                    # Check to see if it's a valid data file.
                    if name.endswith(('.trek')):
                        dir_files[0].append(subdir)
                        dir_files[1].append(name)

    return dir_files, sc_path, out_path


def exp_func(x, e_min=0.015, e_max=0.3):
    '''
    Define exponential function to assign errors.
    '''
    m_min, m_max = min(x), max(x)
    b = (1. / (m_max - m_min)) * np.log(e_max / e_min)
    a = e_max / (np.exp(b * m_max))
    sigma = a * np.exp(b * x)
    # Clip errors.
    sigma[sigma > e_max] = e_max
    return sigma


def error_scatter(sc_path, sub_dir, myfile, max_mag_lim):
    '''
    Reject stars in either the cluster or a field according to the magnitude
    limit set, add errors to the photometric data and scatter in mag and color
    given that errors.
    '''

    # Open cluster/field data file.
    data = np.loadtxt(join(sc_path, sub_dir, myfile), unpack=True)

    if myfile == 'field.plot':
        i_s = 0
    else:
        i_s = 5
    # Store raw MASSCLEAN values.
    x_raw, y_raw, V_raw, UB_raw, BV_raw, VI_raw, JH_raw, HK_raw =\
        data[8 + i_s], data[9 + i_s], data[2 + i_s],\
        data[0 + i_s] - data[1 + i_s], data[1 + i_s] - data[2 + i_s],\
        data[2 + i_s] - data[4 + i_s], data[5 + i_s] - data[6 + i_s],\
        data[6 + i_s] - data[7 + i_s]

    # Filter out stars beyond the maximum magnitude limit.
    x_data = x_raw[V_raw < max_mag_lim]
    y_data = y_raw[V_raw < max_mag_lim]
    V_data = V_raw[V_raw < max_mag_lim]
    UB_data = UB_raw[V_raw < max_mag_lim]
    BV_data = BV_raw[V_raw < max_mag_lim]
    VI_data = VI_raw[V_raw < max_mag_lim]
    JH_data = JH_raw[V_raw < max_mag_lim]
    HK_data = HK_raw[V_raw < max_mag_lim]

    # An ID that starts with a 1 identifies the star as a cluster member.
    # An ID that starts with a 2 identifies the star as a field member.
    k = '2' if myfile == 'field.plot' else '1'
    ids = [int(k + str(i)) for i in range(len(V_data))]

    # Mass data.
    if myfile != 'field.plot':
        mi_raw, ma_raw = data[1], data[2]
        mi_data = mi_raw[V_raw < max_mag_lim]
        ma_data = ma_raw[V_raw < max_mag_lim]
    else:
        mi_data = [0.] * len(V_data)

    # Generate errors. We use a single function for all.
    e_V = exp_func(V_data)

    # Randomly move mag and color through a Gaussian function given
    # their error values.
    V_data = V_data + np.random.normal(0, 1, size=len(e_V)) * e_V
    UB_data = UB_data + np.random.normal(0, 1, size=len(e_V)) * e_V
    BV_data = BV_data + np.random.normal(0, 1, size=len(e_V)) * e_V
    VI_data = VI_data + np.random.normal(0, 1, size=len(e_V)) * e_V
    JH_data = JH_data + np.random.normal(0, 1, size=len(e_V)) * e_V
    HK_data = HK_data + np.random.normal(0, 1, size=len(e_V)) * e_V

    if myfile != 'field.plot':
        x_vals, y_vals, mag_vals, col_vals, mass_vals = [x_raw, x_data], \
            [y_raw, y_data], [V_raw, V_data],\
            [BV_raw, UB_data, BV_data, VI_data, JH_data, HK_data],\
            [mi_raw, ma_raw, mi_data, ma_data]
    else:
        x_vals, y_vals, mag_vals, col_vals, mass_vals = x_data, y_data, \
            V_data, [UB_data, BV_data, VI_data, JH_data, HK_data], mi_data

    return ids, x_vals, y_vals, mag_vals, col_vals, e_V, mass_vals


def write_out_file(out_path_sub, clust_name, region_compl):
    """
    Write synthetic cluster data to file.
    """
    # Unpack lists of data after completeness removal.
    id_cl_fl, x_cl_fl, y_cl_fl, V_cl_fl, e_V_cl_fl, UB_cl_fl,\
        BV_cl_fl, VI_cl_fl, JH_cl_fl, HK_cl_fl, mass_cl_fl = region_compl

    # Check if subdir already exists, if not create it
    if not exists(out_path_sub):
        makedirs(out_path_sub)
    # Store merged cluster+field data to file.
    with open(join(out_path_sub, str(clust_name) + '.DAT'), "w") as f_out:
        f_out.write('#id             x        y        V       eV       UB' +
                    '     eUB        BV      eBV       VI      eVI       ' +
                    'JH      eJH       HK      eHK    m_ini\n')
        for idx, item in enumerate(id_cl_fl):
            f_out.write(
                "{:<8.0f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} "
                "{:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} "
                "{:>8.2f} {:>8.2f} {:>8.2f} {:>8.4f}\n".format(
                    item, x_cl_fl[idx], y_cl_fl[idx], V_cl_fl[idx],
                    e_V_cl_fl[idx], UB_cl_fl[idx], e_V_cl_fl[idx],
                    BV_cl_fl[idx], e_V_cl_fl[idx], VI_cl_fl[idx],
                    e_V_cl_fl[idx], JH_cl_fl[idx], e_V_cl_fl[idx],
                    HK_cl_fl[idx], e_V_cl_fl[idx], mass_cl_fl[idx]))


def make_plots(
        max_mag_lim, metal, age, x_raw, y_raw, V_raw, BV_raw, V_clust,
        BV_cl_fl, dist, vis_abs, mass, mass_vals, id_cl_fl_f, V_cl_fl_f,
        e_V_cl_fl_f, V_cl_fl, e_V_cl_fl, region_full, histos, mass_cl_fl,
        BV_field_c, V_field_c, BV_data_c, V_data_c, x_field_c, y_field_c,
        x_data_c, y_data_c, out_path_sub, clust_name):
    """
    Plot synthetic clusters.
    """
    # Plot field + cluster.
    # figsize(x1, y1), GridSpec(y2, x2) --> To have square plots: x1/x2 =
    # y1/y2 = 2.5
    fig = plt.figure(figsize=(15, 10))  # create the top-level container
    gs = gridspec.GridSpec(4, 6)  # create a GridSpec object

    # x,y finding chart of cluster stars.
    ax1 = plt.subplot(gs[0:2, 0:2])
    # Get max and min values in x,y
    x_delta, y_delta = (max(x_raw) - min(x_raw)) / 10., \
        (max(y_raw) - min(y_raw)) / 10.
    x_min, x_max = min(x_raw) - x_delta, max(x_raw) + x_delta
    y_min, y_max = min(y_raw) - y_delta, max(y_raw) + y_delta
    # Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    # Set axis labels
    plt.xlabel('x (px)', fontsize=12)
    plt.ylabel('y (px)', fontsize=12)
    # Set minor ticks
    ax1.minorticks_on()
    plt.text(0.65, 0.94, '$N_{cluster} = %d$' % len(x_raw),
             transform=ax1.transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=13)
    plt.scatter(x_raw, y_raw, marker='o', c='r',
                s=250. * np.exp(-0.0027 * np.asarray(V_clust) ** 2.5))

    # Cluster's stars CMD.
    ax0 = plt.subplot(gs[0:2, 2:4])
    # Set plot limits
    BV_min, BV_max = min(BV_cl_fl) - 0.2, max(BV_cl_fl) + 0.2
    V_min, V_max = max(max_mag_lim + 0.5, max(V_raw) + 0.2), \
        min(V_raw) - 0.2
    plt.xlim(BV_min, BV_max)
    plt.ylim(V_min, V_max)
    # Set axis labels
    plt.xlabel('$(B-V)$', fontsize=18)
    plt.ylabel('$V$', fontsize=18)
    # Set minor ticks
    ax0.minorticks_on()
    # Set grid
    ax0.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
    text1 = '$N_{cluster} = %d$' '\n' % len(BV_raw)
    text2 = '$V_{min} = %0.1f$' '\n' % max_mag_lim
    text3 = '$\log(age) = %.1f$' '\n' % age
    text4 = '$z = %.4f$' '\n' % metal
    text5 = '$dist = %.1f \,kpc$' '\n' % (dist / 1000.)
    text6 = '$A_V = %.1f \, mag$' '\n' % vis_abs
    text7 = '$M = %d \, M_{\odot}$\n' % mass
    text8 = '$\sum m_i = %d \, M_{\odot}$' % sum(mass_vals[0])
    text = text1 + text2 + text3 + text4 + text5 + text6 + text7 + text8
    plt.text(0.63, 0.54, text, transform=ax0.transAxes,
             bbox=dict(facecolor='white', alpha=0.75), fontsize=13)
    # Plot stars.
    plt.scatter(BV_raw, V_raw, marker='o', c='r', s=15., lw=0.3)
    plt.axhline(y=max_mag_lim, linestyle='-', color='green', lw=1.5)

    # Separate cluster from field stars BEFORE completeness function.
    V_data_f, e_V_data_f, e_BV_data_f, V_field_f, e_V_field_f,\
        e_BV_field_f = [], [], [], [], [], []
    for idx, id_star in enumerate(id_cl_fl_f):
        # Cluster star.
        if str(id_star)[0] == '1':
            V_data_f.append(V_cl_fl_f[idx])
            e_V_data_f.append(e_V_cl_fl_f[idx])
            e_BV_data_f.append(e_V_cl_fl_f[idx])
        # Field star.
        elif str(id_star)[0] == '2':
            V_field_f.append(V_cl_fl_f[idx])
            e_V_field_f.append(e_V_cl_fl_f[idx])
            e_BV_field_f.append(e_V_cl_fl_f[idx])

    # V magnitude error
    ax3 = plt.subplot(gs[0, 4:6])
    # Set plot limits
    plt.xlim(min(V_data_f) - 0.5, max(V_data_f) + 0.5)
    plt.ylim(-0.005, max(e_V_cl_fl) + 0.1)
    # Set axis labels
    plt.ylabel('$\sigma_{V}$', fontsize=18)
    plt.xlabel('$V$', fontsize=18)
    # Set minor ticks
    ax3.minorticks_on()
    plt.text(0.25, 0.85, '$N_{cluster} = %d$' % len(V_data_f),
             transform=ax3.transAxes,
             bbox=dict(facecolor='white', alpha=0.75), fontsize=13)
    # Plot stars errors.
    plt.scatter(V_field_f, e_V_field_f, marker='o', c='k', s=12, lw=0.25,
                facecolors='none')
    plt.scatter(V_data_f, e_V_data_f, marker='o', c='r', lw=0., s=4,
                zorder=3)
    plt.axvline(x=max_mag_lim, linestyle='-', color='green', lw=1.5)

    # (B-V) color error
    ax4 = plt.subplot(gs[1, 4:6])
    # Set plot limits
    plt.xlim(min(V_data_f) - 0.5, max(V_data_f) + 0.5)
    plt.ylim(-0.005, max(e_V_cl_fl) + 0.1)
    # Set axis labels
    plt.ylabel('$\sigma_{(B-V)}$', fontsize=18)
    plt.xlabel('$V$', fontsize=18)
    # Set minor ticks
    ax4.minorticks_on()
    # Plot stars errors.
    plt.scatter(V_field_f, e_BV_field_f, marker='o', c='k', s=12, lw=0.25,
                facecolors='none')
    plt.scatter(V_data_f, e_BV_data_f, marker='o', c='r', lw=0., s=4,
                zorder=3)
    plt.axvline(x=max_mag_lim, linestyle='-', color='green', lw=1.5)

    # Histograms
    ax2 = plt.subplot(gs[2:4, 0:2])
    # Set plot limits
    # Set axis labels
    plt.xlabel('$V$', fontsize=18)
    plt.ylabel('$N$', fontsize=18)
    # Set minor ticks
    ax2.minorticks_on()
    # Backg color.
    ax2.set_axis_bgcolor('#D8D8D8')
    # Set grid
    ax2.grid(b=True, which='major', color='w', linestyle='-', lw=1, zorder=1)
    # Plot stars.
    plt.hist(region_full[3], bins=histos[0], color='blue',
             label='Before removal (N=%d)' % len(region_full[3]),
             histtype='step', zorder=4)
    plt.hist(V_cl_fl, bins=histos[1], color='red',
             label='After removal    (N=%d)' % len(V_cl_fl), histtype='step',
             ls='dashed', hatch="/", zorder=4)
    plt.text(0.05, 0.75, '$\\frac{1}{1 + exp(x - %0.2f)}$' % histos[2],
             transform=ax2.transAxes,
             bbox=dict(facecolor='white', alpha=0.7), fontsize=16)
    # Legends.
    leg2 = plt.legend(fancybox=True, loc='upper left', numpoints=1,
                      fontsize=10)
    # Set the alpha value of the legend.
    leg2.get_frame().set_alpha(0.7)

    # Full region CMD.
    ax5 = plt.subplot(gs[2:4, 2:4])
    # Set plot limits
    BV_min, BV_max = min(BV_cl_fl) - 0.2, max(BV_cl_fl) + 0.2
    V_min, V_max = max(V_cl_fl) + 0.2, min(V_cl_fl) - 0.2
    plt.xlim(BV_min, BV_max)
    plt.ylim(V_min, V_max)
    # Set axis labels
    plt.xlabel('$(B-V)$', fontsize=18)
    plt.ylabel('$V$', fontsize=18)
    # Set minor ticks
    ax5.minorticks_on()
    # Set grid
    ax5.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
             zorder=1)
    # Calculate total number of stars whitin cluster's radius.
    text1 = '$N_{cluster} = %d$\n' % len(BV_data_c)
    text2 = '$\sum m_i = %d \, M_{\odot}$' % sum(mass_cl_fl)
    text = text1 + text2
    plt.text(0.65, 0.87, text, transform=ax5.transAxes,
             bbox=dict(facecolor='white', alpha=0.75), fontsize=13)
    # Plot stars.
    plt.scatter(BV_field_c, V_field_c, marker='o', c='k', s=12, lw=0.25,
                facecolors='none', zorder=3)
    plt.scatter(BV_data_c, V_data_c, marker='o', c='r', s=12., lw=0.5,
                zorder=4)

    # x,y finding chart of full frame
    ax6 = plt.subplot(gs[2:4, 4:6])
    # Get max and min values in x,y
    x_min, x_max = 0., 2048.
    y_min, y_max = 0., 2048.
    # Set plot limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    # Set axis labels
    plt.xlabel('x (px)', fontsize=12)
    plt.ylabel('y (px)', fontsize=12)
    plt.text(0.65, 0.93, '$r_{tidal} = 250 px$', transform=ax6.transAxes,
             bbox=dict(facecolor='white', alpha=0.9), fontsize=13)
    # Set minor ticks
    ax6.minorticks_on()
    circle = plt.Circle((1024., 1024.), 250., color='r', fill=False, lw=2.)
    fig.gca().add_artist(circle)
    plt.scatter(x_field_c, y_field_c, marker='o', c='black',
                s=250. * np.exp(-0.0037 * np.asarray(V_field_c) ** 2.5))
    plt.scatter(x_data_c, y_data_c, marker='o', c='r',
                s=250. * np.exp(-0.0037 * np.asarray(V_data_c) ** 2.5))

    fig.tight_layout()
    # Generate output file for each data file.
    plt.savefig(join(out_path_sub, str(clust_name) + '.png'), dpi=150)

    # Close to release memory.
    plt.clf()
    plt.close()


def main(ver, cmags, max_mag_lim, dir_out):
    """
    Merge cluster and field output files from MASSCLEAN to feed OCAAT.
    Calculates errors and applies Gaussian noise to data.
    """
    # Get names and paths of synthetic clusters.
    dir_files, sc_path, out_path = get_clust_names(ver, dir_out)

    # Iterate through all synthetic clusters inside sub-dirs.
    for f_indx, sub_dir in enumerate(dir_files[0]):

        # Store name of file in 'myfile'.
        cl_file = dir_files[1][f_indx]

        clust_name = cl_file[:-5]
        print(sub_dir, clust_name)

        # Separate mass, distance, visual absorption, metallicity and age
        # values.
        mass, dist, vis_abs = map(float, sub_dir.split('_'))
        metal = float('0.' + clust_name[5:8])
        age = float(clust_name[9:]) / 100.

        # Get raw data and error-scattered photometric data for the cluster.
        id_clust, x_cl, y_cl, V_cl, col_cl, e_V_c, mass_vals = \
            error_scatter(sc_path, sub_dir, cl_file, max_mag_lim)
        x_raw, y_raw, V_raw, BV_raw = x_cl[0], y_cl[0], V_cl[0],\
            col_cl[0]
        x_clust, y_clust, V_clust, col_clust, mass_clust = x_cl[1], y_cl[1],\
            V_cl[1], col_cl[1:], mass_vals[2]

        # Get raw data and error-scattered photometric data for the field.
        id_field, x_field, y_field, V_field, col_field, e_mag_f, mass_f =\
            error_scatter(sc_path, sub_dir, 'field.plot', max_mag_lim)

        # Merge cluster and field.
        id_cl_fl_f = np.concatenate([id_clust, id_field])
        x_cl_fl_f = np.concatenate([x_clust, x_field])
        y_cl_fl_f = np.concatenate([y_clust, y_field])
        V_cl_fl_f = np.concatenate([V_clust, V_field])
        e_V_cl_fl_f = np.concatenate([e_V_c, e_mag_f])
        UB_cl_fl_f = np.concatenate([col_clust[0], col_field[0]])
        BV_cl_fl_f = np.concatenate([col_clust[1], col_field[1]])
        VI_cl_fl_f = np.concatenate([col_clust[2], col_field[2]])
        JH_cl_fl_f = np.concatenate([col_clust[3], col_field[3]])
        HK_cl_fl_f = np.concatenate([col_clust[4], col_field[4]])
        mass_cl_fl_f = np.concatenate([mass_clust, mass_f])

        # Store lists into single list of merged data.
        region_full = [id_cl_fl_f, x_cl_fl_f, y_cl_fl_f, V_cl_fl_f,
                       e_V_cl_fl_f, UB_cl_fl_f, BV_cl_fl_f, VI_cl_fl_f,
                       JH_cl_fl_f, HK_cl_fl_f, mass_cl_fl_f]

        # Call completeness removal function.
        region_compl, histos = compl_removal(region_full, cmags)

        # Write data to file.
        out_path_sub = join(out_path, sub_dir)
        write_out_file(out_path_sub, clust_name, region_compl)

        # Unpack lists of data after completeness removal.
        id_cl_fl, x_cl_fl, y_cl_fl, V_cl_fl, e_V_cl_fl = region_compl[:5]
        BV_cl_fl, mass_cl_fl = region_compl[6], region_compl[-1]

        # Separate cluster from field stars AFTER completeness function.
        x_data_c, y_data_c, V_data_c, e_V_data_c, BV_data_c,\
            e_BV_data_c, x_field_c, y_field_c, V_field_c, e_V_field_c,\
            BV_field_c, e_BV_field_c =\
            [], [], [], [], [], [], [], [], [], [], [], []
        for idx, id_star in enumerate(id_cl_fl):
            # Cluster star.
            if str(id_star)[0] == '1':
                x_data_c.append(x_cl_fl[idx])
                y_data_c.append(y_cl_fl[idx])
                V_data_c.append(V_cl_fl[idx])
                e_V_data_c.append(e_V_cl_fl[idx])
                BV_data_c.append(BV_cl_fl[idx])
                e_BV_data_c.append(e_V_cl_fl[idx])
            # Field star.
            elif str(id_star)[0] == '2':
                x_field_c.append(x_cl_fl[idx])
                y_field_c.append(y_cl_fl[idx])
                V_field_c.append(V_cl_fl[idx])
                e_V_field_c.append(e_V_cl_fl[idx])
                BV_field_c.append(BV_cl_fl[idx])
                e_BV_field_c.append(e_V_cl_fl[idx])

        make_plots(
            max_mag_lim, metal, age, x_raw, y_raw, V_raw, BV_raw, V_clust,
            BV_cl_fl, dist, vis_abs, mass, mass_vals, id_cl_fl_f, V_cl_fl_f,
            e_V_cl_fl_f, V_cl_fl, e_V_cl_fl, region_full, histos, mass_cl_fl,
            BV_field_c, V_field_c, BV_data_c, V_data_c, x_field_c, y_field_c,
            x_data_c, y_data_c, out_path_sub, clust_name)

    print('End.')


if __name__ == "__main__":
    # MASSCLEAN version
    ver = '2.014'
    # Define cut-off minimum magnitude value.
    max_mag_lim = 25.
    # Magnitude value below the minimum magnitude where the completeness
    # removal should start.
    cmags = 2.5
    # Output directory
    dir_out = '/synth_clusters/'

    main(ver, cmags, max_mag_lim, dir_out)
