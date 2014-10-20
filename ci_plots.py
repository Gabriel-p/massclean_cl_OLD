
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 58 15:19:04 2014

@author: gabriel
"""

from os import listdir, walk
from os.path import join
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator

'''
Generate plots of contamination index (CI) for each MASSCLEAN cluster versus
their mass, age, distance (extinction) and metallicity values.

Generate plots of CI versus center and radius deltas.

Generate plots of metallicity, age, distance and extinction deltas.
'''


def get_true_memb_n():
    '''
    Obtaines the true members count from each synthetic cluster.
    '''

    # Location of the MASSCLEAN data files.
    dir_dat_files = '/media/rest/github/massclean_cl/synth_clusters'

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

    clust_memb_num = [[], []]
    # Loop through each file.
    for f_indx, sub_dir in enumerate(dir_files[0]):

        # dir_files[1][f_indx] is the name of the file being processed.
        clust_name = dir_files[1][f_indx][:-4]

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

        # Store name and number for this cluster.
        full_name = str(sub_dir + '/' + clust_name)
        clust_memb_num[0].append(full_name)
        clust_memb_num[1].append(N_T)

    return clust_memb_num


def skip_comments(f):
    '''
    Read lines that DO NOT start with a # symbol.
    '''
    for line in f:
        if not line.strip().startswith('#'):
            yield line


# Obtain the true number of cluster members in each synthetic cluster.
clust_memb_num = get_true_memb_n()


# Read ocaat_output.dat file to obtain the MASSCLEAN clusters actual
# data and parameters.
out_file = '/media/rest/github/asteca/output/ASteCA - paper/massclean/ocaat_output.dat'
f = open(out_file)
names, params = [], []
for line in skip_comments(f):
    names.append(line.split()[0])
    params.append(line.split()[1:])

# Separate into masses, distances, metallicities and ages of MASSCLEAN
# clusters.
mass, dist, extinc, metal, age, memb_num = [], [], [], [], [], []
for clust_str in names:
    # Get index for this cluster in members list.
    cl_indx = clust_memb_num[0].index(clust_str)
    memb_num.append(clust_memb_num[1][cl_indx])
    first, second = clust_str.split('/')
    first_a, first_b, first_c = first.split('_')
    mass.append(float(first_a))
    dist.append(float(first_b) / 1000.)
    extinc.append(float(first_c) / 3.1)
    metal.append(float('0.' + second[5:8]))
    age.append(float(second[9:]) / 100.)

# Read clusters parameters obtained by OCAAT.
cenx, ceny, e_cen, rad, e_rad, ci, memb_n, prob, metal_ocaat, e_met, \
age_ocaat, e_age, dist_ocaat, e_dist, ext_ocaat, e_ext = \
[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []
for par_str in params:
    cenx.append(float(par_str[0]))
    ceny.append(float(par_str[1]))
    e_cen.append(float(par_str[2]))
    rad.append(float(par_str[3]))
    e_rad.append(float(par_str[4]))
    ci.append(float(par_str[9]))
    memb_n.append(float(par_str[10]))
    prob.append(float(par_str[12]))
    metal_ocaat.append(float(par_str[15]))
    e_met.append(float(par_str[16]))
    age_ocaat.append(float(par_str[17]))
    e_age.append(float(par_str[18]))
    ext_ocaat.append(float(par_str[19]))
    e_ext.append(float(par_str[20]))
    dist_ocaat.append(float(par_str[21]))
    e_dist.append(float(par_str[22]))


# Separate between those clusters for which the center was assigned less
# than 250px from the actual center.
rad_i, e_rad_i, metal_i, metal_ocaat_i, age_i, age_ocaat_i, e_age_i, dist_i, \
dist_ocaat_i, e_dist_i, extinc_i, ext_ocaat_i, e_ext_i, ci_i, memb_true_i, \
memb_ocaat_i, cent_i, mass_i, prob_i, e_met_i, e_cent_i, names_i = \
[], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [],\
[], []

mass_o, dist_o, ci_o, prob_o = [], [], [], []

cent_diff = np.sqrt((np.array(cenx) - 1024.) ** 2 +
    (np.array(ceny) - 1024.) ** 2)

for i, cl in enumerate(names):
    #if (cent_diff[i] - 80.) > rad[i]:
    if cent_diff[i] > 90.:
        #print 'out', cl, cent_diff[i], ci[i], prob[i]
        mass_o.append(mass[i])
        dist_o.append(dist[i])
        ci_o.append(ci[i])
        prob_o.append(prob[i])
    else:
        names_i.append(cl)
        cent_i.append(cent_diff[i])
        e_cent_i.append(e_cen[i])
        rad_i.append(rad[i])
        e_rad_i.append(e_rad[i])
        mass_i.append(mass[i])
        metal_i.append(metal[i])
        e_met_i.append(e_met[i])
        # Avoid infinity values ahead.
        if metal_ocaat[i] == 0.:
            metal_ocaat_i.append(0.0005)
        else:
            metal_ocaat_i.append(metal_ocaat[i])
        age_i.append(age[i])
        age_ocaat_i.append(age_ocaat[i])
        e_age_i.append(e_age[i])
        dist_i.append(dist[i])
        dist_ocaat_i.append(dist_ocaat[i])
        e_dist_i.append(e_dist[i])
        extinc_i.append(extinc[i])
        ext_ocaat_i.append(ext_ocaat[i])
        e_ext_i.append(e_ext[i])
        ci_i.append(ci[i])
        memb_true_i.append(memb_num[i])
        memb_ocaat_i.append(memb_n[i])
        prob_i.append(prob[i])

print 'Number of assigned centers out of synth clusters region: %d' % len(ci_o)

#print 'Not distant clusters with delta-cent > 20'
#for i, cl in enumerate(names):
    #if cent_diff[i] > 20. and dist[i] < 0.8:
        #print cl, int(cent_diff[i])

print 'Clusters with prob<0.4 value'
for i, cl in enumerate(names):
    if prob[i] < 0.4:
        print cl, int(cent_diff[i]), prob[i], memb_n[i] / memb_num[i]

# Value that holds 50% of clusters.
val_c = sorted(cent_i)[int(0.5 * len(cent_i))]
print '\n 50% limit:', val_c
print 'Center percentages'
print '<10', float(sum(abs(i) < 10. for i in cent_diff)) / len(cent_diff)
print '<20', float(sum(abs(i) < 20. for i in cent_diff)) / len(cent_diff)
print '<50', float(sum(abs(i) < 50. for i in cent_diff)) / len(cent_diff)
print '<80', float(sum(abs(i) < 80. for i in cent_diff)) / len(cent_diff)
print '<90', float(sum(abs(i) < 90. for i in cent_diff)) / len(cent_diff)

# Radius delta: true - OCAAT.
rad_diff_i = 250 - np.array(rad_i)
val_r = sorted(abs(rad_diff_i))[int(0.5 * len(rad_diff_i))]
print '\n 50% limit:', val_r
val_r_90 = sorted(abs(rad_diff_i))[int(0.9 * len(rad_diff_i))]
print '90% limit:', val_r_90
#rad_diff_io = np.array(rad) - 250.
print 'Radius percentages (%d)' % len(rad_diff_i)
print '<10', float(sum(abs(i) < 10. for i in rad_diff_i)) / len(rad_diff_i)
print '<20', float(sum(abs(i) < 20. for i in rad_diff_i)) / len(rad_diff_i)
print '<50', float(sum(abs(i) < 50. for i in rad_diff_i)) / len(rad_diff_i)
print '<80', float(sum(abs(i) < 80. for i in rad_diff_i)) / len(rad_diff_i)
print '<90', float(sum(abs(i) < 90. for i in rad_diff_i)) / len(rad_diff_i)

# Number of members, relative error.
memb_diff_i = (np.array(memb_ocaat_i) - np.array(memb_true_i)) / \
    np.array(memb_true_i)

# Print names, number of members and CI.
# print '\n'
# for i, nm_t in enumerate(memb_true_i):
#     print names_i[i], memb_ocaat_i[i], nm_t, ci_i[i]
# raw_input()

val_memb = sorted(abs(memb_diff_i))[int(0.5 * len(memb_diff_i))]
print '\n 50% limit:', val_memb
val_memb_90 = sorted(abs(memb_diff_i))[int(0.9 * len(memb_diff_i))]
print '90% limit:', val_memb_90
memb_diff_io = (np.array(memb_n) - np.array(memb_num)) / \
    np.array(memb_num)
print 'Member number percentages'
print '<0.1', float(sum(abs(i) < 0.1 for i in memb_diff_io)) / len(memb_diff_io)
print '<0.2', float(sum(abs(i) < 0.2 for i in memb_diff_io)) / len(memb_diff_io)
print '<0.5', float(sum(abs(i) < 0.5 for i in memb_diff_io)) / len(memb_diff_io)
print '<0.8', float(sum(abs(i) < 0.8 for i in memb_diff_io)) / len(memb_diff_io)
print '<0.9', float(sum(abs(i) < 0.9 for i in memb_diff_io)) / len(memb_diff_io)

# Metallicity delta: true - OCAAT.
delta_met_i = np.log10(np.array(metal_i) / 0.019) - \
np.log10(np.array(metal_ocaat_i) / 0.019)
#delta_met_i = np.array(metal_i) - np.array(metal_ocaat_i)
val_m = sorted(abs(delta_met_i))[int(0.5 * len(delta_met_i))]
print '\n 50% limit:', val_m
#delta_met_io = np.log10(np.array(metal) / 0.019) - \
#np.log10(np.array(metal_ocaat) / 0.019)
print 'Metallicity percentages'
print '<0.1', float(sum(abs(i) < 0.1 for i in delta_met_i)) / len(delta_met_i)
print '<0.2', float(sum(abs(i) < 0.2 for i in delta_met_i)) / len(delta_met_i)
print '<0.5', float(sum(abs(i) < 0.5 for i in delta_met_i)) / len(delta_met_i)
print '<1.', float(sum(abs(i) < 1. for i in delta_met_i)) / len(delta_met_i)
print '<1.5', float(sum(abs(i) < 1.5 for i in delta_met_i)) / len(delta_met_i)
print '<1.8', float(sum(abs(i) < 1.8 for i in delta_met_i)) / len(delta_met_i)
print '<2.5', float(sum(abs(i) < 2.5 for i in delta_met_i)) / len(delta_met_i)
# Transform to [Fe/H] errors.
e_feh_i = (1. / np.log(10.)) * (np.array(e_met_i) / np.array(metal_i))

# Age delta: true - OCAAT.
#delta_age_i = (10 ** (np.array(age_i))) / 1.e09 - \
#(10 ** (np.array(age_ocaat_i))) / 1.e09
delta_age_i = np.array(age_i) - np.array(age_ocaat_i)
val_a = sorted(abs(delta_age_i))[int(0.5 * len(delta_age_i))]
print '\n 50% limit:', val_a
#delta_age_io = np.array(age) - np.array(age_ocaat)
print 'Age percentages'
print '<0.1', float(sum(abs(i) < 0.1 for i in delta_age_i)) / len(delta_age_i)
print '<0.2', float(sum(abs(i) < 0.2 for i in delta_age_i)) / len(delta_age_i)
print '<0.5', float(sum(abs(i) < 0.5 for i in delta_age_i)) / len(delta_age_i)
print '<1.3', float(sum(abs(i) < 1.3 for i in delta_age_i)) / len(delta_age_i)
print '<2', float(sum(abs(i) < 2 for i in delta_age_i)) / len(delta_age_i)
print '<3', float(sum(abs(i) < 3 for i in delta_age_i)) / len(delta_age_i)
print '<4', float(sum(abs(i) < 4 for i in delta_age_i)) / len(delta_age_i)


# Distance delta: true - OCAAT.
d_ocaat_i = (10 ** ((np.array(dist_ocaat_i) + 5) / 5)) / 1000.
e_dist_i_dm = 0.2 * np.log(10.) * d_ocaat_i * np.array(e_dist_i)
delta_dist_i = np.array(dist_i) - np.array(d_ocaat_i)
val_d = sorted(abs(delta_dist_i))[int(0.5 * len(delta_dist_i))]
print '\n 50% limit:', val_d
#d_ocaat_io = (10 ** ((np.array(dist_ocaat) + 5) / 5)) / 1000.
#delta_dist_io = np.array(dist) - np.array(d_ocaat_io)
print 'Dist percentages'
print '<0.1.', float(sum(abs(i) < 0.1 for i in delta_dist_i)) / \
    len(delta_dist_i)
print '<0.5', float(sum(abs(i) < 0.5 for i in delta_dist_i)) / \
len(delta_dist_i)
print '<1.', float(sum(abs(i) < 1. for i in delta_dist_i)) / len(delta_dist_i)
print '<1.5', float(sum(abs(i) < 1.5 for i in delta_dist_i)) / \
len(delta_dist_i)
print '<2.', float(sum(abs(i) < 2. for i in delta_dist_i)) / len(delta_dist_i)
print '<3.', float(sum(abs(i) < 3. for i in delta_dist_i)) / len(delta_dist_i)
print '<4.', float(sum(abs(i) < 4. for i in delta_dist_i)) / len(delta_dist_i)
print '<5.', float(sum(abs(i) < 5. for i in delta_dist_i)) / len(delta_dist_i)

# Extinction delta: true - OCAAT.
delta_ext_i = np.array(extinc_i) - np.array(ext_ocaat_i)
val_e = sorted(abs(delta_ext_i))[int(0.5 * len(delta_ext_i))]
print '\n 50% limit:', val_e
#delta_ext_io = np.array(extinc) - np.array(ext_ocaat)
print 'Extinc percentages'
print '<0.025', float(sum(abs(i) < 0.025 for i in delta_ext_i)) / \
    len(delta_ext_i)
print '<0.05', float(sum(abs(i) < 0.05 for i in delta_ext_i)) / \
len(delta_ext_i)
print '<0.1', float(sum(abs(i) < 0.1 for i in delta_ext_i)) / len(delta_ext_i)
print '<0.2', float(sum(abs(i) < 0.2 for i in delta_ext_i)) / len(delta_ext_i)
print '<0.4', float(sum(abs(i) < 0.4 for i in delta_ext_i)) / len(delta_ext_i)


# Check for correlations.
# http://mathworld.wolfram.com/StatisticalCorrelation.html
# https://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient
# http://surveymethodsaddicts.blogspot.com.ar/2008/09/what-is-difference-
# between-correlation.html
print '\n Correlations:'
data = np.array([delta_met_i, delta_age_i, delta_dist_i, delta_ext_i])
print np.corrcoef(data)
## Check
#d_m, e_m = np.mean(delta_dist_i), np.mean(delta_ext_i)
#d_st, e_st = np.std(delta_dist_i), np.std(delta_ext_i)
#de_cov = 0.
#for i, dist in enumerate(delta_dist_i):
    #de_cov += ((dist - d_m) * (delta_ext_i[i] - e_m))
#print (de_cov / len(delta_dist_i)) / (d_st * e_st)
#raw_input()

# Make plot.
fig = plt.figure(figsize=(14, 25))  # create the top-level container
gs = gridspec.GridSpec(4, 3, width_ratios=[1, 1, 0.05])
cm = plt.cm.get_cmap('RdYlBu_r')
xy_font_s = 21

# Y axis parameter.
ci_param_i = np.log(np.array(ci_i))
ymin, ymax = -4.1, 0.
#ci_param_i = np.array(ci_i)
#ymin, ymax = -0.01, 1.01

# Order.
order_i = np.argsort(-(np.array(mass_i) / 4.))
z1_i = np.take(((np.array(mass_i) / 5.) + 5.), order_i)
z2_i = np.take(dist_i, order_i)
# Define age markers and labels.
mrk = {7.: ('o', '$\log(age)=7.$'), 8.: ('s', '$\log(age)=8.$'),
    9.: ('D', '$\log(age)=9.$')}
z3_i = np.take(age_i, order_i)

order_o = np.argsort(-(np.array(mass_o) / 4.))
z1_o = np.take(((np.array(mass_o) / 4.)), order_o)
z2_o = np.take(dist_o, order_o)


# Prob vs CI.
ax0 = plt.subplot(gs[0])
plt.ylim(0.005, 1.05)
plt.xlim(0., 1.05)
plt.ylabel('$prob$', fontsize=xy_font_s)
plt.xlabel('$CI$', fontsize=xy_font_s)
ax0.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
plt.scatter(ci_o, prob_o, c=z2_o, cmap=cm, s=z1_o, marker='D', lw=0.5)
# Order before plotting.
x = np.take(ci_i, order_i)
y = np.take(prob_i, order_i)
plt.scatter(x, y, c=z2_i, cmap=cm, s=z1_i + 25., lw=0.4)

# Memb num vs CI.
axp = plt.subplot(gs[1])
plt.xlabel('$e_{rel} MN$', fontsize=xy_font_s)
plt.ylabel('$\log(CI)$', fontsize=xy_font_s)
plt.ylim(ymin, ymax)
plt.xlim(-0.5, 0.5)
axp.minorticks_on()
axp.yaxis.set_major_locator(MultipleLocator(1.0))
axp.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5)
plt.axvspan(-val_memb, val_memb, facecolor='grey', alpha=0.5, zorder=1)
# Order before plotting.
x = np.take(memb_diff_i, order_i)
y = np.take(ci_param_i, order_i)
#plt.scatter(x, y, c=z2_i, cmap=cm, s=z1_i, zorder=3)
for key, value in sorted(mrk.items()):
    s1 = (z3_i == key)
    SC = plt.scatter(x[s1], y[s1],
        marker=value[0], label=value[1],
        s=z1_i[s1],
        c=z2_i[s1], cmap=cm, lw=0.4, zorder=3)
# Colorbar
axp2 = plt.subplot(gs[2])
cbar = plt.colorbar(SC, cax=axp2)
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('$dist\,(kpc)$', fontsize=xy_font_s, labelpad=-15, y=0.35)

# Delta center vs log(CI)
ax00 = plt.subplot(gs[3])
plt.xlim(-2., 98.)
plt.ylim(ymin, ymax)
plt.ylabel('$\log(CI)$', fontsize=xy_font_s)
plt.xlabel('$\Delta center\,(px)$', fontsize=xy_font_s)
ax00.minorticks_on()
ax00.yaxis.set_major_locator(MultipleLocator(1.0))
ax00.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5)
plt.axvspan(0., val_c, facecolor='grey', alpha=0.5, zorder=1)
plt.errorbar(cent_i, ci_param_i, xerr=e_cent_i, ls='none', color='grey',
    zorder=1)
# Order before plotting.
x = np.take(cent_i, order_i)
y = np.take(ci_param_i, order_i)
for key, value in sorted(mrk.items()):
    s1 = (z3_i == key)
    plt.scatter(x[s1], y[s1],
        marker=value[0], label=value[1],
        s=z1_i[s1],
        c=z2_i[s1], cmap=cm, lw=0.4, zorder=3)
# Plot legend.
leg = plt.legend(loc="lower right", markerscale=0.7, scatterpoints=1,
    fontsize=17)
for i in range(len(mrk)):
    leg.legendHandles[i].set_color('k')
    leg.get_frame().set_alpha(0.5)

# Delta r_cl vs log(CI)
ax01 = plt.subplot(gs[4])
plt.xlim(-148., 148.)
plt.ylim(ymin, ymax)
plt.xlabel('$\Delta r_{cl}\,(px)$', fontsize=xy_font_s)
# make these tick labels invisible
plt.setp(ax01.get_yticklabels(), visible=False)
ax01.minorticks_on()
ax01.yaxis.set_major_locator(MultipleLocator(1.0))
ax01.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5)
plt.axvspan(-val_r, val_r, facecolor='grey', alpha=0.5, zorder=1)
plt.errorbar(rad_diff_i, ci_param_i, xerr=e_rad_i, ls='none', color='grey',
    zorder=1)
# Order before plotting.
x = np.take(rad_diff_i, order_i)
y = np.take(ci_param_i, order_i)
# Clusters inside rad/cent boundary.
for key, value in sorted(mrk.items()):
    s1 = (z3_i == key)
    SC = plt.scatter(x[s1], y[s1],
        marker=value[0], label=value[1],
        s=z1_i[s1],
        c=z2_i[s1], cmap=cm, lw=0.4, zorder=3)
# Colorbar.
ax01 = plt.subplot(gs[5])
cbar = plt.colorbar(SC, cax=ax01)
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('$dist\,(kpc)$', fontsize=xy_font_s, labelpad=-15, y=0.35)

# Delta metallicity vs loc(CI)
ax1 = plt.subplot(gs[6])
plt.xlim(-1.8, 1.8)
plt.ylim(ymin, ymax)
plt.xlabel('$\Delta [Fe/H]$', fontsize=xy_font_s)
plt.ylabel('$\log(CI)$', fontsize=xy_font_s)
ax1.minorticks_on()
ax1.yaxis.set_major_locator(MultipleLocator(1.0))
ax1.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
x = np.take(delta_met_i, order_i)
y = np.take(ci_param_i, order_i)
plt.errorbar(delta_met_i, ci_param_i, xerr=e_feh_i, ls='none', color='grey',
    zorder=1)
for key, value in sorted(mrk.items()):
    s1 = (z3_i == key)
    plt.scatter(x[s1], y[s1],
        marker=value[0], label=value[1],
        s=z1_i[s1],
        c=z2_i[s1], cmap=cm, lw=0.4, zorder=3)
# Vertical shaded area.
plt.axvspan(-val_m, val_m, facecolor='grey', alpha=0.5, zorder=1)

# Delta log(age) vs log(CI)
ax2 = plt.subplot(gs[7])
plt.xlim(-1.9, 1.9)
plt.ylim(ymin, ymax)
plt.xlabel('$\Delta \log(age)$', fontsize=xy_font_s)
# make these tick labels invisible
plt.setp(ax2.get_yticklabels(), visible=False)
ax2.minorticks_on()
ax2.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5)
ax2.yaxis.set_major_locator(MultipleLocator(1.0))
# Order before plotting.
x = np.take(delta_age_i, order_i)
y = np.take(ci_param_i, order_i)
plt.errorbar(delta_age_i, ci_param_i, xerr=e_age_i, ls='none', color='grey',
    zorder=1)
#SC = plt.scatter(x, y, c=z2_i, cmap=cm, s=z1_i, zorder=3, lw=0.5)
for key, value in sorted(mrk.items()):
    s1 = (z3_i == key)
    SC = plt.scatter(x[s1], y[s1],
        marker=value[0], label=value[1],
        s=z1_i[s1],
        c=z2_i[s1], cmap=cm, lw=0.4, zorder=3)
# Vertical shaded area.
plt.axvspan(-val_a, val_a, facecolor='grey', alpha=0.5, zorder=1)
# Colorbar.
ax22 = plt.subplot(gs[8])
cbar = plt.colorbar(SC, cax=ax22)
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('$dist\,(kpc)$', fontsize=xy_font_s, labelpad=-15, y=0.35)

# Delta dist vs log(CI)
ax3 = plt.subplot(gs[9])
plt.xlim(-3.5, 3.5)
plt.ylim(ymin, ymax)
plt.xlabel('$\Delta dist (kpc)$', fontsize=xy_font_s)
plt.ylabel('$\log(CI)$', fontsize=xy_font_s)
ax3.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5)
ax3.minorticks_on()
ax3.yaxis.set_major_locator(MultipleLocator(1.0))
# Order before plotting.
x = np.take(delta_dist_i, order_i)
y = np.take(ci_param_i, order_i)
plt.errorbar(delta_dist_i, ci_param_i, xerr=e_dist_i_dm, ls='none',
    color='grey', zorder=1)
#plt.scatter(x, y, c=z2_i, cmap=cm, s=z1_i, zorder=3, lw=0.5)
for key, value in sorted(mrk.items()):
    s1 = (z3_i == key)
    plt.scatter(x[s1], y[s1],
        marker=value[0], label=value[1],
        s=z1_i[s1],
        c=z2_i[s1], cmap=cm, lw=0.4, zorder=3)
# Vertical shaded area.
plt.axvspan(-val_d, val_d, facecolor='grey', alpha=0.5, zorder=1)

# Delta extinction vs log(CI)
ax4 = plt.subplot(gs[10])
plt.xlim(-0.39, 0.39)
plt.ylim(ymin, ymax)
plt.xlabel('$\Delta E_{(B-V)}$', fontsize=xy_font_s)
# make these tick labels invisible
plt.setp(ax4.get_yticklabels(), visible=False)
ax4.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5)
ax4.minorticks_on()
ax4.yaxis.set_major_locator(MultipleLocator(1.0))
# Order before plotting.
x = np.take(delta_ext_i, order_i)
y = np.take(ci_param_i, order_i)
plt.errorbar(delta_ext_i, ci_param_i, xerr=e_ext_i, ls='none', color='grey',
    zorder=1)
#SC = plt.scatter(x, y, c=z2_i, cmap=cm, s=z1_i, zorder=3, lw=0.5)
for key, value in sorted(mrk.items()):
    s1 = (z3_i == key)
    SC = plt.scatter(x[s1], y[s1],
        marker=value[0], label=value[1],
        s=z1_i[s1],
        c=z2_i[s1], cmap=cm, lw=0.4, zorder=3)
# Vertical shaded area.
plt.axvspan(-val_e, val_e, facecolor='grey', alpha=0.5, zorder=1)
# Colorbar
ax42 = plt.subplot(gs[11])
cbar = plt.colorbar(SC, cax=ax42)
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('$dist\,(kpc)$', fontsize=xy_font_s, labelpad=-15, y=0.35)


plt.tight_layout()
# Output png file.
plt.savefig('ci_out.png', dpi=150)

print '\nEnd.'
