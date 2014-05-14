
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
out_file = '/media/rest/github/ocaat/output/massclean/ocaat_output.dat'
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
cenx, ceny, e_cen, rad, e_rad, ci, prob, metal_ocaat, e_met, age_ocaat, e_age,\
dist_ocaat, e_dist, ext_ocaat, e_ext = [], [], [], [], [], [], [], [], \
[], [], [], [], [], [], []
for par_str in params:
    cenx.append(float(par_str[0]))
    ceny.append(float(par_str[1]))
    e_cen.append(float(par_str[2]))
    rad.append(float(par_str[3]))
    e_rad.append(float(par_str[4]))
    ci.append(float(par_str[9]))
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
dist_ocaat_i, e_dist_i, extinc_i, ext_ocaat_i, e_ext_i, ci_i, cent_i, mass_i, \
prob_i, e_met_i, e_cent_i = [], [], [], [], [], [], [], [], [], [], [], [], \
[], [], [], [], [], [], []

mass_o, dist_o, ci_o, prob_o = [], [], [], []

cent_diff = np.sqrt((np.array(cenx) - 1024.) ** 2 +
    (np.array(ceny) - 1024.) ** 2)

for i, cl in enumerate(names):
    if (cent_diff[i] - 250.) > rad[i]:
    #if cent_diff[i] > 250.:
        #print 'out', cl, cent_diff[i], ci[i], prob[i]
        #cent_o.append(cent_diff[i])
        #rad_o.append(rad[i])
        mass_o.append(mass[i])
        #metal_o.append(metal[i])
        #e_met_o.append(e_met[i])
        #metal_ocaat_o.append(metal_ocaat[i])
        #age_o.append(age[i])
        #age_ocaat_o.append(age_ocaat[i])
        dist_o.append(dist[i])
        #dist_ocaat_o.append(dist_ocaat[i])
        #extinc_o.append(extinc[i])
        #ext_ocaat_o.append(ext_ocaat[i])
        ci_o.append(ci[i])
        prob_o.append(prob[i])
    else:
        #print 'in', cl, cent_diff[i]
        cent_i.append(cent_diff[i])
        e_cent_i.append(e_cen[i])
        rad_i.append(rad[i])
        e_rad_i.append(e_rad[i])
        mass_i.append(mass[i])
        metal_i.append(metal[i])
        e_met_i.append(e_met[i])
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
        prob_i.append(prob[i])

print 'Number of assigned centers out of synth clusters region: %d' % len(ci_o)

# Value that holds 50% of clusters.
val_c = sorted(cent_i)[int(0.5 * len(cent_i))]
print '\n Center percentages'
print '<10', float(sum(abs(i) < 10. for i in cent_i)) / len(cent_i)
print '<20', float(sum(abs(i) < 20. for i in cent_i)) / len(cent_i)
print '<50', float(sum(abs(i) < 50. for i in cent_i)) / len(cent_i)
print '<80', float(sum(abs(i) < 80. for i in cent_i)) / len(cent_i)

# Radius in/out.
rad_diff_i = np.array(rad_i) - 250.
val_r = sorted(abs(rad_diff_i))[int(0.5 * len(rad_diff_i))]
print '\n Radius percentages'
print '<10', float(sum(abs(i) < 10. for i in rad_diff_i)) / len(rad_diff_i)
print '<20', float(sum(abs(i) < 20. for i in rad_diff_i)) / len(rad_diff_i)
print '<50', float(sum(abs(i) < 50. for i in rad_diff_i)) / len(rad_diff_i)

# Metallicity in/out.
delta_met_i = np.log10(np.array(metal_i) / 0.019) - \
np.log10(np.array(metal_ocaat_i) / 0.019)
val_m = sorted(abs(delta_met_i))[int(0.5 * len(delta_met_i))]
print '\n Metallicity percentages'
print '<0.1', float(sum(abs(i) < 0.1 for i in delta_met_i)) / len(delta_met_i)
print '<0.2', float(sum(abs(i) < 0.2 for i in delta_met_i)) / len(delta_met_i)
print '<0.5', float(sum(abs(i) < 0.5 for i in delta_met_i)) / len(delta_met_i)
print '<1.3', float(sum(abs(i) < 1.3 for i in delta_met_i)) / len(delta_met_i)
# Transform to [Fe/H] errors.
e_feh_i = (1. / np.log(10.)) * (np.array(e_met_i) / np.array(metal_i))

# Age in/out.
#delta_age_i = (10 ** (np.array(age_i))) / 1.e09 - \
#(10 ** (np.array(age_ocaat_i))) / 1.e09
delta_age_i = np.array(age_i) - np.array(age_ocaat_i)
val_a = sorted(abs(delta_age_i))[int(0.5 * len(delta_age_i))]
print '\n Age percentages'
print '<0.1.', float(sum(abs(i) < 0.1 for i in delta_age_i)) / len(delta_age_i)
print '<0.2', float(sum(abs(i) < 0.2 for i in delta_age_i)) / len(delta_age_i)
print '<0.5', float(sum(abs(i) < 0.5 for i in delta_age_i)) / len(delta_age_i)
print '<1.3', float(sum(abs(i) < 1.3 for i in delta_age_i)) / len(delta_age_i)

# Distance in/out.
d_ocaat_i = (10 ** ((np.array(dist_ocaat_i) + 5) / 5)) / 1000.
delta_dist_i = np.array(dist_i) - np.array(d_ocaat_i)
val_d = sorted(abs(delta_dist_i))[int(0.5 * len(delta_dist_i))]
print '\n Dist percentages'
print '<0.1.', float(sum(abs(i) < 0.1 for i in delta_dist_i)) / \
    len(delta_dist_i)
print '<0.5', float(sum(abs(i) < 0.5 for i in delta_dist_i)) / len(delta_dist_i)
print '<1.', float(sum(abs(i) < 1. for i in delta_dist_i)) / len(delta_dist_i)
print '<1.5', float(sum(abs(i) < 1.5 for i in delta_dist_i)) / len(delta_dist_i)

# Extinction in/out.
delta_ext_i = np.array(extinc_i) - np.array(ext_ocaat_i)
val_e = sorted(abs(delta_ext_i))[int(0.5 * len(delta_ext_i))]
print '\n Extinc percentages'
print '<0.025', float(sum(abs(i) < 0.025 for i in delta_ext_i)) / \
    len(delta_ext_i)
print '<0.05', float(sum(abs(i) < 0.05 for i in delta_ext_i)) / len(delta_ext_i)
print '<0.1', float(sum(abs(i) < 0.1 for i in delta_ext_i)) / len(delta_ext_i)
print '<0.2', float(sum(abs(i) < 0.2 for i in delta_ext_i)) / len(delta_ext_i)
print '<0.4', float(sum(abs(i) < 0.4 for i in delta_ext_i)) / len(delta_ext_i)


# Make plot.
fig = plt.figure(figsize=(14, 25))  # create the top-level container
#gs = gridspec.GridSpec(20, 7)  # create a GridSpec object
gs = gridspec.GridSpec(4, 3, width_ratios=[1, 1, 0.05])
cm = plt.cm.get_cmap('RdYlBu_r')

# Y axis parameter.
ci_param_i = np.log(np.array(ci_i))
ymin, ymax = -4.1, 0.
#ci_param_i = np.array(ci_i)
#ymin, ymax = -0.01, 1.01

# Order.
order_i = np.argsort(-(np.array(mass_i) / 4.))
z1_i = np.take(((np.array(mass_i) / 4.)), order_i)
z2_i = np.take(dist_i, order_i)

order_o = np.argsort(-(np.array(mass_o) / 4.))
z1_o = np.take(((np.array(mass_o) / 4.)), order_o)
z2_o = np.take(dist_o, order_o)

ax0 = plt.subplot(gs[0])
plt.ylim(ymin, ymax)
plt.ylabel('$\log(CI)$', fontsize=14)
plt.xlabel('$CI$', fontsize=12)
plt.xlim(-0.05, 1.05)
ax0.minorticks_on()
ax0.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
x = np.take(ci_i, order_i)
y = np.take(ci_param_i, order_i)
plt.scatter(x, y, c=z2_i, cmap=cm, s=z1_i)

#ax5 = plt.subplot(gs[15:18, 0:3])
axp = plt.subplot(gs[1])
plt.ylim(0., 1.05)
plt.xlim(0., 1.0)
plt.ylabel('prob', fontsize=12)
plt.xlabel('$CI$', fontsize=14)
axp.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
x = np.take(ci_i, order_i)
y = np.take(prob_i, order_i)
plt.scatter(ci_o, prob_o, c=z2_o, cmap=cm, s=z1_o, marker='s',
    facecolor='none', lw=1.5)
SC = plt.scatter(x, y, c=z2_i, cmap=cm, s=z1_i)
axp2 = plt.subplot(gs[2])
cbar = plt.colorbar(SC, cax=axp2)
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

#ax00 = plt.subplot(gs[5:8, 0:3])
ax00 = plt.subplot(gs[3])
plt.xlim(-2., 90.)
plt.ylim(ymin, ymax)
plt.ylabel('$\log(CI)$', fontsize=14)
plt.xlabel('$\Delta center\,(px)$', fontsize=14)
ax00.minorticks_on()
ax00.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5)
plt.axvspan(0., val_c, facecolor='grey', alpha=0.5, zorder=1)
plt.errorbar(cent_i, ci_param_i, xerr=e_cent_i, ls='none', color='grey',
    zorder=1)
# Order before plotting.
x = np.take(cent_i, order_i)
y = np.take(ci_param_i, order_i)
plt.scatter(x, y, c=z2_i, cmap=cm, s=z1_i, zorder=3)

#ax01 = plt.subplot(gs[5:8, 3:6])
ax01 = plt.subplot(gs[4])
plt.xlim(-150., 150.)
plt.ylim(ymin, ymax)
plt.xlabel('$\Delta r_{cl}\,(px)$', fontsize=14)
# make these tick labels invisible
plt.setp(ax01.get_yticklabels(), visible=False)
ax01.minorticks_on()
ax01.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5)
plt.axvspan(-val_r, val_r, facecolor='grey', alpha=0.5, zorder=1)
plt.errorbar(rad_diff_i, ci_param_i, xerr=e_rad_i, ls='none', color='grey',
    zorder=1)
# Order before plotting.
x = np.take(rad_diff_i, order_i)
y = np.take(ci_param_i, order_i)
# Clusters inside rad/cent boundary.
SC = plt.scatter(x, y, c=z2_i, cmap=cm, s=z1_i, zorder=3)
ax01 = plt.subplot(gs[5])
cbar = plt.colorbar(SC, cax=ax01)
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

#ax1 = plt.subplot(gs[8:11, 0:3])
ax1 = plt.subplot(gs[6])
plt.xlim(-1.4, 1.4)
plt.ylim(ymin, ymax)
plt.xlabel('$\Delta [Fe/H]$', fontsize=14)
plt.ylabel('$\log(CI))$', fontsize=14)
ax1.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5)
ax1.minorticks_on()
# Order before plotting.
x = np.take(delta_met_i, order_i)
y = np.take(ci_param_i, order_i)
plt.errorbar(delta_met_i, ci_param_i, xerr=e_feh_i, ls='none', color='grey',
    zorder=1)
plt.scatter(x, y, c=z2_i, cmap=cm, s=z1_i, zorder=3, lw=0.5)
plt.axvspan(-val_m, val_m, facecolor='grey', alpha=0.5, zorder=1)

#ax2 = plt.subplot(gs[8:11, 3:6])
ax2 = plt.subplot(gs[7])
plt.xlim(-1.25, 1.25)
plt.ylim(ymin, ymax)
plt.xlabel('$\Delta \log(age/yr)$', fontsize=14)
# make these tick labels invisible
plt.setp(ax2.get_yticklabels(), visible=False)
ax2.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5)
ax2.minorticks_on()
# Order before plotting.
x = np.take(delta_age_i, order_i)
y = np.take(ci_param_i, order_i)
plt.errorbar(delta_age_i, ci_param_i, xerr=e_age_i, ls='none', color='grey',
    zorder=1)
SC = plt.scatter(x, y, c=z2_i, cmap=cm, s=z1_i, zorder=3, lw=0.5)
plt.axvspan(-val_a, val_a, facecolor='grey', alpha=0.5, zorder=1)
ax22 = plt.subplot(gs[8])
cbar = plt.colorbar(SC, cax=ax22)
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

#ax3 = plt.subplot(gs[11:14, 0:3])
ax3 = plt.subplot(gs[9])
plt.xlim(-1.6, 1.6)
plt.ylim(ymin, ymax)
plt.xlabel('$\Delta dist (kpc)$', fontsize=14)
plt.ylabel('$\log(CI)$', fontsize=14)
ax3.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5)
ax3.minorticks_on()
# Order before plotting.
x = np.take(delta_dist_i, order_i)
y = np.take(ci_param_i, order_i)
plt.errorbar(delta_dist_i, ci_param_i, xerr=e_dist_i, ls='none', color='grey',
    zorder=1)
plt.scatter(x, y, c=z2_i, cmap=cm, s=z1_i, zorder=3, lw=0.5)
plt.axvspan(-val_d, val_d, facecolor='grey', alpha=0.5, zorder=1)

#ax4 = plt.subplot(gs[11:14, 3:6])
ax4 = plt.subplot(gs[10])
plt.xlim(-0.4, 0.4)
plt.ylim(ymin, ymax)
plt.xlabel('$\Delta E_{(B-V)}$', fontsize=14)
# make these tick labels invisible
plt.setp(ax4.get_yticklabels(), visible=False)
ax4.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5)
ax4.minorticks_on()
# Order before plotting.
x = np.take(delta_ext_i, order_i)
y = np.take(ci_param_i, order_i)
plt.errorbar(delta_ext_i, ci_param_i, xerr=e_ext_i, ls='none', color='grey',
    zorder=1)
SC = plt.scatter(x, y, c=z2_i, cmap=cm, s=z1_i, zorder=3, lw=0.5)
plt.axvspan(-val_e, val_e, facecolor='grey', alpha=0.5, zorder=1)
ax42 = plt.subplot(gs[11])
cbar = plt.colorbar(SC, cax=ax42)
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')


plt.tight_layout()
# Output png file.
plt.savefig('ci_out.png', dpi=150)

print 'End.'

#ax11 = plt.subplot(gs[5:8, 4:7], aspect=1)
#plt.xlim(-0.01, 0.04)
#plt.ylim(-0.01, 0.04)
#plt.xlabel('$z_{ocaat}$', fontsize=14)
#plt.ylabel('$z_{real}$', fontsize=14)
#ax11.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
## Order before plotting.
#x = np.take(metal_ocaat_i, order_i)
#y = np.take(metal_i, order_i)
#plt.scatter(x, y, c=z2_i, cmap=cm, s=z1_i, lw=0.3)
##plt.errorbar(metal_ocaat, metal, xerr=e_met, ls='none', color='grey',
    ##elinewidth=0.8)
## 1:1 line
#plt.plot([0., 0.01, 0.02, 0.03], [0., 0.01, 0.02, 0.03], 'k-', ls='--')
#cbar = plt.colorbar()
#cbar.set_ticks([0.5, 1., 3., 5.])
#cbar.set_ticklabels([0.5, 1., 3., 5.])
#cbar.set_label('dist (kpc)')

##ax12 = plt.subplot(gs[5:8, 8:11])
##plt.xlim(-5., 1.5)
##plt.xlabel('$e_z^{rel}$', fontsize=16)
##plt.ylabel('$\log(CI)$', fontsize=14)
##ax12.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
### Order before plotting.
##delta_met_rel = delta_met / np.array(metal)
##x = np.take(delta_met_rel, order)
##y = np.take(ci_param, order)
##plt.scatter(x, y, c=z2, cmap=cm, s=z1, zorder=2)
##plt.axvspan(-0.5, 0.5, facecolor='grey', alpha=0.5, zorder=1)
##cbar = plt.colorbar()
##cbar.set_ticks([0.5, 1., 3., 5.])
##cbar.set_ticklabels([0.5, 1., 3., 5.])
##cbar.set_label('dist (kpc)')


#ax21 = plt.subplot(gs[9:12, 4:7], aspect=1)
#plt.xlim(6.5, 9.5)
#plt.ylim(6.5, 9.5)
#plt.xlabel('$log(age/yr)_{ocaat}$', fontsize=14)
#plt.ylabel('$log(age/yr)_{real}$', fontsize=14)
#ax21.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
## Order before plotting.
#x = np.take(age_ocaat_i, order_i)
#y = np.take(age_i, order_i)
#plt.scatter(x, y, c=z2_i, cmap=cm, s=z1_i, lw=0.3)
##plt.errorbar(age_ocaat, age, xerr=e_age, ls='none', color='grey',
    ##elinewidth=0.8)
## 1:1 line
#plt.plot([6., 10.], [6., 10.], 'k-', ls='--')
#cbar = plt.colorbar()
#cbar.set_ticks([0.5, 1., 3., 5.])
#cbar.set_ticklabels([0.5, 1., 3., 5.])
#cbar.set_label('dist (kpc)')

#ax22 = plt.subplot(gs[9:12, 8:11])
#plt.xlim(-2., 2.)
#plt.xlabel('$e_{age/Gyr}^{rel}$', fontsize=16)
#plt.ylabel('$\log(CI)$', fontsize=14)
#ax22.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
## Order before plotting.
#delta_age_rel = delta_age / (10 ** (np.array(age)) / 1.e09)
#x = np.take(delta_age_rel, order)
#y = np.take(ci_param, order)
#plt.scatter(x, y, c=z2, cmap=cm, s=z1, zorder=2)
#plt.axvspan(-0.5, 0.5, facecolor='grey', alpha=0.5, zorder=1)
#cbar = plt.colorbar()
#cbar.set_ticks([0.5, 1., 3., 5.])
#cbar.set_ticklabels([0.5, 1., 3., 5.])
#cbar.set_label('dist (kpc)')


#ax31 = plt.subplot(gs[13:16, 4:7], aspect=1)
#plt.xlim(-0.5, 7.5)
#plt.ylim(-0.5, 7.5)
#plt.xlabel('$DM_{ocaat}$', fontsize=14)
#plt.ylabel('$DM_{real}$', fontsize=14)
#ax31.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
## Order before plotting.
#x = np.take(d_ocaat_i, order_i)
#y = np.take(dist_i, order_i)
#plt.scatter(x, y, c=z2_i, cmap=cm, s=z1_i, lw=0.3)
##plt.errorbar(d_ocaat_i, dist_i, xerr=e_dist, ls='none', color='grey',
    ##elinewidth=0.8)
## 1:1 line
#plt.plot([-0.5, 8.], [-0.5, 8.], 'k-', ls='--')
#cbar = plt.colorbar()
#cbar.set_ticks([0.5, 1., 3., 5.])
#cbar.set_ticklabels([0.5, 1., 3., 5.])
#cbar.set_label('dist (kpc)')

#ax32 = plt.subplot(gs[13:16, 8:11])
##plt.xlim(-1.5, 1.5)
#plt.xlabel('$e_{d(kpc)}^{rel}$', fontsize=16)
#plt.ylabel('$\log(CI)$', fontsize=14)
#ax32.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
## Order before plotting.
#delta_dist_rel = delta_dist / np.array(dist)
#x = np.take(delta_dist_rel, order)
#y = np.take(ci_param, order)
#plt.scatter(x, y, c=z2, cmap=cm, s=z1, zorder=2)
#plt.axvspan(-0.5, 0.5, facecolor='grey', alpha=0.5, zorder=1)
#cbar = plt.colorbar()
#cbar.set_ticks([0.5, 1., 3., 5.])
#cbar.set_ticklabels([0.5, 1., 3., 5.])
#cbar.set_label('dist (kpc)')


#ax41 = plt.subplot(gs[17:20, 4:7], aspect=1)
#plt.xlim(-0.3, 1.2)
#plt.ylim(-0.3, 1.2)
#plt.xlabel('$E(B-V)_{ocaat}$', fontsize=14)
#plt.ylabel('$E(B-V)_{real}$', fontsize=14)
#ax41.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
## Order before plotting.
#x = np.take(ext_ocaat_i, order_i)
#y = np.take(extinc_i, order_i)
#plt.scatter(x, y, c=z2_i, cmap=cm, s=z1_i, lw=0.3)
##plt.errorbar(ext_ocaat, extinc, xerr=e_dist, ls='none', color='grey',
    ##elinewidth=0.8)
## 1:1 line
#plt.plot([-0.3, 0.5, 1.5], [-0.3, 0.5, 1.5], 'k-', ls='--')
#cbar = plt.colorbar()
#cbar.set_ticks([0.5, 1., 3., 5.])
#cbar.set_ticklabels([0.5, 1., 3., 5.])
#cbar.set_label('dist (kpc)')

#ax42 = plt.subplot(gs[17:20, 8:11])
##plt.xlim(-1.5, 1.5)
#plt.xlabel('$e_{E(B-V)}^{rel}$', fontsize=16)
#plt.ylabel('$\log(CI)$', fontsize=14)
#ax42.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
## Order before plotting.
#delta_ext_rel = delta_ext / (1 + np.array(extinc))
#x = np.take(delta_ext_rel, order)
#y = np.take(ci_param, order)
#plt.scatter(x, y, c=z2, cmap=cm, s=z1, zorder=2)
#plt.axvspan(-0.5, 0.5, facecolor='grey', alpha=0.5, zorder=1)
## Colorbar
#cbar = plt.colorbar()
#cbar.set_ticks([0.5, 1., 3., 5.])
#cbar.set_ticklabels([0.5, 1., 3., 5.])
#cbar.set_label('dist (kpc)')
