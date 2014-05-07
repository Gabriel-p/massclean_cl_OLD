
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


fig = plt.figure(figsize=(30, 70))  # create the top-level container
gs = gridspec.GridSpec(28, 12)  # create a GridSpec object
cm = plt.cm.get_cmap('RdYlBu_r')

cent_diff = np.sqrt((np.array(cenx) - 1024.) ** 2 +
    (np.array(ceny) - 1024.) ** 2)
rad_diff = np.array(rad) - 250.
delta_met = np.array(metal) - np.array(metal_ocaat)
delta_age = (10 ** (np.array(age))) / 1.e09 - \
(10 ** (np.array(age_ocaat))) / 1.e09
d_ocaat = (10 ** ((np.array(dist_ocaat) + 5) / 5)) / 1000.
delta_dist = np.array(dist) - np.array(d_ocaat)
delta_ext = np.array(extinc) - np.array(ext_ocaat)

# Y axis parameter.
ci_param = np.log(np.array(ci))

# Order.
order = np.argsort(-((np.array(mass) / 50.) ** 2))
z1 = np.take(((np.array(mass) / 50.) ** 2), order)
z2 = np.take(dist, order)

ax0 = plt.subplot(gs[0:3, 0:3])
plt.ylabel('$\log(CI)$', fontsize=14)
plt.xlabel('CI', fontsize=12)
plt.xlim(-0.05, 1.05)
#ax0.minorticks_on()
ax0.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
x = np.take(ci, order)
y = np.take(ci_param, order)
plt.scatter(x, y, c=z2, cmap=cm, s=z1)
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax00 = plt.subplot(gs[0:2, 4:8])
plt.ylabel('$\Delta center\,(px)$', fontsize=14)
#plt.xlim(0., min(max(ci) + 0.1, 1))
#ax00.minorticks_on()
ax00.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# make these tick labels invisible
plt.setp(ax00.get_xticklabels(), visible=False)
# Order before plotting.
x = np.take(ci_param, order)
y = np.take(cent_diff, order)
plt.scatter(x, y, c=z2, cmap=cm, s=z1)
#plt.scatter(ci, cent_diff, c=dist, cmap=cm, s=((np.array(mass) / 40.) ** 2))
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax01 = plt.subplot(gs[2:4, 4:8])
plt.xlabel('$\log(CI)$', fontsize=14)
plt.ylabel('$\Delta r_{cl}\,(px)$', fontsize=14)
plt.ylim(-100., 100.)
#ax01.minorticks_on()
ax01.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
x = np.take(ci_param, order)
y = np.take(rad_diff, order)
plt.scatter(x, y, c=z2, cmap=cm, s=z1)
#ax01.scatter(ci, rad_diff, c=dist, cmap=cm, s=((np.array(mass) / 40.) ** 2))
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax1 = plt.subplot(gs[5:8, 0:3])
#plt.ylim(-100, 40100)
plt.xlabel('$\Delta z$', fontsize=14)
plt.ylabel('$\log(CI)$', fontsize=14)
ax1.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
x = np.take(delta_met, order)
y = np.take(ci_param, order)
plt.errorbar(x, y, xerr=e_met, ls='none', color='grey', zorder=1)
plt.scatter(x, y, c=z2, cmap=cm, s=z1, zorder=3)
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax11 = plt.subplot(gs[5:8, 4:7], aspect=1)
plt.xlim(-0.01, 0.04)
plt.ylim(-0.01, 0.04)
plt.xlabel('$z_{ocaat}$', fontsize=14)
plt.ylabel('$z_{real}$', fontsize=14)
ax11.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
x = np.take(metal_ocaat, order)
y = np.take(metal, order)
plt.scatter(x, y, c=z2, cmap=cm, s=z1)
plt.errorbar(metal_ocaat, metal, xerr=e_met, ls='none', color='grey',
    elinewidth=0.8)
# 1:1 line
plt.plot([0., 0.01, 0.02, 0.03], [0., 0.01, 0.02, 0.03], 'k-', ls='--')
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax12 = plt.subplot(gs[5:8, 8:11])
plt.xlim(-5., 1.5)
plt.xlabel('$e_z^{rel}$', fontsize=16)
plt.ylabel('$\log(CI)$', fontsize=14)
ax12.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
delta_met_rel = delta_met / np.array(metal)
x = np.take(delta_met_rel, order)
y = np.take(ci_param, order)
plt.scatter(x, y, c=z2, cmap=cm, s=z1, zorder=2)
plt.axvspan(-0.5, 0.5, facecolor='grey', alpha=0.5, zorder=1)
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax2 = plt.subplot(gs[9:12, 0:3])
#plt.ylim(-100, 40100)
plt.xlabel('$\Delta (age/Gyr)$', fontsize=14)
plt.ylabel('$\log(CI)$', fontsize=14)
ax2.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
x = np.take(delta_age, order)
y = np.take(ci_param, order)
plt.errorbar(x, y, xerr=e_age, ls='none', color='grey', zorder=1)
plt.scatter(x, y, c=z2, cmap=cm, s=z1, zorder=3)
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax21 = plt.subplot(gs[9:12, 4:7], aspect=1)
plt.xlim(6.5, 9.5)
plt.ylim(6.5, 9.5)
plt.xlabel('$log(age/yr)_{ocaat}$', fontsize=14)
plt.ylabel('$log(age/yr)_{real}$', fontsize=14)
ax21.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
x = np.take(age_ocaat, order)
y = np.take(age, order)
plt.scatter(x, y, c=z2, cmap=cm, s=z1)
plt.errorbar(age_ocaat, age, xerr=e_age, ls='none', color='grey',
    elinewidth=0.8)
# 1:1 line
plt.plot([6., 10.], [6., 10.], 'k-', ls='--')
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax22 = plt.subplot(gs[9:12, 8:11])
plt.xlim(-2., 2.)
plt.xlabel('$e_{age/Gyr}^{rel}$', fontsize=16)
plt.ylabel('$\log(CI)$', fontsize=14)
ax22.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
delta_age_rel = delta_age / (10 ** (np.array(age)) / 1.e09)
x = np.take(delta_age_rel, order)
y = np.take(ci_param, order)
plt.scatter(x, y, c=z2, cmap=cm, s=z1, zorder=2)
plt.axvspan(-0.5, 0.5, facecolor='grey', alpha=0.5, zorder=1)
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax3 = plt.subplot(gs[13:16, 0:3])
#plt.ylim(-100, 40100)
plt.xlabel('$\Delta dist (kpc)$', fontsize=14)
plt.ylabel('$\log(CI)$', fontsize=14)
ax3.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
x = np.take(delta_dist, order)
y = np.take(ci_param, order)
plt.errorbar(x, y, xerr=e_dist, ls='none', color='grey', zorder=1)
plt.scatter(x, y, c=z2, cmap=cm, s=z1, zorder=3)
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax31 = plt.subplot(gs[13:16, 4:7], aspect=1)
plt.xlim(-0.5, 7.5)
plt.ylim(-0.5, 7.5)
plt.xlabel('$DM_{ocaat}$', fontsize=14)
plt.ylabel('$DM_{real}$', fontsize=14)
ax31.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
x = np.take(d_ocaat, order)
y = np.take(dist, order)
plt.scatter(x, y, c=z2, cmap=cm, s=z1)
plt.errorbar(d_ocaat, dist, xerr=e_dist, ls='none', color='grey',
    elinewidth=0.8)
# 1:1 line
plt.plot([-0.5, 8.], [-0.5, 8.], 'k-', ls='--')
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax32 = plt.subplot(gs[13:16, 8:11])
#plt.xlim(-1.5, 1.5)
plt.xlabel('$e_{d(kpc)}^{rel}$', fontsize=16)
plt.ylabel('$\log(CI)$', fontsize=14)
ax32.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
delta_dist_rel = delta_dist / np.array(dist)
x = np.take(delta_dist_rel, order)
y = np.take(ci_param, order)
plt.scatter(x, y, c=z2, cmap=cm, s=z1, zorder=2)
plt.axvspan(-0.5, 0.5, facecolor='grey', alpha=0.5, zorder=1)
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax4 = plt.subplot(gs[17:20, 0:3])
#plt.ylim(-100, 40100)
plt.xlabel('$\Delta E_{(B-V)}$', fontsize=14)
plt.ylabel('$\log(CI)$', fontsize=14)
ax4.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
x = np.take(delta_ext, order)
y = np.take(ci_param, order)
plt.errorbar(x, y, xerr=e_ext, ls='none', color='grey', zorder=1)
plt.scatter(x, y, c=z2, cmap=cm, s=z1, zorder=3)
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax41 = plt.subplot(gs[17:20, 4:7], aspect=1)
plt.xlim(-0.3, 1.2)
plt.ylim(-0.3, 1.2)
plt.xlabel('$E(B-V)_{ocaat}$', fontsize=14)
plt.ylabel('$E(B-V)_{real}$', fontsize=14)
ax41.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
x = np.take(ext_ocaat, order)
y = np.take(extinc, order)
plt.scatter(x, y, c=z2, cmap=cm, s=z1)
plt.errorbar(ext_ocaat, extinc, xerr=e_dist, ls='none', color='grey',
    elinewidth=0.8)
# 1:1 line
plt.plot([-0.3, 0.5, 1.5], [-0.3, 0.5, 1.5], 'k-', ls='--')
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax42 = plt.subplot(gs[17:20, 8:11])
#plt.xlim(-1.5, 1.5)
plt.xlabel('$e_{E(B-V)}^{rel}$', fontsize=16)
plt.ylabel('$\log(CI)$', fontsize=14)
ax42.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
delta_ext_rel = delta_ext / (1 + np.array(extinc))
x = np.take(delta_ext_rel, order)
y = np.take(ci_param, order)
plt.scatter(x, y, c=z2, cmap=cm, s=z1, zorder=2)
plt.axvspan(-0.5, 0.5, facecolor='grey', alpha=0.5, zorder=1)
# Colorbar
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax5 = plt.subplot(gs[21:24, 0:3])
#plt.xlim(-1.5, 1.5)
plt.xlabel('prob', fontsize=12)
plt.ylabel('CI', fontsize=12)
ax5.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
x = np.take(prob, order)
y = np.take(ci, order)
plt.scatter(x, y, c=z2, cmap=cm, s=z1)
#plt.scatter(prob, ci, c=dist, cmap=cm, s=((np.array(mass) / 40.) ** 2))
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

# Output png file.
plt.savefig('ci_out.png', dpi=150)

print 'End.'