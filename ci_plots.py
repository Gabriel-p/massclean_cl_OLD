
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
delta_age = np.array(age) - np.array(age_ocaat)
dist_mod = np.log10(np.array(dist) * 1000.) * 5 - 5
delta_dist = np.array(dist_mod) - np.array(dist_ocaat)
delta_ext = np.array(extinc) - np.array(ext_ocaat)

# Order.
order = np.argsort(-((np.array(mass) / 50.) ** 2))
z1 = np.take(((np.array(mass) / 50.) ** 2), order)
z2 = np.take(dist, order)

ax0 = plt.subplot(gs[0:3, 0:3])
plt.xlabel('log(age/yr)', fontsize=12)
plt.ylabel('CI', fontsize=12)
#ax0.minorticks_on()
ax0.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
x = np.take(age, order)
y = np.take(ci, order)
plt.scatter(x, y, c=z2, cmap=cm, s=z1)
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax00 = plt.subplot(gs[0:2, 4:8])
plt.ylabel('$\Delta center\,(px)$', fontsize=14)
plt.xlim(0., min(max(ci) + 0.1, 1))
#ax00.minorticks_on()
ax00.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# make these tick labels invisible
plt.setp(ax00.get_xticklabels(), visible=False)
# Order before plotting.
x = np.take(ci, order)
y = np.take(cent_diff, order)
plt.scatter(x, y, c=z2, cmap=cm, s=z1)
#plt.scatter(ci, cent_diff, c=dist, cmap=cm, s=((np.array(mass) / 40.) ** 2))
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax01 = plt.subplot(gs[2:4, 4:8])
plt.xlabel('CI', fontsize=12)
plt.ylabel('$\Delta r_{cl}\,(px)$', fontsize=14)
plt.xlim(0., min(max(ci) + 0.1, 1))
#ax01.minorticks_on()
ax01.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
x = np.take(ci, order)
y = np.take(rad_diff, order)
plt.scatter(x, y, c=z2, cmap=cm, s=z1)
#ax01.scatter(ci, rad_diff, c=dist, cmap=cm, s=((np.array(mass) / 40.) ** 2))
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax1 = plt.subplot(gs[5:8, 0:3])
plt.xlim(-0.06, 0.06)
plt.xlabel('$\Delta metal$', fontsize=14)
plt.ylabel('CI', fontsize=12)
ax1.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
x = np.take(delta_met, order)
y = np.take(ci, order)
plt.scatter(x, y, c=z2, cmap=cm, s=z1)
#plt.scatter(delta_met, ci, c=dist, cmap=cm, s=((np.array(mass) / 40.) ** 2))
plt.errorbar(delta_met, ci, xerr=e_met, ls='none', color='k', elinewidth=0.8)
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax11 = plt.subplot(gs[5:8, 4:7], aspect=1)
plt.xlim(-0.01, 0.04)
plt.ylim(-0.01, 0.04)
plt.xlabel('$metal_{ocaat}$', fontsize=14)
plt.ylabel('$metal_{real}$', fontsize=14)
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

ax2 = plt.subplot(gs[9:12, 0:3])
plt.xlim(-4, 4)
plt.xlabel('$\Delta log(age/yr)$', fontsize=14)
plt.ylabel('CI', fontsize=12)
ax2.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
x = np.take(delta_age, order)
y = np.take(ci, order)
plt.scatter(x, y, c=z2, cmap=cm, s=z1)
#plt.scatter(delta_age, ci, c=dist, cmap=cm, s=((np.array(mass) / 40.) ** 2))
plt.errorbar(delta_age, ci, xerr=e_age, ls='none', color='k', elinewidth=0.8)
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax12 = plt.subplot(gs[9:12, 4:7], aspect=1)
plt.xlim(6., 10.)
plt.ylim(6., 10.)
plt.xlabel('$age_{ocaat}$', fontsize=14)
plt.ylabel('$age_{real}$', fontsize=14)
ax12.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
x = np.take(age_ocaat, order)
y = np.take(age, order)
plt.scatter(x, y, c=z2, cmap=cm, s=z1)
plt.errorbar(age_ocaat, age, xerr=e_age, ls='none', color='grey',
    elinewidth=0.8)
# 1:1 line
plt.plot([6., 7., 8., 9., 10.], [6., 7., 8., 9., 10.], 'k-', ls='--')
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax3 = plt.subplot(gs[13:16, 0:3])
#plt.xlim(-1.5, 1.5)
plt.xlabel('$\Delta dist\,mod$', fontsize=14)
plt.ylabel('CI', fontsize=12)
ax3.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
x = np.take(delta_dist, order)
y = np.take(ci, order)
plt.scatter(x, y, c=z2, cmap=cm, s=z1)
#plt.scatter(delta_dist, ci, c=dist, cmap=cm, s=((np.array(mass) / 40.) ** 2))
plt.errorbar(delta_dist, ci, xerr=e_dist, ls='none', color='k', elinewidth=0.8)
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax13 = plt.subplot(gs[13:16, 4:7], aspect=1)
plt.xlim(6., 16.)
plt.ylim(6., 16.)
plt.xlabel('$DM_{ocaat}$', fontsize=14)
plt.ylabel('$DM_{real}$', fontsize=14)
ax13.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
x = np.take(dist_ocaat, order)
y = np.take(dist_mod, order)
plt.scatter(x, y, c=z2, cmap=cm, s=z1)
plt.errorbar(dist_ocaat, dist_mod, xerr=e_dist, ls='none', color='grey',
    elinewidth=0.8)
# 1:1 line
plt.plot([6., 8., 10., 12., 16.], [6., 8., 10., 12., 16.], 'k-', ls='--')
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax4 = plt.subplot(gs[17:20, 0:3])
plt.xlim(-1.5, 1.5)
plt.xlabel('$\Delta E_{(B-V)}$', fontsize=14)
plt.ylabel('CI', fontsize=12)
ax4.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
# Order before plotting.
x = np.take(delta_ext, order)
y = np.take(ci, order)
plt.scatter(x, y, c=z2, cmap=cm, s=z1)
#plt.scatter(delta_ext, ci, c=dist, cmap=cm, s=((np.array(mass) / 40.) ** 2))
plt.errorbar(delta_ext, ci, xerr=e_ext, ls='none', color='k', elinewidth=0.8)
cbar = plt.colorbar()
cbar.set_ticks([0.5, 1., 3., 5.])
cbar.set_ticklabels([0.5, 1., 3., 5.])
cbar.set_label('dist (kpc)')

ax14 = plt.subplot(gs[17:20, 4:7], aspect=1)
plt.xlim(-0.3, 1.2)
plt.ylim(-0.3, 1.2)
plt.xlabel('$E(B-V)_{ocaat}$', fontsize=14)
plt.ylabel('$E(B-V)_{real}$', fontsize=14)
ax14.grid(b=True, which='both', color='gray', linestyle='--', lw=0.5)
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