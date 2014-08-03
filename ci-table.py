
import numpy as np
import bisect


def mod_to_kpc(d_mod):
    '''
    Convert distance modulus to kpc.
    '''
    return (10 ** ((np.array(d_mod) + 5.) / 5.)) / 1000.


def skip_comments(f):
    '''
    Read lines that DO NOT start with a # symbol.
    '''
    for line in f:
        if not line.strip().startswith('#'):
            yield line


# Read ocaat_output.dat file to obtain the MASSCLEAN clusters actual
# data and parameters.
out_file = '/media/rest/github/ocaat/output/massclean/ocaat_output.dat'
f = open(out_file)
names, params = [], []
for line in skip_comments(f):
    names.append(line.split()[0])
    params.append(line.split()[1:])

# Separate into masses, distances, metallicities and ages of MASSCLEAN
# clusters.
mass, dist, absor, metal, age = [], [], [], [], []
for clust_str in names:
    first, second = clust_str.split('/')
    first_a, first_b, first_c = first.split('_')
    #mass.append(float(first_a))
    dist.append(float(first_b) / 1000.)
    absor.append(float(first_c))
    age.append(float(second[9:]) / 100.)
    metal.append(float('0.' + second[5:8]))

# Read clusters parameters obtained by OCAAT.
ci, metal_ocaat, age_ocaat, dist_ocaat, ext_ocaat = [], [], [], [], []
for par_str in params:
    ci.append(float(par_str[9]))
    dist_ocaat.append(float(par_str[21]))
    #e_dist.append(float(par_str[22]))
    ext_ocaat.append(float(par_str[19]) * 3.1)
    #e_ext.append(float(par_str[20]))
    age_ocaat.append(float(par_str[17]))
    #e_age.append(float(par_str[18]))
    metal_ocaat.append(float(par_str[15]))
    #e_met.append(float(par_str[16]))


'''
Create the table(s) used in the OCAAT article summarizing the results
obtained in the validation process.

Steps:

1- Each SOC is stored in a list according to the value of the parameter
used to generate it, one list for each parameter (dist, ext, age, met)

2- These lists are separated in sub-lists according to the CI values
each one of those SOCs has, for given CI ranges.

3- The mean and standard deviation of the values assigned by OCAAT to all the
SOCs with a given parameter value in a given CI range are obtained.

'''

# Hardcoded values.
# CI ranges.
ci_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
# True parameter values used to generate the MASSCLEAN SOCs.
d_values = [0.5, 1.0, 3., 5.]
e_values = [0.1, 0.5, 1., 3.]
a_values = [7., 8., 9.]
m_values = [0.002, 0.008, 0.019, 0.03]

# Define empty CI list.
ci_lst = [[[[] for _ in range(len(d_values))],
[[] for _ in range(len(e_values))], [[] for _ in range(len(a_values))],
[[] for _ in range(len(m_values))]] for _ in range((len(ci_values) + 1))]

# Store indexes pointing to the appropriate sub-list in ci_lst for each
# ci value.
ci_indxs = [bisect.bisect(ci_values, x) for x in ci]

# Store OCAAT parameters values.
for indx, ci_ind in enumerate(ci_indxs):

    params = [dist_ocaat[indx], ext_ocaat[indx], age_ocaat[indx],
        metal_ocaat[indx]]
    # Skip clusters where the params could not be assigned.
    if not any(i < 0. for i in params):

        # Index of dist sub-list.
        params_ind = [d_values.index(dist[indx]), e_values.index(absor[indx]),
            a_values.index(age[indx]), m_values.index(metal[indx])]

        for indx2, p_ind in enumerate(params_ind):
            # Store OCAAT values in list.
            ci_lst[ci_ind][indx2][p_ind].append(params[indx2])

# Print results.
for ci_rang in ci_lst:

    # Dist.
    d_mean, d_std = [], []
    for d_vals in ci_rang[0]:
        if d_vals:
            d_kpc = mod_to_kpc(d_vals)
            d_kpc_m, d_kpc_st = np.mean(d_kpc), np.std(d_kpc)
            d_mean.append(round(d_kpc_m, 1))
            d_std.append(round(d_kpc_st, 1))
        else:
            d_mean.append('-')
            d_std.append('')
    # Format mean and std dev values in a single list.
    f_lst = [item for sublist in zip(d_mean, d_std) for item in sublist]
    print 'd: ${}\pm{}$ & ${}\pm{}$ & ${}\pm{}$ & ${}\pm{}$ &'.format(*f_lst)

    # Absorption/Extinction
    e_mean, e_std = [], []
    for e_vals in ci_rang[1]:
        if e_vals:
            e_m, e_st = np.mean(e_vals), np.std(e_vals)
            e_mean.append(round(e_m, 1))
            e_std.append(round(e_st, 1))
        else:
            e_mean.append('-')
            e_std.append('')
    # Format mean and std dev values in a single list.
    f_lst = [item for sublist in zip(e_mean, e_std) for item in sublist]
    print 'e: ${}\pm{}$ & ${}\pm{}$ & ${}\pm{}$ & ${}\pm{}$ &'.format(*f_lst)

    # log(Age)
    a_mean, a_std = [], []
    for a_vals in ci_rang[2]:
        if a_vals:
            a_m, a_st = np.mean(a_vals), np.std(a_vals)
            a_mean.append(round(a_m, 1))
            a_std.append(round(a_st, 1))
        else:
            a_mean.append('-')
            a_std.append('')
    # Format mean and std dev values in a single list.
    f_lst = [item for sublist in zip(a_mean, a_std) for item in sublist]
    print 'a: ${}\pm{}$ & ${}\pm{}$ & ${}\pm{}$ &'.format(*f_lst)

    # Metallicity (z)
    m_mean, m_std = [], []
    for m_vals in ci_rang[3]:
        if m_vals:
            m_m, m_st = np.mean(m_vals), np.std(m_vals)
            m_mean.append(round(m_m, 3) * 100.)
            m_std.append(round(m_st, 3) * 100.)
        else:
            m_mean.append('-')
            m_std.append('')
    # Format mean and std dev values in a single list.
    f_lst = [item for sublist in zip(m_mean, m_std) for item in sublist]
    print 'mx1e02: ${}\pm{}$ & ${}\pm{}$ & ${}\pm{}$ & ${}\pm{}$'.format(*f_lst)

    print '\n'