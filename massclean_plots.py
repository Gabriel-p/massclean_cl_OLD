
from os import getcwd
from os.path import join, realpath, dirname
import numpy as np
import matplotlib.pyplot as plt


'''
Generate plot of field and cluster stars to check the files produced.
'''

# Path where the code is running
mypath = realpath(join(getcwd(), dirname(__file__)))


# Cluster data file.
cluster_file = 'is1_p019_0650.trek'
# Open cluster data file and read the data.
data = np.loadtxt(join(mypath, cluster_file), unpack=True)

# Store values.        
mag_data, col1_data = data[7], (data[6]-data[7])

        
# Read data from field stars file.
field_file = '2500_field.plot'
data2 = np.loadtxt(join(mypath, field_file), unpack=True)

# Store values.
mag_field, col1_field = data2[2], (data2[1]-data2[2])


# Merge cluster and field.
mag_cl_fl = mag_data.tolist() + mag_field.tolist()
col1_cl_fl = col1_data.tolist() + col1_field.tolist()


# Plot CMD.
plt.figure(figsize=(10, 15))
plt.ylim(max(mag_cl_fl)+1, min(mag_cl_fl)-1)
plt.scatter(col1_field, mag_field, marker='o', c='k', 
            s=1, zorder=1)        
plt.scatter(col1_data, mag_data, marker='x', c='r', 
            s=15, zorder=2)

plt.savefig(join(mypath, str(cluster_file[:-5])+'.png'), dpi=150)

                
print 'End.'