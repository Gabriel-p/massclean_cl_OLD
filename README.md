MASSCLEAN cluster generator
============

The set of bash/python scripts in this repository are intended to generate synthetic
clusters making use of the [MASSCLEAN](http://www.physics.uc.edu/~bogdan/massclean/)
package and to compare their real parameters (center, radius, metallicity, age,
distance, extinction) with those obtained by the [OCAAT](https://github.com/Gabriel-p/ocaat) code for each one.

1st step
---------
The bash script `massclean_run_all.sh` calls *MASSCLEAN* to generate a number of synthetic open clusters along with fields of stars created with the same parameter values.

Within the script the radius, metallicities, ages, distances, extinction and initial masses of the clusters can be set. The *MASSCLEAN* output files will be stored in the `/massclean2.013/clusters/` folder.

Each sub-folder respects the following naming convention: `(initial mass)/100` for the first two characters, distance in parsecs for the characters after the first underscore and visual absorption after the second underscore.

2nd step
---------
After the synthetic clusters and field plots are generated the script calls the `massclean_cluster_field_merge.py` code to merge cluster and field while adding errors and removing random stars in the frame, mimicking completeness effects.
The final merged `.DAT` files are stored in `/synth_clusters/` where images for each synthetic cluster are also generated.


3rd step
---------
The `ci_plots.py` script compares the real *MASSCLEAN* cluster parameters with those obtained by the *OCAAT* code versus the contamination index (CI) obtained for each synthetic cluster.


4th step
---------
This `massclean_member_index.py` script is used to establish the accuracy of the decontamination algorithm applied on the synthetic clusters vy *OCAAT*. It calculates the members indexes (MI) defined (which estimates how good real cluster members were recovered by the algorithm) and plots the output as a CI vs MI diagram for all synthetic clusters.


Semi-deprecated
---------
`massclean_avrg_MI_algors.py` is SEMI-DEPRECATED. It plots the average member_index values for all the decontamination algorithms used, for each of the member indexes defined.
