MASSCLEAN cluster generator
============

The bash script `massclean_run_all.sh` calls the [MASSCLEAN](http://www.physics.uc.edu/~bogdan/massclean/)
package to generate a number of open clusters along with fields of stars created with the same parameter values.

Within the script the metallicities, ages, distances, extinction and initial masses of the
clusters can be set. The *MASSCLEAN* output files will be stored in the `/massclean2.013/clusters/`
folder.

Each sub-folder respects the following naming convention: `(initial mass)/100` for the first two characters,
distance in parsecs for the characters after the first underscore and visual absortption after the
second underscore.

After the synthetic clusters and field plots are generated the script calls the `massclean_cluster_field_merge.py`
code to merge cluster and field while adding errors and removing random stars in the frame, mimicking
completeness effects.
The final merged `.DAT` files are stored in `/synth_clusters/` where images for each synthetic cluster
are also generated.


The `ci_plots.py` script compares the real *MASSCLEAN* cluster parameters with those obtained by the
OCAAT code versus the contamination index (CI) obtained for each synthetic cluster.

The `massclean_member_index.py` calculates the members indexes (MI) defined and plots the output as a
CI vs MI diagram for each one for all synthetic clusters.

`massclean_avrg_MI_algors.py` is SEMI-DEPRECATED. It plots the average member_index values for all 
the decontamination algorithms used, for each of the member indexes defined.
