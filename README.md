MASSCLEAN cluster generator
============

The bash script `massclean_run_all.sh` calls the [MASSCLEAN](http://www.physics.uc.edu/~bogdan/massclean/)
package to generate a number of open clusters along with fields of stars created with the same parameter values.

Within the script the metallicities, ages, distances, extinction and initial masses of the
clusters can be set. The `MASSCLEAN` output files will be stored in the `/massclean2.013/clusters`
folder.

Each sub-folder respects the following naming convention: (initial mass)/100 for the first two characters
and distance in parsecs for the characters after the underscore.

After the synthetic clusters and field plots are generated the script calls the `massclean_cluster_field_merge.py`
code to merge cluster and field while adding errors and removing stars mimicking completeness in the frame.
The final merged `.DAT` files are stored in `/synth_clusters/` where images for each synthetic cluster
are generated.



3) ** massclean_member_index.py **
Calculates all member_index defined, for a given decontamination algorithm. Plots the output as a
CI vs MI diagram for each member index defined.


4) ** massclean_avrg_MI_algors.py ** SEMI-DEPRECATED
Plots the average member_index values for all the decontamination algorithms used, for each of the member indexes defined.


5) ** massclean_real_vs_fit_params.py **
Plots the real parameter values (z, E(B-V), age, dist) of MASSCLEAN clusters versus the ones obtained with the
automatic isochrone fitting algorithm.



Los parámetros CI y MI son respectivamente el 'contamination index' y el 'member index'. El primero estima la contaminación de estrellas
de campo presente en la región del cúmulo (definida por el centro y radio asignados) e implica un cúmulo más contaminado para valores
más cercanos a 1. El segundo es un estimador de cuán bueno es el comportamiento del algoritmo de decontaminación aplicado, un valor de
1 implica un funcionamiento perfecto del algoritmo (o sea: sólo los miembros reales del cúmulo fueron detectados como tales con un 
peso o probabilidad de 100% asignado a todos ellos). A medida que este parámetro decrece implica un funcionamiento más pobre del algoritmo.
Se aplicaron 3 MI diferentes, cada uno más estricto que el anterior.

Las variables usadas para definir los 3 MI son:

N_T: cantidad total *real* de estrellas miembros del cúmulo sintético (o sea: cantidad de estrellas en el cúmulo generado por MASSCLEAN).
n_m: cantidad de *miembros reales* que el algoritmo asignó como tales.
n_f: cantidad de *estrellas de campo* que el algoritmo asignó como miembros reales.
w_m: peso o probabilidad de pertenecer al cúmulo (entre 0 y 1) asignado a un miembro real detectado como tal.
w_f: peso o probabilidad de pertenecer al cúmulo (entre 0 y 1) asignado a una estrella de campo detectada como miembro real.

En general los 3 MI pueden describirse como sigue:

MI_1: solamente tiene en cuenta la proporción de miembros reales recuperados vs la cantidad total real de miembros. Para este
índice, un valor de 0.5 significa que la mitad de los miembros reales fueron recuperados (sin tener en cuenta los pesos dados
a cada uno ni la cantidad de estrellas de campo asignadas como miembros).
MI_2: tiene en cuenta los pesos asignados a cada estrella y castiga al algoritmo cuando asigna probabilidades altas a estrellas
de campo y bajas a miembros reales.
MI_3: igual a MI_2 con el agregado de que también castiga cuando el número *total* de estrellas asignadas como miembros del cúmulo
(n_m + n_f) se aleja del número real N_T.