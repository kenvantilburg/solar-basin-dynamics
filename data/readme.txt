The files in constraints/ are the experimental constraints
used in code/Summary/constraints.ipynb.

The files in diffusion/ are tabulated values of the 
diffusivity and the probability for a particle to scatter 
between two values of semimajor axis in some time window 
dt in some semimajor axis range d(∆a).  These expressions 
are implicitly defined in the companion paper, in Eqs. 5.6 
and 5.12.

Earth_hits/ contains data for earth hits from forward 
simulations.  To use this data in 
code/Summary/constraints.ipynb, unzip all files in 
Earth_hits and merge them into one file with
"copy hit_record_* hit_record.csv".

The files in paperplots/ are used in the Mathematica
notebook code/PaperPlots.nb--see that file for the usage
of each data file.

hist.npy is an output of code/secular_pt.ipynb: a 3D 
(histogram) array (x,y,z) of the time-averaged solar basin 
density.

Q_DP_solar.csv is the energy loss rate Q for dark photons in 
the sun.  It is tabualated in code/Summary/constraints.ipynb.

struct_b16_agss09.dat contains information about the solar
model used in code/Summary/constraints.ipynb.


