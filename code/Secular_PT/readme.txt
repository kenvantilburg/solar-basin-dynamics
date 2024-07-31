Generate data by running the command
    python3 secular_pt_run.py
many times (in our implementation, more than 10^4 times) with a script, and saving the output data hist_strap_XXXXXXXXX.npy in an appropriate directory.
The notebook secular_pt.ipynb processes these files, and generates:
-- hist.npy: a 3D (histogram) array (x,y,z) of the time-averaged solar basin density
-- hist_strap.npy: a list of 2D (histogram) array (x,y,z=0) of the time-averaged solar basin density of each hist_strap_XXXXXXXXX.npy file
-- dens_sample_strap.npy: a list of 1D histograms of density along Earth's orbit for each hist_strap_XXXXXXXXX.npy file
-- dens_sample_fft_strap.npy: list of discrete Fourier transforms of the former
Among other products, the same notebook then generates the figures dens_DFT_M.pdf and dens_CL_M.pdf of the main text.
