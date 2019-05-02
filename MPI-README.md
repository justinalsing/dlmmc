**Running in parallel with MPI**

It's often necessary to perform regression of a large number time-series (eg., over a grid of observations at different altitudes/latitudes/longitudes) and is advantageous to be able to run these in parallel. If you want to run multiple DLMs in parallel with MPI, you will need to install [openmpi](https://www.open-mpi.org) and [mpi4py](https://mpi4py.readthedocs.io/en/stable/install.html) (again easily done with `conda install`).

The python script `dlm_lat_alt_mpi_run.py` is a template for how to run the DLM code over a grid of time-series at different latitudes/altitudes in parallel using MPI, and save the results to a single netCDF file. This script has the additional dependency [tqdm](https://tqdm.github.io) if you want it to work with a progress bar (`conda install tqdm`). Provided you have MPI working, you can run this script with the following command (using eg. 4 hyperthreaded processes):

`mpirun --use-hwthread-cpus -np 4 python3 dlm_lat_alt_mpi_run.py`

I recommend you run this with a very small number of samples first (eg set `iter=3` and `warmup=1` at the top of the python script `dlm_lat_alt_mpi_run.py`) to check it runs through, before embarking on a long run. Note: this will run very short MCMC chains (only 3 samples!) so stan will throw lots of warnings about lack of convergence, which you can ignore for this test - as long as you get no errors, you're good to go.
