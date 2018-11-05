# dlm

**Dynamical Linear Modelling** (DLM) regression code in python for analysis of time-series data. The code is targeted at atmospheric time-series analysis, with a worked example (and data) included for ozone, but is a fairly general state space model that can be applied to a range of problems.

**Dependencies and installation** 

The code is Python 3. Install the following dependencies (using eg., pip install):

[numpy](http://www.numpy.org)
[scipy](https://www.scipy.org)
[matplotlib](https://matplotlib.org)
[netCDF4](https://pypi.org/project/netcdf/)
[pystan](https://pystan.readthedocs.io/en/latest/)

If you want to run multiple DLMs in parallel with MPI, you will also need [openmpi](https://www.open-mpi.org) and [mpi4py](https://mpi4py.readthedocs.io/en/stable/install.html) (again easily done with pip).

Once you have downloaded the code from this repository and installed the dependencies, run the following script (make sure in python3):

`python compile_stan_models.py`

This compiles the all models on your machine, then you're ready to start DLMing!

**Usage** 

A detailed annotated tutorial walk-through of how to use the code is given in the jupyter notebook `dlm_tutorial.ipynb` -- this tutorial analyses ozone time-series data as a case study. Please work through this notebook to get to grips with the code and contact me at justin.alsing@fysik.su.se if you have any questions.

The python script `dlm_lat_alt_mpi_run.py` is a template for how to run the DLM code over a grid of time-series at different latitudes/altitudes and save the results to a netCDF file, in parallel using MPI. This script has the additional dependency [tqdm](https://tqdm.github.io) if you want it to work with a progress bar. Once you have MPI working, you can run this script with the command (using eg. 4 hyperthreaded processes, again make sure you run with python3)

`mpirun --use-hwthread-cpus -np 4 python dlm_lat_alt_mpi_run.py`

I recommend you run this with a very small number of samples first (eg iter=3, warmup=1) to check it runs through, before embarking on a long run.

**Citing this code**

There is a paper in preparation to accompany the code (appearing soon). Until then, please contact us if you intend to use the code for a major project (justin.alsing@fysik.su.se).

A reasonably close description of the DLM model implemented here can be found in [Laine et al 2014](https://www.atmos-chem-phys.net/14/9707/2014/acp-14-9707-2014.pdf), and the model/code was used for analyzing ozone data in [Ball et al 2017](https://www.research-collection.ethz.ch/handle/20.500.11850/202027) and [Ball et al 2018](https://www.atmos-chem-phys.net/18/1379/2018/acp-18-1379-2018.html). Please consider citing these papers along with the paper accompanying this code (in prep - appearing soon) when you use this code. 

