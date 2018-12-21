# DLM

**Dynamical Linear Modelling** (DLM) regression code in python for analysis of time-series data. The code is targeted at atmospheric time-series analysis, with a detailed worked example (and data) included for stratospheric ozone, but is a fairly general state space model that can be applied or extended to a wide range of problems.

The core of this package is a suite of DLM models (implemented in [stan](https://mc-stan.org)), using a combination of efficient HMC sampling and Kalman filtering to infer the DLM model parameters (trend, seasonal cycle, auto-regressive processes etc) given some data. To make the code as accessible as possible, I provide a detailed easy-to-follow tutorial in python for how to read in your data and run the DLM model(s) implemented in this package on your own data.

### Dependencies and installation

The code is python3 and has the following dependencies (which can be installed using eg., pip install):

[numpy](http://www.numpy.org)<br/>
[scipy](https://www.scipy.org)<br/>
[matplotlib](https://matplotlib.org)<br/>
[netCDF4](https://pypi.org/project/netcdf/)<br/>
[pystan](https://pystan.readthedocs.io/en/latest/)<br/>

If you want to run multiple DLMs in parallel with MPI, you will also need [openmpi](https://www.open-mpi.org) and [mpi4py](https://mpi4py.readthedocs.io/en/stable/install.html) (again easily done with pip).

Once you have downloaded the code from this repository and installed the dependencies, run the following script (make sure in python3):

`python compile_stan_models.py`

This pre-compiles all of the models on your machine, saves them in `models/`, and then you're ready to start DLMing!

### Usage

**Functionality**

A detailed annotated tutorial walk-through of how to use the code is given in the jupyter notebook `dlm_tutorial.ipynb` -- this tutorial analyses stratospheric ozone time-series data as a case study. The notebook takes you step-by-step through the complete functionality of the code: loading in your own data, running the DLM model, and processing and plotting the results.  

**Running in parallel with MPI**

It's often necessary to perform regression of a large number time-series (eg., over a grid of observations at different altirudes/latitudes/longitudes) and is advantageous to be able to run these in parallel. The python script `dlm_lat_alt_mpi_run.py` is a template for how to run the DLM code over a grid of time-series at different latitudes/altitudes in parallel using MPI, and save the results to a netCDF file. This script has the additional dependency [tqdm](https://tqdm.github.io) if you want it to work with a progress bar. Provided you have MPI working, you can run this script with the following command (using eg. 4 hyperthreaded processes, again make sure you run with python3):

`mpirun --use-hwthread-cpus -np 4 python dlm_lat_alt_mpi_run.py`

I recommend you run this with a very small number of samples first (eg iter=3, warmup=1) to check it runs through, before embarking on a long run.

**Model descriptions**

Mathematical descriptions of each of the DLM models implemented in this package can be found in the file `models/model_descriptions.pdf`. This file contains a concise description of the parameters of each model, their physical meanings, and how to refer to them in the code: make sure you have read and understand the model description before running a new model!

### Citing this code

There is a JOSS paper in preparation to accompany the code (appearing soon). Until then, please contact me if you intend to use the code for a major project (justin.alsing@fysik.su.se).

A reasonably close description of the vanilla DLM model implemented here can be found in [Laine et al 2014](https://www.atmos-chem-phys.net/14/9707/2014/acp-14-9707-2014.pdf), and this model/code was used for analyzing ozone data in [Ball et al 2017](https://www.research-collection.ethz.ch/handle/20.500.11850/202027) and [Ball et al 2018](https://www.atmos-chem-phys.net/18/1379/2018/acp-18-1379-2018.html). Please consider citing these papers along with the paper accompanying this code (in prep - appearing soon) when you use this code.

### Contributions, reporting issues and feature requests

If you want to contribute (eg., extended models) to this package, please contact me at justin.alsing@fysik.su.se (I welcome contributors). If you would like to request a feature to be included in the next update, or report an issue, please use the issues channel associated with this Git repository.

