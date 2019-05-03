[![DOI](http://joss.theoj.org/papers/10.21105/joss.01157/status.svg)](https://doi.org/10.21105/joss.01157)<br/>
[![DOI](https://zenodo.org/badge/148547735.svg)](https://zenodo.org/badge/latestdoi/148547735)

# DLMMC

**Dynamical Linear Modelling** (DLM) regression code in python for analysis of time-series data. The code is targeted at atmospheric time-series analysis, with a detailed worked example (and data) included for stratospheric ozone, but is a fairly general suite of state space model that can be applied or extended to a wide range of problems.

The core of this package is a suite of DLM models implemented in [stan](https://mc-stan.org), using a combination of HMC sampling and Kalman filtering to infer the DLM model parameters (trend, seasonal cycle, auto-regressive processes etc) given some time-series data. To make the code as accessible as possible, I provide a step-by-step tutorial in python for how to read in your data, run the DLM model(s), and process the outputs to make nice plots. Once you've worked through this tutorial you should have all the tools you need to apply DLM to your own data!

Note: some basic working knowledge of python and jupyter notebooks is required to make the most of this package (but it's not too onerous).

### Installation

Once you have downloaded the code from this repository you're ready to install dependencies and get set-up.

The code is python3 and has the following dependencies: [numpy](http://www.numpy.org), [scipy](https://www.scipy.org), [matplotlib](https://matplotlib.org), [jupyter](https://jupyter.org/install), [ipython](https://ipython.org/install.html), [netCDF4](https://pypi.org/project/netcdf/), [pystan](https://pystan.readthedocs.io/en/latest/).

**Installation with conda (recommended)**

The most painless way to get set up is using the [Anaconda python distribution](https://www.anaconda.com/distribution/) (recommended), which comes with most of the dependencies as default. The remaining dependencies can then be installed using `conda install` and the DLM models compiled by running:

```
conda install pystan netCDF4
python3 compile_stan_models.py
```

This second line compiles all of the DLM models on your machine, saves them in `models/`, and then you're ready to start DLMing! Jump straight into the jupyter notebook tutorial `dlm_tutorial.ipynb` (see below), or if you prefer you can run a test suite to check that the install worked and all models run smoothly by executing (this will take some minutes to run through):

```
jupyter-nbconvert --to notebook --execute --ExecutePreprocessor.timeout=100000 dlm_validation_tests.ipynb
```

Finally, if you want to see what a successful installation looks like, see `INSTALL.md`.

**Installation with pip (at your own risk)**

Anaconda is not a _requirement_ for installing dlmmc, but is recommended because it works robustly with pystan. If you would rather use a different python distribution and `pip3` for installing dependencies, you are welcome to (at your own risk); see the [pystan readthedocs](https://pystan.readthedocs.io/en/latest/installation_beginner.html) for advice on installing pystan using `pip3` if you run into problems. Note that if you do not use Anaconda you will also have to install the other dependencies listed above, ie., 

```
pip3 install numpy scipy ipython[all] jupyter matplotlib netCDF4 pystan
python3 compile_stan_models.py
```

**Platforms** 

dlmmc has been successfully installed on Mac, Linux and Windows. Note that there are some limitations to the functionality of [pystan on Windows](https://pystan.readthedocs.io/en/latest/windows.html), but these do not restrict the use of the dlmmc package for Windows users.

### Usage

**Functionality**

A detailed annotated tutorial walk-through of how to use the code is given in the jupyter notebook `dlm_tutorial.ipynb` -- this tutorial analyses stratospheric ozone time-series data as a case study. The notebook takes you step-by-step through the complete functionality of the code: loading in your own data, running the DLM model, and processing and plotting the results. The tutorial also serves as a test that the install was successful and the compiled models run smoothly (for a more comprehensive test suite see below).

**Running in parallel with MPI**

It's often necessary to perform regression of a large number time-series (eg., over a grid of observations at different altitudes/latitudes/longitudes) and is advantageous to be able to run these in parallel. Although not a central part of this package, I provide a template example for doing large MPI runs in `dlm_lat_alt_mpi_run.py` - see `MPI-README.md` for a description of getting set up with MPI runs.

**Model descriptions**

Mathematical descriptions of each of the DLM models implemented in this package can be found in the file `models/model_descriptions/model_descriptions.pdf`. This file contains a concise description of the parameters of each model, their physical meanings, and how to refer to them in the code: make sure you have read and understand the model description before running a new model!

**Test suite and code validation**

A more comprehensive test suite is provided in `dlm_validation_tests.ipynb`. In this notebook I run through the suite of DLM models in dlmmc, generating mock data and running the DLM on those mock data, for each model in turn. This acts as both a test suite to check the install has worked robustly (ie., all of the models run to completion without error), and also serves as a set of validation tests demonstrating that the input parameters are recovered correctly, within posterior uncertainties, for each model.

### Citing this code

There is a JOSS paper to accompany this paper, which can be found [here](http://joss.theoj.org/papers/10.21105/joss.01157). Please cite this paper if you use or refer to this code, as:

_Alsing, (2019). dlmmc: Dynamical linear model regression for atmospheric time-series analysis. Journal of Open Source Software, 4(37), 1157, https://doi.org/10.21105/joss.01157_

A close description of the vanilla DLM model implemented here can be found in [Laine et al 2014](https://www.atmos-chem-phys.net/14/9707/2014/acp-14-9707-2014.pdf), and this model/code was used for analyzing ozone data in [Ball et al 2017](https://www.research-collection.ethz.ch/handle/20.500.11850/202027) and [Ball et al 2018](https://www.atmos-chem-phys.net/18/1379/2018/acp-18-1379-2018.html). Please consider citing these papers too when you use this code.

### Contributions, reporting issues and feature requests

If you want to contribute eg., extended models to this package, I provide some instructions for how to develop your own model within the dlmmc framework in the `CONTRIBUTING` file. If you are interested in collaborating on extended models, feel free to contact me at justin.alsing@fysik.su.se. If you would like to request a feature to be included in the next update, or report an issue, please use the issues channel associated with this Git repository.

