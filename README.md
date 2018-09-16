# dlm

**Dynamical Linear Modelling** (DLM) regression code for analysis of time-series data. The code is targeted at atmospheric time-series analysis, with a worked example (and data) included for ozone, but is a fairly general DLM code that can be applied to a range of problems.

There is a paper in preparation to accompany the code (appearing soon). Until then, please contact us if you intend to use the code for a major project.

A close description of the DLM model implemented here can be found in [Laine et al 2014](https://www.atmos-chem-phys.net/14/9707/2014/acp-14-9707-2014.pdf), and the model/code was used for analyzing ozone data in [Ball et al 2017](https://www.research-collection.ethz.ch/handle/20.500.11850/202027) and [Ball et al 2018](https://www.atmos-chem-phys.net/18/1379/2018/acp-18-1379-2018.html). Please consider citing these papers along with the paper accompanying this code (in prep; appearing soon) when you use this code. 

**Dependencies** The code is Python 3. Install the following dependencies (using eg., pip install):

[numpy](http://www.numpy.org)
[scipy](https://www.scipy.org)
[matplotlib](https://matplotlib.org)
[netCDF4](https://pypi.org/project/netcdf/)
[pystan](https://pystan.readthedocs.io/en/latest/)

**Usage** A detailed annotated tutorial walk-through of how to use the code is given in the jupyter notebook `dlm_kalman.ipynb` -- this tutorial analyses ozone time-series data as a case study. Please work through this notebook and contact us if you have any questions.
