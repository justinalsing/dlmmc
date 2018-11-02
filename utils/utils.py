# Import stuff
import netCDF4
import os
import numpy as np

class suppress_stdout_stderr(object):
    '''
    A context manager for doing a "deep suppression" of stdout and stderr in
    Python, i.e. will suppress all print, even if the print originates in a
    compiled C/Fortran sub-function.
       This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).

    '''
    def __init__(self):
        # Open a pair of null files
        self.null_fds = [os.open(os.devnull, os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = (os.dup(1), os.dup(2))

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0], 1)
        os.dup2(self.null_fds[1], 2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0], 1)
        os.dup2(self.save_fds[1], 2)
        # Close the null files
        os.close(self.null_fds[0])
        os.close(self.null_fds[1])


# Create a netCDF file saving all the results
def create_results_netcdf(results_filename, L, P, T, n_samples, nregs):
    
    # Make the netcdf file
    results = netCDF4.Dataset(results_filename, 'w')

    # Create dimensions
    results.createDimension('latitude', len(L))
    results.createDimension('pressure', len(P))
    results.createDimension('time', len(T))
    results.createDimension('nsamples', n_samples)
    results.createDimension('nregs', nregs)

    # Create variables
    LATS = results.createVariable('latitude', int, ('latitude',), zlib=True)
    PRES = results.createVariable('pressure', int, ('pressure',), zlib=True)
    TIME = results.createVariable('time', int, ('time',), zlib=True)

    TREND_SAMPLES = results.createVariable('trend', float, ('nsamples', 'pressure', 'latitude', 'time'), zlib=True)

    TREND_MEAN = results.createVariable('trend_mean', float, ('pressure', 'latitude', 'time'), zlib=True)
    TREND_STD = results.createVariable('trend_std', float, ('pressure', 'latitude', 'time'), zlib=True)

    SEASONAL_MEAN = results.createVariable('seasonal_mean', float, ('pressure', 'latitude', 'time'), zlib=True)
    SEASONAL_STD = results.createVariable('seasonal_std', float, ('pressure', 'latitude', 'time'), zlib=True)

    RESIDUALS_MEAN = results.createVariable('residuals_mean', float, ('pressure', 'latitude', 'time'), zlib=True)
    RESIDUALS_STD = results.createVariable('residuals_std', float, ('pressure', 'latitude', 'time'), zlib=True)

    SLOPE_MEAN = results.createVariable('slope_mean', float, ('pressure', 'latitude', 'time'), zlib=True)
    SLOPE_STD = results.createVariable('slope_std', float, ('pressure', 'latitude', 'time'), zlib=True)

    REGRESSOR_COEFFICIENTS = results.createVariable('regressor_coefficients', float, ('nsamples', 'pressure', 'latitude', 'nregs'), zlib=True)

    SIGMA_TREND = results.createVariable('sigma_trend', float, ('nsamples', 'pressure', 'latitude'), zlib=True)
    SIGMA_SEAS = results.createVariable('sigma_seas', float, ('nsamples', 'pressure', 'latitude'), zlib=True)
    SIGMA_AR = results.createVariable('sigma_AR', float, ('nsamples', 'pressure', 'latitude'), zlib=True)
    RHO_AR = results.createVariable('rho_AR', float, ('nsamples', 'pressure', 'latitude'), zlib=True)

    results.close()

        
# Put the relevant bits into the netCDF...
def add_results_to_netcdf(results_filename, fit, pressure, latitude):

    # Open the netCDF
    results = netCDF4.Dataset(results_filename, 'r+')

    # Trend
    results['trend'][:, pressure, latitude, :] = fit.extract()['trend'][:,:]
    results['trend_mean'][pressure, latitude, :] = np.mean(fit.extract()['trend'][:,:], axis = 0)
    results['trend_std'][pressure, latitude, :] = np.std(fit.extract()['trend'][:,:], axis = 0)
    
    # Slope of the DLM trend
    results['slope_mean'][pressure, latitude, :] = np.mean(fit.extract()['slope'][:,:], axis = 0)
    results['slope_std'][pressure, latitude, :] = np.std(fit.extract()['slope'][:,:], axis = 0)

    # Seasonal cycle
    results['seasonal_mean'][pressure, latitude, :] = np.mean(fit.extract()['seasonal'][:,:], axis = 0)
    results['seasonal_std'][pressure, latitude, :] = np.std(fit.extract()['seasonal'][:,:], axis = 0)
    
    # Residuals
    results['residuals_mean'][pressure, latitude, :] = np.mean(fit.extract()['residuals'][:,:], axis = 0)
    results['residuals_std'][pressure, latitude, :] = np.std(fit.extract()['residuals'][:,:], axis = 0)
    
    # Regressor coefficients
    results['regressor_coefficients'][:, pressure, latitude, :] = fit.extract()['beta'][:,:]

    # DLM hyper parameters
    results['sigma_trend'][:, pressure, latitude] = fit.extract()['sigma_trend']
    results['sigma_seas'][:, pressure, latitude] = fit.extract()['sigma_seas']
    results['sigma_AR'][:, pressure, latitude] = fit.extract()['sigma_AR']
    results['rho_AR'][:, pressure, latitude] = fit.extract()['sigma_AR']
    
    # Close the netCDF file
    results.close()

