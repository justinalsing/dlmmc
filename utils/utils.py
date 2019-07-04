# Import stuff
import netCDF4
import os
import numpy as np
import pickle

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
        self.save_fds = [os.dup(1), os.dup(2)]

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0], 1)
        os.dup2(self.null_fds[1], 2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0], 1)
        os.dup2(self.save_fds[1], 2)
        # Close the null files
        for fd in self.null_fds + self.save_fds:
            os.close(fd)

# Create a netCDF file saving all the results
def create_results_netcdf(results_dir, results_filename, L, P, T, n_samples, nregs):
    
    # Try to remove the netCDF if it already exists
    try:
        os.remove('{}/{}.nc'.format(results_dir, results_filename))
    except:
        pass
    
    # Make the netcdf file
    results = netCDF4.Dataset('{}/{}.nc'.format(results_dir, results_filename), 'w')

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

    SLOPE_MEAN = results.createVariable('slope_mean', float, ('pressure', 'latitude', 'time'), zlib=True)
    SLOPE_STD = results.createVariable('slope_std', float, ('pressure', 'latitude', 'time'), zlib=True)

    REGRESSOR_COEFFICIENTS = results.createVariable('regressor_coefficients', float, ('nsamples', 'pressure', 'latitude', 'nregs'), zlib=True)

    SIGMA_TREND = results.createVariable('sigma_trend', float, ('nsamples', 'pressure', 'latitude'), zlib=True)
    SIGMA_SEAS = results.createVariable('sigma_seas', float, ('nsamples', 'pressure', 'latitude'), zlib=True)
    SIGMA_AR = results.createVariable('sigma_AR', float, ('nsamples', 'pressure', 'latitude'), zlib=True)
    rhoAR1 = results.createVariable('rhoAR1', float, ('nsamples', 'pressure', 'latitude'), zlib=True)

    results.close()

        
# Put the relevant bits into the netCDF...
def add_results_to_netcdf(results_dir, results_filename, fit, pressure, latitude):

    # Open the netCDF
    results = netCDF4.Dataset('{}/{}.nc'.format(results_dir, results_filename), 'r+')

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
    
    # Regressor coefficients
    results['regressor_coefficients'][:, pressure, latitude, :] = fit.extract()['beta'][:,:]

    # DLM hyper parameters
    results['sigma_trend'][:, pressure, latitude] = fit.extract()['sigma_trend']
    results['sigma_seas'][:, pressure, latitude] = fit.extract()['sigma_seas']
    results['sigma_AR'][:, pressure, latitude] = fit.extract()['sigma_AR']
    results['rhoAR1'][:, pressure, latitude] = fit.extract()['rhoAR1']
    
    # Close the netCDF file
    results.close()

# Save the results
def save_results(results_dir, results_filename, fit, pressure, latitude):

    results = {'trend':fit.extract()['trend'][:,:],
                'trend_mean':np.mean(fit.extract()['trend'][:,:], axis = 0),
                'trend_std':np.std(fit.extract()['trend'][:,:], axis = 0),
                'slope_mean':np.mean(fit.extract()['slope'][:,:], axis = 0),
                'slope_std':np.std(fit.extract()['slope'][:,:], axis = 0),
                'seasonal_mean':np.mean(fit.extract()['seasonal'][:,:], axis = 0),
                'seasonal_std':np.std(fit.extract()['seasonal'][:,:], axis = 0),
                'regressor_coefficients':fit.extract()['beta'][:,:],
                'sigma_trend':fit.extract()['sigma_trend'],
                'sigma_seas':fit.extract()['sigma_seas'],
                'sigma_AR':fit.extract()['sigma_AR'],
                'rhoAR1':fit.extract()['rhoAR1']
                }

    with open('{}/{}_pres{}_lat{}.pkl'.format(results_dir, results_filename, pressure, latitude), 'wb') as f:
        pickle.dump(results, f)

# Convert results to netCDF
def convert_to_netcdf(results_dir, results_filename, P, L):
    
    # Open the netCDF
    results = netCDF4.Dataset('{}/{}.nc'.format(results_dir, results_filename), 'r+')

    # Fill the pressures and latitudes one by one
    for pressure in range(len(P)):
        for latitude in range(len(L)):
            try:
                f = open('{}/{}_pres{}_lat{}.pkl'.format(results_dir, results_filename, pressure, latitude), 'rb')
                res = pickle.load(f)
                f.close()
                
                # Trend
                results['trend'][:, pressure, latitude, :] = res['trend']
                results['trend_mean'][pressure, latitude, :] = res['trend_mean']
                results['trend_std'][pressure, latitude, :] = res['trend_std']
                
                # Slope of the DLM trend
                results['slope_mean'][pressure, latitude, :] = res['slope_mean']
                results['slope_std'][pressure, latitude, :] = res['slope_std']
                
                # Seasonal cycle
                results['seasonal_mean'][pressure, latitude, :] = res['seasonal_mean']
                results['seasonal_std'][pressure, latitude, :] = res['seasonal_std']
                
                # Regressor coefficients
                results['regressor_coefficients'][:, pressure, latitude, :] = res['regressor_coefficients']
                
                # DLM hyper parameters
                results['sigma_trend'][:, pressure, latitude] = res['sigma_trend']
                results['sigma_seas'][:, pressure, latitude] = res['sigma_seas']
                results['sigma_AR'][:, pressure, latitude] = res['sigma_AR']
                results['rhoAR1'][:, pressure, latitude] = res['rhoAR1']
                    
                # Now try and delete the file
                try:
                    os.remove('{}/{}_pres{}_lat{}.pkl'.format(results_dir, results_filename, pressure, latitude))
                except:
                    pass
            except:
                pass

    # Close the netCDF file
    results.close()



def sampling_rate(rate):
    if rate == 'monthly':
        return 1.
    if rate == 'daily':
        return 12./365.25
    if rate == 'annual':
        return 12.
    else:
        print("sampling rate not recognised: must be one of `monthly`, `annual` or `daily`. Please choose an appropriate data sampling rate and try again.")

def prepare_missing_data(d, s):
    d_masked = np.zeros(len(d))
    s_masked = np.zeros(len(s))
    
    mask = np.where(np.isnan(d) == True)[0]
    not_mask = np.where(np.isnan(d) == False)[0]
    d_masked[mask] = np.mean(d[not_mask])
    d_masked[not_mask] = d[not_mask]
    s_masked[mask] = 1e20
    s_masked[not_mask] = s[not_mask]

    return d_masked, s_masked
