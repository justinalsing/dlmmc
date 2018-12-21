import numpy as np
import scipy.optimize as opt
import scipy.interpolate as interpolate

# MLR for initialization

def regression_model(coefficients, regressors, T):
    
    # initialize
    y = np.zeros(len(T))
    
    # proxies (including seasonal cycles and trend)
    for i in range(len(regressors)):
        y += coefficients[i]*regressors[i]
    
    return y

def chi2(coefficients, regressors, T, y):
    
    y_model = regression_model(coefficients, regressors, T)
    
    return sum((y - y_model)**2)

# Do the MLR initialization
def mlr_initialization(regressors, nreg, T, d):
    
    nt = len(T)

    # Package the regressors
    mlr_regs = [0]*(nreg+4+2)
    for i in range(0, nreg):
        mlr_regs[i] = regressors[i]

    mlr_regs[nreg] = np.sin(2*np.pi*np.arange(len(T))/12.)
    mlr_regs[nreg+1] = np.cos(2*np.pi*np.arange(len(T))/12.)
    mlr_regs[nreg+2] = np.sin(4*np.pi*np.arange(len(T))/12.)
    mlr_regs[nreg+3] = np.cos(4*np.pi*np.arange(len(T))/12.)

    mlr_regs[nreg+4] = np.arange(len(T))
    mlr_regs[nreg+5] = np.ones(len(T))

    # Minimize chi2
    res = opt.minimize(chi2, x0 = np.ones(len(mlr_regs))*1e-1, args = (mlr_regs, T, d), tol = 1e-10)

    y_model = regression_model(res.x, mlr_regs, T)
    residuals = d - y_model

    # Extract the MLR trend and seasonal cycle
    b_mlr = res.x[0:nreg+4]
    slope_mlr = np.concatenate([np.array([0]), np.ones(len(T)-1)*res.x[-2]])
    sigma_AR_mlr = np.std(residuals)
    rho_AR_mlr = np.corrcoef(residuals[1:nt-1], residuals[2:nt])[0,1]
    seasonal_components_mlr = np.column_stack([mlr_regs[nreg+1]*res.x[nreg+1] + mlr_regs[nreg]*res.x[nreg],
                                           mlr_regs[nreg+1]*res.x[nreg+1] - mlr_regs[nreg]*res.x[nreg],
                                           mlr_regs[nreg+3]*res.x[nreg+3] + mlr_regs[nreg+2]*res.x[nreg+2],
                                           mlr_regs[nreg+3]*res.x[nreg+3] - mlr_regs[nreg+2]*res.x[nreg+2]])

    return residuals, b_mlr, slope_mlr, seasonal_components_mlr, sigma_AR_mlr, rho_AR_mlr
