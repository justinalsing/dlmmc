import numpy as np
import scipy

### DLM_VANILLA_AR1 ###

# Generate realisation of x and y from state space equations given a fixed x0
def simulate_x_y_dlm_vanilla_ar1(x0, sigma_trend, sigma_seas, sigma_AR, r1, regressors, S, V, nt, nx):
    
    # Generate model matrices
    G, W, F, C = model_matrices_dlm_vanilla_ar1(sigma_trend, sigma_seas, sigma_AR, r1, nx, nt, regressors, S)
    
    # Arrays for x and y
    x = np.zeros((nx, nt+1))
    y = np.zeros((nt+1))
    
    # Fix initial x_0
    x[:,0] = x0
    
    # Generate x and y iteratively
    for t in range(1, nt+1):
        x[:,t] = np.dot(G, x[:, t-1]) + np.sqrt(np.diag(W))*np.random.normal(0,1,nx)
        y[t] = np.dot(F[:, t], x[:, t]) + np.sqrt(V[t])*np.random.normal(0, 1)
    
    return x, y


def model_matrices_dlm_vanilla_ar1(sigma_trend, sigma_seas, sigma_AR, r1, nx, nt, regressors, S):
    
    # Fill in the theta values
    nreg = len(regressors.T)
    
    # G matrix
    G = np.zeros((nx, nx))
    
    # Trend
    G[0:2, 0:2] = np.array(([[1, 1],[0, 1]]))
    
    # Seasonal
    G[2:4, 2:4] = np.array([[np.cos(2*np.pi/12), np.sin(2*np.pi/12)],[-np.sin(2*np.pi/12), np.cos(2*np.pi/12)]])
    G[4:6, 4:6] = np.array([[np.cos(2*2*np.pi/12), np.sin(2*2*np.pi/12)],[-np.sin(2*2*np.pi/12), np.cos(2*2*np.pi/12)]])
    
    # Proxies
    G[6:6+nreg, 6:6+nreg] = np.eye(nreg)
    
    # AR3 process
    G[6+nreg, 6+nreg] = r1
    
    # W matrix
    W = np.zeros((nx, nx))
    W[-1, -1] = sigma_AR**2
    W[1, 1] = sigma_trend**2
    W[2,2] = sigma_seas**2
    W[3,3] = sigma_seas**2
    W[4,4] = sigma_seas**2
    W[5,5] = sigma_seas**2
    
    # F matrix
    F = np.zeros((nx, nt+1))
    F[0, :] = 1 # Trend
    F[1, :] = 0 # Trend
    F[2, :] = 1 # Seasonal
    F[3, :] = 0 # Seasonal
    F[4, :] = 1 # Seasonal
    F[5, :] = 0 # Seasonal
    for k in range(0, nreg):
        F[6+k,1:] = regressors.T[k,:] # Regressors
    F[6+nreg, :] = 1 # AR process
    
    # Covariance of prior on x0
    C = np.eye(nx)*S
    C[-1, -1] = sigma_AR**2/(1-r1**2)
    
    return G, W, F, C

### DLM_VANILLA_AR2 ###

# Generate realisation of x and y from state space equations given a fixed x0
def simulate_x_y_dlm_vanilla_ar2(x0, sigma_trend, sigma_seas, sigma_AR, r1, r2, regressors, S, V, nt, nx):
    
    # Generate model matrices
    G, W, F, C = model_matrices_dlm_vanilla_ar2(sigma_trend, sigma_seas, sigma_AR, r1, r2, nx, nt, regressors, S)
    
    # Arrays for x and y
    x = np.zeros((nx, nt+1))
    y = np.zeros((nt+1))
    
    # Fix initial x_0
    x[:,0] = x0
    
    # Generate x and y iteratively
    for t in range(1, nt+1):
        x[:,t] = np.dot(G, x[:, t-1]) + np.sqrt(np.diag(W))*np.random.normal(0,1,nx)
        y[t] = np.dot(F[:, t], x[:, t]) + np.sqrt(V[t])*np.random.normal(0, 1)
    
    return x, y


def model_matrices_dlm_vanilla_ar2(sigma_trend, sigma_seas, sigma_AR, r1, r2, nx, nt, regressors, S):
    
    # Fill in the theta values
    nreg = len(regressors.T)
    
    # G matrix
    G = np.zeros((nx, nx))
    
    # Trend
    G[0:2, 0:2] = np.array(([[1, 1],[0, 1]]))
    
    # Seasonal
    G[2:4, 2:4] = np.array([[np.cos(2*np.pi/12), np.sin(2*np.pi/12)],[-np.sin(2*np.pi/12), np.cos(2*np.pi/12)]])
    G[4:6, 4:6] = np.array([[np.cos(2*2*np.pi/12), np.sin(2*2*np.pi/12)],[-np.sin(2*2*np.pi/12), np.cos(2*2*np.pi/12)]])
    
    # Proxies
    G[6:6+nreg, 6:6+nreg] = np.eye(nreg)
    
    # AR3 process
    G[6+nreg, 6+nreg] = r1
    G[6+nreg, 6+nreg+1] = r2
    G[6+nreg+1, 6+nreg] = 1
    
    # W matrix
    W = np.zeros((nx, nx))
    W[-2, -2] = sigma_AR**2
    W[1, 1] = sigma_trend**2
    W[2,2] = sigma_seas**2
    W[3,3] = sigma_seas**2
    W[4,4] = sigma_seas**2
    W[5,5] = sigma_seas**2
    
    # F matrix
    F = np.zeros((nx, nt+1))
    F[0, :] = 1 # Trend
    F[1, :] = 0 # Trend
    F[2, :] = 1 # Seasonal
    F[3, :] = 0 # Seasonal
    F[4, :] = 1 # Seasonal
    F[5, :] = 0 # Seasonal
    for k in range(0, nreg):
        F[6+k,1:] = regressors.T[k,:] # Regressors
    F[6+nreg, :] = 1 # AR process
    F[6+nreg+1, :] = 0 # AR process

    # Covariance of prior on x0
    C = np.eye(nx)*S
    C[-1, -1] = C[-2, -2] = sigma_AR**2/(1-r1**2 - r2**2 - (2*r2*r1**2)/(1-r2) )
    C[-1, -2] = C[-2, -1] = C[-1, -1]*r1/(1-r2)
    
    return G, W, F, C


### DLM_NOREGS_AR1 ###

# Generate realisation of x and y from state space equations given a fixed x0
def simulate_x_y_dlm_noregs_ar1(x0, sigma_trend, sigma_seas, sigma_AR, r1, S, V, nt, nx):
    
    # Generate model matrices
    G, W, F, C = model_matrices_dlm_noregs_ar1(sigma_trend, sigma_seas, sigma_AR, r1, nx, nt, S)
    
    # Arrays for x and y
    x = np.zeros((nx, nt+1))
    y = np.zeros((nt+1))
    
    # Fix initial x_0
    x[:,0] = x0
    
    # Generate x and y iteratively
    for t in range(1, nt+1):
        x[:,t] = np.dot(G, x[:, t-1]) + np.sqrt(np.diag(W))*np.random.normal(0,1,nx)
        y[t] = np.dot(F[:, t], x[:, t]) + np.sqrt(V[t])*np.random.normal(0, 1)
    
    return x, y


def model_matrices_dlm_noregs_ar1(sigma_trend, sigma_seas, sigma_AR, r1, nx, nt, S):
    
    # G matrix
    G = np.zeros((nx, nx))
    
    # Trend
    G[0:2, 0:2] = np.array(([[1, 1],[0, 1]]))
    
    # Seasonal
    G[2:4, 2:4] = np.array([[np.cos(2*np.pi/12), np.sin(2*np.pi/12)],[-np.sin(2*np.pi/12), np.cos(2*np.pi/12)]])
    G[4:6, 4:6] = np.array([[np.cos(2*2*np.pi/12), np.sin(2*2*np.pi/12)],[-np.sin(2*2*np.pi/12), np.cos(2*2*np.pi/12)]])
    
    # AR1 process
    G[6, 6] = r1
    
    # W matrix
    W = np.zeros((nx, nx))
    W[-1, -1] = sigma_AR**2
    W[1, 1] = sigma_trend**2
    W[2,2] = sigma_seas**2
    W[3,3] = sigma_seas**2
    W[4,4] = sigma_seas**2
    W[5,5] = sigma_seas**2
    
    # F matrix
    F = np.zeros((nx, nt+1))
    F[0, :] = 1 # Trend
    F[1, :] = 0 # Trend
    F[2, :] = 1 # Seasonal
    F[3, :] = 0 # Seasonal
    F[4, :] = 1 # Seasonal
    F[5, :] = 0 # Seasonal
    F[6, :] = 1 # AR process
    
    # Covariance of prior on x0
    C = np.eye(nx)*S
    C[-1, -1] = sigma_AR**2/(1-r1**2)
    
    return G, W, F, C


### DLM_DYNREGS_AR1 ###

# Generate realisation of x and y from state space equations given a fixed x0
def simulate_x_y_dlm_dynregs_ar1(x0, sigma_trend, sigma_seas, sigma_AR, r1, sigma_reg, regressors, S, V, nt, nx):
    
    # Generate model matrices
    G, W, F, C = model_matrices_dlm_dynregs_ar1(sigma_trend, sigma_seas, sigma_AR, r1, sigma_reg, nx, nt, regressors, S)
    
    # Arrays for x and y
    x = np.zeros((nx, nt+1))
    y = np.zeros((nt+1))
    
    # Fix initial x_0
    x[:,0] = x0
    
    # Generate x and y iteratively
    for t in range(1, nt+1):
        x[:,t] = np.dot(G, x[:, t-1]) + np.sqrt(np.diag(W))*np.random.normal(0,1,nx)
        y[t] = np.dot(F[:, t], x[:, t]) + np.sqrt(V[t])*np.random.normal(0, 1)
    
    return x, y


def model_matrices_dlm_dynregs_ar1(sigma_trend, sigma_seas, sigma_AR, r1, sigma_reg, nx, nt, regressors, S):
    
    # Fill in the theta values
    nreg = len(regressors.T)
    
    # G matrix
    G = np.zeros((nx, nx))
    
    # Trend
    G[0:2, 0:2] = np.array(([[1, 1],[0, 1]]))
    
    # Seasonal
    G[2:4, 2:4] = np.array([[np.cos(2*np.pi/12), np.sin(2*np.pi/12)],[-np.sin(2*np.pi/12), np.cos(2*np.pi/12)]])
    G[4:6, 4:6] = np.array([[np.cos(2*2*np.pi/12), np.sin(2*2*np.pi/12)],[-np.sin(2*2*np.pi/12), np.cos(2*2*np.pi/12)]])
    
    # Proxies
    G[6:6+nreg, 6:6+nreg] = np.eye(nreg)
    
    # AR3 process
    G[6+nreg, 6+nreg] = r1
    
    # W matrix
    W = np.zeros((nx, nx))
    W[-1, -1] = sigma_AR**2
    W[1, 1] = sigma_trend**2
    W[2,2] = sigma_seas**2
    W[3,3] = sigma_seas**2
    W[4,4] = sigma_seas**2
    W[5,5] = sigma_seas**2
    for i in range(nreg):
        W[6+i, 6+i] = sigma_reg[i]**2
    
    # F matrix
    F = np.zeros((nx, nt+1))
    F[0, :] = 1 # Trend
    F[1, :] = 0 # Trend
    F[2, :] = 1 # Seasonal
    F[3, :] = 0 # Seasonal
    F[4, :] = 1 # Seasonal
    F[5, :] = 0 # Seasonal
    for k in range(0, nreg):
        F[6+k,1:] = regressors.T[k,:] # Regressors
    F[6+nreg, :] = 1 # AR process
    
    # Covariance of prior on x0
    C = np.eye(nx)*S
    C[-1, -1] = sigma_AR**2/(1-r1**2)
    
    return G, W, F, C
