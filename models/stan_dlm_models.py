code_kalman_ar1 = """

data {
    int N; // number of time steps
    int nreg; // number of regressors
    vector[N] time_series; // n-composites vectors of N-time-steps -- the data
    vector[N] stddev; // corresponding std deviation for each data point
    vector[nreg] regressors[N]; // nreg regressors of N-time-steps
    real<lower=0> sigma_trend_prior;
    real<lower=0> sigma_seas_prior;
    real<lower=0> sigma_AR_prior;
    real<lower=0> S;
    
}

transformed data {

    // Re-scaled data
    real data_mean = mean(time_series);
    real data_range = max(time_series) - min(time_series);
    vector[N] d = (time_series - data_mean)/data_range;
    vector[N] s = stddev/data_range;

    // Set value of nx
    int nx = nreg+7;
    
    // Declare F-vector (observation projection vector)
    vector[nx] F[N+1];
    
    // Initialize
    for(t in 1:N+1){
        F[t] = rep_vector(0, nx);
    }
    
    // First time point
    F[1][1] = 1;
    F[1][2] = 0;
    F[1][3] = 1;
    F[1][4] = 0;
    F[1][5] = 1;
    F[1][6] = 0;
    for(i in 1:nreg){
        F[1][6+i] = 0;
    }
    F[1][6+nreg+1] = 1;
    
    // Loop over next time points
    for(t in 2:N+1){
        F[t][1] = 1;
        F[t][2] = 0;
        F[t][3] = 1;
        F[t][4] = 0;
        F[t][5] = 1;
        F[t][6] = 0;
        for(i in 1:nreg){
            F[t][6+i] = regressors[t-1][i];
        }
        F[t][6+nreg+1] = 1;
    }
}

parameters {
    real<lower=0> sigma_trend;
    real<lower=0> sigma_seas;
    real<lower=0> sigma_AR;
    real<lower=0, upper=1> rhoAR;
}

transformed parameters{

    // Declare intermediate Kalman objects
    vector[nx] x_hat[N+1];
    vector[nx] x_bar[N+1];
    matrix[nx,nx] C_hat[N+1];
    vector[N+1] C_y;
    vector[nx] K[N+1];
    matrix[nx,nx] C_bar[N+1];
    vector[N+1] v;

    // Declare state-space model matrices
    matrix[nx, nx] G; // G-matrix (state transition matrix)
    matrix[nx, nx] W; // W-matrix (state covariance)
    matrix[nx, nx] C; // C-matrix (observation covariance)
    
    // Initialize state-space matrices
    G = rep_matrix(0, nx, nx);
    W = rep_matrix(0, nx, nx);
    C = rep_matrix(0, nx, nx);
    
    // Initialize Kalman filter stuff
    C_y = rep_vector(0, N+1);
    v = rep_vector(0, N+1);
    for(t in 1:N+1){
        x_hat[t] = rep_vector(0, nx);
        x_bar[t] = rep_vector(0, nx);
        K[t] = rep_vector(0, nx);
        C_bar[t] = rep_matrix(0, nx, nx);
        C_hat[t] = rep_matrix(0, nx, nx);
    }
    
    // Assign the G matrix
    
    // Trend part
    G[1,1] = 1;
    G[1,2] = 1;
    G[2,1] = 0;
    G[2,2] = 1;
    
    // 12-month seasonal part
    G[3,3] = cos(2*pi()/12);
    G[3,4] = sin(2*pi()/12);
    G[4,3] = -sin(2*pi()/12);
    G[4,4] = cos(2*pi()/12);
    
    // 6-month seasonal part
    G[5,5] = cos(2*pi()/6);
    G[5,6] = sin(2*pi()/6);
    G[6,5] = -sin(2*pi()/6);
    G[6,6] = cos(2*pi()/6);
    
    // Regressor part
    for(i in 1:nreg){
        G[6+i,6+i] = 1;
    }
    
    // AR process part
    G[6+nreg+1,6+nreg+1] = rhoAR;
    
    // Assign W matrix
    W[1,1] = 1e-10;
    W[2,2] = sigma_trend^2;
    W[3,3] = sigma_seas^2;
    W[4,4] = sigma_seas^2;
    W[5,5] = sigma_seas^2;
    W[6,6] = sigma_seas^2;
    for(i in 1:nreg){
        W[6+i,6+i] = 1e-10;
    }
    W[nx,nx] = sigma_AR^2;
    
    // Assign C matrix
    
    // Prior covariance of regression coefficients (and seasonal cycle and trend amplitudes)
    for(i in 1:(nx-1)){
        C[i,i] = S;
    }
    
    // Prior (stationary) covariance of the AR1 process
    C[nx,nx] = sigma_AR^2/(1-rhoAR^2);
    
    // Kalman filtering
    
    // Initialize C_bar
    C_bar[1] = C;
    
    // Do the forward filtering
    for(t in 2:N+1){
        x_hat[t] = G*x_bar[t-1];
        C_hat[t] = G*(C_bar[t-1]*G') + W;
        C_y[t] = F[t]'*C_hat[t]*F[t] + s[t-1]^2;
        K[t] = (C_hat[t]*F[t])/C_y[t];
        v[t] = d[t-1] - F[t]'*x_hat[t];
        x_bar[t] = x_hat[t] + K[t]*v[t];
        C_bar[t] = C_hat[t] - (K[t]*F[t]')*C_hat[t];
    }
}

model {

    // Priors on hyper parameters
    sigma_trend ~ normal(0, sigma_trend_prior);
    sigma_seas ~ normal(0, sigma_seas_prior);
    sigma_AR ~ normal(0, sigma_AR_prior);
    rhoAR ~ uniform(0, 1);
    
    // Construct target density
    for(t in 2:N+1){
        target += normal_lpdf(v[t] | 0, sqrt(C_y[t]));
    }
}

generated quantities {

    // Declare Kalman smoother quantities
    vector[N] trend;
    vector[N] slope;
    vector[N] seasonal;
    vector[N] ar;
    vector[nreg] beta;
    vector[N] residuals;
    vector[nx] x[N];
    vector[nx] x_realization[N+1];
    vector[nx] x_ubar[N+1];
    vector[nx] x_realization_hat[N+1];
    vector[N+1] y_realization;
    vector[N+1] v_realization;
    vector[nx] r[N+2];
    vector[nx] r_realization[N+2];
    matrix[nx,nx] L[N+2];
    matrix[nx,nx] M[N+2];
    matrix[nx,nx] C_twidle[N+1];
    vector[nx] x_twidle[N+1];
    vector[nx] x_realization_twidle[N+1];
    
    // Initialize matrix quantities and end points of Kalman smoothing vectors
    for(t in 1:N+1){
        L[t] = rep_matrix(0, nx, nx);
        M[t] = rep_matrix(0, nx, nx);
        C_twidle[t] = rep_matrix(0, nx, nx);
    }
    L[N+2] = rep_matrix(0, nx, nx);
    M[N+2] = rep_matrix(0, nx, nx);
    r[N+2] = rep_vector(0, nx);
    r_realization[N+2] = rep_vector(0, nx);
 
    // Generate a draw from the state-space equations
    x_realization[1] = multi_normal_rng(rep_vector(0, nx), C);
    y_realization[1] = F[1]'*x_realization[1] + normal_rng(0, 10*s[1]);
    for(t in 2:N+1){
        x_realization[t] = G*x_realization[t-1] + multi_normal_rng(rep_vector(0, nx), W);
        y_realization[t] = F[t]'*x_realization[t] + normal_rng(0, s[t-1]);
    }
    
    // Generate Kalman filter of realization; x_realization_hat
    x_ubar[1] = rep_vector(0, nx);
    for(t in 2:N+1){
        x_realization_hat[t] = G*x_ubar[t-1];
        v_realization[t] = y_realization[t] - F[t]'*x_realization_hat[t];
        x_ubar[t] = x_realization_hat[t] + K[t]*v_realization[t];
    }
    
    // Generate Kalman smoother of realization and of data
    for(t in 1:N+1){
    
        L[N+2-t] = G - G*(K[N+2-t]*F[N+2-t]');
        r[N+2-t] = F[N+2-t]*v[N+2-t]/C_y[N+2-t] + L[N+2-t]'*r[N+2-t+1];
        r_realization[N+2-t] = F[N+2-t]*v_realization[N+2-t]/C_y[N+2-t] + L[N+2-t]'*r_realization[N+2-t+1];
        M[N+2-t] = (F[N+2-t]*F[N+2-t]')/C_y[N+2-t] + L[N+2-t]'*(M[N+2-t+1]*L[N+2-t]);
        x_twidle[N+2-t] = x_hat[N+2-t] + C_hat[N+2-t]*r[N+2-t];
        x_realization_twidle[N+2-t] = x_realization_hat[N+2-t] + C_hat[N+2-t]*r_realization[N+2-t];
        C_twidle[N+2-t] = C_hat[N+2-t] - C_hat[N+2-t]*(M[N+2-t]*C_hat[N+2-t]);
    }
    
    for(t in 1:N){
        x[t] = x_realization[t+1] - x_realization_twidle[t+1] + x_twidle[t+1];
        trend[t] = x[t][1]*data_range + data_mean;
        slope[t] = x[t][2]*data_range;
        seasonal[t] = (x[t][3] + x[t][5])*data_range;
        ar[t] = x[t][nx]*data_range;
        residuals[t] = (d[t] - (F[t]'*x[t] - ar[t]/data_range))*data_range;
    }
    for(i in 1:nreg){
        beta[i] = x[1][5+i]*data_range;
    }
}

"""


code_kalman_ar2 = """

data {
    int N; // number of time steps
    int nreg; // number of regressors
    int nx; // length of state vector
    vector[N] time_series; // n-composites vectors of N-time-steps -- the data
    vector[N] stddev; // corresponding std deviation for each data point
    vector[nreg] regressors[N]; // nreg regressors of N-time-steps
    real<lower=0> sigma_trend_prior;
    real<lower=0> sigma_seas_prior;
    real<lower=0> sigma_AR_prior;
    real<lower=0> S;
}

transformed data {

    // Re-scaled data
    real data_mean = mean(time_series);
    real data_range = max(time_series) - min(time_series);
    vector[N] d = (time_series - data_mean)/data_range;
    vector[N] s = stddev/data_range;
    
    // Declare F-vector (observation projection vector)
    vector[nx] F[N+1];
    
    // Initialize
    for(t in 1:N+1){
        F[t] = rep_vector(0, nx);
    }
    
    // First time point
    F[1][1] = 1;
    F[1][2] = 0;
    F[1][3] = 1;
    F[1][4] = 0;
    F[1][5] = 1;
    F[1][6] = 0;
    for(i in 1:nreg){
        F[1][6+i] = 0; 
    }  
    F[1][6+nreg+1] = 1;
    
    // Loop over next time points
    for(t in 2:N+1){
        F[t][1] = 1;
        F[t][2] = 0;
        F[t][3] = 1;
        F[t][4] = 0;
        F[t][5] = 1;
        F[t][6] = 0;
        for(i in 1:nreg){
            F[t][6+i] = regressors[t-1][i]; 
        }  
        F[t][6+nreg+1] = 1;
    }
}

parameters {
    real<lower=0> sigma_trend;
    real<lower=0> sigma_seas;
    real<lower=0> sigma_AR;
    real<lower=0, upper=1> rhoAR1;
    real<lower=0, upper=1> rhoAR2;
}

transformed parameters{

    // Declare intermediate Kalman objects
    vector[nx] x_hat[N+1];
    vector[nx] x_bar[N+1];
    matrix[nx,nx] C_hat[N+1];
    vector[N+1] C_y;
    vector[nx] K[N+1];
    matrix[nx,nx] C_bar[N+1];
    vector[N+1] v;

    // Declare state-space model matrices
    matrix[nx, nx] G; // G-matrix (state transition matrix)
    matrix[nx, nx] W; // W-matrix (state covariance)
    matrix[nx, nx] C; // C-matrix (observation covariance)
    
    // Initialize state-space matrices
    G = rep_matrix(0, nx, nx);
    W = rep_matrix(0, nx, nx);
    C = rep_matrix(0, nx, nx);
    
    // Initialize Kalman filter stuff
    C_y = rep_vector(0, N+1);
    v = rep_vector(0, N+1);
    for(t in 1:N+1){
        x_hat[t] = rep_vector(0, nx);
        x_bar[t] = rep_vector(0, nx);
        K[t] = rep_vector(0, nx);
        C_bar[t] = rep_matrix(0, nx, nx);
        C_hat[t] = rep_matrix(0, nx, nx);
    }
    
    // Assign the G matrix
    
    // Trend part
    G[1,1] = 1;
    G[1,2] = 1;
    G[2,1] = 0;
    G[2,2] = 1;
    
    // 12-month seasonal part
    G[3,3] = cos(2*pi()/12);
    G[3,4] = sin(2*pi()/12);
    G[4,3] = -sin(2*pi()/12);
    G[4,4] = cos(2*pi()/12);
    
    // 6-month seasonal part
    G[5,5] = cos(2*pi()/6);
    G[5,6] = sin(2*pi()/6);
    G[6,5] = -sin(2*pi()/6);
    G[6,6] = cos(2*pi()/6);
    
    // Regressor part
    for(i in 1:nreg){
        G[6+i,6+i] = 1;
    }
    
    // AR process part
    G[6+nreg+1,6+nreg+1] = rhoAR1;
    G[6+nreg+1,6+nreg+2] = rhoAR2;
    G[6+nreg+2, 6+nreg+1] = 1;

    // Assign W matrix
    W[1,1] = 1e-10;
    W[2,2] = sigma_trend^2;
    W[3,3] = sigma_seas^2;
    W[4,4] = sigma_seas^2;
    W[5,5] = sigma_seas^2;
    W[6,6] = sigma_seas^2;
    for(i in 1:nreg){
        W[6+i,6+i] = 1e-10;
    }
    W[nx,nx] = sigma_AR^2;
    
    // Assign C matrix
    
    // Prior covariance of regression coefficients (and seasonal cycle and trend amplitudes)
    for(i in 1:(nx-1)){
        C[i,i] = S;
    }
    
    // Prior (stationary) covariance of the AR1 process
    C[nx, nx] = sigma_AR^2/(1-(rhoAR1^2 + rhoAR2^2) - 2*rhoAR2*rhoAR1^2);
    C[nx-1, nx-1] = sigma_AR^2/(1-(rhoAR1^2 + rhoAR2^2) - 2*rhoAR2*rhoAR1^2);
    C[nx, nx-1] = C[nx, nx]*rhoAR1/(1-rhoAR2);
    C[nx-1, nx] = C[nx, nx]*rhoAR1/(1-rhoAR2);

    // Kalman filtering
    
    // Initialize C_bar
    C_bar[1] = C;
        
    // Do the forward filtering
    for(t in 2:N+1){
        x_hat[t] = G*x_bar[t-1];
        C_hat[t] = G*(C_bar[t-1]*G') + W;
        C_y[t] = F[t]'*C_hat[t]*F[t] + s[t-1]^2;
        K[t] = (C_hat[t]*F[t])/C_y[t];
        v[t] = d[t-1] - F[t]'*x_hat[t];
        x_bar[t] = x_hat[t] + K[t]*v[t];
        C_bar[t] = C_hat[t] - (K[t]*F[t]')*C_hat[t];
    }
}

model {
 
    // Priors on hyper parameters
    sigma_trend ~ normal(0, sigma_trend_prior);
    sigma_seas ~ normal(0, sigma_seas_prior);
    sigma_AR ~ normal(0, sigma_AR_prior);
    rhoAR1 ~ uniform(0, 1);
    rhoAR2 ~ uniform(0, sqrt(rhoAR1^4+1-rhoAR1^2)-rhoAR1^2);

    // Construct target density
    for(t in 2:N+1){
        target += normal_lpdf(v[t] | 0, sqrt(C_y[t]));
    } 
}

generated quantities {

    // Declare Kalman smoother quantities
    vector[N] trend;
    vector[N] slope;
    vector[N] seasonal;
    vector[N] ar;
    vector[nreg] beta;
    vector[nx] x[N];
    vector[N] residuals;
    vector[nx] x_realization[N+1];
    vector[nx] x_ubar[N+1];
    vector[nx] x_realization_hat[N+1];
    vector[N+1] y_realization;
    vector[N+1] v_realization;
    vector[nx] r[N+2];
    vector[nx] r_realization[N+2];
    matrix[nx,nx] L[N+2];
    matrix[nx,nx] M[N+2];
    matrix[nx,nx] C_twidle[N+1];
    vector[nx] x_twidle[N+1];
    vector[nx] x_realization_twidle[N+1];
    
    // Initialize matrix quantities and end points of Kalman smoothing vectors
    for(t in 1:N+1){
        L[t] = rep_matrix(0, nx, nx);
        M[t] = rep_matrix(0, nx, nx);
        C_twidle[t] = rep_matrix(0, nx, nx);
    }
    L[N+2] = rep_matrix(0, nx, nx);
    M[N+2] = rep_matrix(0, nx, nx);
    r[N+2] = rep_vector(0, nx);
    r_realization[N+2] = rep_vector(0, nx);
 
    // Generate a draw from the state-space equations
    x_realization[1] = multi_normal_rng(rep_vector(0, nx), C);
    y_realization[1] = F[1]'*x_realization[1] + normal_rng(0, 10*s[1]);
    for(t in 2:N+1){
        x_realization[t] = G*x_realization[t-1] + multi_normal_rng(rep_vector(0, nx), W);
        y_realization[t] = F[t]'*x_realization[t] + normal_rng(0, s[t-1]);
    }
    
    // Generate Kalman filter of realization; x_realization_hat
    x_ubar[1] = rep_vector(0, nx);
    for(t in 2:N+1){
        x_realization_hat[t] = G*x_ubar[t-1];
        v_realization[t] = y_realization[t] - F[t]'*x_realization_hat[t];
        x_ubar[t] = x_realization_hat[t] + K[t]*v_realization[t];
    }
    
    // Generate Kalman smoother of realization and of data
    for(t in 1:N+1){
    
        L[N+2-t] = G - G*(K[N+2-t]*F[N+2-t]');
        r[N+2-t] = F[N+2-t]*v[N+2-t]/C_y[N+2-t] + L[N+2-t]'*r[N+2-t+1];
        r_realization[N+2-t] = F[N+2-t]*v_realization[N+2-t]/C_y[N+2-t] + L[N+2-t]'*r_realization[N+2-t+1];
        M[N+2-t] = (F[N+2-t]*F[N+2-t]')/C_y[N+2-t] + L[N+2-t]'*(M[N+2-t+1]*L[N+2-t]);
        x_twidle[N+2-t] = x_hat[N+2-t] + C_hat[N+2-t]*r[N+2-t];
        x_realization_twidle[N+2-t] = x_realization_hat[N+2-t] + C_hat[N+2-t]*r_realization[N+2-t];
        C_twidle[N+2-t] = C_hat[N+2-t] - C_hat[N+2-t]*(M[N+2-t]*C_hat[N+2-t]);        
    }
    
    for(t in 1:N){
        x[t] = x_realization[t+1] - x_realization_twidle[t+1] + x_twidle[t+1];
        trend[t] = x[t][1]*data_range + data_mean;
        slope[t] = x[t][2]*data_range;
        seasonal[t] = (x[t][3] + x[t][5])*data_range;
        ar[t] = x[t][nx]*data_range;
        residuals[t] = (d[t] - (F[t]'*x[t] - ar[t]/data_range))*data_range;
    }
    for(i in 1:nreg){
        beta[i] = x[1][5+i]*data_range;
    }
}

"""






code_seasonal = """

data {
    int N; // number of time steps
    int nreg; // number of regressors
    vector[N] d; // n-composites vectors of N-time-steps -- the data
    vector[N] s; // corresponding standard deviation for each data point
    vector[nreg] regressors[N]; // nreg regressors of N-time-steps
    matrix[nreg, nreg] C_GP_prior; // prior covariance on regressor coefficients
    matrix[4,4] G; // Seasonal transition matrix
    real<lower=0> sigma_trend_prior;
    real<lower=0> sigma_seas_prior;
    real<lower=0> sigma_AR_prior;
}

parameters {
    vector[nreg] b; // regression coefficients
    vector[N] ar; // ar-process
    vector[4] seasonal_components[N]; // seasonal cycle
    vector[N] slope; // dynamical linear trend [slope]
    real x0;
    real<lower=0> sigma_trend;
    real<lower=0> sigma_seas;
    real<lower=0> sigma_AR;
    real<lower=0, upper=1> rho;
}

transformed parameters{
    vector[N] trend; // dynamical linear trend
    vector[N] seasonal; // seasonal cycle
    vector[N] regression_model; // complete regression model

    // Model for the latent dynamical model
    trend[1] = x0;
    seasonal[1] = seasonal_components[1][1] + seasonal_components[1][3];
    regression_model[1] = dot_product(b, regressors[1]);
    for(t in 2:N){

        trend[t] = trend[t-1] + slope[t];
        seasonal[t] = seasonal_components[t][1] + seasonal_components[t][3];
        regression_model[t] = dot_product(b, regressors[t]);
    }
}

model {

    // Priors on hyper parameters
    sigma_trend ~ normal(0, sigma_trend_prior);
    sigma_seas ~ normal(0, sigma_seas_prior);
    sigma_AR ~ normal(0, sigma_AR_prior);
    seasonal_components[1] ~ normal(0, 10);
    rho ~ uniform(0, 1);
    b ~ multi_normal(rep_vector(0, nreg), C_GP_prior);
    
    // Prior on AR process initial value
    ar[1] ~ normal(0, sqrt(sigma_AR^2/(1-rho^2)));

    // Model for the latent dynamical model
    for(t in 2:N){
    
        // AR-process
        ar[t] ~ normal(ar[t-1]*rho, sigma_AR);
        
        // Seasonal cycle
        seasonal_components[t] ~ normal(G*seasonal_components[t-1], sigma_seas);

        // Stochastic DLM trend
        slope[t] ~ normal(slope[t-1], sigma_trend);
    }

    # Construct target density
    for(t in 1:N){
        target += normal_lpdf(ar[t] + trend[t] + seasonal[t] + regression_model[t] | d[t], s[t]);
    }
}

generated quantities{

    vector[N] residuals;
    residuals = d - regression_model - seasonal - trend;
}

"""
