import pystan
import pickle
from models.stan_dlm_models import *

model_kalman_ar1 = pystan.StanModel(model_code=code_kalman_ar1)
f = open('models/dlm_kalman_ar1.pkl', 'wb')
pickle.dump(model_kalman_ar1, f)
f.close()

model_kalman_ar2 = pystan.StanModel(model_code=code_kalman_ar2)
f = open('models/dlm_kalman_ar2.pkl', 'wb')
pickle.dump(model_kalman_ar2, f)
f.close()

