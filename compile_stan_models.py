import pystan
import pickle
from models.stan_dlm_models import *

model_vanilla_ar1 = pystan.StanModel(model_code=dlm_vanilla_ar1)
f = open('models/dlm_vanilla_ar1.pkl', 'wb')
pickle.dump(model_vanilla_ar1, f)
f.close()

model_vanilla_ar2 = pystan.StanModel(model_code=dlm_vanilla_ar2)
f = open('models/dlm_vanilla_ar2.pkl', 'wb')
pickle.dump(model_vanilla_ar2, f)
f.close()

model_noregs_ar1 = pystan.StanModel(model_code=dlm_noregs_ar1)
f = open('models/dlm_noregs_ar1.pkl', 'wb')
pickle.dump(model_noregs_ar1, f)
f.close()

model_dynregs_ar1 = pystan.StanModel(model_code=dlm_dynregs_ar1)
f = open('models/dlm_dynregs_ar1.pkl', 'wb')
pickle.dump(model_dynregs_ar1, f)
f.close()

model_vanilla_ar1_noseasonal = pystan.StanModel(model_code=dlm_vanilla_ar1_noseasonal)
f = open('models/dlm_vanilla_ar1_noseasonal.pkl', 'wb')
pickle.dump(dlm_vanilla_ar1_noseasonal, f)
f.close()

model_vanilla_ar2_noseasonal = pystan.StanModel(model_code=dlm_vanilla_ar2_noseasonal)
f = open('models/dlm_vanilla_ar2_noseasonal.pkl', 'wb')
pickle.dump(dlm_vanilla_ar2_noseasonal, f)
f.close()


