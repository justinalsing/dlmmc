# Contributing

Here I provide a short summary of how to contribute extended DLM models to this package. 

Note that in order to contribute new models, you will need to be able to write code in [stan](https://mc-stan.org) (click on the link for documentation and access to the stan manual). This might sound like pain for the uninitiated, but stan is a very intuitive language and you will only need to make relatively minor changes to the existing models code in order to implement new models (see below).

## Conding extended models

The suite of DLM models implemented here are completely characterized by the matrices **F**, **G**, **W** and **C** (described in the file `model_descriptions.pdf`). In order to develop a new DLM model within dlmmc, more-or-less all you need to do is work out what the matrices of your model are, copy one of the existing stan models in `models/stan_dlm_models.py` and replace the matrices with your own needs. The likelihood and Kalman filtering steps that makes up the bulk of the code for each model does not need changing for new models, so the modifications required to implement a new model are relatively light.

A more detailed work-flow for implementing a new model is as follows:

1. **Copy-paste one of the existing stan models** inside `stan_dlm_models.py` and re-name it.
2. **Update the `data{}` block** to add any additional data / quantities that are required to build your model matrices.
3. **Update the `transformed_data{}` block**, replacing the model matrix **F** construction with your own requirements and making sure the length of the state vector `nx` is correct for your new model.
4. **Update the `transformed_parameters{}` block** replacing the model matrices **G**, **W** and **C** constructions with your derived model matrices.
5. **Update the `model{}` block** to include any new hyper-parameters that your model has. New hyper-parameters need to be declared along with the other existing hyper-parameters at the top of this code block, and then a prior can be imposed. This is painless if you follow what was done for the other hyper-parameters (it requires just two new lines of code per additional hyper-parameter).
6. **Update the `generated_quantities{}` block** to extract the new bits of your model from the state vector. At the end of this block of code (the last loop in each stan model code), there is a loop which takes the state vector **x** and extracts from it the various components of your model (trend, seasonal cycle, etc). This is useful to make post-processing the output easier. In your new model, you should add lines of code here (following what was done for the other components of the state vector) to extract the various components of your model from **x**.
7. **Compile and run the model.** Once you think you're done coding your new model, you'll need to compile it. This can be done with the following python code:

```import pystan
from models.stan_dlm_models import *

my_new_model = pystan.StanModel(model_code=my_new_model)
```

If it compiles OK, go ahead and run it on your data own (following `dlm_tutorial.ipynb`).

## Once you are happy with your new model

When you are happy that your new model is working and is stable, you should add the following (example) code to `compile_stan_models.py` so that users can automatically compile your model when they download/pull the code:

```my_new_model = pystan.StanModel(model_code=my_new_model)
f = open('models/my_new_model.pkl', 'wb')
pickle.dump(my_new_model, f)
f.close()
```

Finally, you should incude a short description of your extended model into the document `model_descriptions.pdf`. This can be done by modifying `model_descriptions.tex` (and `model_descriptions.bib`) and re-compiling the tex file to update the pdf. You should also include a validation test in the test suite `dlm_validation_tests.ipynb` (following what was done for the other models there; generate a draw from your model as mock data, run your model on the mock data, plot the recovery of the state-vector and hyper-parameters).

## Want a new model but need assistance?

If you desperately want your new model implemented, but are uncomfortable coding in stan, please contact me either directly at justin.alsing@fysik.su.se or via the issues channel associated with this Git repo. I'll be happy to help.
