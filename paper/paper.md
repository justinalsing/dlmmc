---
title: 'dlmmc: Dynamical linear model regression for atmospheric time-series analysis'
tags:
  - python
  - stratospheric ozone
  - atmospheric time-series
  - time-series analysis
authors:
  - name: Justin A. Alsing
    orcid: 0000-0003-4618-3546
    affiliation: "1, 2, 3"
affiliations:
 - name: Oskar Klein Center for Cosmoparticle Physics, Stockholm University, Stockholm
   index: 1
 - name: Center for Computational Astrophysics, Flatiron Institute, New York
   index: 2
 - name: Imperial Centre for Inference and Cosmology, Imperial College London, London
   index: 3
date: 19 December 2018
bibliography: paper.bib
---

# Summary

Regression analysis of atmospheric time-series observations is a key endeavour for identifying long-term trends and studying underlying drivers of variability in the data. Multiple linear regression (MLR) has been a standard tool for climate and atmospheric time-series data analysis. However, MLR has a number of well-known shortcomings that can, in the worst case, lead to biased scientific inferences being drawn from the data. Dynamical linear modeling (DLM) provides a more flexible regression framework that addresses a number of the key issues faced by MLR (reviewed below), providing an attractive and more robust alternative.

``DLMMC`` provides a suite of DLM models for atmospheric time-series analysis, with a user friendly python interface, for use by the community. Models are implemented in stan [@stan], and use a combination of efficient Hamiltonian Monte Carlo (HMC) sampling and Kalman filtering to explore the (Bayesian) posterior distribution of the model parameters given the data.

The basic set of DLM models are based on [@laine2014] and have four main components: a dynamic seasonal cycle (with 6- and 12-month components), a smooth non-linear background trend, forcings from any number of user defined proxy variables (eg., solar activity, quasi-biennial oscillations, the El-Nino southern oscillation, etc), and an auto-regressive (AR) process. These models were originally built with stratospheric research in mind, but can be applied or modified for other climate/atmospheric time-series analysis (or more general time-series analysis problems).

The DLM models in this suite can be readily extended to more complex models with additional components; ``DLMMC`` provides a general-purpose and efficient HMC sampling-Kalman filtering framework for DLM regression analysis.

The code has already been used in a number of scientific publications, including [@ball2017] and [@ball2018]. For an introduction to DLM in the context of stratospheric time-series analysis, see [@laine2014] (on which this code was inspired), or see [@durbinkoopman2012] for a more comprehensive review of time-series analysis with state-space models.

## Advantages of DLM over traditional MLR methods

Dynamical linear models provide a rich framework for time-series regression analysis. Let us briefly review some of the key advantages offered by DLM over commonly used MLR approaches (in the context of stratospheric time-series analysis):

*Flexible non-linear background trends*<br/>
DLM uses a flexible (non-parametric) model for the non-linear background trend, where the degree of non-linearity is included as a free parameter that is fit along with the rest of the regression model. That way, the data are free to choose how rapidly (in time) the background trend is allowed to vary. In contrast, independent-linear trends (ILT) or piecewise-linear trends (PWLT) that are commonplace in MLR trend analyses provide very restrictive models for the background trend. For example, typically the positions of any inflection point(s) in the ILT/PWLT must be fixed in advance even though they are not necessarily known a priori, and the resulting trends are sensitive to these choices. Such restrictive trend-model assumptions such as ILT or PWLT severely hamper our ability to let the data speak for themselves in terms of itentifying long-term background trends, particularly when those trends might be surprising relative to our prior expectations.

*Dynamical seasonal and regressor modulation*<br/>
DLM allows the amplitudes of the seasonal cycle and forcings via proxy variables (such as solar activity, quasi-biennial oscillations, the El-Nino southern oscillation, etc) to vary dynamically with time, whereas MLR has fixed (in time) regression coefficients. This extra flexibility in the model enables DLM to capture more of the variability in the data due to the evolving state of the atmosphere, or evolution in the observing conditions and/or sampling that can lead to changes in the distributional properties of the data with time.

*Principled treatment of auto-regressive processes*<br/>
In many situations, phenomenological regression models (such as DLM and MLR) only capture some fraction of the total variability in the data, ie., there are correlated processes left in the residuals that are not captured by the main features of the regression model. This is often mitigated by including an auto-regressive (AR) process as a surrogate for the un-modelled variability. DLM infers an auto-regressive process simultaneously with the rest of the model parameters, and carefully propagates uncertainties by formally marginalizing over uncertainties in the AR-process when reporting errors on the recovered trend, seasonal cycle etc. In contrast, MLR typically applies a post-hoc correction to the recovered regression coefficients and error bars to account for auto-regressive correlations (eg., the Cochrane–Orcutt correction); this is often approximate, and less principled in its propagation of uncertainties.

*Principled treatment of time-varying uncertainties (heteroscedasticity)*<br/>
Finally, many implementations of MLR in the literature use ordinary least squares (OLS) estimation for the regression coefficients. This approach relies on the assumption that the error distribution on observations is constant in time (homoscedasticity), and can lead to biased parameter estimates when this condition is not met (heteroscedasticity). DLM implements heteroscedastic time-varying uncerainties as standard.

# Acknowledgements

I thank Will Ball, Daniel Mortlock, Marko Laine and Sean Davis for useful discussions and comments on the code.

# References
