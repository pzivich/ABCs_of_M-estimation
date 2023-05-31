####################################################################################################################
# ABC's of M-estimation
#   Code for applied examples (Section 2)
#
# Paul Zivich (2023/05/31)
####################################################################################################################

############################################
# Loading dependencies

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
import scipy as sp
from scipy.optimize import root
from scipy.optimize import approx_fprime
import delicatessen as deli
from delicatessen import MEstimator
from delicatessen.estimating_equations import ee_regression
from delicatessen.utilities import inverse_logit

print("versions")
print('NumPy:        ', np.__version__)
print('Pandas:       ', pd.__version__)
print('SciPy:        ', sp.__version__)
print('Statsmodels:  ', sm.__version__)
print('Delicatessen: ', deli.__version__, '\n')

############################################################################################
# Example 1: Logistic Regression

############################################
# Loading data

d = pd.DataFrame()
d['X'] = [0, 0, 0, 0, 1, 1, 1, 1]
d['W'] = [0, 0, 1, 1, 0, 0, 1, 1]
d['Y'] = [0, 1, 0, 1, 0, 1, 0, 1]
d['n'] = [496, 74, 113, 25, 85, 15, 15, 3]
d['intercept'] = 1
d = pd.DataFrame(np.repeat(d.values, d['n'], axis=0),   # Expanding compact data frame
                 columns=d.columns)                     # ... keeping column names
d = d[['intercept', 'X', 'W', 'Y']].copy()              # Dropping the n column
n = d.shape[0]                                          # Number of observations

############################################
# Regression by MLE

fam = sm.families.Binomial()           # Binomial family
mle = smf.glm("Y ~ X + W",             # GLM with formula
              data=d,                  # ... for loaded data
              family=fam               # ... logistic model
              ).fit()                  # ... then fitting model
beta_mle = np.asarray(mle.params)      # Extracting beta parameters
sigma_mle = np.asarray(mle.bse**2)     # Extracting variance estimates

############################################
# Defining estimating equations


def estimating_functions(theta):
    # Calculating predicted probability of Y
    yhat = inverse_logit(theta[0]*1 + theta[1]*d['X'] + theta[2]*d['W'])

    # Calculating the residual (observed Y minus predicted Y)
    residual = d['Y'] - yhat

    # Calculating the corresponding score functions
    score = [residual*1,         # Score for intercept
             residual*d['X'],    # Score for X
             residual*d['W']]    # Score for W

    # Returning scores as a NumPy array
    return np.asarray(score)


def estimating_equations(theta):
    ef = estimating_functions(theta=theta)
    vals = ()                    # Create empty tuple
    rows = ef.shape[0]           # Determine how many rows / parameters are present
    for i in range(rows):        # Go through each individual theta in the stack
        row = ef[i, :]           # ... extract corresponding row
        vals += (np.sum(row), )  # ... then add the theta sum to the tuple of thetas
    # Return evaluated estimating equations
    return np.asarray(vals)


############################################
# Root-finding

proc = root(estimating_equations,       # Function to find root(s) of
            x0=np.array([0, 0, 0, ]),   # ... starting values for root-finding procedure
            method='lm')                # ... algorithm to use (Levenberg-Marquardt here)
beta_mest = proc.x                      # Extract beta parameters

############################################
# Baking the bread (approximate derivative)

deriv = approx_fprime(xk=beta_mest,             # Approximating partial derivatives at beta
                      f=estimating_equations,   # ... of the estimating equations
                      epsilon=1e-9)             # ... with specified deviation from point
bread = -1*deriv / n                            # Creating bread

############################################
# Cooking the filling (matrix algebra)

filling = np.dot(estimating_functions(beta_mest),
                 estimating_functions(beta_mest).T)
filling = filling / n

############################################
# Assembling the sandwich (matrix algebra)

bread_invert = np.linalg.inv(bread)
sandwich = np.dot(np.dot(bread_invert, filling),
                  bread_invert.T) / n
sigma_mest = np.diag(sandwich)

############################################
# Using delicatessen


def psi(theta):
    # Using built-in regression model functionality from delicatessen
    return ee_regression(theta=theta,
                         y=d['Y'],
                         X=d[['intercept', 'X', 'W']],
                         model='logistic')


mestr = MEstimator(psi, init=[0, 0, 0])
mestr.estimate(solver='lm')

beta_deli = mestr.theta
sigma_deli = np.diag(mestr.variance)

############################################
# Comparing results

print("1: Logistic Regression")
print("Point Estimates")
print("MLE:         ", beta_mle)
print("By-hand:     ", beta_mest)
print("Delicatessen:", beta_deli)

print("Variance Estimates")
print("MLE:         ", sigma_mle ** 0.5)
print("By-hand:     ", sigma_mest ** 0.5)
print("Delicatessen:", sigma_deli ** 0.5)
print("")

############################################################################################
# Example 2a: Standardization by g-computation

# Copies of data with policies applied
d1 = d.copy()
d1['X'] = 1
d0 = d.copy()
d0['X'] = 0

############################################
# Using delicatessen


def psi(theta):
    # Dividing parameters into corresponding parts and labels from slides
    beta = theta[0:3]                     # Logistic model coefficients
    mu1, mu0 = theta[3], theta[4]         # Causal risks
    delta1, delta2 = theta[5], theta[6]   # Causal contrasts

    # Using built-in regression model functionality from delicatessen
    ee_logit = ee_regression(theta=beta,
                             y=d['Y'],
                             X=d[['intercept', 'X', 'W']],
                             model='logistic')

    # Transforming logistic model coefficients into causal parameters
    y1_hat = inverse_logit(np.dot(d1[['intercept', 'X', 'W']], beta))  # Prediction under a=1
    y0_hat = inverse_logit(np.dot(d0[['intercept', 'X', 'W']], beta))  # Prediction under a=0
    # Estimating function for causal risk under a=1
    ee_r1 = y1_hat - mu1
    # Estimating function for causal risk under a=0
    ee_r0 = y0_hat - mu0
    # Estimating function for causal risk difference
    ee_rd = np.ones(d.shape[0])*((mu1 - mu0) - delta1)
    # Estimating function for causal risk ratio
    ee_rr = np.ones(d.shape[0])*(np.log(mu1 / mu0) - delta2)

    # Returning stacked estimating functions in order of parameters
    return np.vstack([ee_logit,   # EF of logistic model
                      ee_r1,      # EF of causal risk a=1
                      ee_r0,      # EF of causal risk a=0
                      ee_rd,      # EF of causal risk difference
                      ee_rr])     # EF of causal log risk ratio


mestr = MEstimator(psi, init=[0, 0, 0, 0.5, 0.5, 0, 0])
mestr.estimate(solver='lm')

print("2a: Causal parameters -- g-computation")
print("Risk 1:         ", np.round(mestr.theta[3], 3))
print("95% CI:         ", np.round(mestr.confidence_intervals()[3, :], 3))
print("Risk 0:         ", np.round(mestr.theta[4], 3))
print("95% CI:         ", np.round(mestr.confidence_intervals()[4, :], 3))
print("Risk Difference:", np.round(mestr.theta[5], 3))
print("95% CI:         ", np.round(mestr.confidence_intervals()[5, :], 3))
print("Risk Ratio:     ", np.round(np.exp(mestr.theta[6]), 3))
print("95% CI:         ", np.round(np.exp(mestr.confidence_intervals()[6, :]), 3))
print("")

############################################################################################
# Example 2b: Standardization by IPW

############################################
# Using delicatessen


def psi(theta):
    # Dividing parameters into corresponding parts and labels from slides
    alpha = theta[0:2]                    # Logistic model coefficients
    mu1, mu0 = theta[2], theta[3]         # Causal risks
    delta1, delta2 = theta[4], theta[5]   # Causal contrasts

    # Using built-in regression model functionality from delicatessen
    ee_logit = ee_regression(theta=alpha,
                             y=d['X'],
                             X=d[['intercept', 'W']],
                             model='logistic')

    # Transforming logistic model coefficients into causal parameters
    pscore = inverse_logit(np.dot(d1[['intercept', 'W']], alpha))  # Propensity score
    wt = d['X']/pscore + (1-d['X'])/(1-pscore)                     # Corresponding weights
    # Estimating function for causal risk under a=1
    ee_r1 = d['X']*d['Y']*wt - mu1
    # Estimating function for causal risk under a=0
    ee_r0 = (1-d['X'])*d['Y']*wt - mu0
    # Estimating function for causal risk difference
    ee_rd = np.ones(d.shape[0])*((mu1 - mu0) - delta1)
    # Estimating function for causal risk ratio
    ee_rr = np.ones(d.shape[0])*(np.log(mu1 / mu0) - delta2)

    # Returning stacked estimating functions in order of parameters
    return np.vstack([ee_logit,   # EF of logistic model
                      ee_r1,      # EF of causal risk a=1
                      ee_r0,      # EF of causal risk a=0
                      ee_rd,      # EF of causal risk difference
                      ee_rr])     # EF of causal log risk ratio


mestr = MEstimator(psi, init=[0, 0, 0.5, 0.5, 0, 0])
mestr.estimate(solver='lm')

print("2b: Causal parameters -- IPW")
print("Risk 1:         ", np.round(mestr.theta[2], 3))
print("95% CI:         ", np.round(mestr.confidence_intervals()[2, :], 3))
print("Risk 0:         ", np.round(mestr.theta[3], 3))
print("95% CI:         ", np.round(mestr.confidence_intervals()[3, :], 3))
print("Risk Difference:", np.round(mestr.theta[4], 3))
print("95% CI:         ", np.round(mestr.confidence_intervals()[4, :], 3))
print("Risk Ratio:     ", np.round(np.exp(mestr.theta[5]), 3))
print("95% CI:         ", np.round(np.exp(mestr.confidence_intervals()[5, :]), 3))
print("")

############################################################################################
# Example 3: Data Fusion

d = pd.DataFrame()
d['R'] = [1, 1, 0, 0, 0, 0]
d['Y'] = [0, 0, 1, 1, 0, 0]
d['W'] = [1, 0, 1, 0, 1, 0]
d['n'] = [680, 270, 204, 38, 18, 71]
d['intercept'] = 1
d = pd.DataFrame(np.repeat(d.values, d['n'], axis=0),   # Expanding compact data frame
                 columns=d.columns)                     # ... keeping column names
d = d[['intercept', 'R', 'W', 'Y']].copy()              # Dropping the n column
n = d.shape[0]                                          # Number of observations

r = np.asarray(d['R'])
w = np.asarray(d['W'])
y = np.asarray(d['Y'])

############################################
# Using delicatessen


def psi(theta):
    # theta[0]: naive mean, theta[1]: sensitivity, theta[2]: specificity, theta[3]: corrected mean
    ee_1 = r*(w - theta[0])                                                                    # EF naive mean
    ee_2 = (1-r) * y * (w - theta[1])                                                          # EF sensitivity
    ee_3 = (1-r) * (1-y) * ((1-w) - theta[2])                                                  # EF specificity
    ee_4 = np.ones(y.shape[0])*theta[3]*(theta[1] + theta[2] - 1) - (theta[0] + theta[2] - 1)  # EF corrected mean

    # Returning stacked estimating functions in order of parameters
    return np.vstack([ee_1,      # EF naive mean
                      ee_2,      # EF sensitivity
                      ee_3,      # EF specificity
                      ee_4])     # EF corrected mean


mestr = MEstimator(psi, init=[0.5, 0.75, 0.75, 0.5])
mestr.estimate(solver='lm')

print("3: Fusion")
print("Uncorrected Mean:", np.round(mestr.theta[0], 3))
print("95% CI:          ", np.round(mestr.confidence_intervals()[0, :], 3))
print("Sensitivity:     ", np.round(mestr.theta[1], 3))
print("95% CI:          ", np.round(mestr.confidence_intervals()[1, :], 3))
print("Specificity:     ", np.round(mestr.theta[2], 3))
print("95% CI:          ", np.round(mestr.confidence_intervals()[2, :], 3))
print("Corrected Mean:  ", np.round(mestr.theta[3], 3))
print("95% CI:          ", np.round(mestr.confidence_intervals()[3, :], 3))
print("")

####################################################################################################################
# EXPECTED OUTPUT (versions may differ):

# versions
# NumPy:         1.22.2
# Pandas:        1.4.1
# SciPy:         1.9.2
# Statsmodels:   0.13.2
# Delicatessen:  1.2
#
# 1: Logistic Regression
# Point Estimates
# MLE:          [-1.89450082  0.11873535  0.36051133]
# By-hand:      [-1.89450082  0.11873535  0.36051133]
# Delicatessen: [-1.89450082  0.11873535  0.36051133]
# Variance Estimates
# MLE:          [0.01496046 0.07764837 0.05660453]
# By-hand:      [0.01484041 0.0777204  0.05652963]
# Delicatessen: [0.01484041 0.0777204  0.05652963]
#
# 2a: Causal parameters -- g-computation
# Risk 1:          0.154
# 95% CI:          [0.089 0.22 ]
# Risk 0:          0.14
# 95% CI:          [0.114 0.165]
# Risk Difference: 0.015
# 95% CI:          [-0.055  0.085]
# Risk Ratio:      1.106
# 95% CI:          [0.697 1.756]
#
# 2b: Causal parameters -- IPW
# Risk 1:          0.153
# 95% CI:          [0.088 0.219]
# Risk 0:          0.14
# 95% CI:          [0.114 0.165]
# Risk Difference: 0.014
# 95% CI:          [-0.057  0.084]
# Risk Ratio:      1.098
# 95% CI:          [0.69  1.747]
#
# 3: Fusion
# Uncorrected Mean: 0.716
# 95% CI:           [0.687 0.744]
# Sensitivity:      0.843
# 95% CI:           [0.797 0.889]
# Specificity:      0.798
# 95% CI:           [0.714 0.881]
# Corrected Mean:   0.801
# 95% CI:           [0.724 0.879]

####################################################################################################################
# END
