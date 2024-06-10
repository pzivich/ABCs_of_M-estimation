####################################################################################################################
# ABC's of M-estimation
#   Code for applied examples (Section 2)
#
# Paul Zivich (2024/06/10)
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
# Example 3:Transport to External Target

# Adding population indicator to ZAPPS
d['zapps'] = 1

# Target Population Data
d_target = pd.DataFrame()
d_target['X'] = [0, 0, 1, 1]
d_target['W'] = [0, 1, 0, 1]
d_target['n'] = [300, 300, 300, 100]
d_target['zapps'] = 0
d_target['intercept'] = 1
d_target = pd.DataFrame(np.repeat(d_target.values, d_target['n'], axis=0),  # Expanding compact data frame
                        columns=d_target.columns)                           # ... keeping column names
d_target = d_target[['intercept', 'X', 'W', 'zapps']].copy()                # Dropping the n column

# Stacking data together
df = pd.concat([d, d_target], ignore_index=True)

############################################################################################
# Example 3a: Standardization by g-computation

# To handle missing, we will fill in a placeholder value, -999. This
#   approach is compatible with the root-finding algorithm and allows
#   us to easily check that only the correct observations are contributing
#   to the nuisance model (it will have trouble converging if these contribute)
y_no_nan = np.nan_to_num(df['Y'], nan=-999)
s = np.asarray(df['zapps'])

############################################
# Using delicatessen


def psi(theta):
    # Dividing parameters into corresponding parts and labels from slides
    beta = theta[0:3]   # Logistic model coefficients
    mu = theta[3]       # Risk

    # Using built-in regression model functionality from delicatessen
    ee_logit = ee_regression(theta=beta,
                             y=y_no_nan,
                             X=df[['intercept', 'X', 'W']],
                             model='logistic')
    ee_logit = ee_logit * s  # Restricting contributions to ZAPPS

    # Transforming logistic model coefficients into causal parameters
    y_hat = inverse_logit(np.dot(df[['intercept', 'X', 'W']], beta))  # Predictions
    # Estimating function for causal risk under a=1
    ee_r1 = (1-s) * (y_hat - mu)

    # Returning stacked estimating functions in order of parameters
    return np.vstack([ee_logit,   # EF of logistic model
                      ee_r1])     # EF of risk in target


mestr = MEstimator(psi, init=[0, 0, 0, 0.5])
mestr.estimate(solver='lm')

print("3a: Transport -- g-computation")
print("Risk:   ", np.round(mestr.theta[-1], 3))
print("95% CI: ", np.round(mestr.confidence_intervals()[-1, :], 3))
print("")

############################################################################################
# Example 3b: Standardization by IOSW

############################################
# Using delicatessen


def psi(theta):
    # Dividing parameters into corresponding parts and labels from slides
    alpha = theta[0:3]   # Logistic model coefficients
    mu = theta[3]       # Risk

    # Using built-in regression model functionality from delicatessen
    ee_logit = ee_regression(theta=alpha,
                             y=s,
                             X=df[['intercept', 'X', 'W']],
                             model='logistic')

    # Computing IOSW
    pi = inverse_logit(np.dot(df[['intercept', 'X', 'W']], alpha))  # Probablility score
    wt = (1-pi)/pi                                                  # Corresponding weights

    # Estimating function for causal risk under a=1
    ee_r1 = s * wt * (y_no_nan - mu)

    # Returning stacked estimating functions in order of parameters
    return np.vstack([ee_logit,   # EF of logistic model
                      ee_r1])     # EF of risk in target


mestr = MEstimator(psi, init=[0, 0, 0, 0.5])
mestr.estimate(solver='lm')

print("3b: Transport -- IOSW")
print("Risk:         ", np.round(mestr.theta[-1], 3))
print("95% CI:         ", np.round(mestr.confidence_intervals()[-1, :], 3))

####################################################################################################################
# EXPECTED OUTPUT (your package versions may differ):

# versions
# NumPy:         1.25.2
# Pandas:        1.4.1
# SciPy:         1.11.2
# Statsmodels:   0.14.1
# Delicatessen:  2.2

# 1: Logistic Regression
# Point Estimates
# MLE:          [-1.89450082  0.11873535  0.36051133]
# By-hand:      [-1.89450082  0.11873535  0.36051133]
# Delicatessen: [-1.89450082  0.11873535  0.36051133]
# Variance Estimates
# MLE:          [0.12231295 0.27865458 0.23791706]
# By-hand:      [0.12182124 0.27878379 0.2377596 ]
# Delicatessen: [0.12182124 0.27878379 0.2377596 ]

# 2a: Causal parameters -- g-computation
# Risk 1:          0.154
# 95% CI:          [0.089 0.22 ]
# Risk 0:          0.14
# 95% CI:          [0.114 0.165]
# Risk Difference: 0.015
# 95% CI:          [-0.055  0.085]
# Risk Ratio:      1.106
# 95% CI:          [0.697 1.756]

# 2b: Causal parameters -- IPW
# Risk 1:          0.153
# 95% CI:          [0.088 0.219]
# Risk 0:          0.14
# 95% CI:          [0.114 0.165]
# Risk Difference: 0.014
# 95% CI:          [-0.057  0.084]
# Risk Ratio:      1.098
# 95% CI:          [0.69  1.747]

# 3a: Transport -- g-computation
# Risk:    0.155
# 95% CI:  [0.121 0.19 ]

# 3b: Transport -- IOSW
# Risk:          0.155
# 95% CI:          [0.115 0.195]

####################################################################################################################
# END
