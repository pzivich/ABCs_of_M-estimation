####################################################################################################################
# ABC's of Estimating Equations
#   Code for M-estimator of the corrected proportion (Part 1)
#
# Paul Zivich (2026/06/16)
####################################################################################################################

############################################
# Loading dependencies

import numpy as np
import pandas as pd
import scipy as sp
from scipy.optimize import root
from scipy.optimize import approx_fprime
import delicatessen as deli
from delicatessen import MEstimator

print("versions")
print('NumPy:        ', np.__version__)
print('Pandas:       ', pd.__version__)
print('SciPy:        ', sp.__version__)
print('Delicatessen: ', deli.__version__, '\n')

############################################
# Data

d = pd.DataFrame()
d['Y*'] = [0, 1, 0, 1]
d['R'] = [1, 1, 0, 0]
d['n'] = [120, 80, 15, 85]
d = pd.DataFrame(np.repeat(d.values, d['n'], axis=0),   # Expanding compact data frame
                 columns=d.columns)                     # ... keeping column names
d = d[['R', 'Y*']].copy()                               # Dropping the n column
n = d.shape[0]                                          # Number of observations

############################################
# By-hand Calculation

mu_star = np.mean(d.loc[d['R'] == 1, 'Y*'])    # Manually calculate the naive proportion
alpha = np.mean(d.loc[d['R'] == 0, 'Y*'])      # Manually calculate the sensitivity
mu = mu_star / alpha                           # Manually calculate the corrected proportion
print("By-Hand")
print("Naive Proportion:    ", mu_star)
print("Sensitivity:         ", alpha)
print("Corrected Proportion:", mu)

############################################
# Defining estimating equation

# Pulling out needed variables from data set
ystar = np.asarray(d['Y*'])
r = np.asarray(d['R'])

def estimating_function(theta):
    mu_tilde, alpha_tilde = theta

    # Parameter-specific estimating functions
    ef_mean = r * (ystar - alpha_tilde*mu_tilde)
    ef_sens = (1-r) * (ystar - alpha_tilde)

    # Stacking the estimating functions together into vectors
    return np.vstack([ef_mean, ef_sens])


def estimating_equation(theta):
    # Summing together the estimating functions into the estimating equation
    estf = np.asarray(estimating_function(theta))  # Return estimating function
    return np.sum(estf, axis=1)                    # Sum over all estimating functions


############################################
# Root-finding

proc = root(estimating_equation,          # Function to find root(s) of
            x0=np.array([0.5, 0.5]),      # ... starting values for root-finding procedure
            method='lm')                  # ... algorithm to use (Levenberg-Marquardt here)
theta_root = proc["x"]                    # Extract the parameter estimates from root-finder object

############################################
# Baking the bread (approximate derivative)

deriv = approx_fprime(xk=proc["x"],            # Array of values to compute derivative at (root of estimating equation)
                      f=estimating_equation,   # ... function to find derivative of
                      epsilon=1e-9)            # ... distance of points for numerical approximation (should be small)
bread = -1*deriv / n

############################################
# Cooking the filling (matrix algebra)

filling = np.dot(estimating_function(theta_root), estimating_function(theta_root).T)
filling = filling / n

############################################
# Assembling the sandwich (matrix algebra)

bread_inv = np.linalg.inv(bread)
sandwich = np.dot(np.dot(bread_inv, filling), bread_inv.T) / n
se = np.sqrt(np.diag(sandwich))

print("Correct Proportion -- EE")
print("Proportion:", theta_root[0])
print("SE:        ", np.round(se[0], 5))
print("95% CI:    ", np.round([theta_root[0] - 1.96*se[0],
                               theta_root[0] + 1.96*se[0]],
                              3))

####################################################################################################################
# Using delicatessen instead of by-hand

mestr = MEstimator(estimating_function, init=[0.5, 0.5])
mestr.estimate(solver='lm')

print("Corrected Proportion -- delicatessen")
print("theta: ", mestr.theta)
print("95% CI:", mestr.confidence_intervals())

####################################################################################################################
# EXPECTED OUTPUT (versions may differ):

# versions
# NumPy:         2.3.5
# Pandas:        2.3.3
# SciPy:         1.16.3
# Delicatessen:  4.3
#
# By-Hand
# Naive Proportion:     0.4
# Sensitivity:          0.85
# Corrected Proportion: 0.4705882352941177
# Correct Proprtion -- EE
# Proportion: 0.4705882352941177
# SE:         0.0453
# 95% CI:     [0.382 0.559]
# Corrected Proportion -- delicatessen
# theta:  [0.47058824 0.85      ]
# 95% CI: [[0.38181031 0.55936616]
#  [0.78001529 0.91998471]]

####################################################################################################################
# END
