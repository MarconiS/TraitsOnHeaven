import sklearn
import csv
import sys
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import scipy 

#find a way to import them from R
spectra = pd.read_csv('/Users/sergiomarconi/OneDrive/Projects/OSBS/inputs/Bootstrap_2/onePix1Crown_N_pct94.csv')
traits = pd.read_csv('/Users/sergiomarconi/OneDrive/Projects/OSBS/inputs/Spectra/CrownTraits.csv')
test_spectra
train_spectra
test_traits
train_traits
specificTrait

def rsquared(x, y):
    """ Return R^2 where x and y are array-like."""

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    return r_value**2
	
#traits = traits.drop(traits.index[len(traits)-1]) 
average_spectra = spectra.groupby(['pixel_crownID'], as_index = False).mean()
sd_spectra = spectra.groupby(['pixel_crownID'], as_index = False).std()
sd_spectra['pixel_crownID'] = average_spectra['pixel_crownID'] 
scaled_spectra = average_spectra - average_spectra.mean()
scaled_spectra['pixel_crownID'] = average_spectra['pixel_crownID'] 
scaled_spectra = scaled_spectra.dropna(axis=1)
from sklearn import gaussian_process
from sklearn.gaussian_process.kernels import Matern, WhiteKernel, ConstantKernel

#find a way to import it from R
kernel = ConstantKernel() + Matern(length_scale=2, nu=3/2) + WhiteKernel(noise_level=1)

X = train_spectra.drop("pixel_crownID", 1)
Y = train_traits['N_pct']
gp = gaussian_process.GaussianProcessRegressor(kernel=kernel)
gp.fit(X,Y)
gp.kernel_
R2_train  = gp.score(X,Y)

x_test = test_spectra.drop("pixel_crownID", 1)
y_pred, sigma = gp.predict(x_test, return_std=True)
R2(rsquared(y_pred,test_traits['N_pct']))


