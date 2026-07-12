# Python code for function 
import numpy as np
from scipy.optimize import *

def f(x, mu, sigma,a):
    x=np.array(x)
    
    res_1 = ((a*np.exp(-((np.log(x[:,0])-  mu)**2)/(2*sigma**2))))
    #res_1 = np.log(1+np.exp(res_1))
    return res_1

def generate_Params(x,y):
    x = np.array(x)
    y = np.array(y)
    #DC_W_UL =  np.ones(x.shape[1]-1) * np.inf
    #indices = np.where(np.sum(x[:,1:],axis=0) == 0)[0] # find the empty drugs concentration
    #DC_W_UL[indices]=0
    params, extras = curve_fit(f, x, y,bounds=([0,0,0],[5.82,np.inf,60]), maxfev=100000)
    
    return params

def f_log(x, mu, sigma,a) :
    return ((a)*np.exp(-((np.log(x)- 
   mu)**2)/(2*sigma**2))) 