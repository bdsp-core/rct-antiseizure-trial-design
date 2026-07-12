# Python code for function 
import numpy as np
from scipy.optimize import *

def f(x, mu, sigma,a,drug_1,drug_2,drug_3,drug_4,drug_5,drug_6,drug_7) :
    x=np.array(x)
    drug_coeff = [drug_1,drug_2,drug_3,drug_4,drug_5,drug_6,drug_7]
    res_1 = (((a*np.exp(-((np.log(x[:,0])-  mu)**2)/(2*sigma**2))))-np.sum(drug_coeff*x[:,1:],axis=1))
    res_1 = np.log(1+np.exp(res_1))
    return res_1

def generate_Params(x,y):
    x = np.array(x)
    y = np.array(y)
    #DC_W_UL =  np.ones(x.shape[1]-1) * np.inf
    indices = np.where(np.sum(x[:,1:],axis=0) == 0)[0] # find the empty drugs concentration
    #DC_W_UL[indices]=0
    params, extras = curve_fit(f, x, y,bounds=([0,0,0,0,0,0,0,0,0,0],[5.82,np.inf,60,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf]), maxfev=100000)
    params[indices+3]=0 # drug indices starts from postion 3 with 0 be mean, 1 be sigma, 2 be amplitude
    return params

def f_log(x, mu, sigma,a) :
    return ((a)*np.exp(-((np.log(x)- 
   mu)**2)/(2*sigma**2))) 