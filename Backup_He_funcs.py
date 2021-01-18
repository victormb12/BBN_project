#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 15:34:42 2020

@author: victormaura12
"""

import numpy as np
from scipy.integrate import quad_vec,solve_ivp
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import interp1d
from matplotlib import rcParams
from scipy.special import zeta
import time
from scipy.optimize import brentq 
import geff_with_neff as ge
import Helium_functions as He

rcParams['lines.linewidth'] = 1.5
rcParams['font.size']=11
plt.style.use('seaborn-ticks')


def func_from_pkl(filename):
    with open(filename+'_func.pickle', 'rb') as handle:
        cases=pickle.load(handle)
    return cases

def func_to_pkl(dictionary,filename):
     with open(filename+'_func.pickle', 'wb') as handle:
        pickle.dump(dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)
        

"""

-------------------------------------------------------------------------------------------------------------------------------------
Calculation of the weak interaction rates accroding to https://arxiv.org/pdf/1812.05605.pdf 
-------------------------------------------------------------------------------------------------------------------------------------

"""

def Mev2_to_cm3s(a, b=0):
    
    conv=((1.97e-11)**3)/(6.58e-28)
    if b==0:
        return a*conv
    elif b==1:
        return a/conv

def MeV_to_1_over_s(T,b=1):
    """

    Parameters
    ----------
    T : Temperature in MeV to be turned into Hz=1/s or vice versa
    b : if b ==1, then the temperature gets turned into Hertz. Else, Hertz gets turned into MeV

    Returns
    -------
    float or array-like
        the Temperature or frequency converted into the right units

    """
    a=6.5821e-22
    if b==1:
        return T/a
    else:
        return T*a





def n_to_p_int(epsilon,z_gamma,sign=1,xi_nu=0):
    T_over_nu=func_from_pkl('T_over_T_nu')
    T=0.511/z_gamma
    K=MeV_to_1_over_s(1/(1.636*880.2),b=0)
    q=1.2933/0.511
    func1=((epsilon-sign*q)**2*np.sqrt((epsilon**2-1))*epsilon)/((1+np.exp(-epsilon*z_gamma))*(1+np.exp((epsilon-sign*q)*z_gamma*T_over_nu(T)-sign*xi_nu)))
    func2=((epsilon+sign*q)**2*np.sqrt((epsilon**2-1))*epsilon)/((1+np.exp(epsilon*z_gamma))*(1+np.exp(-(epsilon+sign*q)*z_gamma*T_over_nu(T)-sign*xi_nu)))
    return K*(func1+func2)


def gamma_p_to_n(z_gamma):
    integral=quad_vec(lambda epsilon: n_to_p_int(epsilon,z_gamma,sign=-1),1,np.inf)
    return integral[0]
    
def gamma_n_to_p(z_gamma):
    integral=quad_vec(lambda epsilon: n_to_p_int(epsilon,z_gamma,sign=1),1,np.inf)
    return integral[0]

def gamma_pn_spline(filename,z_range=(-4,3)):
    z=np.logspace(z_range[0],z_range[1],10000)
    func=gamma_p_to_n(z)
    f2 = interp1d(z, func, kind='cubic')
    func_to_pkl(f2,filename)
    
def gamma_np_spline(filename,z_range=(-4,3)):
    z=np.logspace(z_range[0],z_range[1],10000)
    func=gamma_n_to_p(z)
    f2 = interp1d(z, func, kind='cubic')
    func_to_pkl(f2,filename)
    
def gamma_H_z_pn(filename,z_range,Neff=3.045,xi_nu=0):
    g_eff= lambda T: ge.g_eff_tot(T,ge.N_eff_SMC(T,N_nu=Neff,xi_nu=xi_nu))
    zi=np.logspace(z_range[0],z_range[1],1000)
    func=[gamma_p_to_n(z,xi_nu)/(1.66*np.sqrt(g_eff(0.511/z))*(0.511)**2/1.22e22)*z for z in zi]
    f2 = interp1d(zi, func, kind='cubic')
    func_to_pkl(f2,filename)
    
def gamma_H_z_np(filename,z_range,Neff=3.045,xi_nu=0):
    g_eff= lambda T: ge.g_eff_tot(T,ge.N_eff_SMC(T,N_nu=Neff,xi_nu=xi_nu))
    gamma_np=func_from_pkl('gamma_np')
    zi=np.logspace(z_range[0],z_range[1],1000)
    func=[gamma_n_to_p(z,xi_nu)/(1.66*np.sqrt(g_eff(0.511/z))*(0.511)**2/1.22e22)*z for z in zi]
    f2 = interp1d(zi, func, kind='cubic')
    func_to_pkl(f2,filename)

"""

-------------------------------------------------------------------------------------------------------------------------------------
ODE solver class from my Bachelor's thesis
-------------------------------------------------------------------------------------------------------------------------------------

"""


class ODE(object):
    
    def __init__(self, function, y0, t_span, t_eval=None, params=None, method='Radau',rtol=1e-3, atol=1e-6,dense_output=False, verbose = False,vectorized=False):
        self.function=function
        self.params=params
        self.y0=y0
        self.t_span=t_span
        self.t_eval=t_eval
        self.method=method
        self.rtol=rtol
        self.atol=atol
        self.dense_output=dense_output
        self.verbose=verbose
        self.vectorized=vectorized
        self.func_sol=self.solve()
        self.Y_inf=self.func_sol.y[0][-1]
        
        
    def solve(self):
        
        solved_eq=solve_ivp(self.function,self.t_span,self.y0, t_eval=self.t_eval, args=self.params, dense_output=self.dense_output, method=self.method, rtol=self.rtol, atol=self.rtol,vectorized=self.vectorized)
        if self.verbose:
            print(solved_eq.message)
        return solved_eq
    
    def plot(self,label=None,ax=None,color=None, linestyle='-'):
        
        t=self.func_sol.t
        y=self.func_sol.y
        self.Y_inf=y[0][-1]
        
        if ax==None:
            plt.plot(t,y.T,label=label,color=color,linestyle=linestyle)
        else:
            ax.plot(t,y.T,label=label,color=color)

"""

-------------------------------------------------------------------------------------------------------------------------------------
RHS of the Boltzmann equation. The ratio between gamma and H is calculated before hand and stored as a spline so 
that the values do not span as many orders of magnitude. 
-------------------------------------------------------------------------------------------------------------------------------------

"""
g_eff= lambda T: ge.g_eff_tot(T,ge.N_eff_SMC(T,N_nu=3.045,xi_nu=0))

# g_eff= lambda T: ge.g_eff_tot(T,4)
            
# def gamma_H_z_pn(filename,z_range):
#     gamma_pn=func_from_pkl('gamma_pn')
#     zi=np.logspace(z_range[0],z_range[1],1000)
#     func=[gamma_pn(z)/(1.66*np.sqrt(g_eff(0.511/z))*(0.511)**2/1.22e22)*z for z in zi]
#     f2 = interp1d(zi, func, kind='cubic')
#     func_to_pkl(f2,filename)
    
# def gamma_H_z_np(filename,z_range):
#     gamma_np=func_from_pkl('gamma_np')
#     zi=np.logspace(z_range[0],z_range[1],1000)
#     func=[gamma_np(z)/(1.66*np.sqrt(g_eff(0.511/z))*(0.511)**2/1.22e22)*z for z in zi]
#     f2 = interp1d(zi, func, kind='cubic')
#     func_to_pkl(f2,filename)
            
def Boltzmann_for_n(z_gamma,X_n):
    gammaHz_pn=func_from_pkl('gamma_H_z_pn')
    gammaHz_np=func_from_pkl('gamma_H_z_np')
    func=gammaHz_pn(z_gamma)*(1-X_n)-(gammaHz_np(z_gamma))*X_n#(m_Pl/(1.66*np.sqrt(g)*m**2))*z_gamma*(gamma_pn(z_gamma)*(1-X_n)-gamma_np(z_gamma)*X_n)
    return func    
            
def X_n_eq(z):
    return 1/(1+np.exp(1.293/0.511*z))


"""

-------------------------------------------------------------------------------------------------------------------------------------
Functions to calculate the abundance from Xn, namely the temperature of deuterium production 
-------------------------------------------------------------------------------------------------------------------------------------

"""
def BBN_func(T_inv,eta):
    B=2.23
    m_n=939
    c=(12*zeta(3)*eta)/(m_n**1.5*np.sqrt(np.pi))#Is it 24 or 12? Because it does indeed change the numbers
    f=c*T_inv**(-1.5)*np.exp(B*T_inv)-1
    return f
    
def BBN_func2(zc,eta,Xn):
    z=0.511/2.23*zc
    f=(1.18e13*eta**2/np.sqrt(g_eff(2.23/zc))*(zc)**(-17/6)*np.exp(-1.44*(zc)**(1/3))*np.exp(zc)*Xn(z)*(1-Xn(z)))-1
    return f


def T_BBN(eta,Xn,a=2.23,b=223):
    vec_res = np.array([brentq(BBN_func2, a , b , args=(eta_i,Xn,)) for eta_i in eta])
    return 2.23/vec_res


def T_BBN_2(eta,XD=1e-2):
    return 0.064/(1-0.029*np.log((eta*1e10)))
 
def Y_He(eta,solved_Boltz,Neff=3.045,verbose=False,td=2):
    Xn=solved_Boltz.sol
    if td==1:
        T_d=T_BBN(eta,Xn)
    else:
        T_d=T_BBN_2(eta,Xn)
    m=0.511
    Xn=solved_Boltz.sol(m/T_d)
    if verbose: 
        print("The Helium abundance is Yp = {:.5g}. It shoud be: {:g}\\n The effective number of neutrinos is:{:.5g}".format(2*Xn[0],Y_He_control(eta),3))#ge.N_eff_SMC(float(T_d))))
    return 2*Xn,T_d

def Y_He_control(eta,tau_n=879,delta_n=0,xi_nu=0):
    # Y=0.2462+0.01*np.log(eta/5e-10)*(eta/5e-10)**-0.2
    Y=0.24462+0.012*delta_n+0.00021*(tau_n-887)+0.009*np.log(eta/5e-10)-0.25*xi_nu
    return Y


