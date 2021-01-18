#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 13:50:34 2021

@author: victormaura12
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import geff_with_neff as ge
from scipy.interpolate import InterpolatedUnivariateSpline as interp1d
from scipy.optimize import brentq
import Helium_functions as He
pd.set_option('display.float_format', lambda x: '%.5g' % x)
sigmas=pd.read_csv('deuterium_reaction_rates.txt',delimiter=' ')
print(sigmas.head())
def GK_to_Mev(T,invert=False):
    const=8.617e-2
    if invert:
        return 1/const*T
    return const*T

def Mev2_to_cm3s(a, b=0):
    
    conv=((1.9733e-11)**3)/(6.5821e-22)
    if b==0:
        return a*conv
    elif b==1:
        return a/conv

sigmas['T_in_MeV']=sigmas['T_in_GK'].apply(GK_to_Mev)
sigmas=sigmas[['T_in_MeV','3He','3H']]
Na=6.022e23
print('\n\n\n',sigmas.head())

func=lambda sigma: 1/Na*Mev2_to_cm3s(sigma,b=1)
sigmas[['3He','3H']]=sigmas[['3He','3H']].apply(func)
print('\n\n\n',sigmas.head())
func2=lambda T: 2.23/T
sigmas['zc']=sigmas['T_in_MeV'].apply(func2)
sigmas['3He+3H']=sigmas['3H']+sigmas['3He']
interp_df=sigmas[['zc','3He+3H']]
interp_df=interp_df.set_index('zc')
interp_df=interp_df.sort_index()
interp=interp1d(interp_df.index.to_numpy(),interp_df['3He+3H'].to_numpy(),k=3)


def BBN_func2(zc,eta,Xn,g_eff):
    B=2.23
    mp=938.3
    z=0.511/2.23*zc
    T=0.511/z
    geffs= ge.g_eff_s_tot(T, return_spline=True)
    dgeffs= lambda T: geffs.derivative(1)(T)
    const=3.76e18*eta**2
    dzdt=1/(np.sqrt(g_eff(T)))*(1+1/3*T/geffs(T)*dgeffs(T))
    func=const*dzdt*zc**-3.5*np.exp(zc)*Xn(z)*(1-Xn(z))*interp(zc)*0.03
    return func-1

def T_BBN(eta,Xn,dN_eff=0,xi_nu=0):
    a=2.23
    b=223
    g_eff=lambda T: ge.g_eff_e_star(T,dN_eff=dN_eff,xi_nu=xi_nu)
    vec_res = np.array([brentq(BBN_func2, a , b , args=(eta_i,Xn,g_eff)) for eta_i in eta])
    return 2.23/vec_res
eta=np.array([6.1e-10])
Xn=He.solve_Boltz()
print(T_BBN(eta,Xn.func_sol.sol))



# def rate_func(zc,eta,dN_eff=0,xi_nu=0):
#     B=2.23
#     T=B/zc
#     g_eff= lambda t: ge.g_eff_e_star(t,dN_eff=dN_eff,xi_nu=xi_nu)
#     geffs=ge.g_eff_s_tot(T,return_spline=True)
#     dgeffs= lambda T: geffs.derivative(1)(T)
#     mult=1/(1.66*np.sqrt(g_eff(T))*(B)**2/1.221e22)*(1+1/3*T/geffs(T)*dgeffs(T))*eta*2*1.202/np.pi**2*B**3*zc**-2
#     return mult
# eta=6.1e-10
# sigmas['mult']=sigmas['zc'].apply(lambda zc: rate_func(zc,eta))
# sigmas['3He_rate']=sigmas['3He']*sigmas['mult']
# sigmas['3H_rate']=sigmas['3H']*sigmas['mult']
# sigmas['Rdd']=sigmas['3He_rate']+sigmas['3H_rate']
# print('\n\n\n',sigmas.head())
# sigmas.plot(x='zc',y='Rdd')
# zc=sigmas['zc'].to_numpy()
# plt.plot(zc,2*0.3/0.87*2.4e7*zc**(-4/3)*np.exp(-1.44*zc**(1/3)))
# plt.loglog()
# # rates=sigmas.apply()