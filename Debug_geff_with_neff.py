#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  6 18:54:13 2020

@author: victormaura12
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.integrate import quad
from scipy.interpolate import interp1d
import time
import pandas as pd

def func_from_pkl(filename):
    with open(filename+'_func.pickle', 'rb') as handle:
        cases=pickle.load(handle)
    return cases

def func_to_pkl(dictionary,filename):
     with open(filename+'_func.pickle', 'wb') as handle:
        pickle.dump(dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)
        


def geff_part_e(T,m,gi,alpha=1):
    const=gi*15/np.pi**4
    z=m/T
    func=lambda u: (u**2*np.sqrt(u**2-z**2))/(np.exp(u)+alpha)
    integral=quad(func,z,np.inf)
    return const*integral[0]

def geff_part_s(T,m,gi,alpha=1):
    const=gi*45/(4*np.pi**4)
    z=m/T
    func=lambda u: (u*np.sqrt(u**2-z**2))*(4*u**2-z**2)/(3*u*(np.exp(u)+alpha))
    integral=quad(func,z,np.inf)
    return const*integral[0]

def g_eff_s(T):
    g_e=2+geff_part_s(T,0.511,4)
    return g_e

def g_eff_e(T):
    vecg=np.vectorize(geff_part_e)
    g_e=2+vecg(T,0.511,4)
    return g_e

def g_eff_e_spl(filename,T_range=(-5,2)):
    T=np.logspace(T_range[0]-1,T_range[1]+1,1000)
    func=g_eff_e(T)
    start_time=time.time()
    func2=interp1d(T,func,kind='cubic',assume_sorted=True)
    print("--- %s seconds ---" % (time.time() - start_time))
    func_to_pkl(func2,filename)

T_over_nu=func_from_pkl('T_over_T_nu')

def N_eff_SMC(T,N_nu=3,xi_nu=0):
    n=N_nu*(11/4)**(4/3)*(1/T_over_nu(T))**4*(1+30/(7*np.pi**2)*xi_nu**2+15/(7*np.pi**4)*xi_nu**4)
    return n

    
    

def T_over_T_nu(T_range,spline=False,step=1000,filename='T_over_T_nu'):
    if spline:
            T=np.logspace(T_range[0]-1,T_range[1]+1,step)
            Tnu=[np.cbrt(5.5/g_eff_s(Ti)) for Ti in T]
            func=interp1d(T,Tnu,kind='cubic')
            func_to_pkl(func,filename)
    try:
        func2=func_from_pkl(filename)
    
    except FileNotFoundError: 
        T=np.logspace(T_range[0]-1,T_range[1]+1,step)
        Tnu=[np.cbrt(g_eff_s(Ti)/5.5) for Ti in T]
        func=interp1d(T,Tnu,kind='cubic')
        func_to_pkl(func,filename)
        func2=func_from_pkl(filename)
    
    return func2
        
    


def g_eff_tot(T,N_eff):
    filename='g_eff_e'
    geffspl=func_from_pkl(filename)
    return geffspl(T)+2*7/8*N_eff*(4/11)**(4/3)



T_range=(-3,1)
T=np.logspace(T_range[0]-1,T_range[1]+1,1000)
# # # # filename='g_eff_e'
# # # # g_eff_e_spl(filename,T_range=T_range)
# # # # geffspl=func_from_pkl(filename)

# # # # plt.plot(T,g_eff_e(T),T,geffspl(T))
# # # # plt.semilogx()

# # # # df=pd.read_csv('Neff_Miguel_SM.csv',names=['x','y'], delimiter=' ', decimal=',')
# # # # df=df.set_index('x')
# # # # df.plot()
    

df=pd.read_csv('SMgstar_2019_Laine.dat',skiprows=tuple(range(6)),delimiter='    ')
df.set_index('T[MeV]',inplace=True)
df=df[['g_eff','h_eff','i_eff']]
# df.plot(y='g_eff')
df2=df[df.index<10]
T=df2.index.to_numpy()
gstar=interp1d(T,df2['g_eff'].to_numpy())

# # # # m=0.511
# # # # g_i=4
# # # # T_range=(1e-4,1000)
# # # # step=1000
# # # log_T=(np.log10(T_range[0]),np.log10(T_range[1]))
# # # # T=np.logspace(log_T[0],log_T[1],step)
g_eff= lambda T: g_eff_tot(T,N_eff_SMC(T,N_nu=3,xi_nu=0))
plt.plot(T,100*abs(g_eff(T)-gstar(T))/g_eff(T))
# # # plt.plot(T,T_over_nu(T))
# # # Neff=[N_eff_SMC(Ti) for Ti in T]
# # # plt.plot(T,Neff
# # plt.plot(T,1/T_over_nu(T))
# # plt.plot(T,T_over_T_nu((T_range[0],T_range[1]),spline=False,step=1000)(T))
# # # # g=[(g_eff_s(Ti)/5.5)**(-1/3) for Ti in T]
# # # # # plt.plot(T,g)
# # # # # # # # # geff=[geff_part_e(Ti,m,4) for Ti in T]
# # # # # # # # # geff_s=[geff_part_s(Ti,m,4) for Ti in T]

# # # # # # # # # plt.plot(T,geff,label='ENERGY')
# # # # # # # # # plt.plot(T,geff_s,label='ENTROPY')
# # plt.xlim(T[-1],T[0])
# # plt.semilogx()
# # plt.legend()
# # plt.tight_layout()
# # plt.grid(True)