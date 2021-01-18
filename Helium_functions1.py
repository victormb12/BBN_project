#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 01:36:25 2020

@author: victormaura12
"""

import numpy as np
from scipy.integrate import quad_vec,solve_ivp
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import InterpolatedUnivariateSpline as interp1d
from matplotlib import rcParams
import time
from scipy.optimize import brentq 
import geff_with_neff_1 as ge
from cycler import cycler
import itertools
import pandas as pd

rcParams['lines.linewidth'] = 1.5
rcParams['font.size']=11
rcParams['axes.prop_cycle']= cycler('color',list('grkby'))
plt.style.use('seaborn-ticks')

def dict_from_pkl(filename):
    with open(filename+'_dict.pickle', 'rb') as handle:
        cases=pickle.load(handle)
    return cases

def dict_to_pkl(dictionary,filename):
     with open(filename+'_dict.pickle', 'wb') as handle:
        pickle.dump(dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)

def func_from_pkl(filename):
    with open(filename+'_func.pickle', 'rb') as handle:
        cases=pickle.load(handle)
    return cases

def func_to_pkl(function,filename):
     with open(filename+'_func.pickle', 'wb') as handle:
        pickle.dump(function, handle, protocol=pickle.HIGHEST_PROTOCOL)
        

"""

-------------------------------------------------------------------------------------------------------------------------------------
Calculation of the weak interaction rates accroding to https://arxiv.org/pdf/1812.05605.pdf 
-------------------------------------------------------------------------------------------------------------------------------------

"""

def Mev2_to_cm3s(a, b=0):
    
    conv=((1.9733e-11)**3)/(6.5821e-28)
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
    # print(T_over_nu(T))
    K=MeV_to_1_over_s(1/(1.636*880.2),b=0)
    q=1.2933/0.511
    func1=((epsilon-sign*q)**2*np.sqrt((epsilon**2-1))*epsilon)/((1+np.exp(-epsilon*z_gamma))*(1+np.exp((epsilon-sign*q)*T_over_nu(T)*z_gamma-sign*xi_nu)))
    func2=((epsilon+sign*q)**2*np.sqrt((epsilon**2-1))*epsilon)/((1+np.exp(epsilon*z_gamma))*(1+np.exp(-(epsilon+sign*q)*T_over_nu(T)*z_gamma-sign*xi_nu)))
    return K*(func1+func2)


def gamma_p_to_n(z_gamma,xi_nu):
    integral=quad_vec(lambda epsilon: n_to_p_int(epsilon,z_gamma,sign=-1,xi_nu=xi_nu),1,np.inf)
    return integral[0]
    
def gamma_n_to_p(z_gamma,xi_nu):
    integral=quad_vec(lambda epsilon: n_to_p_int(epsilon,z_gamma,sign=1,xi_nu=xi_nu),1,np.inf)
    return integral[0]

def gamma_pn_spline(filename,z_range=(-4,3),xi_nu=0,return_spline=False):
    z=np.logspace(z_range[0],z_range[1],1000)
    func=gamma_p_to_n(z)
    f2 = interp1d(z, func,k=3)
    if return_spline:
        return f2
    func_to_pkl(f2,filename)


def gamma_np_spline(filename,z_range=(-4,3),xi_nu=0,return_spline=False):
    z=np.logspace(z_range[0],z_range[1],1000)
    func=gamma_n_to_p(z)
    f2 = interp1d(z, func, k=3)
    if return_spline:
        return f2
    func_to_pkl(f2,filename)




    
def rate_spline(z_range,Neff=3,xi_nu=0):
    filename_Hz_np='gamma_H_z_np'
    filename_Hz_pn='gamma_H_z_pn'
    gamma_H_z_np(filename_Hz_np,z_range=z_range,Neff=Neff,xi_nu=xi_nu)
    gamma_H_z_pn(filename_Hz_pn,z_range=z_range,Neff=Neff,xi_nu=xi_nu)

def gamma_H_z_pn(filename,z_range,dN_eff=0,xi_nu=0,return_spline=False):
    zi=np.logspace(np.log10(z_range[0])-1,np.log10(z_range[1])+1,1000)
    m=0.511
    T=m/zi
    """
    In order to use the time temperature relation by Laine et al just uncomment these two lines
    """
    # dT_dt=ge.dT_dt(T)
    # func=-m/zi**2*1/dT_dt(T)*gamma_p_to_n(zi,xi_nu)
    g_eff= lambda t: ge.g_eff_e_star(t,dN_eff=dN_eff,xi_nu=xi_nu)
    geffs=ge.g_eff_s_tot(T,return_spline=True)
    dgeffs= lambda T: geffs.derivative(1)(T)
    func=gamma_p_to_n(zi,xi_nu)/(1.66*np.sqrt(g_eff(0.511/zi))*(0.511)**2/1.221e22)*zi*(1+1/3*T/geffs(T)*dgeffs(T))
    f2 = interp1d(zi, func, k=3)
    if return_spline:
        return f2
    func_to_pkl(f2,filename)
    
def gamma_H_z_np(filename,z_range,dN_eff=0,xi_nu=0,return_spline=False):
    zi=np.logspace(np.log10(z_range[0])-1,np.log10(z_range[1])+1,1000)
    m=0.511
    T=m/zi
    """
    In order to use the time temperature relation by Laine et al just uncomment these two lines
    """
    # dT_dt=ge.dT_dt(T)
    # func=-m/zi**2*1/dT_dt(T)*gamma_n_to_p(zi,xi_nu)
    g_eff= lambda t: ge.g_eff_e_star(t,dN_eff=dN_eff,xi_nu=xi_nu)
    geffs=ge.g_eff_s_tot(T,return_spline=True)
    dgeffs= lambda T: geffs.derivative(1)(T)
    func=gamma_p_to_n(zi,xi_nu)/(1.66*np.sqrt(g_eff(0.511/zi))*(0.511)**2/1.221e22)*zi*(1+1/3*T/geffs(T)*dgeffs(T))
    func=gamma_n_to_p(zi,xi_nu)/(1.66*np.sqrt(g_eff(0.511/zi))*(0.511)**2/1.221e22)*zi*(1+1/3*T/geffs(T)*dgeffs(T))
    
    f2 = interp1d(zi, func, k=3)
    if return_spline:
        return f2
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



def Boltzmann_for_n(z_gamma,X_n,pn_rate=func_from_pkl('gamma_H_z_pn'),np_rate=func_from_pkl('gamma_H_z_np')):
    func=pn_rate(z_gamma)*(1-X_n)-(np_rate(z_gamma))*X_n#(m_Pl/(1.66*np.sqrt(g)*m**2))*z_gamma*(gamma_pn(z_gamma)*(1-X_n)-gamma_np(z_gamma)*X_n)
    return func    
            
def X_n_eq(z,xi_nu=0):
    return 1/(1+np.exp(1.293/0.511*z+xi_nu))


"""

-------------------------------------------------------------------------------------------------------------------------------------
Functions to calculate the abundance from Xn, namely the temperature of deuterium production 
-------------------------------------------------------------------------------------------------------------------------------------

"""


def BBN_func2(zc,eta,Xn,g_eff):
    z=0.511/2.23*zc
    T=0.511/z
    geffs= ge.g_eff_s_tot(T, return_spline=True)
    dgeffs= lambda T: geffs.derivative(1)(T)
    f=(1.18e13*eta**2/np.sqrt(g_eff(2.23/zc))*(zc)**(-17/6)*np.exp(-1.44*(zc)**(1/3))*np.exp(zc)*Xn(z)*(1-Xn(z))*(1+1/3*T/geffs(T)*dgeffs(T)))-1
    return f

def T_BBN(eta,Xn,N_eff=3,xi_nu=0):
    a=2.23
    b=223
    g_eff=lambda T: ge.g_eff_tot(T,ge.N_eff_SMC(T,N_nu=N_eff,xi_nu=xi_nu))
    vec_res = np.array([brentq(BBN_func2, a , b , args=(eta_i,Xn,g_eff)) for eta_i in eta])
    return 2.23/vec_res

def T_BBN_2(eta,XD=1e-2):
    return 0.064/(1-0.029*np.log((eta*1e10)))

def Y_He(eta,solved_Boltz,dN_eff=0,xi_nu=0,verbose=False,Td=1):
    Xn=solved_Boltz.sol
    if Td==1:
        T_d=T_BBN(eta,Xn,dN_eff,xi_nu)
    elif Td==2:
        T_d=T_BBN_2(eta,Xn)
    else:
        T_d=Td#taken from BBN Miguel
    m=0.511
    Xn=solved_Boltz.sol(m/T_d)
    if verbose: 
        try:
            print("The Helium abundance is Yp = {:.5g}. It shoud be: {:g}\\n The effective number of neutrinos is:{:.5g}".format(2*Xn[0],Y_He_control(eta),3))#ge.N_eff_SMC(float(T_d))))
        except Exception as e: 
            print(e)
    return 2*Xn

def Y_He_control(eta,tau_n=879,dN_eff=0,xi_nu=0,steigman=False,Td=1):
    if steigman:
        df=pd.read_csv('n_over_p_Steigman.csv',names=['x','y'], delimiter=' ', decimal=',')
        df['x']=df['x'].div(1000)               
        df['x']=df['x'].apply(lambda row: 0.511/row)
        df=df.set_index('x')
        np=interp1d(df.index.to_numpy(),df['y'].to_numpy(),k=3)
        Xn=lambda z: 1/(1+1/np(z))
        if Td==1:
            T_d=T_BBN(eta,Xn,dN_eff,xi_nu)
        elif Td==2:
            T_d=T_BBN_2(eta,Xn)
        else:
            T_d=Td
        Xn=1/(1+1/np(0.511/T_d))
        return 2*Xn
    
    df=pd.read_csv('Y_He_Pdg.csv',names=['eta','Y'], delimiter=' ', decimal=',')
    df=df.set_index('eta')
    Y_spl=interp1d(df.index.to_numpy(),df['Y'].to_numpy())
    return Y_spl(eta)


def solve_Boltz(T_range=(2,1e-2),dN_eff=0,xi_nu=0,create_spline=False,H_z_pn_rate=None,H_z_np_rate=None ,rtol=1e-6,atol=1e-20,y0=None):
    g_eff= lambda T: ge.g_eff_e_star(T,dN_eff=dN_eff,xi_nu=xi_nu)
    start_time=time.time()
    m=0.511
    z_range=(m/T_range[0],m/T_range[1])
    rtol=rtol
    atol=atol
    
    if y0 is None:
        y0= X_n_eq(z_range[0],xi_nu=xi_nu)
        if type(y0)!=type(np.array([0])):
            y0=np.array([y0])
    if create_spline:
        filename_np='gamma_np'
        filename_pn='gamma_pn'
        filename_Hz_np='gamma_H_z_np'
        filename_Hz_pn='gamma_H_z_pn'
        z_range_2=(np.log10(z_range[0])-1,np.log10(z_range[1])+1)
        gamma_pn_spline(filename_pn,z_range=z_range_2)
        gamma_np_spline(filename_np,z_range=z_range_2)
        gamma_H_z_np(filename_Hz_np,z_range_2,g_eff)
        gamma_H_z_pn(filename_Hz_pn,z_range_2,g_eff)
    if not (H_z_pn_rate is None or H_z_np_rate is None):
        Boltzmann2=ODE(Boltzmann_for_n,y0,z_range,rtol=rtol,atol=atol,params=(H_z_pn_rate,H_z_np_rate),dense_output=True,verbose = True, vectorized=True)
    else:
        Boltzmann2=ODE(Boltzmann_for_n,y0,z_range,rtol=rtol,atol=atol,dense_output=True,verbose = True, vectorized=True)
    print("--- %s seconds ---" % (time.time() - start_time))
    return Boltzmann2

def plot_Boltz(T_range=(2,1e-2),dN_eff=0,xi_nu=0,vline=False,step=1000,n_over_p=False):
    m=0.511
    return_spline=True
    z_range=(m/T_range[0],m/T_range[1])
    z=np.logspace(np.log10(z_range[0]),np.log10(z_range[1]),step)
    np_rate=gamma_H_z_np('',z_range,dN_eff=dN_eff,xi_nu=xi_nu, return_spline=return_spline)
    pn_rate=gamma_H_z_pn('',z_range,dN_eff=dN_eff,xi_nu=xi_nu, return_spline=return_spline)
    Boltzmann2=solve_Boltz(T_range,dN_eff=dN_eff,xi_nu=xi_nu,H_z_pn_rate=pn_rate,H_z_np_rate=np_rate)
    if n_over_p:
        fig,ax=plt.subplots(2,1,sharex=True)
        secax = ax[0].secondary_xaxis('top', functions=(lambda z:0.511e3/z, lambda z: 0.511e3/z))
        df2=pd.read_csv('Free_Neutron_decay.csv',names=['x','y'], delimiter=' ', decimal=',')
        df2['x']=df2['x'].div(1000)               
        df2['x']=df2['x'].apply(lambda row: 0.511/row)
        df2=df2.set_index('x')
        df2.plot(ax=ax[0])
        
        df=pd.read_csv('n_over_p_Steigman.csv',names=['x','y'], delimiter=' ', decimal=',')
        df['x']=df['x'].div(1000)               
        df['x']=df['x'].apply(lambda row: 0.511/row)
        df=df.set_index('x')
        df.plot(ax=ax[0])
        Xn=Boltzmann2.func_sol.sol(z).T
        n_over_p=X_n_eq(z,xi_nu)/(1-X_n_eq(z,xi_nu))
        
        ax[0].plot(z,n_over_p,label=r'$X_n^{eq}$')  
        ax[0].plot(z,Xn/(1-Xn))
        
        
        z_range=(df.index[0],df.index[-1])
        z=np.logspace(np.log10(z_range[0]),np.log10(z_range[1]),step)
        Xn=Boltzmann2.func_sol.sol(z)[0]
        n_over_p=Xn/(1-Xn)
        n_over_p2=interp1d(df.index.to_numpy(),df['y'].to_numpy(),k=3)
        ax[0].set_ylim((0,1))

        dev=100*abs((n_over_p-n_over_p2(z))/n_over_p2(z))
        ax[1].plot(z,dev)
        ax[1].set_ylim((0,100))
        ax[0].semilogx()
        ax[0].grid(True)
        ax[1].grid(True)
        ax[0].set_ylabel(r'$X_n$')
        ax[0].set_xlabel(r'$z=\frac{m_e}{T}$')
        ax[0].legend()
        ax[1].set_ylabel('deviation')
        ax[1].set_xlabel(r'$z=\frac{m_e}{T}$')
        
    
    else:
        fig,ax=plt.subplots(1,1)
        secax = ax.secondary_xaxis('top', functions=(lambda z:0.511e3/z, lambda z: 0.511e3/z))
        ax.plot(z,X_n_eq(z,xi_nu),label=r'$X_n^{eq}$')  
        Boltzmann2.plot(label=r'$\Delta Neff={:.3g},\xi_\nu={:.3g}$'.format(dN_eff[0],xi_nu[0]))
        plt.ylim((0,0.5))
        ax.set_ylabel(r'$X_n$')
        ax.set_xlabel(r'$z=\frac{m_e}{T}$')
        ax.grid(True,which='both')
        ax.semilogx()
        ax.legend()
    
    if vline:
        plt.plot(z,[0.15]*len(z))
        ax.vlines(0.511/0.066,0,1,linestyle='--')
    plt.tight_layout()
        
def plot_He(eta,Neff=3,xi_nu=0,deviation=False,compare=True,T_d=1):
    if type(Neff)!=type(np.array([1])):
        Neff=np.array([Neff])
    if type(xi_nu)!=type(np.array([1])):
        xi_nu=np.array([xi_nu])
    Y_dict=Y_He_fin(eta,Neff,xi_nu,T_d=T_d)
    # if deviation:
    #     fig,ax=plt.subplots(2,1)
    #     ax[1].plot(eta,abs((Y_dict[r'$Y_{He}']-Y_He(eta,Boltzmann2.func_sol,verbose=False,td=1)[0][0])/Y_He_control(eta)), label= 'relative deviation')
    #     ax[0].plot(eta,Y_He(eta,Boltzmann2.func_sol,verbose=False,td=1)[0], label= 'Naive')    
    #     ax[0].plot(eta,Y_He_control(eta),label='Bernstein et al.')
    #     ax[0].plot(eta,[0.245]*len((eta)))
    # else:
    fig,ax=plt.subplots(1,1)
    Y_dict[(Neff[0],xi_nu[0])].plot(ax=ax,label=r'$Neff={:.3g}, \xi_\nu={:.3g}$'.format(Neff[0],xi_nu[0]))    
    ax.set_xlabel(r'$\eta$')
    ax.set_ylabel(r'$Y_{He}$')
    if compare:
        ax.plot(eta,Y_He_control(eta,xi_nu=xi_nu,dN_eff=Neff),label='Bernstein et al.')
    ax.plot(eta,[0.245]*len((eta)))
    plt.semilogx()
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

def plot_Td(eta,xi_nu=0,Neff=3,compare=False):
    return_spline=True
    m=0.511
    T_range=(2,1e-2)#in MeV
    z_range=(m/T_range[0],m/T_range[1])
    np_rate=gamma_H_z_np('',z_range,Neff=Neff,xi_nu=xi_nu, return_spline=return_spline)
    pn_rate=gamma_H_z_pn('',z_range,Neff=Neff,xi_nu=xi_nu, return_spline=return_spline)
    Xn=solve_Boltz(T_range,Neff=Neff,xi_nu=xi_nu,H_z_pn_rate=pn_rate,H_z_np_rate=np_rate)
    fig,ax=plt.subplots(1,1)
    ax.plot(eta,T_BBN(eta,Xn.func_sol.sol),label='My calculation')
    if compare:
        plt.plot(eta,T_BBN_2(eta,Xn), label='Dolgov')     
    ax.semilogx()
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    
def Y_He_fin(etas=np.array([3.251e-10]),dN_effs=np.array([0]),xi_nus=np.array([0]),verbose=False,T_d=1):
    start_time=time.time()
   
    if type(etas)!=type(np.array([1])):
            etas=np.array([etas])
    if type(xi_nus)!=type(np.array([1])):
            xi_nus=np.array([xi_nus])
    if type(dN_effs)!=type(np.array([1])):
            dN_effs=np.array([dN_effs])
    Y_dict={}
    return_spline=True
    m=0.511
    T_range=(2,1e-2)#in MeV
    z_range=(m/T_range[0],m/T_range[1])
    for dN_eff,xi_nu in itertools.product(dN_effs,xi_nus):
        print(dN_eff,xi_nu)
        np_rate=gamma_H_z_np('',z_range,dN_eff=dN_eff,xi_nu=xi_nu, return_spline=return_spline)
        pn_rate=gamma_H_z_pn('',z_range,dN_eff=dN_eff,xi_nu=xi_nu, return_spline=return_spline)
        Xn=solve_Boltz(T_range,dN_eff=dN_eff,xi_nu=xi_nu,H_z_pn_rate=pn_rate,H_z_np_rate=np_rate)
        # boltzmann_dict[(1,xi_nu)]=Xn
        Y_dict[(dN_eff,xi_nu)]=pd.DataFrame(Y_He(etas,Xn.func_sol,dN_eff=dN_eff,xi_nu=xi_nu,Td=T_d).T,index=etas, columns=[r'$Y_{He}$'])
    print("--- %s seconds ---" % (time.time() - start_time))
    return Y_dict
    

def N_v_xi(eta,Y,Neffs,xi_nus,step=100,epsilon=0.01,create_dict=False,plot=True):
    filename='N_vs_xi'
    if create_dict:    
        if type(eta)!=type(np.array([1])):
            eta=np.array([eta])
        Y_dict=Y_He_fin(eta,Neffs,xi_nus)
        dict_to_pkl(Y_dict,filename)
    else:
        Y_dict=dict_from_pkl(filename)
    plot_xi=[]
    plot_Neff=[]
    for Neff,xi_nu in Y_dict.keys():
        if Y-epsilon <= Y_dict[(Neff,xi_nu)].iloc[0][0]<=Y+epsilon:
            plot_xi.append(xi_nu)
            plot_Neff.append(Neff)
    if plot:
        fig,ax=plt.subplots(1,1)
        ax.set_xlabel(r'$\xi_\nu$')
        ax.set_ylabel(r'$N_{eff}$')
        ax.plot(plot_xi,plot_Neff,label=r'$\eta={0:.3g},Y={1:.3g}$'.format(eta[0],Y))
        ax.legend()
        ax.grid(True)
        plt.tight_layout()
    return plot_Neff,plot_xi
        

def Y_v_N(eta,xi_nu,Neffs,create_dict=False,plot=True):
    filename='Y_vs_N'
    if create_dict:    
        if type(eta)!=type(np.array([1])):
            eta=np.array([eta])
        if type(xi_nu)!=type(np.array([1])):
            xi_nu=np.array([xi_nu])
        Y_dict=Y_He_fin(eta,Neffs,xi_nu)
        dict_to_pkl(Y_dict,filename)
    else:
        Y_dict=dict_from_pkl(filename)
    plot_Y=[]
    plot_Neff=[]
    for Neff,xi_nu in Y_dict.keys():
        plot_Y.append(Y_dict[(Neff,xi_nu)].iloc[0][0])
        plot_Neff.append(Neff)
    if plot:
        fig,ax=plt.subplots(1,1)
        ax.set_xlabel(r'$N_{eff}$')
        ax.set_ylabel(r'$Y_{He}(N_{eff})$')
        ax.plot(plot_Neff,plot_Y,label=r'$\eta={0:.3g},\xi_\nu={1:.3g}$'.format(eta[0],xi_nu))
        ax.legend()
        ax.grid(True)
        plt.tight_layout()
    return plot_Neff,plot_Y
    
def Y_v_xi(eta,Neff,xi_nu,create_dict=False,plot=True):
    filename='Y_vs_xi'
    if create_dict:    
        if type(eta)!=type(np.array([1])):
            eta=np.array([eta])
        if type(Neff)!=type(np.array([1])):
            Neff=np.array([Neff])
        Y_dict=Y_He_fin(eta,Neff,xi_nu)
        dict_to_pkl(Y_dict,filename)
    else:
        Y_dict=dict_from_pkl(filename)
    plot_xi=[]
    plot_Y=[]
    for Neff,xi_nu in Y_dict.keys():
        plot_Y.append(Y_dict[(Neff,xi_nu)].iloc[0][0])
        plot_xi.append(xi_nu)
    if plot:
        fig,ax=plt.subplots(1,1)
        ax.set_xlabel(r'$\xi_{\nu}$')
        ax.set_ylabel(r'$Y_{He}(\xi_\nu)$')
        ax.plot(plot_xi,plot_Y,label=r'$\eta={0:.3g},Neff={1:.3g}$'.format(eta[0],Neff))
        ax.legend()
        ax.grid(True)
        plt.tight_layout()
    return plot_xi,plot_Y

def eta_v_xi(Y,Neff,xi_nus,etas,epsilon=0.001,create_dict=False,plot=True):
    filename='eta_vs_xi'
    if create_dict:    
        if type(Neff)!=type(np.array([1])):
            Neff=np.array([Neff])
        Y_dict=Y_He_fin(etas,Neff,xi_nus)
        dict_to_pkl(Y_dict,filename)
    else:
        Y_dict=dict_from_pkl(filename)
    plot_eta=[]
    plot_xi=[]
    for Neff,xi_nu in Y_dict.keys():
            df=Y_dict[(Neff,xi_nu)]
            eta_mean= np.mean(list(df[(df[r'$Y_{He}$'] >= Y-epsilon) & (df[r'$Y_{He}$'] <= Y+epsilon)].index))
            plot_xi.append(xi_nu)
            plot_eta.append(eta_mean)
    if plot:
        fig,ax=plt.subplots(1,1)
        ax.set_xlabel(r'$\eta$')
        ax.set_ylabel(r'$\xi_{\nu}$')
        ax.plot(plot_eta,plot_xi,label=r'$Y={0:.3g},Neff={1:.3g}$'.format(Y,Neff))
        ax.legend()
        ax.grid(True)
        ax.semilogx()
        plt.tight_layout()
    return plot_eta,plot_xi



# df=pd.read_csv('n_over_p_Steigman.csv',names=['x','y'], delimiter=' ', decimal=',')
# df=df.set_index('x')
# df.plot()





# create_dict=False

# etas=np.array([5e-10])#np.logspace(-10,-9,100)
# xi_nus=np.linspace(0,0.1,10)
# Neffs=np.linspace(2,4,10)

# plot_N,plot_xi=N_v_xi(etas,0.245,Neffs=Neffs,xi_nus=xi_nus,epsilon=0.001,create_dict=create_dict)

# # plt.plot(plot_xi,plot_N)

# # xi_nus=[0.2]
# # for xi_nu in xi_nus:
# #     plot_N,plot_Y,Y_dict=Y_v_N(etas,xi_nu,N_range=(Neffs[0],Neffs[-1]),step=10,create_dict=create_dict)
# #     plt.plot(plot_N,plot_Y)

# create_dict=False
# Neff=3.045
# xi_nus=np.linspace(0,1,100)
# etas=np.logspace(-12,-7,100)
# plot_xi,plot_Y=Y_v_xi(etas, Neff, xi_nus,create_dict=create_dict)



# plot_xi,plot_eta=eta_v_xi(0.245,Neff,xi_nus,etas,epsilon=0.001,create_dict=create_dict)
# # plt.plot(plot_eta,plot_xi)
# # plt.semilogx()


# # plt.grid(True)




