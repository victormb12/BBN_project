#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 00:44:10 2020

@author: victormaura12
"""

import numpy as np
import matplotlib.pyplot as plt

# def p_to_n_int_1(epsilon,z_gamma,sign=1,xi_nu=0):
#     """
#     f_anti_nu*(1-f_anti_e), E_nu=E_e+dm
#     integrated from 1 to infinity

#     Parameters
#     ----------
#     epsilon : TYPE
#         DESCRIPTION.
#     z_gamma : TYPE
#         DESCRIPTION.
#     sign : TYPE, optional
#         DESCRIPTION. The default is 1.
#     xi_nu : TYPE, optional
#         DESCRIPTION. The default is 0.

#     Returns
#     -------
#     TYPE
#         DESCRIPTION.

#     """
#     T_over_nu=func_from_pkl('T_over_T_nu')
#     T=0.511/z_gamma
#     K=MeV_to_1_over_s(1/(1.636*880.2),b=0)
#     q=1.2933/0.511
#     func1=((epsilon+q)**2*np.sqrt((epsilon**2-1))*epsilon)/((1+np.exp(-epsilon*z_gamma))*(1+np.exp((epsilon+q)*z_gamma*T_over_nu(T)+xi_nu)))
#     return K*func1

# def p_to_n_int_2(epsilon,z_gamma,sign=1,xi_nu=0):
#     """
#     f_e*(1-f_nu), E_e=E_nu+dm
#     integrated from q to infinity
#     Parameters
#     ----------
#     epsilon : Electron energy over electron mass, parameter to be integrated.
#     z_gamma : TYPE
#         DESCRIPTION.
#     sign : TYPE, optional
#         DESCRIPTION. The default is 1.
#     xi_nu : TYPE, optional
#         Neutrino Chemical potential. The default is 0.

#     Returns
#     -------
#     TYPE
#         DESCRIPTION.

#     """
#     T_over_nu=func_from_pkl('T_over_T_nu')
#     T=0.511/z_gamma
#     K=MeV_to_1_over_s(1/(1.636*880.2),b=0)
#     q=1.2933/0.511
#     func2=((epsilon-q)**2*np.sqrt((epsilon**2-1))*epsilon)/((1+np.exp(epsilon*z_gamma))*(1+np.exp(-(epsilon-q)*z_gamma*T_over_nu(T)+xi_nu)))
#     return K*func2

# def p_to_n_int_3(epsilon,z_gamma,sign=1,xi_nu=0):
#     """
#     f_anti_nu*f_e, E_e+E_nu=dm
#     integrated from 1 to q

#     Parameters
#     ----------
#     epsilon : TYPE
#         DESCRIPTION.
#     z_gamma : TYPE
#         DESCRIPTION.
#     sign : TYPE, optional
#         DESCRIPTION. The default is 1.
#     xi_nu : TYPE, optional
#         DESCRIPTION. The default is 0.

#     Returns
#     -------
#     TYPE
#         DESCRIPTION.

#     """
#     T_over_nu=func_from_pkl('T_over_T_nu')
#     T=0.511/z_gamma
#     K=MeV_to_1_over_s(1/(1.636*880.2),b=0)
#     q=1.2933/0.511
#     func3=((epsilon-q)**2*np.sqrt((epsilon**2-1))*epsilon)/((1+np.exp(epsilon*z_gamma))*(1+np.exp(-(epsilon-q)*z_gamma*T_over_nu(T)+xi_nu)))
#     return K*func3



# def n_to_p_int_1(epsilon,z_gamma,sign=1,xi_nu=0):
#     """
#     f_nu*(1-f_e), E_e=E_nu+dm
#     integrated from q to infinity

#     Parameters
#     ----------
#     epsilon : TYPE
#         DESCRIPTION.
#     z_gamma : TYPE
#         DESCRIPTION.
#     sign : TYPE, optional
#         DESCRIPTION. The default is 1.
#     xi_nu : TYPE, optional
#         DESCRIPTION. The default is 0.

#     Returns
#     -------
#     TYPE
#         DESCRIPTION.

#     """
#     T_over_nu=func_from_pkl('T_over_T_nu')
#     T=0.511/z_gamma
#     K=MeV_to_1_over_s(1/(1.636*880.2),b=0)
#     q=1.2933/0.511
#     func1=((epsilon-q)**2*np.sqrt((epsilon**2-1))*epsilon)/((1+np.exp(-epsilon*z_gamma))*(1+np.exp((epsilon-q)*z_gamma*T_over_nu(T)-xi_nu)))
#     return K*func1

# def n_to_p_int_2(epsilon,z_gamma,sign=1,xi_nu=0):
#     """
#     f_e*(1-f_nu), E_e=E_nu+dm
#     integrated from q to infinity
#     Parameters
#     ----------
#     epsilon : Electron energy over electron mass, parameter to be integrated.
#     z_gamma : TYPE
#         DESCRIPTION.
#     sign : TYPE, optional
#         DESCRIPTION. The default is 1.
#     xi_nu : TYPE, optional
#         Neutrino Chemical potential. The default is 0.

#     Returns
#     -------
#     TYPE
#         DESCRIPTION.

#     """
#     T_over_nu=func_from_pkl('T_over_T_nu')
#     T=0.511/z_gamma
#     K=MeV_to_1_over_s(1/(1.636*880.2),b=0)
#     q=1.2933/0.511
#     func2=((epsilon-q)**2*np.sqrt((epsilon**2-1))*epsilon)/((1+np.exp(epsilon*z_gamma))*(1+np.exp(-(epsilon-q)*z_gamma*T_over_nu(T)+xi_nu)))
#     return K*func2

# def n_to_p_int_3(epsilon,z_gamma,sign=1,xi_nu=0):
#     """
#     f_anti_nu*f_e, E_e+E_nu=dm
#     integrated from q to infinity

#     Parameters
#     ----------
#     epsilon : TYPE
#         DESCRIPTION.
#     z_gamma : TYPE
#         DESCRIPTION.
#     sign : TYPE, optional
#         DESCRIPTION. The default is 1.
#     xi_nu : TYPE, optional
#         DESCRIPTION. The default is 0.

#     Returns
#     -------
#     TYPE
#         DESCRIPTION.

#     """
#     T_over_nu=func_from_pkl('T_over_T_nu')
#     T=0.511/z_gamma
#     K=MeV_to_1_over_s(1/(1.636*880.2),b=0)
#     q=1.2933/0.511
#     func3=((epsilon-q)**2*np.sqrt((epsilon**2-1))*epsilon)/((1+np.exp(epsilon*z_gamma))*(1+np.exp(-(epsilon-q)*z_gamma*T_over_nu(T)-xi_nu)))
#     return K*func3

# def gamma_p_to_n(z_gamma):
#     q=1.2933/0.511
#     integral1=quad_vec(lambda epsilon: n_to_p_int_1(epsilon,z_gamma,sign=-1),q,np.inf)
#     integral2=quad_vec(lambda epsilon: n_to_p_int_2(epsilon,z_gamma,sign=-1),1,np.inf)
#     return integral1[0]+integral2[0]




# def T_BBN(eta):#problem is in the vectorization somehow!!!!!!
#     B=2.23
#     try:
#         guess=[B/0.08]*len(eta)
#     except TypeError:
#         guess=B/0.08
#     sol=newton(BBN_func2,guess,args=(eta,),full_output=False)
#     return B/sol
    
    # B=2.23
    # m_n=939
    # c=m_n**1.5*np.sqrt(np.pi)/(12*zeta(3)*eta)#Is it 24 or 12? Because it does indeed change the numbers
    # T_inv = Symbol('T_inv')
    # T_n=np.array(solve(T_inv**-1.5*sp.exp(B*T_inv)-c, T_inv))
    # return (T_n**-1)[-1]