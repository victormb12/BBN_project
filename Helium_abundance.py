 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 17:25:21 2020

@author: victormaura12
"""
import numpy as np
import matplotlib.pyplot as plt
import Helium_functions as He
import time
import pandas as pd
start_time=time.time()



# """

# -------------------------------------------------------------------------------------------------------------------------------------
# Boltzmann equation solution
# -------------------------------------------------------------------------------------------------------------------------------------

# """

# T_range=(4,1e-2)
# dN_eff=np.array([0])
# xi_nu=np.array([0])
# He.plot_Boltz(T_range=T_range,dN_eff=dN_eff,xi_nu=xi_nu,n_over_p=True)


# """

# -------------------------------------------------------------------------------------------------------------------------------------
# Helium abundnace as a function of eta for given Neff, xi_nu
# -------------------------------------------------------------------------------------------------------------------------------------

# """
# order=3
# compare=False
# etas=np.logspace(-10,-9,20)      
# eta0=6.1e-10
# detas=etas-eta0
# dN_eff=np.array([0])
# xi_nu=np.array([0])
# plot_Y1=He.plot_He(etas,xi_nu=xi_nu,dN_eff=dN_eff,T_d=1,plot=not compare)
# # plt.gca().set_xlim(1e-10,1e-9)

# if compare: 
#     Y0=0.245144
#     fig,ax=plt.subplots(1,1)
#     ax.plot(detas,(plot_Y1-Y0)/Y0,label='My calculation')
#     ax.plot(detas,He.Delta_Y_P_control(xi_nu,dN_eff,detas,0,order=order,return_absolute=False),label='Pytrou et al.')
#     coeffs=np.polyfit(detas,(plot_Y1-Y0)/Y0,order)
#     print(coeffs)
#     ax.set_xlabel(r'$\Delta \eta_B$')
#     ax.set_ylabel(r'$\frac{\Delta Y_{He}}{Y_{He}}$')
#     ax.set_title(r'Relative Deviation of $Y_{He}$ as a Function of $\Delta \eta_B$ ')
#     # ax.plot(detas,np.polyval(coeffs,detas),label='My calculation')
#     # ax.plot(dN_effs2,He.Delta_Y_P_control(xi_nu2,dN_effs2,0,0,order=order),label='Pytrou et al.')
#     ax.legend()
#     ax.grid(True)
#     plt.tight_layout()

"""

-------------------------------------------------------------------------------------------------------------------------------------
Helium abundnace for given eta, Neff, xi_nu
-------------------------------------------------------------------------------------------------------------------------------------

"""

eta=np.array([6e-10])
dN_eff=np.array([0])#np.linspace(-1,1,2)#
xi_nu=np.array([0])#np.linspace(0,0.012,3)#
Y=He.Y_He_fin(eta, dN_eff, xi_nu,T_d=1,verbose=True)
# df=pd.DataFrame(Y)
print('The Helium Abundance is: Y_He={:.7g}'.format(Y[(dN_eff[0],xi_nu[0])]['$Y_{He}$'].iloc[0]))


# """

# -------------------------------------------------------------------------------------------------------------------------------------
# dN_eff as  function of xi for given Y, eta
# -------------------------------------------------------------------------------------------------------------------------------------

# """

# create_dict1=False
# order=3
# eta1=np.array([6.1e-10])
# Y1=0.245
# xi_nus1=np.linspace(0,0.12,20)
# dN_effs1=np.linspace(-1,1,20)
# fig,ax=plt.subplots(1,1)
# plot_N1,plot_xi1=He.N_v_xi(eta1,Y1,dN_effs=dN_effs1,xi_nus=xi_nus1,epsilon=0.003,create_dict=create_dict1,ax=ax)
# # plot_N1,plot_xi1=He.N_v_xi_control(eta1,Y1,dN_effs=dN_effs1,xi_nus=xi_nus1,epsilon=0.003,order=order)
# plot_Nd,plot_xi2=He.N_v_xi_deuterium(eta1,2.547,dN_effs=dN_effs1,xi_nus=xi_nus1,epsilon=0.25,order=order,ax=ax)


# """

# -------------------------------------------------------------------------------------------------------------------------------------
# Y as a function of Neff for given xi, eta
# -------------------------------------------------------------------------------------------------------------------------------------

# """

# order=3
# compare=True
# create_dict2=False
# eta2=np.array([6.1e-10])
# xi_nu2=np.array([0])
# dN_effs2=np.linspace(-1,1,200)
# fig,ax=plt.subplots(1,1,sharex=True)
# plot_N2,plot_Y2=He.Y_v_N(eta2,xi_nu2,dN_effs2,create_dict=create_dict2,plot=not compare,ax=ax)
# # plot_N2,plot_Y2=He.Y_v_N_control(eta2,xi_nu2,dN_effs2,plot=not compare,order=order,ax=ax)
# plt.tight_layout()

# if compare: 
#     Y0=0.245144
#     fig,ax=plt.subplots(1,1)
#     ax.plot(plot_N2,(plot_Y2-Y0)/Y0,label='My calculation')
#     ax.plot(dN_effs2,He.Delta_Y_P_control(xi_nu2,dN_effs2,0,0,order=order,return_absolute=False),label='Pytrou et al.')
#     coeffs=np.polyfit(plot_N2,(plot_Y2-Y0)/Y0,order)
#     print(coeffs)
#     ax.set_xlabel(r'$\Delta N_{eff}$')
#     ax.set_ylabel(r'$\frac{\Delta Y_{He}}{Y_{He}}$')
#     ax.set_title(r'Relative Deviation of $Y_{He}$ as a Function of $\Delta N_{eff}$ ')
#     # ax.plot(dN_effs2,np.polyval(coeffs,dN_effs2),label='My calculation')
#     # ax.plot(dN_effs2,He.Delta_Y_P_control(xi_nu2,dN_effs2,0,0,order=order),label='Pytrou et al.')
#     ax.legend()
#     ax.grid(True)
#     plt.tight_layout()

# """

# -------------------------------------------------------------------------------------------------------------------------------------
# Y as a function of xi for given Neff,eta
# -------------------------------------------------------------------------------------------------------------------------------------

# """
# compare=True
# create_dict3=False
# dN_eff3=0
# eta3=np.array([6.1e-10])
# xi_nus3=np.linspace(0,0.12,200)
# order=1

# plot_xi3,plot_Y3=He.Y_v_xi(eta3, dN_eff3, xi_nus3,create_dict=create_dict3,plot=not compare)
# # plot_xi3,plot_Y3=He.Y_v_xi_control(eta3, dN_eff3, xi_nus3,plot=not compare,order=order)


# if compare:
#     Y0=0.245144
#     fig,ax=plt.subplots(1,1)
#     # fig,ax=plt.subplots(1,1)
#     # ax.plot(plot_xi3,plot_Y3,label='My calculation')
#     # ax.plot(plot_xi3,He.Delta_Y_P_control(xi_nus3,dN_eff3,0,0,order=order,return_absolute=True),label='Pytrou et al.')
#     ax.plot(plot_xi3,(plot_Y3-Y0)/Y0,label='My calculation')
#     coeffs=np.polyfit(plot_xi3,(plot_Y3-0.244)/0.244,order)
#     print(coeffs)
#     # ax.plot(xi_nus3,np.polyval(coeffs,xi_nus3),label='My calculation')
#     ax.plot(xi_nus3,He.Delta_Y_P_control(xi_nus3,dN_eff=dN_eff3,deta=0,dTau_n=0),label='Pytrou et al.')
#     ax.set_xlabel(r'$\xi _{\nu}$')
#     ax.set_ylabel(r'$\frac{\Delta Y_{He}}{Y_{He}}$')
#     ax.set_title(r'Relative Deviation of $Y_{He}$ as a Function of $\xi _{\nu}$ ')
#     ax.legend()
#     ax.grid(True)
#     plt.tight_layout()
    
# """

# -------------------------------------------------------------------------------------------------------------------------------------
# xi as a function of eta for given Neff,Y
# -------------------------------------------------------------------------------------------------------------------------------------

# """

# create_dict4=False
# dN_eff4=0
# Y4=0.245
# etas4=np.logspace(-12,-7,10)
# xi_nus4=np.linspace(0,0.12,10)

# plot_xi,plot_eta=He.eta_v_xi(Y4,dN_eff4,xi_nus4,etas4,epsilon=0.01,create_dict=create_dict4)


# print("--- %s seconds ---" % (time.time() - start_time))
