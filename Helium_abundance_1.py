 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 17:25:21 2020

@author: victormaura12
"""
import numpy as np
import matplotlib.pyplot as plt
import Helium_functions1 as He
import time
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


"""

-------------------------------------------------------------------------------------------------------------------------------------
Helium abundnace as a function of eta for given Neff, xi_nu
-------------------------------------------------------------------------------------------------------------------------------------

"""

etas=np.logspace(-9.9,-9.1,100)
Neff=np.array([3])
xi_nu=np.array([0])
He.plot_He(etas,xi_nu=xi_nu,Neff=Neff,T_d=2)



# """

# -------------------------------------------------------------------------------------------------------------------------------------
# Helium abundnace for given eta, Neff, xi_nu
# -------------------------------------------------------------------------------------------------------------------------------------

# """

# eta=np.array([3.25e-10])
# dN_eff=np.array([0])
# xi_nu=np.array([0])
# Y=He.Y_He_fin(eta, dN_eff, xi_nu,T_d=0.073,verbose=True)
# print('The Helium Abundance is: Y_He={:.3g}'.format(Y[(dN_eff[0],xi_nu[0])]['$Y_{He}$'].iloc[0]))


# """

# -------------------------------------------------------------------------------------------------------------------------------------
# Neff as  function of xi for given Y, eta
# -------------------------------------------------------------------------------------------------------------------------------------

# """

# create_dict1=False

# eta1=np.array([3.25e-10])#np.logspace(-10,-9,100)
# Y1=0.245
# xi_nus1=np.linspace(0,0.12,20)
# Neffs1=np.linspace(2,4,40)
# plot_N1,plot_xi1=He.N_v_xi(eta1,Y1,Neffs=Neffs1,xi_nus=xi_nus1,epsilon=0.003,create_dict=create_dict1)

# """

# -------------------------------------------------------------------------------------------------------------------------------------
# Y as a function of Neff for given xi, eta
# -------------------------------------------------------------------------------------------------------------------------------------

# """
# create_dict2=False
# eta2=np.array([5e-10])
# xi_nu2=np.array([0.1])
# Neffs2=np.linspace(2,4,10)
# plot_N2,plot_Y2=He.Y_v_N(eta2,xi_nu2,Neffs2,create_dict=create_dict2,td=2)
# """

# -------------------------------------------------------------------------------------------------------------------------------------
# Y as a function of xi for given Neff,eta
# -------------------------------------------------------------------------------------------------------------------------------------

# """
# create_dict3=False
# Neff3=3.045
# eta3=np.array([5e-10])
# xi_nus3=np.linspace(0,1,100)
# plot_xi3,plot_Y3=He.Y_v_xi(eta3, Neff3, xi_nus3,create_dict=create_dict3)

# """

# -------------------------------------------------------------------------------------------------------------------------------------
# xi as a function of eta for given Neff,Y
# -------------------------------------------------------------------------------------------------------------------------------------

# """
# create_dict4=False
# Neff4=3.045
# Y4=0.245
# etas4=np.logspace(-12,-7,100)
# xi_nus4=np.linspace(0,1,100)

# plot_xi,plot_eta=He.eta_v_xi(Y4,Neff4,xi_nus4,etas4,epsilon=0.001,create_dict=create_dict4)


print("--- %s seconds ---" % (time.time() - start_time))
