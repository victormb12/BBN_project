#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 19:14:05 2021

@author: victormaura12
"""

import numpy as np
import matplotlib.pyplot as plt
import Helium_functions as He
import pandas as pd
import itertools
from matplotlib import rcParams
from cycler import cycler

rcParams['lines.linewidth'] = 1.5
rcParams['font.size']=11
rcParams['axes.prop_cycle']= cycler('color',list('grkby'))
plt.style.use('seaborn-ticks')

"""

-------------------------------------------------------------------------------------------------------------------------------------
Helium abundnace for given eta, Neff, xi_nu
-------------------------------------------------------------------------------------------------------------------------------------

"""

eta=np.array([6e-10])
dN_eff=np.linspace(-1,1,20)#np.array([0])
xi_nu=np.linspace(0,0.12,20)#,np.array([0])
Y=He.dict_from_pkl('test_dict')#He.Y_He_fin(eta, dN_eff, xi_nu,T_d=1,verbose=True)#
array=np.zeros((len(xi_nu),len(dN_eff)))

xixi,dNeff=np.meshgrid(xi_nu,dN_eff)
for i,j in itertools.product(range(len(xi_nu)),range(len(dN_eff))):
    array[i,j]=Y[(dN_eff[j],xi_nu[i])]['$Y_{He}$'].iloc[0]
Z=He.Delta_D_H_control(xixi,dNeff,0,0,order=3,return_absolute=True)
z2=array.T+Z
plt.contourf(xixi,dNeff,array.T-0.245,[-0.003,0.003],colors='green',alpha=0.75)
plt.contourf(xixi,dNeff,Z-2.547,[-0.025,0.025],colors='red',alpha=0.75)
plt.contour(xixi,dNeff,array.T-0.245,[-0.003,0.003],colors='green',alpha=1,linestyles='solid')
plt.contour(xixi,dNeff,Z-2.547,[-0.025,0.025],colors='red',alpha=1,linestyles='solid')
plt.grid(True)
plt.plot([0],[0],color='green',label='Helium')
plt.plot([0],[0],color='red',label='Deuterium')
plt.legend(loc='lower right')
plt.gca().set_xlabel(r'$\xi_{\nu}$')
plt.gca().set_ylabel(r'$\Delta N_{eff}$')
plt.gca().set_title(r'Parameter Space of allowed $\Delta N_{eff}$ as a Function of $\xi_\nu$')
plt.gca().set_xlim((0,0.06))
plt.gca().set_ylim((-0.25,1))
plt.tight_layout()
