#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 16:56:31 2020

@author: victormaura12
"""

from sympy.solvers import solve
from sympy import Symbol
import sympy as sp
import numpy as np
from scipy.special import zeta
import matplotlib.pyplot as plt




def T_BBN(eta):
    B=2.23
    m_n=939
    c=m_n**1.5*np.sqrt(np.pi)/(24*zeta(3)*eta)
    T_inv = Symbol('T_inv')
    T_n=np.array(solve(T_inv**-1.5*sp.exp(B*T_inv)-c, T_inv))
    return (T_n**-1)[-1]
    


eta=np.logspace(-12,-6,10)
T=[T_BBN(e) for e in eta]
plt.plot(eta,T)
plt.semilogx()

