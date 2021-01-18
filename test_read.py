#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 12:49:13 2021

@author: victormaura12
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

df=pd.read_csv('SMgstar_2019_Laine.dat',skiprows=tuple(range(6)),delimiter='    ')
df.set_index('T[MeV]',inplace=True)
df=df[['g_eff','h_eff','i_eff']]
df.plot(y='g_eff')
print(df.head())