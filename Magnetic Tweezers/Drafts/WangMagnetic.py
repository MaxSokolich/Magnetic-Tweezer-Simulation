#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 15:24:23 2021

@author: bizzaro
"""

import numpy as np

milli = (10**(-3))
micro = (10**(-6))
nano = (10**(-9))
pico = (10**(-12))

Pole_Pole_Distance = (1* milli)

bead_diameter = 10*micro

perm0 = ((4*np.pi)*10**(-7))   
perm_bead = 40 


km = perm0/(4*np.pi)
Bead_V = ((4*np.pi)/3)*((bead_diameter/2)**3)

L=(Pole_Pole_Distance/2)






kq = ((3*Bead_V*(km**2)*(perm_bead-perm0)))/(2*perm0*(L**4)*(perm_bead+(2*perm0)))
print(kq)
'''
uj = np.array([0,0,0])
uk = np.array([0,0,0])

M = [[u11,u12,u13,u14,u15,u16],
    [u21,u22,u23,u24,u25,u26],
    [u31,u32,u33,u34,u35,u36],
    [u41,u42,u43,u44,u45,u46],
    [u51,u52,u53,u54,u55,u56],
    [u61,u62,u63,u64,u65,u66]]
'''




