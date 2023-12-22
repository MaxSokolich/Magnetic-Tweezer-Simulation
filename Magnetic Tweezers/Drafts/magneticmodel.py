#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 12:07:40 2021

@author: bizzaro
"""
#imports
import magpylib as mag3
s = mag3.magnet.Cylinder(magnetization=(0,0,350), dimension=(4,5)) 
observer_pos = (4,4,4)
print(s.getB(observer_pos))