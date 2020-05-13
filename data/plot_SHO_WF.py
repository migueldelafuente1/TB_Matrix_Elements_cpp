#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 23:43:57 2018

@author: miguel
"""

from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":  
    file  = 'SHO.txt'
    data = np.loadtxt(file,delimiter = ' ',dtype = str)

    r = map(float, data[:,0])
    psi = map(float, data[:,1])
    psi_cuadrado = map(float, data[:,2])
    
    # Integral en el radio 
    suma = psi_cuadrado[0]
    dr = r[1] - r[0]
    for i in range(1,len(psi)-1):
        suma += 2*(psi_cuadrado[i])*(r[i]**2)
    
    suma += psi_cuadrado[len(psi)-1]*(r[len(psi)-1]**2)
    suma *= 0.5 * dr
    print ('Norma=',suma)
    
    fig = plt.figure()
    plt.title('Comprobacion de las funciones de onda del oscilador armonico')
    plt.ylabel(r'$\psi_{n,l}$')
    plt.xlabel('r (fm)')
    plt.plot(r,psi)