# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 13:44:28 2018

@author: Miguel
"""
from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt

class Q_Numbers:
#    def __init__(self, n1,l1,j1,n2,l2,j2,J,M):
    def __init__(self, array):
        aux = array.split(' ')
        self.n1 = int(aux[0])
        self.l1 = int(aux[1])
        self.j1 = aux[2]
        self.n2 = int(aux[3])
        self.l2 = int(aux[4])
        self.j2 = aux[5]
        self.J = int(aux[6])
        self.M = int(aux[7])
        
#        self.n1 = n1
#        self.l1 = l1
#        self.j1 = j1 
#        self.n2 = n2
#        self.l2 = l2
#        self.j2 = j2
#        self.J = J
#        self.M = M
        
#    def paso_argumentos(self, array):
#        aux = array.split(' ')
#        self.n1 = int(aux[0])
#        self.l1 = int(aux[1])
#        self.j1 = float(aux[2]) 
#        self.n2 = int(aux[3])
#        self.l2 = int(aux[4])
#        self.j2 = float(aux[5])
#        self.J = int(aux[6])
#        self.M = int(aux[7])
        
        
    
if __name__ == "__main__":   
#%% TALMI RESULS    
    file  = 'ME_Talmi_for_Python.txt'
    DATA = np.loadtxt(file,delimiter = '   ',dtype = str)
    
    data_WF, data = [],[]
    
    for i in range(len(DATA)):
        
        data.append([DATA[i][1],DATA[i][2],DATA[i][3],DATA[i][4],DATA[i][5],DATA[i][6],DATA[i][7]])
        data_WF.append( DATA[i][0].split(' / / '))  
    # [0] for BRA  [1] for KET
     
#   
    Len,M_Central,M_LS,M_Tensor,Art_LS,Art_Tensor = [],[],[],[],[],[]
    diff_Tensor, diff_LS= [],[]
#    n1,l1,j1,n2,l2,j2,J,M = [],[],[],[],[],[],[],[]
#    n1_q,l1_q,j1_q,n2_q,l2_q,j2_q,J_q,M_q = [],[],[],[],[],[],[],[]
    
    data_number1,data_number2 = [],[]
    

    for i in range(len(data)):
        
        data_number1.append(Q_Numbers(data_WF[i][0]))
        data_number2.append(Q_Numbers(data_WF[i][1]))
        
        
        
        Len.append(int (data[i][0]))
        M_Central.append(float(data[i][1]))
        M_LS.append(float(data[i][2]))
        M_Tensor.append(float(data[i][3]))
        Art_LS.append(float(data[i][4]))
        Art_Tensor.append(float(data[i][5]))
        
        diff_LS.append(Art_LS[-1] - M_LS[-1])
        diff_Tensor.append(Art_Tensor[-1] - M_Tensor[-1])
    
        
    fig = plt.figure()
    plt.title('Matrix element for LS interaction')
    plt.ylabel('Frecuecia')
    plt.xlabel('difference (MeV)')
#    plt.plot(Len[0:],diff_LS[0:])
    plt.hist(diff_LS[0:],100)
    plt.xlim(-0.05,0.05)
    
    fig = plt.figure()
    plt.title('Matrix element for LS interaction')
    plt.ylabel('difference (MeV)')
    plt.plot(Len[0:],diff_LS[0:],'.')

    fig = plt.figure()
    plt.title('Matrix element for Tensor interaction')
    plt.ylabel('difference (MeV)')
    plt.plot(Len[0:],diff_Tensor[0:],'.')

    
#%% MULTIPOLAR DESCOMPOSITION
    file  = 'ME_Multipolar.txt'
    dat = np.loadtxt(file,delimiter = ' ',dtype = str)
    
    Order_aux = map(int,dat[:,0])
    J_aux = map(int,dat[:,1])
    numeric_aux= map(float,dat[:,2]) 
    
    J,Order = [],[]
    numeric = []
    for i in range(len(Order_aux)):
        if J_aux[i] not in J:
            J.append(J_aux[i])
            numeric.append([])
        
        if Order_aux[i] not in Order:
            Order.append(Order_aux[i])
    
    for i in range(len(Order_aux)):
        for j in range(len(J)):
            if J_aux[i] == J[j]:
                numeric[j].append(numeric_aux[i])
    
    fig = plt.figure()
    for j in range(len(J)):
        plt.plot(Order,numeric[j],'.-',label = J[j])
    
    plt.xlabel('order')
    plt.legend()
    

    