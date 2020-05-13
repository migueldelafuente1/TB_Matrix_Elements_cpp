# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 10:41:32 2018

@author: Miguel
"""
from __future__ import division, print_function
from sympy.physics.wigner import racah, clebsch_gordan, wigner_9j
from sympy import S
import numpy as np
import re

if __name__ == "__main__":
#    
#%% Checking Racah Coeficients
    file  = 'racah_for_python.txt'
    data = np.loadtxt(file,delimiter = ' ',dtype = str)

    A = map(int, data[:,0])
    B = map(int, data[:,1])
    C = map(int, data[:,2])
    D = map(int, data[:,3]) 
    E = map(int, data[:,4])
    F = map(int, data[:,5])
#    W_coeff = map(float, data[:,6])
#    errors,non_zero = 0,0
#    tol = 1e-5
#    Erroneos,No_erroneos = [],[]
#    
#    for i in range(len(A)):
#        W,a,b,c,d,e,f = W_coeff[i],A[i],B[i],C[i],D[i],E[i],F[i]
#        w_analitic = racah(a,b,c,d,e,f)
#        if abs(W-w_analitic)>tol:
#            errors += 1
#            non_zero +=1
##            print(a,b,c,d,e,f,"W=",W," // w_an=",w_analitic," //Error")
#            Erroneos.append([[a,b,c,d,e,f],W,w_analitic]);
##            raw_input("press any key")
#        else:
##            print(a,b,c,d,e,f,"W=",W," // w_an=",w_analitic)
#            if(abs(W)>tol):
#                No_erroneos.append([[a,b,c,d,e,f],W,w_analitic]);
#                non_zero += 1
#            pass
#    print("Racah Coefficients:: ",errors ," errores de ",len(A)," /Ceros No triviales" )
##    

#%% Cheching Half integer Racah coefficients  
    ## far easier than the ccgg process
#    file  = 'racah_hf_for_python.txt'
#    data = np.loadtxt(file,delimiter = ' ',dtype = str)
#    J1,J2,J,M1,M2,M = [],[],[],[],[],[]
#    aux = []
#    Data = []
#    
#    for i in range(len(data)): 
#        Data_aux = []
#        # Data aux = [j1,j2,j,m1,m2,m]
#        for j in range(len(data[i])-1):
#            if '/' in data[i][j]:
#                aux = re.split('/',data[i][j])
#                aux[0] = int(aux[0])
#                aux[1] = int(aux[1])
#            else:
#                aux = [int(data[i][j]),1]
#            Data_aux.append(aux)
#        Data_aux.append(float(data[i][6]))
#        Data.append(Data_aux)
#        
#    errors,non_zero = 0,0
#    tol = 1e-5
#    Erroneos,No_erroneos = [],[]
#    
#    for i in range(len(Data)):
#        if Data[i][0][1]==2:
#            a = S(Data[i][0][0])/2
#        else:
#            a = Data[i][0][0]
#        
#        if Data[i][1][1]==2:
#            b = S(Data[i][1][0])/2
#        else:
#            b = Data[i][1][0]
#        
#        if Data[i][2][1]==2:
#            c = S(Data[i][2][0])/2
#        else:
#            c = Data[i][2][0]
#        
#        if Data[i][3][1]==2:
#            d = S(Data[i][3][0])/2
#        else:
#            d = Data[i][3][0]
#        
#        if Data[i][4][1]==2:
#            e = S(Data[i][4][0])/2
#        else:
#            e = Data[i][4][0]
#        
#        if Data[i][5][1]==2:
#            f = S(Data[i][5][0])/2
#        else:
#            f = Data[i][5][0]
#        
#
#        
#        w_numeric = Data[i][6]
#
#        try:
#            w_analitic = racah(a,b,c,d,e,f)
#        except ValueError:
#            w_analitic = 0
#        if abs(w_analitic-w_numeric)>tol:
#            Erroneos.append([Data[i],w_analitic])
#            errors+=1
#    print("Racah Coeficientes:: ",errors ," Errores de ",len(Data))

#%% Checking CG
#    file  = 'ccgg_for_python.txt'
#    data = np.loadtxt(file,delimiter = ' ',dtype = str)
#
#    J1 = map(int, data[:,0])
#    J2 = map(int, data[:,1])
#    J = map(int, data[:,2])
#    M1 = map(int, data[:,3]) 
#    M2 = map(int, data[:,4])
#    M = map(int, data[:,5])
#    CG = map(float, data[:,6])
#    errors,non_zero = 0,0
#    tol = 1e-5
#    Erroneos,No_erroneos = [],[]
#    
#    for i in range(len(J1)):
#        cg,j1,j2,j,m1,m2,m = CG[i],J1[i],J2[i],J[i],M1[i],M2[i],M[i]
#        cg_analitic = clebsch_gordan(j1, j2, j, m1, m2, m)
#        if abs(cg-cg_analitic)>tol:
#            errors += 1
#            non_zero +=1
##            print(j1, j2, j, m1, m2, m,"cg=",cg," // cg_an=",cg_analitic," //Error")
#            Erroneos.append([[j1,j2,j,m1,m2,m],cg,cg_analitic]);
##            raw_input("press any key")
##        else:
##            print(j1, j2, j, m1, m2, m,"cg=",cg," // cg_an=",cg_analitic)
###            if(abs(cg)>tol):
##                No_erroneos.append([[j1,j2,j,m1,m2,m],cg,cg_analitic]);
##                non_zero += 1
##            pass
#    print("Clebsh-Gordan:: ",errors ," Errores de ",len(J1))

#%% Checking Half Integer CG
#
#    file  = 'ccgg_hf_for_python.txt'
#    data = np.loadtxt(file,delimiter = ' ',dtype = str)
#    J1,J2,J,M1,M2,M = [],[],[],[],[],[]
#    aux = []
#    Data = []
#    
#    for i in range(len(data)): 
#        Data_aux = []
#        # Data aux = [j1,j2,j,m1,m2,m]
#        for j in range(len(data[i])-1):
#            if '/' in data[i][j]:
#                aux = re.split('/',data[i][j])
#                aux[0] = int(aux[0])
#                aux[1] = int(aux[1])
#            else:
#                aux = [int(data[i][j]),1]
#            Data_aux.append(aux)
#        Data_aux.append(float(data[i][6]))
#        Data.append(Data_aux)
#    
#
#
#    errors,non_zero = 0,0
#    tol = 1e-5
#    Erroneos,No_erroneos = [],[]
#    
#    for DA in Data:
#        if DA[0][1] == 1:
#            if DA[1][1] == 1:
#                if DA[2][1] == 1:
#                    # EEE
#                    cg_analitic = clebsch_gordan(DA[0][0],DA[1][0],DA[2][0],\
#                                                 DA[3][0],DA[4][0],DA[5][0])
#                    if abs(DA[6]-cg_analitic)>tol:
#                        errors += 1
#                        Erroneos.append([DA[0:6],DA[6],cg_analitic]);
##                        print(DA[0:5],"cg=",DA[6]," // cg_an=",cg_analitic," //Error")
##                        raw_input("press any key")
#                    else:
#                        pass
#                        # no almacena el correcto
#                    
#            else:
#                if DA[2][1] == 2:
#                    # ESS
#                    cg_analitic = clebsch_gordan(DA[0][0],S(DA[1][0])/DA[1][1],\
#                                                 S(DA[2][0])/DA[2][1],DA[3][0],\
#                                                 S(DA[4][0])/DA[4][1],\
#                                                 S(DA[5][0])/DA[5][1])
#                    if abs(DA[6]-cg_analitic)>tol:
#                        errors += 1
#                        Erroneos.append([DA[0:6],DA[6],cg_analitic]);
##                        print(DA[0:6],"cg=",DA[6]," // cg_an=",cg_analitic," //Error")
##                        raw_input("press any key")
#                    else:
#                        pass
#                        # no almacena el correcto
#                else:
#                    pass
#        else:
#            if DA[1][1] == 1:
#                if DA[2][1] == 2:
#                    # SES
#                    cg_analitic = clebsch_gordan(S(DA[0][0])/DA[0][1],DA[1][0],\
#                                                 S(DA[2][0])/DA[2][1],\
#                                                 S(DA[3][0])/DA[3][1],DA[4][0],\
#                                                 S(DA[5][0])/DA[5][1])
#                    if abs(DA[6]-cg_analitic)>tol:
#                        errors += 1
#                        Erroneos.append([DA[0:6],DA[6],cg_analitic]);
##                        print(DA[0:5],"cg=",DA[6]," // cg_an=",cg_analitic," //Error")
##                        raw_input("press any key")
#                    else:
#                        pass
#                        # no almacena el correcto
#                else:
#                    pass
#            else:
#                if DA[2][1] == 1:
#                    # SSE
#                    cg_analitic = clebsch_gordan(S(DA[0][0])/DA[0][1],\
#                                                 S(DA[1][0])/DA[1][1],DA[2][0],\
#                                                 S(DA[3][0])/DA[3][1],\
#                                                 S(DA[4][0])/DA[4][1],DA[5][0])
#                    if abs(DA[6]-cg_analitic)>tol:
#                        errors += 1
#                        Erroneos.append([DA[0:6],DA[6],cg_analitic]);
##                        print(DA[0:5],"cg=",DA[6]," // cg_an=",cg_analitic," //Error")
##                        raw_input("press any key")
#                    else:
#                        pass
#                        # no almacena el correcto
#                else:
#                    pass
#                    
#    print("Clebsh-Gordan:: ",errors ," Errores de ",len(Data))
##    
#%% Checking 9j Symbols
#    file  = '9j_for_python.txt'
#    data = np.loadtxt(file,delimiter = ' ',dtype = str)
#
#    j1 = map(int, data[:,0])
#    j2 = map(int, data[:,1])
#    j3 = map(int, data[:,2])
#    j4 = map(int, data[:,3]) 
#    j5 = map(int, data[:,4])
#    j6 = map(int, data[:,5])
#    j7 = map(int, data[:,6])
#    j8 = map(int, data[:,7])
#    j9 = map(int, data[:,8])
#    W_9j = map(float, data[:,9])
#    errors,non_zero = 0,0
#    tol = 1e-5
#    Erroneos,No_erroneos = [],[]
#    
#    for i in range(len(j1)):
#        w_9j,j_1, j_2, j_3, j_4, j_5, j_6, j_7, j_8, j_9 = W_9j[i],j1[i], j2[i], j3[i], j4[i], j5[i], j6[i], j7[i], j8[i], j9[i]
#        w9j_analitic = wigner_9j(j_1, j_2, j_3, j_4, j_5, j_6, j_7, j_8, j_9, prec=None)
#        if abs(w_9j-w9j_analitic)>tol:
#            errors += 1
#            non_zero +=1
##            print(j1, j2, j, m1, m2, m,"cg=",cg," // cg_an=",cg_analitic," //Error")
#            Erroneos.append([[j_1, j_2, j_3, j_4, j_5, j_6, j_7, j_8, j_9],w_9j,w9j_analitic]);
##            raw_input("press any key")
##        else:
##            print(j1, j2, j, m1, m2, m,"cg=",cg," // cg_an=",cg_analitic)
###            if(abs(cg)>tol):
##                No_erroneos.append([[j1,j2,j,m1,m2,m],cg,cg_analitic]);
##                non_zero += 1
##            pass
#    print("9j Coefficients:: ",errors ," Errores de ",len(j1))
#    raw_input("pulse una tecla para cerrar")
    
#%% Checking Half Integer 9j Symbols
    file  = '9j_for_python.txt'
    data = np.loadtxt(file,delimiter = ' ',dtype = str)

    j1 = map(float, data[:,0])
    j2 = map(float, data[:,1])
    j3 = map(float, data[:,2])
    j4 = map(float, data[:,3]) 
    j5 = map(float, data[:,4])
    j6 = map(float, data[:,5])
    j7 = map(float, data[:,6])
    j8 = map(float, data[:,7])
    j9 = map(float, data[:,8])
    W_9j = map(float, data[:,9])
    errors,non_zero = 0,0
    tol = 1e-5
    Erroneos,No_erroneos = [],[]
    
    for i in range(len(j1)):
        w_9j, j_1, j_2, j_3, j_4, j_5, j_6, j_7, j_8, j_9 = \
        W_9j[i], j1[i], j2[i], j3[i], j4[i], j5[i], j6[i], j7[i], j8[i], j9[i]
        try:
            w9j_analitic = wigner_9j(j_1, j_2, j_3,
                                 j_4, j_5, j_6, j_7, j_8, j_9, prec=None)
        except ValueError:
            w9j_analitic = 0
        if abs(w_9j-w9j_analitic)>tol:
            errors += 1
            non_zero +=1
#            print(j1, j2, j, m1, m2, m,"cg=",cg," // cg_an=",cg_analitic," //Error")
            Erroneos.append([[j_1, j_2, j_3, j_4, j_5, j_6, j_7, j_8, j_9],w_9j,w9j_analitic]);
            raw_input("press any key")
#        else:
#            print(j1, j2, j, m1, m2, m,"cg=",cg," // cg_an=",cg_analitic)
##            if(abs(cg)>tol):
#                No_erroneos.append([[j1,j2,j,m1,m2,m],cg,cg_analitic]);
#                non_zero += 1
#            pass
    print("9j Coefficients:: ",errors ," Errores de ",len(j1))
    raw_input("pulse una tecla para cerrar")