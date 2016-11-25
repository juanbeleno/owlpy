'''--------------------------- Core Script ---------------------------------'''
'''
    Description: This library is based on the algorithms described in 
    [1] Chin-Chia Michael Yeh, Yan Zhu, Liudmila Ulanova, Nurjahan Begum, 
        Yifei Ding, Hoang Anh Dau, Diego Furtado Silva, Abdullah Mueen, 
        Eamonn Keogh (2016). Matrix Profile I: All Pairs Similarity Joins 
        for Time Series: A Unifying View that Includes Motifs, Discords and 
        Shapelets. IEEE ICDM 2016.
'''
import numpy as np
import pandas as pd
import time

import matplotlib.pyplot as plt

def slidingDotProduct(Q, T):
    n = len(T)
    m = len(Q)
    
    Ta = list(T) # Copy the T without modifying the original
    [Ta.append(0) for i in range(0, n)]
    
    Qr = Q[::-1]
    
    Qra = list(Qr) # Copy the T without modifying the original
    [Qra.append(0) for i in range(0, 2*n - m)]
    
    Qraf = np.fft.fft(Qra)
    Taf = np.fft.fft(Ta)
    
    QT = np.fft.ifft(np.multiply(Qraf, Taf))
    
    return QT[m:n]
    
def calculateDistanceProfile(Q, T, QT, mean_Q, std_Q, MT, ET):
    n = len(T)
    m = len(Q)
    D = []
    
    for i in range(0, n - m):
        a = QT[i] - m*mean_Q*MT[i]
        b = m*std_Q*ET[i]
        c = a/b
        
        di = np.sqrt(2*m*(1 - c))
        D.append(di)
        
    return D

# The code below takes O(m) for each subsequence
# you should replace it for MASS
def computeMeanStd(Q, T):
    n = len(T)
    m = len(Q)
    
    mean_Q = np.mean(Q)
    std_Q = np.std(Q)
    MT = [np.mean(T[i : i + m]) for i in range(0, n - m)]
    ET = [np.std(T[i : i + m]) for i in range(0, n - m)]
    
    return mean_Q, std_Q, MT, ET
    
# MUEEN’S ALGORITHM FOR SIMILARITY SEARCH (MASS)
def mass(Q, T):
    QT = slidingDotProduct(Q, T)
    mean_Q, std_Q, MT, ET = computeMeanStd(Q, T)
    D = calculateDistanceProfile(Q, T, QT, mean_Q, std_Q, MT, ET)
    
    return D
    
def elementWiseMin(Pab, Iab, D, idx):
    Pab_copy = [1000000 for i in range(0, len(Pab))]
    flag = False
    
    for i in range(0, len(D)):
        if(Pab[i] == float('Inf')):
            Pab_copy[i] = D[i]
            Iab[i] = idx
            flag = True
        elif(D[i] < Pab[i]):
            Pab[i] = D[i]
            Iab[i] = idx
            
    if flag:
        return Pab_copy, Iab
    else:
        return Pab, Iab
    
    

def stamp(Ta, Tb, m):
    nb = len(Tb)
    na = len(Ta)
    Pab = [float('Inf') for i in range(0, na - m)]
    Iab = [0 for i in range(0, na - m)]
    idxes = range(0, nb - m)
    
    for idx in idxes:
        print(idx)
        D = mass(Tb[idx : idx + m], Ta)
        Pab, Iab = elementWiseMin(Pab, Iab, D, idx)
        
    return Pab, Iab

# Quick Test
def test_stamp(Ta, Tb, m):
    start_time = time.time()
    
    df_serie_a = pd.read_csv(Ta, header= 0)
    df_serie_b = pd.read_csv(Tb, header= 0)
    
    ts_a = np.ravel(df_serie_a['value']).tolist()
    ts_b = np.ravel(df_serie_b['value']).tolist()
    
    Pab, Iab = stamp(ts_a, ts_b, m)
    print("--- %s seconds ---" % (time.time() - start_time))
    
    plot_graphics(Pab)
    
def plot_graphics(values):
    fig_width = 16
    fig_height = 8
    fig_dpi = 100
    plt.figure(figsize=(fig_width, fig_height), dpi=fig_dpi)

    plt.plot(range(0,len(values)), values, '#ff5722')
    plt.xlabel('Index')
    plt.ylabel('Value')
    plt.show()


