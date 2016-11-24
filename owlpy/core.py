'''--------------------------- Core Script ---------------------------------'''
'''
    Description: This library is based on the algorithms described in 
    [1] Chin-Chia Michael Yeh, Yan Zhu, Liudmila Ulanova, Nurjahan Begum, 
        Yifei Ding, Hoang Anh Dau, Diego Furtado Silva, Abdullah Mueen, 
        Eamonn Keogh (2016). Matrix Profile I: All Pairs Similarity Joins 
        for Time Series: A Unifying View that Includes Motifs, Discords and 
        Shapelets. IEEE ICDM 2016.
'''

import numpy

def slidingDotProduct(Q, T):
    n = len(T)
    m = len(Q)
    
    Ta = T
    [Ta.append(0) for i in range(0, n)]
    
    Qr = reversed(Q)
    
    Qra = Qr
    [Qra.append(0) for i in range(0, 2*n - m)]
    
    Qraf = numpy.fft.fft(Qra)
    Taf = numpy.fft.fft(Ta)
    
    QT = numpy.fft.ifft(numpy.multiply(Qraf, Taf))
    
    return QT
    
def calculateDistanceProfile(Q, T, QT, mean_Q, std_Q, MT, ET):
    n = len(T)
    m = len(Q)
    D = []
    
    for i in range(0, n - m):
        a = QT[i] - m*mean_Q*MT[i]
        b = m*std_Q*ET[i]
        c = a/b
        
        di = numpy.sqrt(2*m*(1 - c))
        D.append(di)
        
    return D
    
    
# MUEEN’S ALGORITHM FOR SIMILARITY SEARCH (MASS)
def mass(Q, T):
    n = len(T)
    m = len(Q)
    
    QT = slidingDotProduct(Q, T)
    
    mean_Q = numpy.mean(Q)
    std_Q = numpy.std(Q)
    MT = []
    ET = []
    
    # The code below takes O(m) for each subsequence
    # you should replace it for MASS
    [MT.append(numpy.mean(T[range(i, i + m)])) for i in range(0, n - m)]
    [ET.append(numpy.std(T[range(i, i + m)])) for i in range(0, n - m)]
    
    D = calculateDistanceProfile(Q, T, QT, mean_Q, std_Q, MT, ET)
    
    return D
    
def stamp(Ta, Tb, m):
    nb = len(Tb)
    Pab = 0
    
    
    