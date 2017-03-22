# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 17:52:37 2017

@author: 93137
"""
a = -0.1312312312312312312312312312312
print int(a/abs(a))

def float2fraction(num='numerator',den='denominator',tol='tolerance'):
    if (abs(num) <= 1E-15)&(abs(den) <= 1E-15):
        return [0,0]
    elif abs(num) <= 1E-15:
        return [0,int(den/abs(den))]
    elif abs(num) <= 1E-15:
        return [int(num/abs(num)),0]
    else:
        intden = int(round(1 / tol))
        intnum = int((float(num) / float(den)) * intden) 
        gcd_numden = gcd(intden,intnum)
        return [intnum/gcd_numden,intden/gcd_numden]
    
def gcd(a,b): #calculate the Greatest Common Divisor of a and b
    while (b != 0) :
        c = a % b
        a = b
        b = c
    return a
    
print float2fraction(0,-23.23,2e-5)