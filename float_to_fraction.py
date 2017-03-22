# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 17:52:37 2017

@author: 93137
"""
a = -0.1312312312312312312312312312312
print int(a/abs(a))

def ratioconvert(numerator,denominator,tolerance):
    '''Convert a ratio of two floats into the ratio of two integers.
    Tolerance means the tolerance of the bias bewteen these two ratios.  
    '''
    intden = int(round(1 / tolerance))
    intnum = int((float(numerator) / float(denominator)) * intden) 
    if (abs(numerator) <= 1E-15)&(abs(denominator) <= 1E-15):
        return [0,0]
    elif abs(numerator) <= 1E-15:
        return [0,intden]
    elif abs(denominator) <= 1E-15:
        intnum = int(float(numerator) * intden)
        return [intnum,0]
    else:
        return [intnum,intden]    
    
print ratioconvert(1.2323,1e-16,1e-5)
print ratioconvert(4.2323,1e-16,1e-5)