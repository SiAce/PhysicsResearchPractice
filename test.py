# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 18:44:03 2017

@author: 93137
"""
import ase
import ase.io.vasp as vp
import numpy as np

tolerance = 1e-5 #the tolerance of 

abc = ase.Atoms.get_cell(vp.read_vasp('POSCAR'))

def gcd(a,b): 
    '''Calculate the Greatest Common Divisor of a and b
    '''
    while (b != 0) :
        c = a % b
        a = b
        b = c
    return a
    
def stridx2intidx(stridx):
    '''Convert the index like '1,m3,1' to '1,-3,1' 
    '''
    strarray = stridx.split(',')
    intidx = [0,0,0]
    for i in range(3):
        if strarray[i][0] == 'm':
            intidx[i] = -1 * int(strarray[i][1:])
        else:
            intidx[i] = int(strarray[i])
    return intidx
    
def ratioconvert(numerator,denominator,tolerance):
    '''Convert a ratio of two floats into the ratio of two integers.
    Tolerance means the tolerance of the bias bewteen these two ratios.  
    '''
    prec = 1e-3
    multiple = int(round(np.log10(1 / tolerance)))
    scale = int(round(np.log10(abs(float(numerator) / float(denominator)))))
    scale_num = int(round(np.log10(abs(float(numerator)))))
    intden = 10**multiple  
    if (abs(numerator) <= prec)&(abs(denominator) <= prec):
        return [0,0]
    elif abs(numerator) <= prec:
        return [0,intden]
    elif abs(denominator) <= prec:
        intnum = int(round(float(numerator),multiple - scale_num) * intden)
        return [intnum,0]
    else:
        intnum = int(round(float(numerator) / float(denominator),multiple - scale) * intden) 
        return [intnum,intden]
            
def absvec(abc):
    '''Make the coordinates of vector positive
    '''
    n = 0
    for i in range(3):
        if (abc[i] < 0):
            n = n + 1
    if n >= 2:
        for i in range(3):
            abc[i] = (-1)*abc[i]
    return abc

def hkl2uvw(hkl,abc,tolerance):
    '''Convert the miller index h,k,l into its related orientation index u,v,w
    '''
    h,k,l = hkl[0],hkl[1],hkl[2]
    a,b,c = abc[0],abc[1],abc[2]
    aa,bb,cc,ab,ac,bc = np.dot(a,a),np.dot(b,b),np.dot(c,c),np.dot(a,b),np.dot(a,c),np.dot(b,c)
    u1 = (l*bc-k*cc)*(k*ab-h*bb)-(k*ac-h*bc)*(l*bb-k*bc)
    v1 = (l*bc-k*cc)*(h*ab-k*aa)-(k*ac-h*bc)*(k*ac-l*ab)
    w1 = (k*ab-h*bb)*(k*ac-l*ab)-(l*bb-k*bc)*(h*ab-k*aa)
    #We have already solve the value of v/u and w/u
    #Code below aims at turn u,v,w into their smallest integers
    u2 = ratioconvert(u1,v1,tolerance)[0]
    v2 = ratioconvert(u1,v1,tolerance)[1]
    w2 = ratioconvert(w1,v1,tolerance)[0]
    gcd_uvw = gcd(gcd(abs(u2),abs(v2)),abs(w2))
    u = u2 / gcd_uvw
    v = v2 / gcd_uvw
    w = w2 / gcd_uvw
    return absvec([u,v,w])

def uvw2hkl(uvw,abc,tolerance):
    '''Convert the orientation index h,k,l into its related miller index u,v,w
    '''
    u,v,w = uvw[0],uvw[1],uvw[2]
    a,b,c = abc[0],abc[1],abc[2]
    aa,bb,cc,ab,ac,bc = np.dot(a,a),np.dot(b,b),np.dot(c,c),np.dot(a,b),np.dot(a,c),np.dot(b,c)
    h1 = u*aa + v*ab + w*ac
    k1 = u*ab + v*bb + w*bc
    l1 = u*ac + v*bc + w*cc 
    #We have already solve the value of h/k and l/k
    #Code below aims at turn h,k,l into their smallest integers
    h2 = ratioconvert(h1,k1,tolerance)[0]
    k2 = ratioconvert(h1,k1,tolerance)[1]
    l2 = ratioconvert(l1,k1,tolerance)[0]
    gcd_hkl = gcd(gcd(abs(h2),abs(k2)),abs(l2))
    h = h2 / gcd_hkl
    k = k2 / gcd_hkl
    l = l2 / gcd_hkl
    return absvec([h,k,l])

print hkl2uvw([1,2,3],abc,tolerance)
