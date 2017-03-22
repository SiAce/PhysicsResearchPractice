# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 18:44:03 2017

@author: 93137
"""
import ase
import ase.io.vasp as vp
import numpy as np
import argparse

tolerance = 1e-5 #the tolerance of 


parser = argparse.ArgumentParser(description='Convert miller index.')
parser.add_argument('index', metavar='INDEX', type=str,
                    help='the index can either be miller index or orientation index')
parser.add_argument('poscar', metavar='POSCAR', type=str,
                    help='poscar file where basis vectors come from')
parser.add_argument('-m', dest='miller', action='store_true',                    
                    help='convert miller index to orientation index')
parser.add_argument('-i', dest='orientation', action='store_true',
                    help='convert orientation index to miller index')
args = parser.parse_args()
abc = ase.Atoms.get_cell(vp.read_vasp(args.poscar))

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

if args.miller :
    print hkl2uvw(stridx2intidx(args.index),abc,tolerance)
if args.orientation :
    print tuple(uvw2hkl(stridx2intidx(args.index),abc,tolerance))
