# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 18:44:03 2017

@author: 93137

Edited on Thu Mar 16 23:15:51 2017

@author: 93137
"""
import ase.io.vasp as vp
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Convert miller index.')
parser.add_argument('integers', metavar='N', type=int, nargs='+',
                    help='an integer for the accumulator')
parser.add_argument('-m', dest='miller', action='store_const',
                    const='hkl2uvw()', default=max,
                    help='convert miller index to orientation index')
parser.add_argument('-i', dest='orientation', action='store_const',
                    const=sum, default=max,
                    help='convert orientation index to miller index')
parser.add_argument('integers', metavar='N', type=int, nargs='+',
                    help='an integer for the accumulator')
                    
args = parser.parse_args()                    
            
def gcd(a,b): #calculate the Greatest Common Divisor of a and b
    while (b != 0) :
        c = a % b
        a = b
        b = c
    return a

def lcm(a,b): #calculate the Lowest Common Multiple of a and b
    return a*b/gcd(a,b)
    
def absvec(abc): #make the coordinates of vector positive
    n = 0
    for i in range(3):
        if (abc[i] < 0):
            n = n + 1
    if n >= 2:
        for i in range(3):
            abc[i] = (-1)*abc[i]
    return abc

def hkl2uvw(hkl,a,b,c):
    h,k,l = hkl[0],hkl[1],hkl[2]
    aa,bb,cc,ab,ac,bc = np.dot(a,a),np.dot(b,b),np.dot(c,c),np.dot(a,b),np.dot(a,c),np.dot(b,c)
    v1 = (l*bc-k*cc)*(h*ab-k*aa)-(k*ac-h*bc)*(k*ac-l*ab)
    u1 = (l*bc-k*cc)*(k*ab-h*bb)-(k*ac-h*bc)*(l*bb-k*bc)
    w1 = (k*ab-h*bb)*(k*ac-l*ab)-(l*bb-k*bc)*(h*ab-k*aa) #we have already solve the value of v/u and w/u
    gcdOfuvw = gcd(gcd(abs(u1),abs(v1)),abs(w1))#code below aims at turn u,v,w into their smallest integers
    u = u1 / gcdOfuvw
    v = v1 / gcdOfuvw
    w = w1 / gcdOfuvw
    return absvec([u,v,w])

def uvw2hkl(uvw,a,b,c):
    u,v,w = uvw[0],uvw[1],uvw[2]
    aa,bb,cc,ab,ac,bc = np.dot(a,a),np.dot(b,b),np.dot(c,c),np.dot(a,b),np.dot(a,c),np.dot(b,c)
    h1 = u*aa + v*ab + w*ac
    k1 = u*ab + v*bb + w*bc
    l1 = u*ac + v*bc + w*cc #we have already solve the value of h/k and l/k
    gcdOfhkl = gcd(gcd(abs(h1),abs(k1)),abs(l1)) #code below aims at turn h,k,l into their smallest integers
    h = h1 / gcdOfhkl
    k = k1 / gcdOfhkl
    l = l1 / gcdOfhkl
    return absvec([h,k,l])



a = np.array([4,2,1])
b = np.array([5,3,2])
c = np.array([3,1,5])
uvw = hkl2uvw([4,2,3],a,b,c)
print uvw
hkl = uvw2hkl(uvw,a,b,c)
print hkl
