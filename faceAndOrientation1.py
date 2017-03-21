# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 18:44:03 2017

@author: 93137
"""
import numpy as np

def gcd(a,b): #calculate the greatest common divisor of a and b
    while (b != 0) :
        c = a % b
        a = b
        b = c
    return a

def lcm(a,b): #calculate the lowest common multiple of a and b
    return a*b/gcd(a,b)

def hkl2uvw(h,k,l):
    a,b,c = [0,0,0],[0,0,0],[0,0,0]
    aa,bb,cc,ab,ac,bc = 0,0,0,0,0,0
    u,v,w,u1,v1,w1,u2,u3,v3,w3 = 0,0,0,0,0,0,0,0,0,0
    for i in range(3):
        print "a[%d]:" %i
        a[i] = input()
    for i in range(3):
        print "b[%d]:" %i
        b[i] = input()
    for i in range(3):
        print "c[%d]:" %i
        c[i] = input()
    for i in range(3):
        aa += a[i]**2
    for i in range(3):
        bb += b[i]**2
    for i in range(3):
        cc += c[i]**2
    for i in range(3):
        ab += a[i]*b[i]
    for i in range(3):
        ac += a[i]*c[i]
    for i in range(3):
        bc += b[i]*c[i]
    v1 = (l*bc-k*cc)*(h*ab-k*aa)-(k*ac-h*bc)*(k*ac-l*ab)
    u1 = (l*bc-k*cc)*(k*ab-h*bb)-(k*ac-h*bc)*(l*bb-k*bc)
    w1 = (l*bb-k*bc)*(h*ab-k*aa)-(k*ab-h*bb)*(k*ac-l*ab)
    u2 = (l*bb-k*bc)*(k*ab-h*bc)-(k*ab-h*bb)*(l*bc-k*cc)
    u3 = lcm(u1,u2)
    v3 = v1 * u3 / u1
    w3 = w1 * u3 / u1
    u = u3 / gcd(gcd(u3,v3),w3)
    v = v3 / gcd(gcd(u3,v3),w3)
    print u,v,w

def uvw2hkl(u,v,w):
    a,b,c = [0,0,0],[0,0,0],[0,0,0]
    aa,bb,cc,ab,ac,bc = 0,0,0,0,0,0  
    h1,k1,l1,h,k,l = 0,0,0,0,0,0  
    for i in range(3):
        print "a[%d]:" %i
        a[i] = input()
    for i in range(3):
        print "b[%d]:" %i
        b[i] = input()
    for i in range(3):
        print "c[%d]:" %i
        c[i] = input()
    for i in range(3):
        aa += a[i]**2
    for i in range(3):
        bb += b[i]**2
    for i in range(3):
        cc += c[i]**2
    for i in range(3):
        ab += a[i]*b[i]
    for i in range(3):
        ac += a[i]*c[i]
    for i in range(3):
        bc += b[i]*c[i]
    h1 = u*aa + v*ab + w*ac
    k1 = u*ab + v*bb + w*bc
    l1 = u*ac + v*bc + w*cc
    h = h1 / gcd(gcd(h1,k1),l1)
    k = k1 / gcd(gcd(h1,k1),l1)
    l = l1 / gcd(gcd(h1,k1),l1)
    print h,k,l

hkl2uvw(1,1,1)