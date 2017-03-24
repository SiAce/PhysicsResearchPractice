#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
Created on Wed Mar 15 18:44:03 2017

@author: 93137
"""
import ase
import ase.io.vasp as vp
import numpy as np
import argparse, sys 
from fractions import Fraction
from numpy.linalg import inv


parser = argparse.ArgumentParser(description='Convert miller index.')
parser.add_argument('poscar', metavar='POSCAR', type=str,
                    help='poscar file where basis vectors come from')
parser.add_argument('-m', dest='miller', type=str, 
                    help='Miller index, to be converted to orientation index')
parser.add_argument('-i', dest='orientation', type=str, 
                    help='orientation index, to be converted to Miller index')
parser.add_argument('-d', dest='digits', type=int, default=8, 
                    help='max digits of denominator in h,k,l or u,v,w')
args = parser.parse_args()
atoms = vp.read_vasp(args.poscar)
abc = atoms.get_cell()
reciprocal_cell = atoms.get_reciprocal_cell()

def gcd(m, n): 
    ''' greatest common divisor '''
    m, n = abs(m), abs(n)
    prec = 1E-12
    if m < prec and n < prec: 
        return 1
    if m < prec or n < prec: 
        return max(m, n)
    while m % n != 0: 
        m, n = n, m % n

    return n
    
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


def sign(a): 
    if a > 0: 
        return 1
    else: 
        return -1

def ratioconvert(numerator, denominator, digits=8): 
    '''Convert a ratio of two float into the ratio of two simpliest integer.
    digits: the returned denominator will be within 10**digits 
    '''

    prec = 1E-12
    
    a = float(numerator)
    b = float(denominator)
    if abs(a) < prec and abs(b) < prec: 
        return 0,0
    if abs(a) < prec: 
        return 0, sign(b)
    if abs(b) < prec: 
        return sign(a), 0

    result = Fraction(a/b).limit_denominator(10**digits)

    return sign(a)*abs(result.numerator),sign(b)*abs(result.denominator)
    
    
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

def hkl2uvw(hkl,abc,digits):
    '''Convert the miller index h,k,l into its related orientation index u,v,w
    '''
    prec = 1E-12 

    fake_atoms = ase.Atoms(cell=np.array(abc), symbols=['C'], scaled_positions=[[0,0,0]])
    reciprocal_cell = fake_atoms.get_reciprocal_cell()
    h,k,l = hkl[0],hkl[1],hkl[2]
    a,b,c = abc[0],abc[1],abc[2]

    # NOTE: the direction determined by hkl in reciprocal space, is the same as 
    # the direction determined by uvw in real space, REF: 
    # https://en.wikipedia.org/wiki/Miller_index
    direction = np.dot(hkl, reciprocal_cell)
    uvw_raw = direction.dot(inv(abc)) # non-integer u,v,w
    ind = np.lexsort([abs(uvw_raw)]) # obtain the index of sorting 
    rind = np.lexsort([ind]) # use this to sort-back the sorted array later
    sorted_uvw_raw = uvw_raw[ind]

    if abs(sorted_uvw_raw[2]) < prec: 
        print >> sys.stderr, "the result does not make sense..." 
        return [0,0,0]

    if abs(sorted_uvw_raw[1]) < prec: 
        result_sorted = np.array([0,0,1])
        return list(result_sorted[rind])

    if abs(sorted_uvw_raw[0]) < prec: 
        vw_rationalized = ratioconvert(sorted_uvw_raw[1], sorted_uvw_raw[2], digits) 
        result_sorted = np.array([0, vw_rationalized[0], vw_rationalized[1]])
        return list(result_sorted[rind])

    uw_rationalized = ratioconvert(sorted_uvw_raw[0], sorted_uvw_raw[2], digits) 
    vw_rationalized = ratioconvert(sorted_uvw_raw[1], sorted_uvw_raw[2], digits) 

    uvw_rationalized_sorted = [uw_rationalized[0] * vw_rationalized[1], 
            vw_rationalized[0] * uw_rationalized[1], 
            uw_rationalized[1] * vw_rationalized[1]]
    uvw_rationalized_sorted = np.array(uvw_rationalized_sorted) / gcd(
            gcd(uvw_rationalized_sorted[0], uvw_rationalized_sorted[1]), 
            uvw_rationalized_sorted[2])

    return list(uvw_rationalized_sorted[rind])

def uvw2hkl(uvw,abc,digits):
    '''Convert the orientation index h,k,l into its related miller index u,v,w
    '''
    prec = 1E-12 

    fake_atoms = ase.Atoms(cell=np.array(abc), symbols=['C'], scaled_positions=[[0,0,0]])
    reciprocal_cell = fake_atoms.get_reciprocal_cell()
    u,v,w = uvw[0],uvw[1],uvw[2]
    a,b,c = abc[0],abc[1],abc[2]

    direction = np.dot(uvw, abc)
    hkl_raw = direction.dot(inv(reciprocal_cell)) # non-integer h,k,l
    ind = np.lexsort([abs(hkl_raw)]) 
    rind = np.lexsort([ind]) # use this to sort-back the sorted array 
    sorted_hkl_raw = hkl_raw[ind]

    if abs(sorted_hkl_raw[2]) < prec: 
        print >> sys.stderr, "the result does not make sense..." 
        return [0,0,0]

    if abs(sorted_hkl_raw[1]) < prec: 
        result_sorted = np.array([0,0,1])
        return list(result_sorted[rind])

    if abs(sorted_hkl_raw[0]) < prec: 
        kl_rationalized = ratioconvert(sorted_hkl_raw[1], sorted_hkl_raw[2], digits) 
        result_sorted = np.array([0, kl_rationalized[0], kl_rationalized[1]], digits)
        return list(result_sorted[rind])

    hl_rationalized = ratioconvert(sorted_hkl_raw[0], sorted_hkl_raw[2], digits) 
    kl_rationalized = ratioconvert(sorted_hkl_raw[1], sorted_hkl_raw[2], digits) 

    hkl_rationalized_sorted = [hl_rationalized[0] * kl_rationalized[1], 
            kl_rationalized[0] * hl_rationalized[1], 
            hl_rationalized[1] * kl_rationalized[1]]
    hkl_rationalized_sorted = np.array(hkl_rationalized_sorted) / gcd(
            gcd(hkl_rationalized_sorted[0], hkl_rationalized_sorted[1]), 
            hkl_rationalized_sorted[2])

    return list(hkl_rationalized_sorted[rind])

if args.miller :
    uvw = hkl2uvw(stridx2intidx(args.miller), abc, args.digits)
    print tuple(stridx2intidx(args.miller)), '->', uvw
    #print np.dot(stridx2intidx(args.miller), reciprocal_cell) / np.dot(uvw, abc)
elif args.orientation :
    hkl = tuple(uvw2hkl(stridx2intidx(args.orientation), abc, args.digits))
    print stridx2intidx(args.orientation), '->', tuple(hkl)
    #print np.dot(stridx2intidx(args.orientation), abc) / np.dot(hkl, reciprocal_cell)
else: 
    print >>sys.stderr, 'err: please specify either -m or -i'
