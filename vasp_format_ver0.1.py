#!/usr/bin/env python
# -*- coding:utf-8 -*-
#keywords:
#Program name: vasp_format_ver0.1.py 
#Author: Libin Wen, Created: 2017-01-18, Last Modified: 2017-01-18
#Dependence: 
#Usage: 
#


from ase.io import vasp
import sys, os, argparse

parser = argparse.ArgumentParser(description='format a poscar using VASP5 format')

parser.add_argument(dest='poscar', type=str, metavar='POSCAR', nargs=1,
        help='the POSCAR file to be translated')
parser.add_argument('-o', '--out', dest='poscar_out', type=str, metavar='POSCAR_formatted', 
        default='POSCAR_formatted', help='output formatted poscar')
parser.add_argument('-l', '--label', dest='label', type=str, metavar='label', 
        default='', help='the label in the first line of output POSCAR file')

args = parser.parse_args()
args.poscar = args.poscar[0]

poscar = vasp.read_vasp(args.poscar)
vasp.write_vasp(args.poscar_out, poscar, direct=True, vasp5=True, label=args.label)
