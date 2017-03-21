# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 16:21:48 2017

@author: 93137
"""

import ase.io.vasp as vp
import ase
import numpy as np

atoms = vp.read_vasp('POSCAR.Diamond_cubic.1x1x1.vasp')
abc = ase.Atoms.get_cell(atoms)
print abc[1]
print int(3.6)
a = np.array([4,2,1])
print a
b = [3,2,4]
asdsda