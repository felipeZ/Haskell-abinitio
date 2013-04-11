#!/usr/bin/python

import numpy as np

N = 100
b = np.random.random_integers(-2000,2000,size=(N,N))
b_symm = (b + b.T)/2

c = np.array(b_symm,dtype=float)

np.savetxt("symm.txt",b_symm)


