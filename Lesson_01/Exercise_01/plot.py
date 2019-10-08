#!/usr/local/python3

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math

X = np.genfromtxt('output01.1.dat').astype(np.dtype('float32'))

x,sum_prog,err_prog=np.hsplit(X, 3)

plt.errorbar(x,sum_prog,yerr=err_prog)
plt.xlabel('#throws')
plt.ylabel('<r>-1/2')
plt.grid(True)
plt.show()
