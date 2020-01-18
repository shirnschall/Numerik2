# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 18:53:53 2020

@author: rafael
"""

import numpy as np
import matplotlib.pylab as plt

def f(k, i):
    return ((np.sqrt(k)-1)/(np.sqrt(k)+1))**i

k = np.arange(1,100,1)

fig = plt.figure()
ax = plt.subplot(111)

ax.plot(k,f(k,1), label='t=1')
ax.plot(k,f(k,2), 'r--', label='t=2')
ax.plot(k,f(k,3), '-.', label='t=3')
ax.plot(k,f(k,4), ':', label='t=4')

ax.grid(which='major', linewidth=0.5)
ax.set_xlabel('K', size='large')
ax.set_ylabel('Bruch', size='large')
ax.legend()
plt.show