# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 18:53:53 2020

@author: rafael
"""

import numpy as np
import matplotlib.pylab as plt
from matplotlib import rc

#rc('text', usetex = True)

def f(k, i):
    return ((np.sqrt(k)-1)/(np.sqrt(k)+1))**i

t = np.arange(1,50,0.01)

fig = plt.figure()
ax = plt.subplot(111)

ax.loglog(t,f(2,t), label='$\kappa = 2$')
ax.loglog(t,f(10,t), 'r--', label='$\kappa = 10$')
ax.loglog(t,f(100,t), '-.', label='$\kappa = 100$')
#ax.loglog(t,f(50,t), ':', label='$\kappa = 50$')

ax.grid(which='major', linewidth=0.5)
ax.set_xlabel('$t$', size='large')
ax.set_ylabel('$\omega^t$', size='large')
plt.title(r'$\omega^{t} := \left(\frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1}\right)^{t}$',fontsize='large')
ax.legend()
plt.show