#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize


#------------------------------------------------------
# data (global)
t_obs=np.array([0., 100., 198., 306.])
y_obs=np.array([1. , 5., 48., 505. ])

#------------------------------------------------------

param_init=np.array([0])

def model(t, p0, tau):
    retour = p0*np.exp(+t/tau)
    print(retour)
    return retour

param_opt, pcov = optimize.curve_fit(model, t_obs, y_obs,(0,50))

print(param_opt)

fig = plt.figure()
ax = fig.add_subplot(111) #, aspect='equal', autoscale_on=True
ax.tick_params(axis='both', labelsize=20)
ax.set_xlabel(r'$t $', fontsize=20)
ax.set_ylabel(r"$y$", fontsize=20)
ax.set_yscale('log')

ax.plot(t_obs, model(t_obs, param_opt[0]), linestyle='-', linewidth=1.0, marker='o', ms=6, color='black')
ax.plot(t_obs, y_obs, linestyle="None", linewidth=1.0, marker='+', ms=7, color='red')


#plt.show()
# plt.savefig("model.pdf")















