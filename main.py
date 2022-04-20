from urllib.request import AbstractDigestAuthHandler
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

dP_coarse = 14.11
dP_medium = 14.25
dP_refined = 14.32

grid_coarse = 45097
grid_medium = 177092
grid_refined = 432621

vol_coarse = 0.187133
vol_medium = 0.187980
vol_refined = 0.18801

h_coarse = (vol_coarse/grid_coarse)**(1/3)
h_medium = (vol_medium/grid_medium)**(1/3)
h_refined = (vol_refined/grid_refined)**(1/3)

# Refined = 1 - Medium = 2 - Coarse = 3

r21 = h_medium/h_refined
r32 = h_coarse/h_medium
eps21 = dP_medium - dP_refined
eps32 = dP_coarse - dP_medium

s = 1*np.sign(eps32/eps21)

def qp_func(p):
    return np.log((r21**p-s)/(r32**p-s))

def p_func(p):
    return (1/np.log(r21))*(np.log(np.abs(eps32/eps21))+qp_func(p)) - p

p = optimize.newton(p_func, 1.5)

dp_ext = ((r21**p)*dP_refined - dP_medium)/((r21**p)-1)

ea21 = np.abs((dP_refined-dP_medium)/dP_refined)*100
GCI_fine_21 = 1.25*ea21/((r21**p)-1)

eext1 = np.abs((dp_ext-dP_refined)/dp_ext)*100
eext2 = np.abs((dp_ext-dP_medium)/dp_ext)*100
eext3 = np.abs((dp_ext-dP_coarse)/dp_ext)*100

# plt.plot([grid_coarse,grid_medium,grid_refined],
#          [dP_coarse,dP_medium,dP_refined])

# p = np.log((dP_coarse-dP_medium)/(dP_medium-dP_refined))/np.log(1.45)

# dP_exact = dP_medium + ((dP_refined-dP_medium)*1.45**p)/(1.45**p-1)

# eps = (dP_medium-dP_refined)/dP_refined

# GCIfine = 1.25*np.abs(eps)/(1.45**p-1)
# GCIcoarse = GCIfine*1.45**p

fig, ax = plt.subplots()
ax.set_ylabel('Perda de carga [Pa]')
ax.set_xlabel('Espaçamento da malha (m)')
ax.axis([-0.0001, 0.02, 12, 15])
# ax.xaxis.set_major_locator(plt.MultipleLocator(0.2))
# ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
# ax.yaxis.set_major_locator(plt.MultipleLocator(1.0))
# ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
ax.plot(
    [h_coarse,h_medium,h_refined],
    [dP_coarse,dP_medium,dP_refined],
    ls='-',
    color='k',
    marker='o'
)
ax.errorbar(
    [h_refined,0],
    [dP_refined,dp_ext],
    yerr=[0,GCI_fine_21/100*dp_ext/2],
    marker='o',
    color='black'
)
ax.grid()
# ax.legend(facecolor="white", framealpha=1, frameon=1)
fig.tight_layout(pad=0.15)
# fig.savefig('velocity.png',dpi=1200,format='png')
