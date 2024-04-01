#!/usr/bin/env python3
########################

# In this code, for radiative elastic events (e + p -> e + p + gamma),
# plot the true Q2 range as a function of electron x. Each curve represents
# a fixed value of electron Q2. Note that for these radiative elastic events,
# the true x is equal to 1.

from cProfile import label
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#Output pdf files
pp1 = PdfPages('Q2_range.pdf')

#Set parameters 
Mp = 0.93827

#Gamma-squared
def Gammasq(Q2,x):
    return 4*Mp*Mp*x*x/Q2

#Minimum true Q2
def Q2_min(Q2,x,xtrue):

    xr = x/xtrue
    gamma2 = Gammasq(Q2,x)

    Q2_true = Q2 * ( 2. * (1.-xr) * (1.-np.sqrt(1.+gamma2)) + gamma2 )
    Q2_true = Q2_true / ( gamma2 + 4.*xr*(1.-xr) )

    return Q2_true

#Maximum true Q2
def Q2_max(Q2,x,xtrue):

    xr = x/xtrue
    gamma2 = Gammasq(Q2,x)

    Q2_true = Q2 * ( 2. * (1.-xr) * (1.+np.sqrt(1.+gamma2)) + gamma2 )
    Q2_true = Q2_true / ( gamma2 + 4.*xr*(1.-xr) )

    return Q2_true

#Set true x
x_true = 1.0

#Set electron values
Q2_electron_range = np.array([1.0,10.0,40.0])

#Plot1
fig,ax = plt.subplots()
plt.subplots_adjust(left=0.125, bottom=0.125, right=0.875, top=0.875)

ax.set_xlabel(r'$x_{elec.}$',fontsize=18)
ax.set_ylabel(r'Limits of $Q^{2}_{True}~[GeV^{2}]$',fontsize=18)
ax.set_title(r'$e + p \rightarrow e + p + \gamma_{Rad.}$',fontsize=18)

ax.tick_params(top=True, right=False)
ax.tick_params(axis="x",direction="in")
ax.tick_params(labelsize=14)
ax.spines["left"].set_linewidth(2)
ax.spines["right"].set_linewidth(2)
ax.spines["bottom"].set_linewidth(2)
ax.spines["top"].set_linewidth(2)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([6e-5,2])
ax.set_ylim([1e-9,1e6])

#Set colors
mycolors = ['blue','green','red']

#Loop over true x and plot curve (minimum true Q2 vs. electron x) for each
for index,Q2_elec in enumerate(Q2_electron_range):
    x_elec = np.logspace(-4,np.log10(x_true),1000)
    Q2_true_min = Q2_min(Q2_elec,x_elec,x_true)
    Q2_true_max = Q2_max(Q2_elec,x_elec,x_true)

    leg_label = (r'$Q^{2}_{elec} = %.0f ~GeV^{2}$' %(Q2_elec))
    plt.plot(x_elec,Q2_true_min,label=leg_label,color=mycolors[index])
    plt.plot(x_elec,Q2_true_max,color=mycolors[index])

plt.legend(loc="lower right", prop={'size': 15})

#Figure Output
fig = plt.gcf()
fig.set_size_inches(18.5/2, 10.5/2)
fig.savefig(pp1, format='pdf')

#Finished plots
#Close output pdfs
pp1.close()

# end of script