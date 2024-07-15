import numpy as np
import matplotlib.pyplot as plt
import spectroFuncs as sf
import sys
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


#Overplot selected line of all dataset with different color for each rotational cycle 

usage = 'HOW TO USE:\n- Run python plotProfileSup.py lineFileName\
        \n Additional keywords:\
        \n- save (default 0): if 1, save the plot in Out/Figures\
        \n- showPlot (default 1): if 0, the plot is not shown\
        \n Ex: python plotProfileSup.py ha_res.out save=0 showPlot=1'

if len(sys.argv)<2 or len(sys.argv)>4:
    print(usage)
    exit()

save,showPlot=0,1

for i, arg in enumerate(sys.argv[1:]):
    if arg[:4]=="save":
        save=int(arg[5:])
    elif arg[:8]=="showPlot":
        showPlot=int(arg[9:])

if save!=0 and save!=1:
    print('ERROR: invalid value for "save" keyword')
    print(usage)
    exit()
if showPlot!=0 and showPlot!=1:
    print('ERROR: invalid value for "showPlot" keyword')
    print(usage)
    exit()

#Informations on star and template
sp = sf.readStellarParams('stellarParams.dat')
lp = sf.lineProf(sys.argv[1])
#Compute phase
phase = (sp.julianDates-sp.oriTime)/sp.period
phase_int = phase.astype(int)
phase = phase-phase_int


############# Line plot ############

#Define continuum level to plot
if lp.stype=='res':
    continuum = 0
    ylabel = 'Residual flux'
else:
    continuum = 1
    ylabel = 'Normalised flux'



rows, column = np.shape(lp.flux)

#Plot
plt.figure(1, figsize=[8.27, 5.85])
for i in range(column):
    plt.plot(lp.vel, lp.flux[:, i], 'C'+str(i), label=str(np.round(lp.days[i]-2450000, decimals=2)))
    plt.axhline(y=continuum, c='k', ls=':')

#########
#Plot parameters
plt.title(r'$\lambda$'+str(np.round(lp.wlLine, decimals=1))+' nm', fontsize=16)
plt.ylabel(ylabel, fontsize=16)
plt.xlabel('v (km/s)', fontsize=16)
# x ticks settings
Xticks = np.array([lp.vel.min(), lp.vel.min()/2,0,lp.vel.max()/2, lp.vel.max()], dtype=int)
Xticks_label = Xticks.astype(str)
ax = plt.gca()
ax.set_xticks(Xticks)
ax.set_xticklabels(Xticks_label)
ax.xaxis.set_minor_locator(MultipleLocator(50))
plt.tick_params(labelsize=16)
plt.legend(fontsize=13)


if save==1:
    plt.savefig('Out/Figures/'+sys.argv[1][:-4]+'Sup.pdf', bbox_inches='tight')

if showPlot==1:
    plt.axvline(x=0, ls=':', c='k')
    plt.tight_layout()
    plt.show()
else:
    plt.close()




