import numpy as np
import matplotlib.pyplot as plt
import spectroFuncs as sf
import sys
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

#Plot selected line of all dataset with different color for each rotational cycle in column


usage = 'HOW TO USE:\n- Run python plotProfileCol.py lineFileName\
        \n Additional keywords:\
        \n- colNB (default 2): number of columns in figure\
        \n- profNB (default 10): number of profiles per column\
        \n- save (default 0): if 1, save the plot in Out/Figures\
        \n- showPlot (default 1): if 0, the plot is not shown\
        \n Ex: python plotProfileCol.py ha_res.out colNB=2 profNB=10 save=0 showPlot=1'

if len(sys.argv)<2 or len(sys.argv)>6:
    print(usage)
    exit()

#Informations on star and template
sp = sf.readStellarParams('stellarParams.dat')
lp = sf.lineProf(sys.argv[1])

save,showPlot=0,1

colNB = int(np.ceil(sp.nbobs/5))
profNB = int(np.ceil(sp.nbobs/colNB))

for i, arg in enumerate(sys.argv[1:]):
    if arg[:4]=="save":
        save=int(arg[5:])
    elif arg[:8]=="showPlot":
        showPlot=int(arg[9:])
    elif arg[:5]=="colNB":
        colNB=int(arg[6:])
    elif arg[:6]=="profNB":
        profNB=int(arg[7:])

if save!=0 and save!=1:
    print('ERROR: invalid value for "save" keyword')
    print(usage)
    exit()
if showPlot!=0 and showPlot!=1:
    print('ERROR: invalid value for "showPlot" keyword')
    print(usage)
    exit()


#Compute phase
phase = (sp.julianDates-sp.oriTime)/sp.period
phase_int = phase.astype(int)
phase = phase-phase_int
phase_int -= phase_int.min()


############# Line plot ############

#Define continuum level to plot
if lp.stype=='res':
    continuum = 0
    ylabel = 'Residual flux'
else:
    continuum = 1
    ylabel = 'Normalised flux'

rows, column = np.shape(lp.flux)

cycle = (sp.julianDates - sp.oriTime)/ sp.period
phase = cycle - cycle.astype(int)
phaseInd = np.argsort(phase)
phase = phase[phaseInd]

lp.flux = lp.flux[:,phaseInd]
sp.vrad = sp.vrad[phaseInd]
sp.vradB = sp.vradB[phaseInd]

if (lp.flux.max()-continuum)>(continuum-lp.flux.min()):
    offset = (lp.flux.max()-continuum)
else:
    offset = (continuum-lp.flux.min())


#Plot
plt.figure(1, figsize=[colNB*3, 8])

for j in range(colNB):
    plt.subplot(1, colNB, j+1)
    plt.axvline(0, ls=':', c='k')
    Xticks = np.array([lp.vel.min(), lp.vel.min()/2,0,lp.vel.max()/2, lp.vel.max()], dtype=int)
    Xticks_label = Xticks.astype(str)
    ax = plt.gca()
    ax.set_xticks(Xticks)
    ax.set_xticklabels(Xticks_label)
    ax.xaxis.set_minor_locator(MultipleLocator(50))
    plt.xlabel('v (km/s)', fontsize=16)
    if j==0:
        plt.ylabel(r'$\lambda$'+str(np.round(lp.wlLine, decimals=1))+' nm', fontsize=16)
    plt.yticks([])

    for i in range(profNB):
        if i+profNB*j<column:
            plt.plot(lp.vel, lp.flux[:, i+profNB*j]- i*offset, c='k')
            plt.plot(lp.vel, np.ones(len(lp.vel))*continuum- i*offset, ':k')
            if sp.SB2=='y':
                if lp.stype=='res':
                    plt.plot(sp.vradB[i+profNB*j]-sp.vrad[i+profNB*j], continuum+ offset/2- i*offset, '|r')
                else:
                    plt.plot(sp.vrad[i+profNB*j], continuum+ offset/2- i*offset, '|b')
                    plt.plot(sp.vradB[i+profNB*j], continuum+ offset/2- i*offset, '|r')

            plt.annotate(str(np.round(sp.julianDates[phaseInd[i+profNB*j]]-2450000, decimals=2)), xy=[150, continuum+offset/3-i*offset])
            plt.annotate(r'$\phi$='+str(np.round(phase[i+profNB*j], decimals=2)), xy=[lp.vel.min(), continuum+offset/3-i*offset])

    plt.ylim(continuum - (i+.5)*offset, continuum+ 1.5*offset)


#Saving procedure
if save==1:
    plt.savefig('Out/Figures/'+sys.argv[1][:-4]+'Col.pdf', bbox_inches='tight')

if showPlot==1:
    plt.axvline(x=0, ls=':', c='k')
    plt.tight_layout()
    plt.show()
else:
    plt.close()
#######



