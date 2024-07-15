import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import spectroFuncs as sf
import sys
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


#Function computing the periodogram and its False Alarm Probability map and returning the grids
#Inputs: arrays produced by the "readSpec" function
def periodo(vel, flux, date):
    #Compute mean velocity grid from all observations
    mean_vel = vel
    #Number of velocity channel on which the periodogram is computed
    columns = len(mean_vel)
    #Compute a first periodogram and FAP to obtain the grid size
    LS = LombScargle(date, flux[0,:])
    frequency, power = LS.autopower(minimum_frequency=0.01, maximum_frequency=1)
    FAP = LS.false_alarm_probability(power, minimum_frequency=0.01, maximum_frequency=1)
    #Prepare array to fill with each periodogram power and FAP for each velocity channel
    periodoPow = np.zeros((len(frequency), columns))
    periodoFAP = np.zeros((len(frequency), columns))
    #Fill 1st column
    periodoPow[:,0] = power
    periodoFAP[:,0] = np.flip(FAP)
    #Computing all periodograms and fill the grid
    for i in range(1,columns):
        #Compute periodogram
        LS = LombScargle(date, flux[i,:])
        frequency, power = LS.autopower(minimum_frequency=0.01, maximum_frequency=1)
        periodoPow[:,i] = power
        #Compute FAP
        FAP = LS.false_alarm_probability(power, minimum_frequency=0.01, maximum_frequency=1)
        periodoFAP[:,i] = np.flip(FAP)

    #Grid x and y axis
    periodo_vel = mean_vel
    periodo_freq = frequency

    #Return X,Y,Z readable by 2D plotting function
    X, Y = np.meshgrid(periodo_vel, periodo_freq)
    Z = periodoPow
    Zfap = periodoFAP
    
    return X,Y,Z,Zfap

usage = 'HOW TO USE:\n- Run python periodo2d.py lineFileName\
        \n Additional keywords:\
        \n- save (default 0): if 1, save the plot in Out/Figures\
        \n- showPlot (default 1): if 0, the plot is not shown\
        \n Ex: python periodo2d.py ha_res.out save=0 showPlot=1'

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

#Computing the periodogram
X,Y,Z, Zfap = periodo(lp.vel, lp.flux, lp.days)

#Compute array for plotting lower panel with line+variation
meanV = lp.vel
meanF = np.mean(lp.flux, axis=1)
maxF = np.mean(lp.flux, axis=1) + np.std(lp.flux, axis=1)
minF = np.mean(lp.flux, axis=1) - np.std(lp.flux, axis=1)

#Define continuum level
if lp.stype == 'res':
    cont = 0
else:
    cont = 1

#Produce line name in LaTeX format for the figure
# label = sf.getName(line)

#Define square figure
fig = plt.figure(1, figsize=[6.5, 6.5])
#Define plots grid abd different subplots
gs = gridspec.GridSpec(6, 1)
ax1 = plt.subplot(gs[:5, 0])
ax2 = plt.subplot(gs[5, 0])
plt.subplots_adjust(hspace=.01, wspace=.01)

#Plotting periodogram
plot = ax1.pcolormesh(X,Y,Z,shading="gouraud")

ax1.axhline(y=1/sp.period, color='w', ls=':')

ax_divider = make_axes_locatable(ax1)
cax = ax_divider.append_axes("top", size="7%", pad="2%")
cbar = fig.colorbar(plot, cax=cax, orientation='horizontal')
cbar.set_label("Power", fontsize=12)
ax1.set_xlim(lp.vel.min(),lp.vel.max())
cax.xaxis.set_ticks_position("top")
cax.xaxis.set_label_position("top")
ax1.set_ylabel(r"Frequency (d$^{-1}$)", fontsize=12)
ax1.tick_params(labelsize=12)
#Define ticks
# velLim = np.array([np.round(meanV.min()/10, decimals=0)*10, np.round(meanV.max()/10, decimals=0)*10])
# xticks = np.arange(velLim[0], velLim[1], 50)

Xticks = np.array([lp.vel.min(), lp.vel.min()/2,0,lp.vel.max()/2, lp.vel.max()], dtype=int)
Xticks_label = Xticks.astype(str)

# ax1.set_xticks(xticks)
ax1.set_xticklabels('')
ax1.xaxis.set_minor_locator(MultipleLocator(50))


#Plotting mean+variance profiles below the periodogram
ax2.plot(meanV, meanF, 'k')
ax2.plot(meanV, np.ones(len(meanV))*cont, ':k')
ax2.fill_between(meanV, minF, maxF)

ax2.set_xlim(lp.vel.min(),lp.vel.max())
ax2.set_xticks(Xticks)
ax2.set_xticklabels(Xticks_label)
ax2.xaxis.set_minor_locator(MultipleLocator(50))
#Define the position to annote line name in the plot
if cont-np.mean(meanF) > 0: #absorption line
    ax2.annotate(r'$\lambda$'+str(np.round(lp.wlLine, decimals=1))+' nm', xy=(meanV.min(), meanF.min()), fontsize=12)
else: #emission line
    ax2.annotate(r'$\lambda$'+str(np.round(lp.wlLine, decimals=1))+' nm', xy=(meanV.min(), 0.9*meanF.max()), fontsize=12)

ax2.set_xlabel("v (km/s)", fontsize=12)
ax2.tick_params(labelsize=12)

#Save the figure
if save==1:
    plt.savefig('Out/Figures/'+sys.argv[1][:-4]+'Periodo.jpg', bbox_inches='tight')
if showPlot==1:
    fig = plt.figure()
    plt.title('False Alarm Probability')
    plot = plt.imshow(Zfap, interpolation='nearest', aspect='auto', extent=[X.min(), X.max(), Y.min(), Y.max()], vmax=0.3)
    plt.axhline(y=1/sp.period, color='k', ls=':')
    cbar = fig.colorbar(plot)
    plt.xlabel("v (km/s)")
    plt.ylabel(r"Frequency (d$^{-1}$)")
    plt.tight_layout()
    plt.show()
else:
    plt.close()


