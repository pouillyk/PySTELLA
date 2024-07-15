import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import matplotlib.gridspec as gridspec
import spectroFuncs as sf
import sys
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

#Function to compute the correlation matrix
def correlation_matrix(vel1, flux1, vel2, flux2):
    #Use velocity arrays's shape to define grid's size
    rows1, column1 = np.shape(flux1)  #line on the y-axis
    rows2, column2 = np.shape(flux2)  #line on the x-axis
    #Define grid array to fill
    corr_mat = np.zeros((rows1,rows2))
    #Fill the grid
    for i in range(rows1):
        for j in range(rows2):
            #Compute linear correlation coefficient between each velocity channel
            r, p_value = stats.pearsonr(flux2[j,:],flux1[i,:])
            corr_mat[i,j] = r

    return corr_mat


usage = 'HOW TO USE:\n- Run python corrMatrix.py lineFileNameX lineFileNameY\
        \n Additional keywords:\
        \n- save (default 0): if 1, save the plot in Out/Figures\
        \n- showPlot (default 1): if 0, the plot is not shown\
        \n Ex: python corrMatrix.py ha_res.out hb_res.out save=0 showPlot=1'

if len(sys.argv)<3 or len(sys.argv)>5:
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
lpX = sf.lineProf(sys.argv[1])  #line on x-axis
lpY = sf.lineProf(sys.argv[2])  #line on y-axis

#Compute the correlation matrix
print('Computing correlation matrix . . .')           
corr_mat = correlation_matrix(lpY.vel, lpY.flux, lpX.vel, lpX.flux)

print('Correlation matrix OK, plotting figure.')

#Plot the figures
#Define grid variable
X, Y = np.meshgrid(lpX.vel, lpY.vel)
Z = corr_mat

#First figure with correlation matrix and corresponding profiles on each axis (for publication)
fig = plt.figure(1, figsize=[5, 5])
gs = gridspec.GridSpec(100, 100)
ax1 = plt.subplot(gs[:85, 15:90])
ax2 = plt.subplot(gs[87:, 15:90])
ax3 = plt.subplot(gs[17:85, :13])
plt.subplots_adjust(hspace=.1, wspace=.1)
#Compute mean and variance profiles to plot on each axis
meanF2 = np.mean(lpX.flux, axis=1)
meanV2 = lpX.vel
maxF2 = meanF2 + np.std(lpX.flux, axis=1)
minF2 = meanF2 - np.std(lpX.flux, axis=1)

meanF1 = np.mean(lpY.flux, axis=1)
meanV1 = lpY.vel
maxF1 = meanF1 + np.std(lpY.flux, axis=1)
minF1 = meanF1 - np.std(lpY.flux, axis=1)

#Plot the correlation matrix
plot = ax1.pcolormesh(X, Y, Z, cmap='coolwarm', shading="gouraud", vmin=-1, vmax=1)
cbar = fig.colorbar(plot, location='top', ax=ax1, ticks=[-1,-0.5,0,0.5,1], pad=0.05)
cbar.set_label("Pearson coefficient", fontsize=15)
cbar.ax.tick_params(labelsize=15)
#Automatically define the velocity limits to set the ticks on each axis
Xticks = np.array([lpX.vel.min(), lpX.vel.min()/2,0,lpX.vel.max()/2, lpX.vel.max()], dtype=int)
Xticks_label = Xticks.astype(str)

Yticks = np.array([lpY.vel.min(), lpY.vel.min()/2,0,lpY.vel.max()/2, lpY.vel.max()], dtype=int)
Yticks_label = Yticks.astype(str)

ax1.set_xticklabels('')
ax1.xaxis.set_minor_locator(MultipleLocator(50))
ax1.set_yticklabels('')
ax1.yaxis.set_minor_locator(MultipleLocator(50))
ax1.set_xlim(meanV2.min(), meanV2.max())
ax1.set_ylim(meanV1.min(), meanV1.max())

#Plot mean and variance profile on x axis
ax2.plot(meanV2, meanF2,'k')
ax2.fill_between(meanV2, minF2, maxF2)
ax2.set_xlabel(r"v (km/s) - $\lambda$"+str(np.round(lpX.wlLine, decimals=1))+' nm', fontsize=15)
ax2.set_yticks([])

ax2.set_xticks(Xticks)
ax2.set_xticklabels(Xticks_label, fontsize=15)
ax2.xaxis.set_minor_locator(MultipleLocator(50))
ax2.set_xlim(Xticks.min(), Xticks.max())
for tick in ax2.get_xticklabels():
    tick.set_horizontalalignment('center')

#Plot mean and variance profile on y axis
ax3.plot(-meanF1, meanV1, 'k')
ax3.fill_betweenx(meanV1, -minF1, -maxF1)
ax3.set_ylabel(r"v (km/s) - $\lambda$"+str(np.round(lpY.wlLine, decimals=1))+' nm', fontsize=15)
ax3.set_xticks([])
#Define ticks

ax3.set_yticks(Yticks)
ax3.set_yticklabels(Yticks_label, fontsize=15)
ax3.yaxis.set_minor_locator(MultipleLocator(50))
ax3.set_ylim(Yticks.min(), Yticks.max())
for tick in ax3.get_yticklabels():
    tick.set_rotation(90)
    tick.set_verticalalignment('center')

plt.tight_layout()
#Save the figure
if save==1:
    plt.savefig('Out/Figures/'+sys.argv[2][:-4]+'X'+sys.argv[1][:-4]+'CM.jpg', bbox_inches='tight', transparent=True)
if showPlot==0:
    plt.close()

#Plot a second figure with only the correlation matrix (for analysis, unsaved, shown if save != 'yes at line 28)
fig = plt.figure(2, figsize=[8.5, 6.5])
fig.subplots_adjust(right=0.94, left=0.2)
gs = gridspec.GridSpec(1, 1)
ax1 = plt.subplot(gs[0, 0])

plot = ax1.imshow(Z, origin='lower',extent=[X.min(), X.max(), Y.min(), Y.max()], cmap='coolwarm')
cbar = fig.colorbar(plot)
plt.contour(X,Y,Z, levels=[0.96,0.97,0.98])
ax1.set_xlabel(r"v (km/s) - $\lambda$"+str(np.round(lpX.wlLine, decimals=1))+' nm', fontsize=16)
ax1.set_ylabel(r"v (km/s) - $\lambda$"+str(np.round(lpY.wlLine, decimals=1))+' nm', fontsize=16)
ax1.tick_params(labelsize=16)
ax1.set_xlim(lpX.vel.min(), lpX.vel.max())
ax1.set_ylim(lpY.vel.min(), lpY.vel.max())
plt.tight_layout()

if showPlot==1:
    plt.show()
else:
    plt.close()

