import numpy as np
import matplotlib.pyplot as plt
import spectroFuncs as sf
import sys
from matplotlib.backend_bases import MouseButton


#Algorithm to renormalise line's continuum



usage = 'HOW TO USE:\n- Run python normLine.py lineFileName\
        \n Additional keywords:\
        \n- deg (default 3): degree of the polynomial fit\
        \n Ex: python normLine.py ha_res.out deg=3'
        # \n- showPlot (default 1): if 0, the plot is not shown\
        # \n Ex: python plotProfileSup.py ha_res.out save=0 showPlot=1'

if len(sys.argv)<2 or len(sys.argv)>3:
    print(usage)
    exit()
print('Double click on the plot to select continuum point, right click on one\
      of them to delete it')
deg = 3
for i, arg in enumerate(sys.argv[1:]):
    if arg[:3]=="deg":
        deg=int(arg[4:])

#Informations on star and template
sp = sf.readStellarParams('stellarParams.dat')
lp = sf.lineProf(sys.argv[1])

#Define continuum level to plot
if lp.stype=='res':
    lp.flux +=1


for i in range(sp.nbobs):
    xc = np.array([])
    yc = np.array([])
    fig = plt.figure()
    plt.title('Double click to select, right click to delete')
    plt.plot(lp.vel, lp.flux[:,i], 'k')
    plt.plot(lp.vel, np.ones(len(lp.vel)), '0.3')
    plt.ylim(lp.flux[:,i].min(), lp.flux[:,i].max())
    def onclick(event):
        if event.dblclick:
            global ix
            ix = event.xdata
            iy = event.ydata
            global xc
            xc = np.append(xc,ix)
            global yc
            yc = np.append(yc,iy)
            plt.clf()
            plt.title('Double click to select, right click to delete')
            plt.plot(lp.vel, lp.flux[:,i], 'k')
            plt.plot(lp.vel, np.ones(len(lp.vel)), '0.3')
            plt.plot(xc, yc, 'xr', ms=10)
            if len(xc)>=4:
                coeff = np.polyfit(xc,yc,deg)
                poly = np.poly1d(coeff)
                plt.plot(lp.vel, poly(lp.vel), 'r')
            plt.ylim(lp.flux[:,i].min(), lp.flux[:,i].max())
        
        if event.button is MouseButton.RIGHT:
            global jx
            jx = event.xdata
            if len(xc)!=0:
                ind = np.where(np.abs(xc-jx)==np.abs(xc-jx).min())[0]
                xc = np.delete(xc, ind)
                yc = np.delete(yc,ind)
                plt.clf()
                plt.title('Double click to select, right click to delete')
                plt.plot(lp.vel, lp.flux[:,i], 'k')
                plt.plot(lp.vel, np.ones(len(lp.vel)), '0.3')
                plt.plot(xc, yc, 'xr', ms=10)
                if len(xc)>=4:
                    coeff = np.polyfit(xc,yc,deg)
                    poly = np.poly1d(coeff)
                    plt.plot(lp.vel, poly(lp.vel), 'r')
                plt.ylim(lp.flux[:,i].min(), lp.flux[:,i].max())

        fig.canvas.draw()

    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.tight_layout()
    plt.show()
    if len(xc)!=0:
        coeff = np.polyfit(xc, yc, deg)
        poly = np.poly1d(coeff)
        # ask = str(input('Save file ? (y/n)'))
        # if ask=='y':
        lp.flux[:,i] = lp.flux[:,i]/poly(lp.vel)

if lp.stype=='res':
    lp.flux -=1
rows,column = np.shape(lp.flux)

with open('Out/Data/norm_'+sys.argv[1], 'w') as fOut:
    fOut.write('# Line '+str(lp.wlLine)+' nm file, re-normalised spectra '+sp.name+'\n')
    fOut.write('# Format: Flux-1st-day  Delta-Flux-1st-day  Flux-2nd-day  Delta-Flux-2nd-day  ect ...\n')
    fOut.write('# Velocity range and step (km/s)\n')
    fOut.write(str(lp.velMin)+'   '+str(int(lp.velMax))+'   '+str(lp.velStep)+'\n')
    fOut.write('# Julian date of observations: \n')
    for i in range(sp.nbobs):
        fOut.write(str(sp.julianDates[i])+'   ')
    fOut.write('\n')
    fOut.write('\n')
    fOut.write('\n')
    for i in range(rows):
        for j in range(column):
            fOut.write('  {:f}  '.format(lp.flux[i,j]))
            fOut.write('  {:f}  '.format(lp.dflux[i,j]))
        fOut.write('\n')
