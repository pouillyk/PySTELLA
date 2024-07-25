import warnings
warnings.filterwarnings("ignore")
import numpy as np
from scipy.interpolate import interp1d
from PyAstronomy import pyasl
import spectroFuncs as sf
import sys
import matplotlib.pyplot as plt


par = sf.readStellarParams('stellarParams.dat')
iS = sf.readSpecShape('infoSpec.dat')

#Produce line file with radial velocity correction
def takeBackLine(outFile, wlLine, velMin, velMax, velStep): 
    print('Make profile file of', wlLine, 'nm line')
    print('Velocity range:', velMin, 'to', velMax, 'km/s')
    print('Velocity step:', velStep, 'km/s') 
    print('### Preparing spectra ... ###')
    velInterp = np.arange(velMin, velMax+velStep, velStep)
    
    #Prepare arrays to write file
    allFlux = np.zeros((len(velInterp), par.nbobs))
    allDflux = np.zeros((len(velInterp), par.nbobs))
        
    #Loop to fill arrays with data
    for i in range(par.nbobs):
        #Load file
        spec = np.loadtxt(par.specPath[i], skiprows=iS.objHead)
        #Check units
        if iS.objUnit == 'A':
            spec[:,iS.objWL] /= 10
        #Correct radial velocity
        spec[:,iS.objWL] = spec[:,iS.objWL] - (par.vrad[i]/3e5*spec[:,iS.objWL])
        #Convert to velocity space
        vel = 3e5*(spec[:,iS.objWL]-wlLine)/spec[:,iS.objWL]
        flux = spec[:,iS.objI]
        deltaFlux = spec[:,iS.objSigI]
        #Sort velocities, flux and uncertainties (usefull for order reconnection)
        ind = np.argsort(vel)
        vel = vel[ind]
        flux = flux[ind]
        deltaFlux = deltaFlux[ind]
        #Interpolation on new velocity grid
        f = interp1d(vel, flux)
        df = interp1d(vel, deltaFlux)
        flux = f(velInterp)
        deltaFlux = df(velInterp)
       
        #Fill the arrays
        allFlux[:,i] = flux 
        allDflux[:,i] = deltaFlux    
        
    rows, columns = np.shape(allFlux)
    print('Spectra OK')
    
    #Write the file
    print('Writting masterfile')
    with open('Out/Data/'+outFile+'.out', 'w') as fOut:
        fOut.write('# Line '+str(wlLine)+' nm file, reduced spectra of '+par.name+', corrected from radial velocity\n')
        fOut.write('# Format: Flux-1st-day  Delta-Flux-1st-day  Flux-2nd-day  Delta-Flux-2nd-day  ect ...\n')
        fOut.write('# Velocity range and step (km/s)\n')
        fOut.write(str(velMin)+'   '+str(velMax)+'   '+str(velStep)+'\n')
        fOut.write('# Julian date of observations: \n')
        for i in range(par.nbobs):
            fOut.write(str(par.julianDates[i])+'   ')
        fOut.write('\n')
        fOut.write('\n')
        fOut.write('\n')
        for i in range(rows):
            for j in range(columns):
                fOut.write('  {:f}  '.format(allFlux[i,j]))
                fOut.write('  {:f}  '.format(allDflux[i,j]))
            fOut.write('\n')

    print('File written: Out/Data/'+outFile+'.out')

#Produce residual line file
def takeBackResLine(outFile, wlLine, velMin, velMax, velStep, inFile):
    print('Make residual profile file of', wlLine, 'nm line')
    print('Velocity range:', velMin, 'to', velMax, 'km/s')
    print('Velocity step:', velStep, 'km/s') 
    print('### Preparing template ... ###')
    #Prepare the template array
    #Load template
    template = np.loadtxt(par.tempPath, skiprows=iS.tempHead)
    #Check unit
    if iS.tempUnit=="A":
        template[:,iS.tempWL]/=10
    ind = np.where((template[:,iS.tempWL]>wlLine-50)&(template[:,iS.tempWL]<wlLine+50))[0]
    template = template[ind, :]
    #Correct from radial velocity shift
    template[:,iS.tempWL] = template[:,iS.tempWL] - par.tempVrad/3e5*template[:,iS.tempWL]
    #Interpolation needed for rotational broadenning (rotBroad function - linear array of velocity asked)
    f_t = interp1d(template[:,iS.tempWL],template[:,iS.tempI])
    f_dt = interp1d(template[:,iS.tempWL],template[:,iS.tempSigI])
    if par.SB2=='y' and par.tempNameB!=par.tempPathB:
        #Load template
        templateB = np.loadtxt(par.tempPathB, skiprows=iS.tempBHead)
        #Check unit
        if iS.tempBUnit=="A":
            templateB[:,iS.tempBWL]/=10
        ind = np.where((templateB[:,iS.tempBWL]>wlLine-50)&(templateB[:,iS.tempBWL]<wlLine+50))[0]
        templateB = templateB[ind, :]
        #Correct from radial velocity shift
        templateB[:,iS.tempBWL] = templateB[:,iS.tempBWL] - par.tempVradB/3e5*templateB[:,iS.tempBWL]
        #Interpolation needed for rotational broadenning (rotBroad function - linear array of velocity asked)
        f_tB = interp1d(templateB[:,iS.tempBWL],templateB[:,iS.tempBI])
        f_dtB = interp1d(templateB[:,iS.tempBWL],template[:,iS.tempBSigI])    
    #Differential rojected velocity needed for rotational broadenning
    broad_vsini = np.sqrt(par.vsini**2-par.tempVsini**2)
    #Define linear array of velocity to use rotBroad function
    template_x = np.linspace(template[0,iS.tempWL],template[len(template)-1,iS.tempWL], len(template))
    template_y = f_t(template_x)
    template_dy = f_dt(template_x)
    #Rotational broadening from Gray function
    template_y = pyasl.rotBroad(template_x, template_y, 0.74, broad_vsini)
    #Convert to velocity space
    template_x = 3e5*(template_x - wlLine) / template_x
    #Keep only a range a bit larger than input line range
    lim = np.where((template_x>=velMin-150)&(template_x<=velMax+150))[0]
    template_x = template_x[lim]
    template_y = template_y[lim]
    template_dy = template_dy[lim]
    #Interpolation for future substraction with spectrum
    fA_t = interp1d(template_x, template_y)
    fA_dt = interp1d(template_x, template_dy)

    if par.SB2 == 'y':
        print('Template A OK')
        #Differential rojected velocity needed for rotational broadenning
        broad_vsini = np.sqrt(par.vsiniB**2-par.tempVsini**2)
        #Define linear array of velocity to use rotBroad function
        if par.tempPath==par.tempPathB:
            template_x = np.linspace(template[0,iS.tempWL],template[len(template)-1,iS.tempWL], len(template))
            template_y = f_t(template_x)
            template_dy = f_dt(template_x)
        else:
            template_x = np.linspace(templateB[0,iS.tempBWL],templateB[len(templateB)-1,iS.tempBWL], len(templateB))
            template_y = f_tB(template_x)
            template_dy = f_dtB(template_x)
        #Rotational broadening from Gray function
        template_y = pyasl.rotBroad(template_x, template_y, 0.74, broad_vsini)
        #Convert to velocity space
        template_x = 3e5*(template_x - wlLine) / template_x
        #Keep only a range a bit larger than input line range
        lim = np.where((template_x>=velMin-150)&(template_x<=velMax+150))[0]
        template_x = template_x[lim]
        template_y = template_y[lim]
        template_dy = template_dy[lim]
        #Interpolation for future substraction with spectrum
        fB_t = interp1d(template_x, template_y)
        fB_dt = interp1d(template_x, template_dy)
        print('Template B OK')
    else:
        print('Template OK')

    print('### Preparing spectra ... ###')

    velInterp = np.arange(velMin, velMax+velStep, velStep)


    #Prepare arrays to write file
    allFlux = np.zeros((len(velInterp), par.nbobs))
    allDflux = np.zeros((len(velInterp), par.nbobs))

    #Loop to fill arrays with data
    if inFile!=None:
        lp = sf.lineProf(inFile)
        if lp.stype!='classic':
            print('ERROR: "inFile" argument only available for radial\
                  velocity-corrected line profile')
            exit()
        if par.SB2=='y':
            velB = np.zeros((len(lp.vel),par.nbobs))
            for i in range(par.nbobs):
                velB[:,i] = lp.vel + par.vrad[i]-par.vradB[i]
            photoCorr = ((fA_t(lp.vel)[:,None]+par.veil)/(1+par.veil)/(1+1/par.LR) + 1/(1+par.LR)) + ((fB_t(velB)+par.veilB)/(1+par.veilB)/(1+par.LR) + par.LR/(1+par.LR)) -1
            dPhotoCorr = fA_dt(lp.vel)[:,None]/(1+1/par.LR) + fB_dt(velB)/(1+par.LR)
        else:
            photoCorr = (fA_t(lp.vel)[:,None]+par.veil)/(1+par.veil)
            dPhotoCorr = fA_dt(lp.vel)[:,None]
        print('Please check radial velocity, vsini, veiling, LR')
        plt.figure()
        plt.title(r'Please check normalisation, V$_r$, vsini, veiling, LR')
        plt.plot(lp.vel, lp.flux[:,0], 'k', label='Observation')
        plt.plot(lp.vel, photoCorr[:,0], 'r', label='Template')
        if par.SB2=='y':
            plt.plot(velInterp, ((fA_t(lp.vel)+par.veil)/(1+par.veil)/(1+1/par.LR) + 1/(1+par.LR)), c='b', label='A')
            plt.plot(velInterp, ((fB_t(velB)[:,0]+par.veilB)/(1+par.veilB)/(1+par.LR) + par.LR/(1+par.LR)), c='g', label='B')
        
        plt.legend()
        plt.show()
        ask = str(input('Everything OK (y/n)? '))
        if ask=='n':
            print('Mission aborted, check parameters')
            exit()
        allFlux=lp.flux - photoCorr
        allDflux = lp.dflux + dPhotoCorr
    
    else:
        for i in range(par.nbobs):
            #Load spectrum
            spec = np.loadtxt(par.specPath[i], skiprows=iS.objHead)
            #Check units
            if iS.objUnit=='A':
                spec[:,iS.objWL] /= 10
            #Correct radial velocity shift
            spec[:,iS.objWL] = spec[:,iS.objWL] - (par.vrad[i]/3e5*spec[:,iS.objWL])
            #Convert to velocity space
            vel = 3e5*(spec[:,iS.objWL]-wlLine)/spec[:,iS.objWL]
            flux = spec[:,iS.objI]
            deltaFlux = spec[:,iS.objSigI]
            #Sort spectrum (treat order reconnection)
            ind = np.argsort(vel)
            vel = vel[ind]
            flux = flux[ind]
            deltaFlux = deltaFlux[ind]
            f = interp1d(vel, flux)
            df = interp1d(vel, deltaFlux)
            flux = f(velInterp)
            deltaFlux = df(velInterp)
            
            #Subtracte veiling corrected photospheric template from the observations
            #if SB2, substract secondary's template
            if par.SB2=='y':
                wlB = spec[:,iS.objWL] - ((par.vradB[i]-par.vrad[i])/3e5*spec[:,iS.objWL])
                velB = 3e5*(wlB-wlLine)/wlB
                velBInterp = np.linspace(np.min(velB), np.max(velB), len(velInterp))
                velBInterp = velInterp - par.vradB[i]+par.vrad[i]
                photoCorr = ((fA_t(velInterp)+par.veil)/(1+par.veil)/(1+1/par.LR) + 1/(1+par.LR)) + ((fB_t(velBInterp)+par.veilB)/(1+par.veilB)/(1+par.LR) + par.LR/(1+par.LR)) -1
                dPhotoCorr = fA_dt(velInterp)/(1+1/par.LR) + fB_dt(velBInterp)/(1+par.LR)
            else:
                photoCorr = (fA_t(velInterp)+par.veil)/(1+par.veil)
                dPhotoCorr = fA_dt(velInterp)

            if i==0:
                print('Please check radial velocity, vsini, veiling, LR')
                # print('Current values:', par.vrad[i], par.vsini, par.veil)
                plt.figure()
                plt.title(r'Please check V$_r$, vsini, veiling, LR')
                plt.axhline(y=1, ls=':', c='0.5', alpha=0.5)
                plt.plot(velInterp, flux, 'k', label='Observation')
                plt.plot(velInterp, photoCorr, 'r', label='Template')
                if par.SB2=='y':
                    plt.plot(velInterp, ((fA_t(velInterp)+par.veil)/(1+par.veil)/(1+1/par.LR) + 1/(1+par.LR)), c='b', label='A')
                    plt.plot(velInterp, ((fB_t(velBInterp)+par.veilB)/(1+par.veilB)/(1+par.LR) + par.LR/(1+par.LR)), c='g', label='B')
                plt.legend()
                plt.show()
                ask = str(input('Everything OK (y/n)? '))
                if ask=='n':
                    print('Mission aborted, check parameters')
                    exit()
            flux = flux - photoCorr 
            deltaFlux = deltaFlux + dPhotoCorr
        
            #Fill the arrays
            allFlux[:,i] = flux 
            allDflux[:,i] = deltaFlux  

    rows, columns = np.shape(allFlux)
    print('Spectra OK')
    print('Writting masterfile')

    with open('Out/Data/'+outFile+'_res.out', 'w') as fOut:
        fOut.write('# Line '+str(wlLine)+' nm file, residual spectra of '+par.name+', photospheric template '+par.tempName+', corrected from radial velocity\n')
        fOut.write('# Format: Flux-1st-day  Delta-Flux-1st-day  Flux-2nd-day  Delta-Flux-2nd-day  ect ...\n')
        fOut.write('# Velocity range and step (km/s)\n')
        fOut.write(str(velMin)+'   '+str(velMax)+'   '+str(velStep)+'\n')
        fOut.write('# Julian date of observations: \n')
        for i in range(par.nbobs):
            fOut.write(str(par.julianDates[i])+'   ')
        fOut.write('\n')
        fOut.write('\n')
        fOut.write('\n')
        for i in range(rows):
            for j in range(columns):
                fOut.write('  {:f}  '.format(allFlux[i,j]))
                fOut.write('  {:f}  '.format(allDflux[i,j]))
            fOut.write('\n')
    
    print('File written: Out/Data/'+outFile+'_res.out')


#Produce line file without radial velocity correction
def takeBackLineUncorrVel(outFile, wlLine, velMin, velMax, velStep):
    print('Make velocity-uncorrected profile file of', wlLine, 'nm line')
    print('Velocity range:', velMin, 'to', velMax, 'km/s')
    print('Velocity step:', velStep, 'km/s')     
    print('### Preparing spectra ... ###') 
    velInterp = np.arange(velMin, velMax+velStep, velStep)
    
    #Prepare arrays to write file
    allFlux = np.zeros((len(velInterp), par.nbobs))
    allDflux = np.zeros((len(velInterp), par.nbobs))
        
    #Loop to fill arrays with data
    for i in range(par.nbobs):
        #Load file
        spec = np.loadtxt(par.specPath[i], skiprows=iS.objHead)
        #Check units
        if iS.objUnit == 'A':
            spec[:,iS.objWL] /= 10
        #Convert to velocity space
        vel = 3e5*(spec[:,iS.objWL]-wlLine)/spec[:,iS.objWL]
        flux = spec[:,iS.objI]
        deltaFlux = spec[:,iS.objSigI]
        #Sort velocities, flux and uncertainties (usefull for order reconnection)
        ind = np.argsort(vel)
        vel = vel[ind]
        flux = flux[ind]
        deltaFlux = deltaFlux[ind]
        #Interpolation on new velocity grid
        f = interp1d(vel, flux)
        df = interp1d(vel, deltaFlux)
        flux = f(velInterp)
        deltaFlux = df(velInterp)
       
        #Fill the arrays
        allFlux[:,i] = flux 
        allDflux[:,i] = deltaFlux    
        
    rows, columns = np.shape(allFlux)
    print('Spectra OK')
    
    #Write the file
    print('Writting masterfile')
    with open('Out/Data/'+outFile+'_unCorrVel.out', 'w') as fOut:
        fOut.write('# Line '+str(wlLine)+' nm file, reduced spectra of '+par.name+', not corrected from radial velocity\n')
        fOut.write('# Format: Flux-1st-day  Delta-Flux-1st-day  Flux-2nd-day  Delta-Flux-2nd-day  ect ...\n')
        fOut.write('# Velocity range and step (km/s)\n')
        fOut.write(str(velMin)+'   '+str(velMax)+'   '+str(velStep)+'\n')
        fOut.write('# Julian date of observations: \n')
        for i in range(par.nbobs):
            fOut.write(str(par.julianDates[i])+'   ')
        fOut.write('\n')
        fOut.write('\n')
        fOut.write('\n')
        for i in range(rows):
            for j in range(columns):
                fOut.write('  {:f}  '.format(allFlux[i,j]))
                fOut.write('  {:f}  '.format(allDflux[i,j]))
            fOut.write('\n')

    print('File written: Out/Data/'+outFile+'_unCorrVel.out')

#####################################################               
        
usage = 'HOW TO USE:\n- Run python makeProfile.py line_wavelength(nm)\
    outfilename(without extension)\
    \n Additional keywords:\
    \n- velMin/velMax (default: -500/500): min/max velocity of the profile (km/s)\
    \n- velStep (default: 1.5): step of the velocity grid \
    \n- res (defaut:0): if 1, the photospheric template is removed from the profile.\
    "_res" is automatically added to the file name.\
    \n- corrVel (default:1): if 0 the radial velocity is not corrected,\
    "_unCorrVel" is automatically added to the file name.\
    \n if res=1, the velocity is corrected anyway. If SB2 the radial velocity\
    corrected is the value set for the primary\
    \n- (Only for residual profile) inFile (default:None) : use a previously\
    computed line profile file using radial velocity correction (a bit quicker and usefull if a\
    re-normalisation was needed)\
    \n Ex: python makeProfile.py 656.279 ha velMin=-500 velMax=500 velStep=1.5 res=0 corrVel=1 inFile=None'

if len(sys.argv)<3 or len(sys.argv)>9:
    print('ERROR: wrong keywords number')
    print(usage)
    exit()

velMin,velMax,velStep,res,corrVel,inFile = -500,500,1.5,0,1,None

for i, arg in enumerate(sys.argv[1:]):
    if arg[:6]=="velMin":
        velMin=float(arg[7:])
    elif arg[:6]=="velMax":
        velMax=float(arg[7:])
    elif arg[:7]=="velStep":
        velStep=float(arg[8:])
    elif arg[:3]=="res":
        res=int(arg[4:])
    elif arg[:7]=="corrVel":
        corrVel=int(arg[8:])
    elif arg[:6]=="inFile":
        inFile=str(arg[7:])

if res!=0 and res!=1:
    print('ERROR: invalid value for "res" keyword')
    print(usage)
    exit()
if corrVel!=0 and corrVel!=1:
    print('ERROR: invalid value for "corrVel" keyword')
    print(usage)
    exit()
if corrVel==0 and res==1:
    print('WARNING: radial velocity correction is mandatory to\
    compute residual profiles. corrVel was set to 1.')
    corrVel=1
if inFile!=None and res!=1:
    print('ERROR: "inFile" argument only available for residual \
          profile option')
    exit()

if res==0 and corrVel==1:
    takeBackLine(sys.argv[2], float(sys.argv[1]), velMin, velMax, velStep)
if res==0 and corrVel==0:
    takeBackLineUncorrVel(sys.argv[2], float(sys.argv[1]), velMin, velMax, velStep)
if res==1:
    takeBackResLine(sys.argv[2], float(sys.argv[1]), velMin, velMax, velStep, inFile)