import numpy as np

#_______________________________________

#readStellarParams: class containing information on the studied object and correspondig template

class readStellarParams:
    def __init__(self, inStellarParams):
        with open("In/"+inStellarParams, 'r') as fIn:
            self.specPath = np.array([])
            self.julianDates = np.array([])
            self.vrad = np.array([])
            self.vradB = np.array([])
            i = 0
            for inLine in fIn:
                if inLine.strip() == '':
                    continue
                if inLine.strip()[0] != '#':
                    if i == 0:
                        self.name = inLine.split()[0]
                    
                    elif i == 1:
                        self.mass = float(inLine.split()[0])
                        self.radius = float(inLine.split()[1])
                        self.teff = float(inLine.split()[2])
                        self.period = float(inLine.split()[3])
                        self.vsini = float(inLine.split()[4])
                        self.inc = float(inLine.split()[5])
                        self.veil = float(inLine.split()[6])
                        self.dist = float(inLine.split()[7])
                    elif i == 2:
                        nb = int(inLine.split()[0])
                        self.nbobs = int(inLine.split()[0])
                    elif i >= 3 and i < 3+nb: 
                        self.specPath = np.append(self.specPath, inLine.split()[0])
                        self.julianDates = np.append(self.julianDates, float(inLine.split()[1]))
                        self.vrad = np.append(self.vrad, float(inLine.split()[2]))
                        self.vradB = np.append(self.vradB, float(inLine.split()[3]))
                    elif i == 3+nb: 
                        self.oriTime = float(inLine.split()[0]) 
                    elif i == 4+nb:
                        self.tempName = inLine.split()[0]
                    elif i == 5+nb:
                        self.tempVrad = float(inLine.split()[0])
                        self.tempVsini = float(inLine.split()[1])
                    elif i == 6+nb:
                        self.tempPath = inLine.split()[0]
                    
                    i += 1

#__________________________________

#readSpecShape: class containing information form 'infoSpec.dat' file about the shape of spectra   
     
class readSpecShape:
    def __init__(self, inSpecShape):
        with open("In/"+inSpecShape, 'r') as fIn:
            i = 0
            for inLine in fIn:
                if inLine.strip() == '':
                    continue
                if inLine.strip()[0] != '#':
                    if i == 0:
                        self.objHead = int(inLine.split()[0])
                    elif i == 1:
                        self.objUnit = inLine.split()[0]
                    elif i == 2:
                        self.objWL = int(inLine.split()[0])
                        self.objI = int(inLine.split()[1])
                        self.objSigI = int(inLine.split()[2])
                        self.objV = int(inLine.split()[3])
                        self.objSigV = int(inLine.split()[4])
                        self.objNull = int(inLine.split()[5])
                        self.objSigNull = int(inLine.split()[6])
                    elif i == 3:
                        self.tempHead = int(inLine.split()[0])
                    elif i == 4:
                        self.tempUnit = inLine.split()[0]
                    elif i == 5 : 
                        self.tempWL = int(inLine.split()[0])
                        self.tempI = int(inLine.split()[1])
                        self.tempSigI = int(inLine.split()[2])
                        self.tempV = int(inLine.split()[3])
                        self.tempSigV = int(inLine.split()[4])
                        self.tempNull = int(inLine.split()[5])
                        self.tempSigNull = int(inLine.split()[6])
                    
                    i += 1

#____________________________________

#lineProf: class containing individual lines information and flux produced with makeProfile.py

class lineProf:
    def __init__(self, lineFile):
        path = 'Out/Data/'+lineFile
        spec = np.loadtxt(path, skiprows=6)
        self.days = np.loadtxt(path, skiprows=5, max_rows=1)
        velPar = np.loadtxt(path, skiprows=3, max_rows=1)
        self.wlLine = float(np.loadtxt(path, comments=None,max_rows=1, dtype=str)[2])
        self.velMin=velPar[0]
        self.velMax=velPar[1]
        self.velStep=velPar[2]
        self.vel= np.arange(velPar[0], velPar[1]+velPar[2], velPar[2])
        self.flux = spec[:, ::2]
        self.dflux = spec[:, 1::2]
        if lineFile[-7:-4]=='res':
            self.stype='res'
        elif lineFile[-13:-4]=='unCorrVel':
            self.stype='unCorrVel'
        else:
            self.stype='classic'

#inGauss: class containing all the input parameters of the multi gaussian fit      
class inGauss:
    def __init__(self, inGaussFile):
        with open(inGaussFile, 'r') as fIn:
            i=0
            self.amp = np.array([])
            self.mean = np.array([])
            self.sigma = np.array([])
            self.boundMin = np.array([])
            self.boundMax = np.array([])
            for inLine in fIn:
                if inLine.strip() == '':
                    continue
                if inLine.strip()[0] != '#':
                    if i==0:
                        self.gaussNB = int(inLine.split()[0])
                    elif i>0 and i%3==1:
                        self.amp = np.append(self.amp, float(inLine.split()[0]))
                        self.mean = np.append(self.mean, float(inLine.split()[1]))
                        self.sigma = np.append(self.sigma, float(inLine.split()[2]))
                    elif i>0 and i%3==2:
                        self.boundMin = np.append(self.boundMin, float(inLine.split()))
                    elif i>0 and i%3==0:
                        self.boundMax = np.append(self.boundMax, float(inLine.split()))

                    i+=1