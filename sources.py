import numpy as np
import matplotlib.pylab as plt


class SphericalSource(object):
    
    def __init__(self, diameter=1e-6):
        self.diameter=diameter
        self.radius = self.diameter/2.0
        self.spectrum = None
        self.spectrum_desc = None
        self.minE = None
        self.maxE = None
        
        
    def loadSpectrum(self, filepath):
        """ This should return an Nx2 array of energies (col0) and intensities (col1) """
        if filepath[:-3] == "ppd":
            print 'true'
        else:
            pass
        x,y = np.loadtxt(filepath,unpack=True)
        #print sfile.shape
        #spect = np.zeros((sfile.shape[0],2))
        #spect[:,0] = sfile[:,0]
        #spect[:,1] = 12398.0/sfile[:,0]  # convert eV to Angstroms
        #spect[:,1] = sfile[:,1]
     
        self.spectrum=np.array([x,y])

        self.spectrum_desc = "Loaded from file...%s" % (filepath)
        
    # define a gaussian spectral feature. Energy and fwhm in eV. 
    # returns two column array of energy and intensity
    def setGaussianSpectrum(self, energy, fwhm):
    
        c = fwhm/2.35482
        b = energy
        
        x = np.arange(b-b/2,b+b/2,.01)
        
        spect = np.exp( -(x-b)*(x-b)/(2*c*c))
        self.spectrum = np.array([x,spect])
        self.minE = x.min()
        self.maxE = x.max()
    
    def _gauss(self,x,a,b,c):
        
        g = a*np.exp( -(x-b)*(x-b)/(2*c*c))
        return g
        

    def setMultiGaussianSpectrum(self, energy_list, intensity_list, fwhm_list):
        
        self.minE = np.array(energy_list).min()
        self.maxE = np.array(energy_list).max()
        
        x = np.arange(self.minE - self.minE*.1,self.maxE + self.maxE*.1,.1)
        
        totspect = np.zeros(x.shape)
        
        for i in range(len(energy_list)):
            
            s = self._gauss(x,intensity_list[i],energy_list[i],fwhm_list[i])
            totspect = totspect + s
        
        self.spectrum = np.array([x,totspect])
        
    
    def setBandSpectrum(self, min_energy, max_energy):
        
        pass
    
 
    def getRandomVolumePoints(self,N):
        
        rand1 = np.random.random((N,))
        rand2 = np.random.random((N,))
        rand3 = np.random.random((N,))
        
        randTheta = 2.0*np.pi*(1.0-rand1)
        randPhi = np.pi*(1.0-rand2)
        randRad = self.radius*(1.0-rand3)
        
        x = randRad*np.cos(randTheta)*np.sin(randPhi)
        y = randRad*np.sin(randTheta)*np.sin(randPhi)
        z = randRad*np.cos(randPhi) 
        
        p = np.array([x,y,z]).T
        return p
    
    def getRandomEnergy(self, N):

        if self.spectrum is not None:
            ranE = (self.maxE - self.minE)*np.random.random(N) + self.minE
            return ranE
        else:
            return np.zeros(N)
     
        
    def getMeshRep(self):
        
        h = np.linspace(0,2.0*np.pi)
        v = np.linspace(0,np.pi)
        
        H,V = np.meshgrid(h,v)
        radius = self.radius
        
        x = np.cos(H)*np.sin(V)
        y = np.sin(H)*np.sin(V)
        z = np.cos(V) 
        
        shape = np.array([x,y,z])
        return shape
        
if __name__ == "__main__":
    
    src = SphericalSource(diameter=1e-6)
    src.setGaussianSpectrum(energy=4750.,fwhm=60.)
    plt.plot(src.spectrum[0],src.spectrum[1])
    plt.show()
    #print Src.getRandVolPoints(5)
        
        