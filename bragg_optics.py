#==============================================================================
# This contains classes for flat, cylindrical and spherical x-ray optics.
# All these inherit from the object 'BraggOptic' which contains methods for 
# setting the reflection profile and whether or not the crystal is mosaic.
# Flat optics are initialized centered at (0,0,0) symmetric about the x-axis.
# Spherical and cylindrical optics are displaced in the positive x direction by
# an amount equal to the horizontal radius of curvature.
# The active surface always faces the negative x direction.
#==============================================================================

import numpy as np


class BraggOptic(object):
    
    def __init__(self, twoD, thick=100e-6):
        
        # 2d plane spacing
        self.twoD = twoD
       
        # thickness of optic
        self.thick = thick
        
        # reflection curve
        #self.reflectProfile=None
        self.fwhm = None
        self.peakReflectivity = None
        self.setGaussianReflectProfile()
        
        # is optic mosaic or not?
        self.mosaic=False
        
    # use a generic Gaussian reflectivity profile, fwhm in radians
    def setGaussianReflectProfile(self, fwhm=1e-4, peakReflectivity=1.0):
        
        self.fwhm = fwhm
        self.peakReflectivity = peakReflectivity
             
    # return the Gaussian profile for convolution
    def getGaussianReflectProfile(self,x):
        
        m = np.median(x)
        c = self.fwhm/2.35482
        g = np.exp(-(x-m)**2/2.0/c/c)
        g =g/g.sum()
        return g
        
    # make the optic a mosaic crystal, mosaic_spread in radians, thickness in microns 
    def setMosaic(self, mosaic_spread=1e-3, thickness = 100):
        
        self.mosaic=True
        
    # load a reflection profile from outside.  typically calculated with XOP...
    def loadReflectProfile(self, fname):
        pass
  
  
class FlatBraggOptic(BraggOptic):
    
    def __init__(self,width, height,twoD = 19.95e-10):
    
        self.width = width
        self.height = height
        self.twoD = twoD
        self.aspect = self.width/self.height
        
        self.refCentroid = np.array([0.0,0.0,0.0])
        self.uMax = self.width/2.0
        self.uMin = -self.width/2.0
        self.vMax = self.height/2.0
        self.vMin = -self.height/2.0
        
        BraggOptic.__init__(self,twoD)
        
    def getPoint(self, u, v):
    
        x = u*0.0
        y = u
        z = v  
        p = np.array([x,y,z])
        return p    
        
    def getNorm(self, y, z):
        
        n = np.array([-1.0,0.0,0.0])
        return n
        
    def getRandomPointNorm(self, n):
        
        rand1 = np.random.random((n,))
        rand2 = np.random.random((n,))
        randU = self.uMin + rand1*(self.uMax-self.uMin)
        randV = self.vMin + rand2*(self.vMax-self.vMin)
        
        randPoints = self.getPoint(randU,randV)         
        randNorms = self.getNorm(randU,randV)
        
        return randPoints, randNorms        
        
    def getOrderedPointNorm(self, hstepsize=0.01, vstepsize=0.01):
        """Return points and normals from an ordered-space sampling"""

        umin = self.uMin
        vmin = self.vMin
        umax = self.uMax
        vmax = self.vMax
      
        pstorage = []
        nstorage = []

        ustep = umin
        vstep = vmin
        while ustep < umax:
            ustep = ustep + hstepsize
            vstep = vmin
            while vstep < vmax:
                newpoint = self.getPoint(ustep,vstep)
                newnorm = self.getNorm(ustep,vstep)
                
                pstorage.append(newpoint)
                nstorage.append(newnorm)
                vstep = vstep + vstepsize
        
                
        parray = np.asarray(pstorage).T
        narray = np.asarray(nstorage).T
        
        return parray,narray
        
        
        
        
        
        
#        u = np.linspace(self.uMin,self.uMax,n)
#        v = np.linspace(self.vMin,self.vMax,n)  
#         
#        p = self.getPoint(u,v)
#        n = self.getNorm(u,v)   
#        
#        return p,n    
        
    def getDispersionPointNorm(self,n):
        
        u = np.linspace(self.uMin,self.uMax,n)
        v = np.zeros(n)
        p = self.getPoint(u,v)
        n = self.getNorm(u,v)
        
        return p,n

    def getSurfaceEnergyProfile(self,distance,normAngle,order,step=5,zoffset=0.0):
        """Returns the X,Y,Z of the mesh with resolution determined by step and 
            the values of the Bragg reflected wavelength at each mesh point to be
            used as the scalar represenation in mayavi plot."""
            
        # reflection order 
        m = order    
        # convert alignment angle to rad
        normAngleRad = np.deg2rad(normAngle)     
        # meshgrid
        u = np.linspace(self.uMin,self.uMax,step)
        v = np.linspace(self.vMin,self.vMax,step)
        U,V = np.meshgrid(u,v)     
        # surface points
        P = self.getPoint(U,V)
        px = P[0,:,:].reshape(-1)
        py = P[1,:,:].reshape(-1)
        pz = P[2,:,:].reshape(-1)
        p = np.array([px,py,pz])   
        print p.shape
        # corresponding normals
        N = self.getNorm(U,V)
        nx = N[0,:,:].reshape(-1)
        ny = N[1,:,:].reshape(-1)
        nz = N[2,:,:].reshape(-1)
        n = np.array([nx,ny,nz]).T        
        # source relative to optic
        srcRelOptic = np.array([distance*np.cos(normAngleRad),
                                -distance*np.sin(normAngleRad),
                                    zoffset])                                   
        # source relative to global origin                   
        srcRelOrigin = self.refCentroid-srcRelOptic   
        # K-vectors from source
        kvecs = p.T - srcRelOrigin      
        # make k-vectors and normals on surface unit length
        kvecs /= np.sqrt((kvecs*kvecs).sum(-1))[...,np.newaxis]
        n /= np.sqrt((n*n).sum(-1))[...,np.newaxis] 
        # angle between k-vector and normal
        knAngle = np.arccos((-kvecs*n).sum(-1))[...,np.newaxis]
        grazAngle = np.pi/2.0 - knAngle       
        # Bragg reflection equation
        lFunc = lambda t: self.twoD*np.sin(t)/m    
        lams = lFunc(grazAngle).reshape(step,step)
    
        return P, lams 
        
    def getSolidAngle(self, point):
        
        sigma = []
        u = np.linspace(self.uMin,self.uMax,100)
        v = np.linspace(self.vMin,self.vMax,100)    
        for i in range(len(u)-1):
            for j in range(len(v)-1):
                
                p0 = self.getPoint(u[i],v[j])
                p1 = self.getPoint(u[i+1],v[j])
                p2 = self.getPoint(u[i],v[j+1])
                p1 = p1 - p0
                p2 = p2 - p0
                N = np.cross(p1,p2)             
                R = point - p0
                r = np.linalg.norm(R)             
                sigma.append(np.dot(R,N)/(r*r*r))
                
        Sigma = np.sum(sigma)
        return Sigma
        
    def getMeshRep(self):
        
        u = np.linspace(self.uMin,self.uMax,5)
        v = np.linspace(self.vMin,self.vMax,5)   
        U,V = np.meshgrid(u,v)
        
        p = self.getSurfPoint(U,V)
        n = self.getSurfNorm(U,V)
        
        return p,n   
        
        
""" 
Spherically-bent Bragg optic class.  The horizontal (meridional) plane 
is x-y and the vertical(sagittal) is z.  
The optic is initially defined to be symmetric about the x axis. 
"""           
class SphericalBraggOptic(FlatBraggOptic,BraggOptic):
    
    def __init__(self, width, height, radius, twoD = 19.95e-10):
        
        self.width = width
        self.height = height
        self.radius = radius
        self.twoD = twoD
        
        # define limits on horizontal and vertical angles based on 
        # optic width, height and radius of curvature     
        
        self.uMax = np.arcsin(self.width/2.0/self.radius)
        self.uMin = -self.uMax
        
        # define phi from meridional plane
        self.vMin = np.pi/2.0 - np.arcsin(self.height/2.0/self.radius)
        self.vMax = np.pi/2.0 + np.arcsin(self.height/2.0/self.radius)
        
        # reference point
        self.refCentroid  = self.getPoint(0.0,np.pi/2.0)
        #self.refCentroidNorm = self.getSurfNorm(0.0,np.pi/2.0)  

        BraggOptic.__init__(self,twoD)
        
    def getPoint(self, h, v):
        
        x = self.radius*np.cos(h)*np.sin(v)
        y = self.radius*np.sin(h)*np.sin(v)
        z = self.radius*np.cos(v)  
        p = np.array([x,y,z])
        return p
        
    def getNorm(self, h, v):
        """ Returns the normal vector at the surface point.
            This vector is NOT unit length. """
        n = self.getPoint(h,v)
        return -n
        
    def getDispersionPointNorm(self,n):
        
        u = np.linspace(self.uMin,self.uMax,n)
        v = np.ones(n)*np.pi/2.0
        p = self.getPoint(u,v)
        n = self.getNorm(u,v)
        
        return p,n
        

""" 
Cylindrically-bent Bragg optic class.  The horizontal (meridional) plane 
is x-y and the vertical(sagittal) is z.  
The optic is initially defined to be symmetric about the x axis. 
"""          
class CylindricalBraggOptic(FlatBraggOptic,BraggOptic):
    
    def __init__(self, width, height, radius, twoD = 19.95e-19):
        
        self.width = width
        self.height = height
        self.radius = radius
        self.twoD = twoD
        
        # define limits on horizontal and vertical angles based on 
        # optic width, height and radius of curvature     
        self.uMax = np.arcsin(self.width/2.0/self.radius)
        self.uMin = -self.uMax
        
        # define phi from meridional plane
        self.vMax = self.height/2.0
        self.vMin = -self.vMax
        
        # reference point
        self.refCentroid  = self.getPoint(0.0,0.0)
        #self.refCentroidNorm = self.getSurfNorm(0.0,np.pi/2.0)  
        
        BraggOptic.__init__(self,twoD)
        
    def getPoint(self, h, v):
        
        x = self.radius*np.cos(h)
        y = self.radius*np.sin(h)
        z = v 
        p = np.array([x,y,z])
        return p
        
    def getNorm(self, h, v):
        """ Returns the normal vector at the surface point.
            This vector is NOT unit length. """
        n = self.getPoint(h,v)
        
        n[2] = 0.0
       
        return -n
        
       
if __name__ == "__main__":

    S = SphericalBraggOptic(50e-3,20e-3,180.0e-3,twoD=19.95e-10)

    from mayavi import mlab
    randP,randN = S.getRandomPointNorm(100)
   
    
    pro,s = S.getSurfaceEnergyProfile(distance=117.39e-3,normAngle=49.0,order=5.0,step=100)
    mlab.mesh(pro[0],pro[1],pro[2],scalars=s)
    mlab.quiver3d(randP[0],randP[1],randP[2],randN[0],randN[1],randN[2])
    mlab.scalarbar()
    mlab.show()

    import matplotlib.pylab as plt 
    plt.plot(S.reflectProfile)
    plt.show()        
  
    
    
    
    
  
          
            