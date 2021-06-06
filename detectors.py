import numpy as np

class Detector(object):
    
    def __init__(self, width, height, pixel_size=20e-6):
        
        self.width=width
        self.height=height
        self.pixelSize = pixel_size
        #self.resPixPerMM = 1.0/(1e-3*pixelSize)

        
    def setFilter(self):
        pass
    
    def getMeshRep(self):
        
        u = np.linspace(-self.width/2.0,self.width/2.0)
        v = np.linspace(-self.height/2.0,self.height/2.0)
        U,V = np.meshgrid(u,v)
        y = U
        z = V
        x = np.zeros(np.shape(y))  
        shape = np.array([x,y,z])
        return shape        
    
if __name__ == "__main__":


    from mayavi import mlab
   
    d = Detector(50e-3,10e-3)
    pro = d.getMeshRep()
    mlab.mesh(pro[0],pro[1],pro[2])

    mlab.show()
