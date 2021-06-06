# class for conical bragg optic

import numpy as np

class ConicalBraggOptic(object):
    
    def __init__(self, r1, r2, height, material=None, mode=None):
        
        self.rmin = r1
        self.rmax= r2
        self.height = height
        self.material = material
        self.mode = mode
        
        # define theta variable range
        self.theta_max= np.arcsin(self.height/2.0/self.rmin)
        
        # define a prime normal vector used to establish initial position
        # and rotation angle relative to an initial k from the source
        _prime_theta = 0.0
        _prime_rad = self.rmin + (self.rmax - self.rmin)/2.0
        
        
    def surface_point(self, r, h):
        
        x = r
        y = r*np.sin(h)
        z = r*np.cos(h)
        
        p = np.array([x,y,z])
        
        return p
        
    def surface_normal(self, r, h):
        
        xn = r
        yn = -r*np.sin(h-np.pi/2.0)
        zn = -r*np.cos(h-np.pi/2.0)
                
        N = np.array([xn,yn,zn])
        n = N/np.sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2])
       
        return n
        
    def get_parametric_rep(self):
        
        # define a vertical and horizontal range
        r = np.linspace(self.rmin,self.rmax,4)
        phi = np.linspace(-self.theta_max,self.theta_max,4)      
        
        # create 2d arrays
        H,V = np.meshgrid(r,phi)    
        
        P = self.surface_point(H,V)
        
        N = self.surface_normal(H,V)
        
        return P,N


if __name__ == "__main__":

    S = ConicalBraggOptic(100.,200.,20.0)

    from mayavi import mlab
    
 
    P,N = S.get_parametric_rep()
    
    # plot the optical surface
    mlab.mesh(P[0],P[1],P[2],color=(1,1,1))
    mlab.quiver3d(P[0],P[1],P[2],N[0],N[1],N[2])
    mlab.points3d(0.,0.,0.,color=(1,0,0))
    
    mlab.show()       
        