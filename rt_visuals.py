# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 12:06:32 2013

@author: cad
"""
import numpy as np
import matplotlib.pylab as plt
from scipy.misc import imsave
from mayavi import mlab
import os

""" Functions for Sombrero """

""" This function draws the system and traces a small subset of rays throughout """

def mayavi_visual(system, showRowland=True, scalarBar=True, units="Angstroms"):
    # find correct directory
    curr_dir = os.getcwd()
    run_dir = curr_dir+"/"+system.sysName
    
    if os.path.exists(run_dir):
        os.chdir(run_dir)
        
    # get the points for visual...there may be millions so pick a subset of these
    #sourcePoints = system.sourcePoints[::1e4]
    sourcePoints = np.loadtxt("source_points.dat")
    #opticPoints = system.opticPoints[::1e4]
    opticPoints = np.loadtxt("random_crystal_intercept_points.dat")
    #finPoints = system.detectorIntersectPoints[::1e4]
    finPoints = np.loadtxt("detector_intercept_points.dat")
    
    X = np.vstack((sourcePoints[:,0],opticPoints[:,0],finPoints[:,0])).T.reshape(1,-1)
    Y = np.vstack((sourcePoints[:,1],opticPoints[:,1],finPoints[:,1])).T.reshape(1,-1)
    Z = np.vstack((sourcePoints[:,2],opticPoints[:,2],finPoints[:,2])).T.reshape(1,-1)
   
    # make connection array
    s1 = np.arange(0,X.size-1)
    s2 = np.arange(1,X.size)
    s = np.vstack((s1,s2)).T
    sout = np.arange(2,X.size,3)
    connections = np.delete(s,sout,0)
  
    # Create the points
    src = mlab.pipeline.scalar_scatter(X,Y,Z)

    # Connect them
    src.mlab_source.dataset.lines = connections
    
    # The stripper filter cleans up connected lines
    lines = mlab.pipeline.stripper(src)
    
    # Finally, display the set of lines
    mlab.pipeline.surface(lines, colormap='Accent', line_width=1, opacity=.1)
    #mlab.outline()
    
    # get the surface rep of the optic
    order = system.diffractionOrder
    dist = system.sourceToOpticDist
    angle = system.initBraggNormAngle
    
    opticSurf,scalars = system.optic.getSurfaceEnergyProfile(dist,angle,order=order,step=20)
    
    sourceMesh = system.source.getMeshRep()
    
    for i in range(3):
        sourceMesh[i] += system.sourceVec[i]
    
    # decide units to show
    fmt = "%2.5f"
    if units == "keV":
        scalars = 12.398/scalars
        title = 'keV'
        
    elif units == "eV":
        scalars = 12398.0/scalars
        title = 'eV'
    else:
        title= "Angstroms"
        
    mlab.mesh(opticSurf[0],opticSurf[1],opticSurf[2],scalars=scalars)
    mlab.mesh(sourceMesh[0],sourceMesh[1],sourceMesh[2],color=(1,0,0))
    
    # decide if scalarbar is shown or not
    if scalarBar:
     
        mlab.scalarbar(nb_labels=2,label_fmt=fmt,title=title)
    
    # decide if rowland circle is shown
    if showRowland:
        rowlandRadius = system.optic.radius/2.0
        rowlandCenter = np.array([rowlandRadius, 0.0, 0.0])
        phiVals = np.linspace(0,2.0*np.pi,100)
        xRow = rowlandCenter[0] + rowlandRadius*np.cos(phiVals)
        yRow = rowlandCenter[1] + rowlandRadius*np.sin(phiVals)
        zRow = np.zeros(100)      
        mlab.plot3d(xRow,yRow,zRow,tube_radius=None,color=(1,0,0),opacity=0.5)
    

    # draw the detector
    
    # rotate and then translate 
    rotangle = np.arccos(np.dot(system.detectorNorm,np.array([-1.0,0,0])))
    trans = system.detectorVec
    
    detMesh = system.detector.getMeshRep()
    xnew = detMesh[0]*np.cos(rotangle) - detMesh[1]*np.sin(rotangle) + detMesh[2]*0.0 + trans[0]
    ynew = detMesh[0]*np.sin(rotangle) + detMesh[1]*np.cos(rotangle) + detMesh[2]*0.0 + trans[1]
    znew = detMesh[2] + trans[2]
    
    mlab.mesh(xnew,ynew,znew,color=(1,0,0))
    
    mlab.show()    
    
    
    
""" This function creates an image of the spectra as it would appear on the detector
    and saves both the image and the dispersion information """
def spectra_image(system, pixelSize=25.0, bgNoise=True, bgLevel=50.0):
    
    # get vectors in plane of detector
    points = system.detectorIntersectPoints - system.detectorVec
    
    yaxis = np.array([0.0,1.0,0.0])
    rotang = np.arccos(np.dot(yaxis,-system.detectorNorm))

    # rotation matrix about z
    Rz = np.matrix([ [np.cos(rotang) , -np.sin(rotang), 0.0],
                     [np.sin(rotang), np.cos(rotang), 0.0],
                      [0.0, 0.0, 1.0] ])
    
    # operate on the points with the rotation                  
    rotPoints = Rz*points.T

    # the non-zero entries become the new x,y coordinates for the image
    # x --> x , z --> y
    x = rotPoints[0] 
    y = rotPoints[2]
    
    # 1/2 detector dimensions
    hw = 0.5*system.detector.width
    hh = 0.5*system.detector.height
    
    # masked arrays to find points outside detector
    xm = np.ma.masked_outside(x,-hw,hw)
    ym = np.ma.masked_outside(y,-hh,hh)
    
    # convert points outside detector ranges to 0.0 value...
    # this places them along the edge of the detector frame
    xm-= hw
    ym-= hh
    xm=xm.filled(0.0)
    ym=ym.filled(0.0)
    
    # convert in order to use as indices for the image array
    
    #system resolution
    rez = 1.0/(1e-3*pixelSize)

    #rez = system.detector.resPixPerMM # system resolution

    xm = np.abs(np.around(xm*rez))
    ym = np.abs(np.around(ym*rez))
  
    cols = xm.astype(int)
    rows = ym.astype(int)
  
    # empty image array with dimensions of detector
    imRowMax = np.abs(np.around(rez*2.0*hh))+1
    imColMax = np.abs(np.around(rez*2.0*hw))+1
    im = np.zeros((imRowMax,imColMax),dtype=np.uint8)
        
        
    vals = system.ewi[:,2]

    # set the image pixels
    im[rows,cols] = np.around(vals*255)
    
    # save the linear dispersion information
    # I only need about 100 pts
    savepoints = np.around(system.numRays/100.0)
    
    # the xvalues converted to mm
    xvals = cols.T[::savepoints]/rez
  
    # the energies
    nrgs = system.ewi[:,0][::savepoints].reshape(xvals.shape)

    # save path is in the current working directory
    cwd = os.getcwd()
    savedata = np.hstack((nrgs,xvals))
    np.savetxt(cwd+"/detector_dispersion.txt",savedata)
    
    # add noise if True
    if bgNoise:
        randNoise = np.abs(np.around(np.random.standard_normal(size=im.shape)))*bgLevel

        im += randNoise 
 
    plt.subplot(111)
    #cmap = "Greys" for black and white
    
    # save the image 

    imsave(cwd+"/test_image.png",im)
    plt.imshow(im,cmap="Greys",extent=[0,system.detector.width,0,system.detector.height],interpolation="nearest")

    plt.show()    