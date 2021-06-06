# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 16:43:58 2013

@author: cad
"""
import numpy as np
from scipy.interpolate import interp1d
import os
import time
from bragg_optics import *

                                      
#==============================================================================
# Improved trace functions
#==============================================================================
def grid_trace_system(system, hpoints = 1000, vpoints = 100):

    nrays = hpoints*vpoints
    system.setTraceType("ordered")
    system.setNumRays(nrays)
    
    u = np.linspace(system.optic.uMin, system.optic.uMax, hpoints)
    v = np.linspace(system.optic.vMin, system.optic.vMax, vpoints)
    U,V = np.meshgrid(u,v)
    
    # save the array shape for easy mlab plotting later...
    system.setMeshGridShape(U.shape)
    
    # Points and normals from optic class
    P = system.optic.getPoint(U,V)
    N = system.optic.getNorm(U,V)
    
    # reshape for future use 
    px = P[0,:,:].reshape(-1)
    py = P[1,:,:].reshape(-1)
    pz = P[2,:,:].reshape(-1)
    Points = np.array([px,py,pz]) 
    nx = N[0,:,:].reshape(-1)
    ny = N[1,:,:].reshape(-1)
    nz = N[2,:,:].reshape(-1)
    Norms = np.array([nx,ny,nz])
   
    # get source points and energy list
    system.setNumRays(hpoints*vpoints)
    randSourcePoints = system.source.getRandomVolumePoints(hpoints*vpoints)
        
    # get location of source    
    sourceVec = system.sourceVec
    
    # get optic info and diffraction order
    twoD = system.optic.twoD
    order = system.diffractionOrder

    
    # create empty arrays to hold the source points, optic normals, kvectors from source, and
    # crystal reflection points
    norms = np.zeros((nrays,3),dtype=np.float32)
    Kvecs = np.zeros((nrays,3),dtype=np.float32)
    opticPoints = np.zeros((nrays,3),dtype=np.float32)
    sourcePoints = np.zeros((nrays,3),dtype=np.float32)
    
    # fill empty arrays
    for i in range(3):
        norms[:,i]= Norms[i].reshape((-1,1)).squeeze()
        opticPoints[:,i] = Points[i].reshape((-1,1)).squeeze()
        #Kvecs[:,i] = opticPoints[:,i] - sourceVec[i]    
        sourcePoints[:,i] = sourceVec[i] + randSourcePoints[:,i]  
        Kvecs[:,i] = opticPoints[:,i] - sourcePoints[:,i]
        
    # ensure optic normals and kvecs are unit length    
    norms /= np.sqrt((norms*norms).sum(-1))[...,np.newaxis]
    kvecs = Kvecs/np.sqrt( (Kvecs*Kvecs).sum(-1))[...,np.newaxis]

    # the reflected k vector and unit normal reflected k vector
    Rvecs = Kvecs - 2.0*((norms*Kvecs).sum(-1)[...,np.newaxis])*norms
    rvecs = Rvecs/np.sqrt( (Rvecs*Rvecs).sum(-1))[...,np.newaxis]
 
    # get all normal incidence angles
    normAngles = np.arccos((-kvecs*norms).sum(-1))[...,np.newaxis]    
    
    # convert to grazing angles
    grazAngles = np.pi/2.0 - normAngles
    
    # inline function to determine wavelength from angle
    wavelenFunc = lambda t: twoD*np.cos(t)/order   
    lams = wavelenFunc(normAngles)
    enrgs = 12398.0/lams/1e10   

    
    # setup file path for saving data
    curr_dir = os.getcwd()
    run_dir = curr_dir+"/"+system.sysName
    
    # check if directory exists, if not then make it
    if not os.path.exists(run_dir):
        print "Directory does not exist...creating new\n"
        os.makedirs(run_dir)
        
    # save general setup information
    print "Saving grid trace optical setup information...\n"
    fsetup = open(run_dir+"/grid_trace_optical_setup_info.txt", "w")
    fsetup.write("Sombrero run '"+system.sysName+"' on "+time.strftime("%d/%m/%Y")+" at "+time.strftime("%H:%M")+"\n\n")
    fsetup.write("##-##-"*10+"\n\n")
    
    fsetup.write("Crystal dimensions (mm): %2.1f x %2.1f \n\n" % (system.optic.width*1e3,system.optic.height*1e3))
      
    optic_type = type(system.optic)            
    if optic_type is SphericalBraggOptic:
        fsetup.write("Crystal is spherically-bent...%2.2fmm radius of curvature\n\n" % (system.optic.radius*1e3))
    elif optic_type is CylindricalBraggOptic:
        fsetup.write("Crystal is cylindrically-bent...%2.2fmm radius of curvature\n\n" % (system.optic.radius*1e3))
    elif optic_type is FlatBraggOptic:
        fsetup.write("Crystal is flat...\n\n")
    else:
        fsetup.write("Crystal shape not specified??\n\n")
 
    fsetup.write("Normal incidence alignment angle (degrees): %2.2f \n\n" % (system.initBraggNormAngle))
    fsetup.write("##-##-"*10+"\n\n")
    
    fsetup.write("Source size (microns): %2.2f \n\n" % (system.source.diameter*1e6))
    
    fsetup.write("Source spectrum: %s \n\n" % (system.source.spectrum_desc))
    
    if system.sourceOnRowland:
        fsetup.write("SOURCE ON ROWLAND CIRCLE!  Source to optic distance (cm): %2.2f \n\n" % (system.sourceToOpticDist*1e2))
    else:
        fsetup.write("Source to optic distance (cm): %2.2f \n\n" % (system.sourceToOpticDist*1e2))
    fsetup.write("##-##-"*10+"\n\n")
    
    fsetup.write("Detector dimensions (mm): %2.1f x %2.1f \n\n" % (system.detector.width*1e3, system.detector.height*1e3))
    if system.detectorOnRowland:
        fsetup.write("DETECTOR ON ROWLAND CIRCLE!  Optic to detector distance (cm): %2.2f \n\n" % (system.opticToDetectorDist*1e2))
    else:
        fsetup.write("Optic to detector distance (cm): %2.2f \n\n" % (system.opticToDetectorDist*1e2))
        
    fsetup.write("System was traced with %d %s rays" % (nrays, system.traceType))
    fsetup.close()
    
    # source location points
    print "Saving source data...\n"
    np.savetxt(run_dir+"/source_points.dat", sourcePoints)

    # crystal reflection points
    print "Saving crystal intercept data...\n"
    np.savetxt(run_dir+"/grid_crystal_intercept_points.dat", opticPoints)
    
    # energy, wavelength, angle, and k-vector data
    print "Saving energy, wavelength, angle and k-vector data...\n"
    np.savetxt(run_dir+"/grid_energies.dat",enrgs)
    np.savetxt(run_dir+"/grid_wavelengths.dat",lams)
    np.savetxt(run_dir+"/grid_normal_incidence_angles.dat",normAngles)
    np.savetxt(run_dir+"/grid_grazing_incidence_angles.dat",grazAngles)
    np.savetxt(run_dir+"/grid_incident_kvectors.dat",kvecs)
    np.savetxt(run_dir+"/grid_reflected_kvectors.dat",rvecs)
    
def random_trace_system(system, nrays=1e5):
    nrays = nrays
    system.setTraceType("random")
    system.setNumRays(nrays)
    
    # get random points and normals on surface of optic
    Points,Norms = system.optic.getRandomPointNorm(nrays)
    # get source points 
    system.setNumRays(nrays)
    randSourcePoints = system.source.getRandomVolumePoints(nrays)
        
    # get location of source    
    sourceVec = system.sourceVec
    
    # get optic info and diffraction order
    twoD = system.optic.twoD
    order = system.diffractionOrder
 
    # create empty arrays to hold the source points, optic normals, kvectors from source, and
    # crystal reflection points
    norms = np.zeros((nrays,3),dtype=np.float32)
    Kvecs = np.zeros((nrays,3),dtype=np.float32)
    opticPoints = np.zeros((nrays,3),dtype=np.float32)
    sourcePoints = np.zeros((nrays,3),dtype=np.float32)
    
    # fill empty arrays
    for i in range(3):
        norms[:,i]= Norms[i].reshape((-1,1)).squeeze()
        opticPoints[:,i] = Points[i].reshape((-1,1)).squeeze()
        #Kvecs[:,i] = opticPoints[:,i] - sourceVec[i]    
        sourcePoints[:,i] = sourceVec[i] + randSourcePoints[:,i]  
        Kvecs[:,i] = opticPoints[:,i] - sourcePoints[:,i]
        
    # ensure optic normals and kvecs are unit length    
    norms /= np.sqrt((norms*norms).sum(-1))[...,np.newaxis]
    kvecs = Kvecs/np.sqrt( (Kvecs*Kvecs).sum(-1))[...,np.newaxis]

    # the reflected k vector and unit normal reflected k vector
    Rvecs = Kvecs - 2.0*((norms*Kvecs).sum(-1)[...,np.newaxis])*norms
    rvecs = Rvecs/np.sqrt( (Rvecs*Rvecs).sum(-1))[...,np.newaxis]
 

    # get all normal incidence angles
    normAngles = np.arccos((-kvecs*norms).sum(-1))[...,np.newaxis]    
    
    # convert to grazing angles
    grazAngles = np.pi/2.0 - normAngles
    
    # inline function to determine wavelength from angle
    wavelenFunc = lambda t: twoD*np.cos(t)/order   
    lams = wavelenFunc(normAngles)
    enrgs = 12398.0/lams/1e10   

    
    # setup file path for saving data
    curr_dir = os.getcwd()
    run_dir = curr_dir+"/"+system.sysName
    
    # check if directory exists, if not then make it
    if not os.path.exists(run_dir):
        print "Directory does not exist...creating new\n"
        os.makedirs(run_dir)
        
    # save general setup information
    print "Saving random trace optical setup information...\n"
    fsetup = open(run_dir+"/random_trace_optical_setup_info.txt", "w")
    fsetup.write("Sombrero run '"+system.sysName+"' on "+time.strftime("%d/%m/%Y")+" at "+time.strftime("%H:%M")+"\n\n")
    fsetup.write("##-##-"*10+"\n\n")
    
    fsetup.write("Crystal dimensions (mm): %2.1f x %2.1f \n\n" % (system.optic.width*1e3,system.optic.height*1e3))
      
    optic_type = type(system.optic)            
    if optic_type is SphericalBraggOptic:
        fsetup.write("Crystal is spherically-bent...%2.2fmm radius of curvature\n\n" % (system.optic.radius*1e3))
    elif optic_type is CylindricalBraggOptic:
        fsetup.write("Crystal is cylindrically-bent...%2.2fmm radius of curvature\n\n" % (system.optic.radius*1e3))
    elif optic_type is FlatBraggOptic:
        fsetup.write("Crystal is flat...\n\n")
    else:
        fsetup.write("Crystal shape not specified??\n\n")
 
    fsetup.write("Normal incidence alignment angle (degrees): %2.2f \n\n" % (system.initBraggNormAngle))
    fsetup.write("##-##-"*10+"\n\n")
    
    fsetup.write("Source size (microns): %2.2f \n\n" % (system.source.diameter*1e6))
    
    fsetup.write("Source spectrum: %s \n\n" % (system.source.spectrum_desc))
    
    if system.sourceOnRowland:
        fsetup.write("SOURCE ON ROWLAND CIRCLE!  Source to optic distance (cm): %2.2f \n\n" % (system.sourceToOpticDist*1e2))
    else:
        fsetup.write("Source to optic distance (cm): %2.2f \n\n" % (system.sourceToOpticDist*1e2))
    fsetup.write("##-##-"*10+"\n\n")
    
    fsetup.write("Detector dimensions (mm): %2.1f x %2.1f \n\n" % (system.detector.width*1e3, system.detector.height*1e3))
    if system.detectorOnRowland:
        fsetup.write("DETECTOR ON ROWLAND CIRCLE!  Optic to detector distance (cm): %2.2f \n\n" % (system.opticToDetectorDist*1e2))
    else:
        fsetup.write("Optic to detector distance (cm): %2.2f \n\n" % (system.opticToDetectorDist*1e2))
        
    fsetup.write("System was traced with %d %s rays" % (nrays, system.traceType))
    fsetup.close()
    
    # source location points
    print "Saving source data...\n"
    np.savetxt(run_dir+"/source_points.dat", sourcePoints)

    # crystal reflection points
    print "Saving crystal intercept data...\n"
    np.savetxt(run_dir+"/random_crystal_intercept_points.dat", opticPoints)
    
    # energy, wavelength, angle, and k-vector data
    print "Saving energy, wavelength, angle and k-vector data...\n"
    np.savetxt(run_dir+"/random_energies.dat",enrgs)
    np.savetxt(run_dir+"/random_wavelengths.dat",lams)
    np.savetxt(run_dir+"/random_normal_incidence_angles.dat",normAngles)
    np.savetxt(run_dir+"/random_grazing_incidence_angles.dat",grazAngles)
    np.savetxt(run_dir+"/random_incident_kvectors.dat",kvecs)
    np.savetxt(run_dir+"/random_reflected_kvectors.dat",rvecs)    
    
#==============================================================================
#  Optic reflection details.  This is best used with the "ordered" trace results.
#==============================================================================
def grid_optic_intercept_info(system):
    
    # find correct directory
    curr_dir = os.getcwd()
    run_dir = curr_dir+"/"+system.sysName
    
    if os.path.exists(run_dir):
        os.chdir(run_dir)
        
    spectrum_func = interp1d(system.source.spectrum[0],system.source.spectrum[1])  
    # load the energy list and modify the spectrum according to the distribution   
    ref_nrgs = np.loadtxt(run_dir+"/grid_energies.dat") 
    init_intensities = np.ones(ref_nrgs.shape)
    src_intensities = spectrum_func(ref_nrgs)
    fin_intensities = src_intensities*init_intensities
                
        
    # load optic reflection points
    optic_points = np.loadtxt(run_dir+"/grid_crystal_intercept_points.dat")
      
    x = optic_points[:,0].reshape(system.meshShape)
    y = optic_points[:,1].reshape(system.meshShape)
    z = optic_points[:,2].reshape(system.meshShape)
    s = fin_intensities.reshape(system.meshShape)

    
    from mayavi import mlab
    mlab.mesh(x,y,z,scalars=s)
    mlab.view(azimuth=180, elevation=90,distance=.09)
    #mlab.scalarbar(nb_labels=5)
    mlab.show()
    
#==============================================================================
#  Detector image details    
#==============================================================================
def detector_intercept_info(system, use_random=True, dispersion_curve=False):
    # consider using mayavi.mlab.imshow for this image...
    # compute intercept points and energies, dispersion curves, etc...

    # find correct directory
    curr_dir = os.getcwd()
    run_dir = curr_dir+"/"+system.sysName
    
    if os.path.exists(run_dir) and curr_dir is not run_dir:
        print "Changing to directory with ray trace data...\n"
        os.chdir(run_dir)   
        
    # get optic reflection points and energies
    if use_random:
        
        optic_points = np.loadtxt("random_crystal_intercept_points.dat")
        rvecs = np.loadtxt("random_reflected_kvectors.dat")
        
        # convolve with rocking curve
        rc = system.optic.getGaussianReflectProfile(system.source.spectrum[0])
        print "Convolving spectrum with rocking curve...\n"
        conv_spect = np.convolve(system.source.spectrum[1],rc,mode='same')
        print "Convolution done...\n"
        conv_spect = conv_spect/conv_spect.max()
        print "Interpolating spectrum...\n"
        spectrum_func = interp1d(system.source.spectrum[0],conv_spect,bounds_error=False,fill_value=0.0) 
        print "Interpolation done...\n"
        #spectrum_func = interp1d(system.source.spectrum[0],system.source.spectrum[1]) 
        
        # load the energy list and modify the spectrum according to the distribution   
        ref_nrgs = np.loadtxt("random_energies.dat") 
        init_intensities = np.ones(ref_nrgs.shape)
        src_intensities = spectrum_func(ref_nrgs)
        fin_intensities = src_intensities*init_intensities
      
    else:
        
        optic_points = np.loadtxt("grid_crystal_intercept_points.dat")
        rvecs = np.loadtxt("grid_reflected_kvectors.dat")
        spectrum_func = interp1d(system.source.spectrum[0],system.source.spectrum[1]) 
        
        # load the energy list and modify the spectrum according to the distribution   
        ref_nrgs = np.loadtxt("grid_energies.dat") 
        init_intensities = np.ones(ref_nrgs.shape)
        src_intensities = spectrum_func(ref_nrgs)
        fin_intensities = src_intensities*init_intensities
        
    
    # detector location
    p = system.detectorVec
    
    # detector normal
    n = system.detectorNorm
    pR = p-optic_points
    
    # intersection distance
    d = (pR*n).sum(-1)[...,np.newaxis]/(rvecs*n).sum(-1)[...,np.newaxis]
    
    # define the end points
    finPoints = optic_points+rvecs*d
    
    # the translated intersect points
    points = finPoints - p
    
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
    rez = 1.0/(system.detector.pixelSize)

    #rez = system.detector.resPixPerMM # system resolution

    xm = np.abs(np.around(xm*rez))
    ym = np.abs(np.around(ym*rez))
  
    cols = xm.astype(int)
    rows = ym.astype(int)
  
    # empty image array with dimensions of detector
    imRowMax = np.abs(np.around(rez*2.0*hh))+1
    imColMax = np.abs(np.around(rez*2.0*hw))+1
    im = np.zeros((imRowMax,imColMax),dtype=np.uint32)
    noise = np.random.random((imRowMax,imColMax))
    noise = noise*.01
        
    vals = fin_intensities

    # set the image pixels
    im[rows,cols] = vals*255
    noise = noise*255
    im = im + noise
 

    import matplotlib.pylab as plt
    
    # lineout average
    avglin = np.average(im[400:800,:],axis=0)

    fig = plt.figure()
    axim =fig.add_subplot(221)
    axim.imshow(im,cmap='binary')
    
    axdat = fig.add_subplot(222)
    axdat.plot(system.source.spectrum[0],conv_spect)
    
    axlin = fig.add_subplot(223)
    axlin.plot(avglin)
    plt.show()
    #plt.imshow(im,cmap='binary')
    #plt.plot(system.source.spectrum[0],conv_spect)
    #from mayavi import mlab
   # mlab.imshow(im,colormap='gray',interpolate=True)
    #mlab.view(0,0)
   # mlab.roll(-90)
   # mlab.show()
    


    
    
    