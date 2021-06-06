# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 10:20:48 2013

@author: cad
"""

import numpy as np

""" Class that contains information on the optical setup. """

class RayPathSystem(object):
    
    def __init__(self, source, optic, detector):
        
        self.source = source
        self.optic = optic
        self.detector = detector
        
        self.initBraggNormAngle = 0.0
        self.initBraggAngle = 90.0 - self.initBraggNormAngle
        self.sourceToOpticDist = 0.0
        self.opticToDetectorDist = 0.0
 
        self.diffractionOrder = 1.0
        self.sourceToOpticVec = np.zeros((3,1))
        self.sourceVec = np.zeros((3,1))
        self.detectorVec = np.zeros((3,1))
        self.detectorToOpticVec = np.zeros((3,1))
        self.detectorToSourceVec = np.zeros((3,1))
        
        # normal to detector surface
        self.detectorNorm = np.zeros((3,1))
        
        # generic system name for storage
        self.sysName = "opt_sys_001"
    
        self.traceType = None
        self.meshShape = None
        
        # is source or detector on Rowland circle?
        self.sourceOnRowland = False
        self.detectorOnRowland = False
        
        # store all the relevant points for visualization and other post-processing
        self.sourcePoints = None
        self.opticPoints = None
        self.apertureIntersectPoints = None
        self.detectorIntersectPoints = None
        self.ewi = None
        self.numRays = None
     
    def setSystemName(self,newName):
        self.sysName = newName
     
    def setTraceType(self, new):
        self.traceType = new
        
    def setNumRays(self,nrays):
        self.numRays = nrays
                
    def setDiffractionOrder(self, order):
        self.diffractionOrder = order
     
    def setMeshGridShape(self,shape):
        self.meshShape = shape
    # the source to optic vector
    def setSourcePosition(self, d, angle, zoffset=0.0):
        
        normAngleRad = np.deg2rad(angle)
        self.initBraggNormAngle = angle
        self.sourceToOpticDist = d
        self.sourceToOpticVec = np.array([d*np.cos(normAngleRad),
                                -d*np.sin(normAngleRad),
                                    zoffset])
                                    
        # source relative to global origin                   
        self.sourceVec = self.optic.refCentroid-self.sourceToOpticVec    
    
    def setSourcePositionOnRowland(self, angle):
        normAngleRad = np.deg2rad(angle)
        dRow = self.optic.radius*np.cos(normAngleRad)
        self.setSourcePosition(dRow,angle)
        self.sourceOnRowland=True
        
    def setDetectorPosition(self, d, angle, zoffset=0.0):
        
        normAngleRad = np.deg2rad(angle)
        self.opticToDetectorDist = d
        self.detectorToOpticVec = np.array([d*np.cos(normAngleRad),
                                   d*np.sin(normAngleRad),
                                    zoffset])
                                    
        # define the location of the detector centroid                            
        self.detectorVec = self.optic.refCentroid - self.detectorToOpticVec    
                                
        # define the detector normal vector
        self.detectorToSourceVec = self.sourceVec - self.detectorVec
        z = np.array([0.0,0.0,1.0])
        self.detectorNorm = -np.cross(self.detectorToSourceVec,z)
        self.detectorNorm /= np.linalg.norm(self.detectorNorm)
        
    def setDetectorPositionOnRowland(self, angle):
        normAngleRad = np.deg2rad(angle)
        dRow = self.optic.radius*np.cos(normAngleRad)
        self.setDetectorPosition(dRow,angle)
        self.detectorOnRowland = True
        
    def setDetectorAtSagFocus(self, zoffset=0.0):
        """ Set the detector distance at the saggital focus point """
        """need to fix this so only spherically-bent can use """
        
        # sagittal focal distance from optic
        r = self.optic.radius
        p = self.sourceToOpticDist
        normAngleRad = np.deg2rad(self.initBraggNormAngle)
        theta = np.pi/2.0 - normAngleRad
        q = 1.0/((2*np.sin(theta)/r)-(1.0/p))
        

        self.opticToDetectorDist = q
        
        self.detectorToOpticVec = np.array([q*np.cos(normAngleRad),
                                   q*np.sin(normAngleRad),
                                    zoffset])
                                    
        # define the location of the detector centroid                            
        self.detectorVec = self.optic.refCentroid - self.detectorToOpticVec    
                                
        # define the detector normal vector
        self.detectorToSourceVec = self.sourceVec - self.detectorVec
        z = np.array([0.0,0.0,1.0])
        self.detectorNorm = -np.cross(self.detectorToSourceVec,z)
        self.detectorNorm /= np.linalg.norm(self.detectorNorm)    
 

# deprecated       
    def setPostProcData(self, source, optic, detector, ewi=None, aperture=None):
        """ Store points for visualization and post processing """
        
        self.sourcePoints = source
        self.opticPoints = optic
        self.apertureIntersectPoints = aperture
        self.detectorIntersectPoints = detector
        self.ewi = ewi
        self.setNumRays(self.sourcePoints.size/3.0)
        
        

        
        
        
        
        
        
        


        
        
        
        
        
    
        