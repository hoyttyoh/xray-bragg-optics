# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 12:48:49 2013

@author: cad
"""

""" Setup for cylindrical KAP (2d = 26.6A) on XP for DEF/IP comparison study """

from BraggOpticsClasses import *
from DetectorClasses import *

from SourceClasses import *
from SetupClasses import *
from VisualizationFunctions import *
from RaytraceFunctions import *

source = SphericalSource(diameter=1e-6)
source.loadSpectrum("source_spectra/Cu_800_5e22.txt")

optic = CylindricalBraggOptic(101.,20.,350.0,twoD=26.6)
optic.setGaussianReflectProfile(fwhm=1e-4)
detector = Detector(33.0,14.4)

setup = RayPathSystem(source=source,optic=optic,detector=detector)
setup.setDiffractionOrder(2)
setup.setSourcePosition(620.0,68.0)
#setup.setSourcePositionOnRowland(22.64) # 22.0 is good for Ti He-alpha
#setup.setDetectorPositionOnRowland(70.0)
setup.setDetectorPosition(150.0, 68.0)

#setup.setDetectorPosition(250.0,22.9)
#setup.setDetectorAtSagFocus()

TraceSystem(setup,nrays=1e6,useReflectivityCurve=False)
MayaviVisual(setup,scalarBar=False,units="Angstrom",showRowland=False)
SpectraImage(setup,bgNoise=False)