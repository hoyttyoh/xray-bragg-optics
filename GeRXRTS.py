# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 13:23:47 2013

@author: cad
"""

""" Refocused XRTS with Ge using Ti He-alpha line """

from bragg_optics import *
from detectors import *
from sources import *
from optical_system import *
from rt_visuals import *
from rt_functions import *

source = SphericalSource(diameter=50e-6)
#source.loadSpectrum("C:/Users/cad/thesis_research/raytracer/sombrero/source_spectra/Ti_900_1e23.txt")
optic = SphericalBraggOptic(50.,20.,180.0,twoD=2.828)
optic.setGaussianReflectProfile(fwhm=1e-4)
detector = Detector(20,10)

setup = RayPathSystem(source=source,optic=optic,detector=detector)
setup.setDiffractionOrder(1)
setup.setSourcePosition(690.0,22.0)
#setup.setSourcePositionOnRowland(22.64) # 22.0 is good for Ti He-alpha
#setup.setDetectorPositionOnRowland(22.64)
#setup.setDetectorPosition(200.0, 22.0)

setup.setDetectorPosition(93.0,22.0)
#setup.setDetectorAtSagFocus()

trace_system_path(setup,nrays=1e5)
#DispersionCurve(setup)
mayavi_visual(setup,scalarBar=True,units="keV")
spectra_image(setup,bgNoise=False)
