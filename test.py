# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 13:21:54 2013

@author: cad
"""


from bragg_optics import *
from detectors import *

from sources import *
from optical_system import *

from rt_functions import *
from rt_visuals import *


# define components
source = SphericalSource(diameter=1e-3)
#source.loadSpectrum("/media/cad/Research/thesis/xrts_analysis/2999/lineouts/fixed_Ti_spectrum")
#source.setMultiGaussianSpectrum([4750.0,4720.0,4690.0],[1.0,.3,.7],[25.0,10,10.0])
source.loadSpectrum("/media/cad/Research/thesis/xrts_analysis/2999/Al_cold_scatter_exo/run2_scatter.ppd")
optic = SphericalBraggOptic(55e-3,18e-3,180.0e-3,twoD=6.7e-10)
optic.setGaussianReflectProfile(fwhm=20) # fwhm in eV , 10eV ~ 2e-3 rad ~ 0.15 deg

detector = Detector(5e-2,2.5e-2)

# define the optical path
setup = RayPathSystem(source=source,optic=optic,detector=detector)
setup.setDiffractionOrder(2.0)
setup.setSourcePosition(120.0e-3,38.8)
setup.setDetectorPosition(120e-3,38.8)
#setup.setDetectorPositionOnRowland(22.64)
#setup.setDetectorAtSagFocus()

# trace thru the system
#grid_trace_system(setup,hpoints=500,vpoints=500)
random_trace_system(setup,nrays=1e5)
#grid_optic_intercept_info(setup)

detector_intercept_info(setup)
#mayavi_visual(setup)



