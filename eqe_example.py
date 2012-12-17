#!/usr/bin/env python

## -*-Pyth-*-
 # ########################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # Author: Jonathan Guyer <guyer@nist.gov>
 #   mail: NIST
 #    www: <http://www.ctcms.nist.gov/fipy/>
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 # 
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are 
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 # 
 # ########################################################################
 ##

from fipy import Grid1D
from fipy.tools import numerix
from fipy.viewers import Viewer

import solar.device
reload(solar.device)
from solar.device import Device
from solar.impurities import Donor, Acceptor
from solar.materials import CdTe, CdSe, Vacuum
from solar.contacts import OhmicContact
from solar.traps import ShockleyReadHallTrap
from solar.lights import AM1_5
from solar.parsimoniousIterator import ParsimoniousIterator

n_thickness = 1e-6 # m
p_thickness = 2e-6 # m
grid_resolution = 2.0e-8 # m

# Creates series of spacings (dx) that decrease rapidly according to a power law.
compression_factor = 0.8
compression_count = 60

compressed_dx = grid_resolution * compression_factor**numerix.arange(compression_count)
compressed_length = compressed_dx.sum()
compressed_dx = list(compressed_dx)

# Creates the spacings for the n-side of the pn junction.  Starts dense (short dx) near the contact, gets longer in the middle of the n region, then gets shorter again at the junction.  The compressed_dx spacing is used to describe the spacings (dx) as you get closer to the interesting regions, which are the contacts and the junction.
n_uncompressed_thickness = n_thickness - 2 * compressed_length
n_dx = n_uncompressed_thickness / round(n_uncompressed_thickness / grid_resolution)
n_dx = compressed_dx[::-1] + [n_dx] * int(n_uncompressed_thickness / n_dx) + compressed_dx

# Same thing is done for the p side.
p_uncompressed_thickness = p_thickness - 2 * compressed_length
p_dx = p_uncompressed_thickness / round(p_uncompressed_thickness / grid_resolution)
p_dx = compressed_dx[::-1] + [p_dx] * int(p_uncompressed_thickness / p_dx + 0.5) + compressed_dx

mesh = Grid1D(dx=n_dx + p_dx)

x = mesh.cellCenters[0]

# To view the node density:
#from pylab import plot
#plot (x, ones(len(x)), 'o')

ntype = x < n_thickness # print ntype to see True and False values
ptype = x >= n_thickness
vacuum = ~(ntype | ptype)

# Assigning material properties to each Cell
cdse_cdte = (CdSe() * ntype + CdTe() * ptype + Vacuum() * vacuum)

nFaces = mesh.facesLeft
pFaces = mesh.facesRight
illuminatedFaces =mesh.facesLeft

nContact = OhmicContact(faces=nFaces)
pContact = OhmicContact(faces=pFaces)

diode = Device(mesh=mesh, 
               temperature=300, # K
               material=cdse_cdte, 
               dopants=(Donor(concentration=1e24 * ntype, # m**-3"
                              ionizationEnergy=-1.), # eV
                        Acceptor(concentration=1e21 * ptype, # m**-3
                                 ionizationEnergy=-1.)), # eV
               traps=(ShockleyReadHallTrap(recombinationLevel=0, # eV
                                           electronMinimumLifetime=3e-11, # s
                                           holeMinimumLifetime=3e-11),), # s
               contacts=(pContact, nContact))

diode.contacts[1].bias.value = 0. # V
diode.contacts[0].bias.value = 0. # V

# solve in dark and then illuminate
diode.solve(solver=None, sweeps=10, outer_tol=1e4)

viewer = Viewer(vars=(diode.Ec, diode.Efn, diode.Efp, diode.Ev))
npViewer = Viewer(vars=(diode.n, diode.p), log=True)

#raw_input("dark done.")

light = AM1_5(orientation=[[1]], 
              faces=illuminatedFaces)
              
# sample the light source from 300 nm to 1 um at 10 nm intervals
diode.illuminate(light(wavelength=numerix.arange(300e-9, 1000e-9, 10e-9))) # m

def view():
    viewer.plot()
    npViewer.plot()
 
view()

#raw_input("light done.")

eqe = diode.EQE(light, numerix.arange(320e-9,1000e-9, 10e-9), path=None, adapt=True)
import pylab
pylab.figure()
wavelength, eqe = zip(*eqe)
pylab.plot(wavelength,eqe)
header = 'afdsa\nafsdafsd\n'
with open('eqe1.csv', 'wb') as f:
    f.write(header)
    for i in zip(wavelength,eqe):
        wavelength_point, eqe_point = i
        f.write('{wavelength}\t{eqe}\n'.format(wavelength=wavelength_point, eqe=eqe_point))

#    
#JV = diode.JV(contact=nContact, 
#              biases=ParsimoniousIterator(start=0., stop=0.9, num=20), # V
#              path="JV.txt",
#              sweeps=11, outer_tol=1e-4, currentContinuity=1.,
#              callback=view)

#raw_input("Voc")

#Voc = diode.Voc(contact=nContact, JV=JV)

#print "Voc = ", Voc

#raw_input("Vmax")
#                            
#Vmax = diode.Vmax(contact=nContact, JV=JV, Vtol=0.01, Ptol=0.01)

#print "Vmax = ", Vmax

#raw_input("done")
                            
