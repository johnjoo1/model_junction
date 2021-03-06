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
from solar.materials import ZnO, Cu2O, Vacuum
from solar.contacts import OhmicContact
from solar.traps import ShockleyReadHallTrap
from solar.lights import AM1_5

import pylab

class EqeDoping(object):
    def __init__(self, Na, Nd):
        n_thickness = 2e-6 # m
        p_thickness = 4e-6 # m
        grid_resolution = 2.0e-8 # m
        self.n_thickness = n_thickness
        self.p_thickness = p_thickness
        self.grid_resolution = grid_resolution

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
        self.mesh = mesh

        x = mesh.cellCenters[0]
        self.x = x

        # To view the node density:
        #from pylab import plot
        #plot (x, ones(len(x)), 'o')

        ntype = x < n_thickness # print ntype to see True and False values
        ptype = x >= n_thickness
        vacuum = ~(ntype | ptype)

        # Assigning material properties to each Cell
        zno_cu2o = (ZnO() * ntype + Cu2O() * ptype + Vacuum() * vacuum)
        self.zno_cu2o = zno_cu2o
        self.ntype = ntype
        self.ptype = ptype

        nFaces = mesh.facesLeft
        pFaces = mesh.facesRight
        illuminatedFaces =mesh.facesLeft
        self.illuminatedFaces = illuminatedFaces

        nContact = OhmicContact(faces=nFaces)
        pContact = OhmicContact(faces=pFaces)
        self.nContact = nContact
        self.pContact = pContact

        self.Na=Na
        self.Nd=Nd

        self.diode = Device(mesh=mesh, 
                       temperature=300, # K
                       material=zno_cu2o, 
                       dopants=(Donor(concentration=Nd * ntype, # m**-3"
                                      ionizationEnergy=-1.), # eV
                                Acceptor(concentration=Na * ptype, # m**-3
                                         ionizationEnergy=-1.)), # eV
                       traps=(ShockleyReadHallTrap(recombinationLevel=0, # eV
                                                   electronMinimumLifetime=3e-11, # s
                                                   holeMinimumLifetime=3e-11),), # s
                       contacts=(pContact, nContact))

        self.diode.contacts[1].bias.value = 0. # V
        self.diode.contacts[0].bias.value = 0. # V

    def change_doping(self, Na):
        self.Na=Na

        self.diode = Device(mesh=self.mesh, 
                       temperature=300, # K
                       material=self.zno_cu2o, 
                       dopants=(Donor(concentration=self.Nd * self.ntype, # m**-3"
                                      ionizationEnergy=-1.), # eV
                                Acceptor(concentration=Na * self.ptype, # m**-3
                                         ionizationEnergy=-1.)), # eV
                       traps=(ShockleyReadHallTrap(recombinationLevel=0, # eV
                                                   electronMinimumLifetime=3e-11, # s
                                                   holeMinimumLifetime=3e-11),), # s
                       contacts=(self.pContact, self.nContact))
                       
        self.diode.contacts[1].bias.value = 0. # V
        self.diode.contacts[0].bias.value = 0. # V
        
    def nt(self, t, A=1, c = 0):
        return A*numerix.log(self.diode.generationRate.value*t) + c
    
    def solve_dark(self, save = False):
        # solve in dark and then illuminate
        self.diode.darken()
        self.diode.solve(solver=None, sweeps=10, outer_tol=1e4)
        if save == True:
            self.diode.save('dark{Na}.band'.format(Na=self.Na))
    
    def solve_light(self,  save = False):
        self.light = AM1_5(orientation=[[1]], 
                      faces=self.illuminatedFaces)                    
        # sample the light source from 300 nm to 1 um at 10 nm intervals
        self.diode.illuminate(self.light(wavelength=numerix.arange(300e-9, 1000e-9, 10e-9))) # m
        self.diode.solve(solver=None, sweeps=10, outer_tol=1e4)
        if save == True:
            self.diode.save('light{Na}.band'.format(Na=self.Na))


    def solve_eqe(self, view=False, save = False):
        self.eqe_spectrum = self.diode.EQE(self.light, numerix.arange(320e-9,700e-9, 1e-9), path=None, adapt=True)
        if view ==True:
            import pylab
            pylab.figure()
            wavelength, eqe = zip(*self.eqe_spectrum)
            pylab.plot(wavelength,eqe)
        if save == True:
            self.save_eqe('eqe{Na}.eqe'.format(Na = self.Na))

    def save_eqe(self,fname):
        header = """
        ZnO _doping:{Nd}
        Cu2O doping: {Na}
        ZnO thickness: {n_thick}
        Cu2O thickness: {p_thick}
        grid resolution: {grid_res}
""".format(Na=self.Na, Nd=self.Nd, n_thick = self.n_thickness, p_thick = self.p_thickness, grid_res = self.grid_resolution)

        with open(fname, 'wb') as f:
            f.write(header)
            for i in self.eqe_spectrum:
                wavelength_point, eqe_point = i
                f.write('{wavelength}\t{eqe}\n'.format(wavelength=wavelength_point, eqe=eqe_point))
                
    def view(self, legend = False):
        '''
        The fipy.viewer method doesn't seem to update after a single loop.
        '''
#        if hasattr(self, 'viewer'):
#            self.viewer.plot()
#        else:
#            self.viewer = Viewer(vars=(self.diode.Ec, self.diode.Efn, self.diode.Efp, self.diode.Ev))
        if hasattr(self, 'plot_band'):
            self.plot_band.plot(self.x.value, self.diode.Ec.value, label = '$E_c$')
            self.plot_band.plot(self.x.value, self.diode.Ev.value, label = '$E_v$')
            self.plot_band.plot(self.x.value, self.diode.Efn.value, label = '$E_fn$')
            self.plot_band.plot(self.x.value, self.diode.Efp.value, label = '$E_fp$')
            pylab.draw()
        else:
            fig = pylab.figure() 
            self.plot_band = fig.add_subplot(111)
            self.plot_band.plot(self.x.value, self.diode.Ec.value, label = '$E_c$')
            self.plot_band.plot(self.x.value, self.diode.Ev.value, label = '$E_v$')
            self.plot_band.plot(self.x.value, self.diode.Efn.value, label = '$E_fn$')
            self.plot_band.plot(self.x.value, self.diode.Efp.value, label = '$E_fp$')
            pylab.draw()
            if legend == True:
                pylab.legend()

    def equilibrium_check(self):
        if [round(x,2) for x in self.diode.Efp.value] == [round(x,2) for x in self.diode.Efn.value]:
            return True #is in equilibrium
        else:
            return False #
            
if __name__ == "__main__":  
    time = 10**numerix.arange(0,100000)
    for i,t in enumerate(time):
        if i == 0:
            Na = 1e20
            eqedoping = EqeDoping(Na,1e24)
        else:
            Na = eqedoping.nt(t)*(eqedoping.x.value>=eqedoping.n_thickness)
            eqedoping.change_doping(Na)
        eqedoping.solve_dark(save=False)
        
        print 'dark'
        if eqedoping.equilibrium_check():
            print 'Efn = Efp'
        else:
            print 'not equilibrium!!!'

        eqedoping.view()
#        raw_input('wait.')
        
        eqedoping.solve_light( save=False)
        print 'illuminated'
        
#        eqedoping.solve_eqe(view=True, save=True)

        if eqedoping.equilibrium_check():
            print 'Efn = Efp'
        else:
            print 'obviously, not equilibrium.  under illumination.'
#        raw_input('wait.')
        
    
                            
