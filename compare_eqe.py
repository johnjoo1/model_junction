import pylab
import glob
from scipy.interpolate import UnivariateSpline

eqe_files = glob.glob('*.eqe')

doping = []
eqe_at_doping = []
for fname in eqe_files:
    
    eqe_spec = pylab.loadtxt(fname, skiprows = 6)
    wavelength = eqe_spec[:,0]
    eqe = eqe_spec[:,1]
    current_doping = float(fname.split('eqe')[1][:-1])
    
    pylab.plot(wavelength,eqe, label=current_doping)
    
    spline = UnivariateSpline(x=wavelength,
                                       y=eqe,
                                       k=1,
                                       s=0.1)
    wavelength_spline = 390e-9
    eqe_spline = spline(wavelength_spline)
    
    doping.append(current_doping)
    eqe_at_doping.append(eqe_spline)

pylab.figure()
pylab.plot(doping, eqe_at_doping, 'o')    
pylab.show()
pylab.legend()



