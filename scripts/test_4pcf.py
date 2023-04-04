#Preamble
import numpy as np
import astropy.io.fits as pyf
import os
import sarabande


hdulist = pyf.open('/Users/toriwilliamson/Desktop/LearningTheISM/dens_t800.fits.gz')
data = np.log(hdulist[0].data.astype(np.float64))

#FORMAT THE DATA LIKE THIS
data = (data-np.mean(data))/(np.std(data))
# data = data[:128, :128, :128]

#string to directory to save data into
save_dir = '/Users/toriwilliamson/Desktop/LearningTheISM/scripts/output/'

#create measure_obj
_4PCF = sarabande.measure(nPCF=4, projected=False,
                          density_field_data = data, save_dir=save_dir,
                          save_name='example', nbins=5, ell_max=0,
                          physical_boxsize=256, rmin=1e-14, rmax=128, normalize=True)

sarabande.calc_zeta(_4PCF, verbose_flag=True, parallelized=False)

np.save("test_4pcf_file.npy", _4PCF.zeta)
