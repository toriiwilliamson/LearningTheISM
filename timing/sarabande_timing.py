#!/usr/bin/env python

import numpy as np
import astropy.io.fits as pyf
import os
import sarabande
import time
from tqdm import tqdm

start = time.time()

hdulist = pyf.open('/blue/zslepian/williamson.v/data/Ms_0p7__Ma_0p7/dens_t500.fits.gz')
data = np.log(hdulist[0].data.astype(np.float64))
data = (data-np.mean(data))/(np.std(data))

save_dir = '/blue/zslepian/williamson.v/output/'

_4PCF = sarabande.measure(nPCF=4, projected=False,
                          density_field_data = data, save_dir=save_dir,
                          save_name='example', nbins=15, ell_max=2,
                          physical_boxsize=256, rmin=1e-14, rmax=128,
			  normalize=True, calc_bin_overlaps=True, calc_odd_modes=True)

sarabande.calc_zeta(_4PCF, verbose_flag=False, parallelized=False)

end = time.time()

print(f"it took {end - start:0.3f} seconds to finish")

np.save("timing_4pcf_file.npy", _4PCF.zeta)

# In[ ]:




