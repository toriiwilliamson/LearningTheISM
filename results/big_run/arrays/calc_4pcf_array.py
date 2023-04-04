import numpy as np
import astropy.io.fits as pyf
import os
import sarabande
import time

idx = int(os.environ["SLURM_ARRAY_TASK_ID"])

print(f'the id of this job is {idx}')
# execute the rest of the script using myparam

#os.mkdir(f'/blue/zslepian/jsunseri/output/output_{idx}')


hdulist = pyf.open(f'/blue/zslepian/williamson.v/data/Ms_0p7__Ma_0p7/dens_t{idx}.fits.gz')
data = np.log(hdulist[0].data.astype(np.float64))

#FORMAT THE DATA LIKE THIS
data = (data-np.mean(data))/(np.std(data))
# data = data[:128, :128, :128]

#string to directory to save data into
save_dir = f'/blue/zslepian/williamson.v/output/output_{idx}/'

#create measure_obj
_4PCF = sarabande.measure(nPCF=4, projected=False,
                          density_field_data = data, save_dir=save_dir,
                          save_name='_Ms_0p7__Ma_0p7_', nbins=32, ell_max=2,
                          physical_boxsize=256, rmin=1e-14, rmax=128, normalize=True)

sarabande.calc_zeta(_4PCF, verbose_flag=False, parallelized=True)

np.save(f"results_4pcf_{idx}.npy", _4PCF.zeta)

