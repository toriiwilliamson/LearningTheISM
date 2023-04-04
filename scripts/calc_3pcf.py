# %%
import numpy as np
import astropy.io.fits as pyf
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1 as axgrid
from mpl_toolkits.axes_grid1 import make_axes_locatable
from itertools import product
from tqdm import tqdm
import time
import h5py
import sarabande

#-----------------------------------------------------
#        Mach Number and Alfvénic Mach Number
#-----------------------------------------------------

import argparse

parser = argparse.ArgumentParser(description='3PCF Computation Script')

parser.add_argument('-M_s', type=float, required=True,
                    help='Mach Number.')

parser.add_argument('-M_A', type=float, required=True,
                   help='Alfvénic Mach Number.')                   

results = parser.parse_args()


M_s, M_A = results.M_s, results.M_A
nbins, m_max = 5,0

#-----------------------------------------------------
#                   Data Formatting
#-----------------------------------------------------

def format_dir_string(M_s, M_A):
    """
    Format Directory String based on Mach Numbers
    """
    if M_s == 0.7:
        M_s_str = '0p7'
    elif M_s == 2.0:
        M_s_str = '2p0'
    else:
        M_s_str = None

    if M_A == 0.7:
        M_A_str = '0p7'
    elif M_A == 2.0:
        M_A_str = '2p0'
    else:
        M_A_str = None

    save_dir = f'/blue/zslepian/williamson.v/data/Ms_{M_s_str}__Ma_{M_A_str}'
    
    return save_dir

def compile_data(save_dir):
    """
    A function to compile all the data for each parameter set
    """
    
    all_data = {}
    
        
    for i in range(500, 901, 50):
        save_name = f'/dens_t{i}'
        file = save_dir + save_name+'.fits.gz'
        hdulist = pyf.open(file)
        data = hdulist[0].data.astype(np.float64)
        data = np.abs((data - np.mean(data))/np.std(data))
        
        
        all_data[f't_{i}'] = data
        
    return all_data

def compute_NPCF_Coefficients(nbins, m_max):
    """
    Measure NPCF Coefficients.
    """
    rmax = 128
    proj_3PCF_coefficients = []

    for i in tqdm(range(len(MHD_data))):
        for slice_ in range(0,256, 5):
            t = 500 + 50*i
            
            #Along x-axis
            density_field = MHD_data[f"t_{t}"]
            proj_3pcf = sarabande.measure(nPCF=3, projected=False, m_max=m_max, 
                                          density_field_data=density_field, nbins=nbins, 
                                          physical_boxsize=256, rmin=1, rmax=rmax, normalize=True,
                                          save_name='test', save_dir='')

            sarabande.calc_zeta(proj_3pcf,verbose_flag=False)
            proj_3PCF_coefficients.append(proj_3pcf.zeta)
            
    
    proj_3PCF_coefficients_real, proj_3PCF_coefficients_imag = [], []
    for i in range(len(proj_3PCF_coefficients)):
        proj_3PCF_coefficients_real.append(proj_3PCF_coefficients[i].real)
        proj_3PCF_coefficients_imag.append(proj_3PCF_coefficients[i].imag)
        
    return proj_3PCF_coefficients_real, proj_3PCF_coefficients_imag

def save_coefficients(save_dir, proj_3PCF_coefficients_real, proj_3PCF_coefficients_imag):
    hf = h5py.File(f'{save_dir}/Coefficients_3PCF_real.h5', 'w')
    hf.create_dataset('data', data=proj_3PCF_coefficients_real)
    hf.close()
    
    hf = h5py.File(f'{save_dir}/Coefficients_3PCF_imag.h5', 'w')
    hf.create_dataset('data', data=proj_3PCF_coefficients_imag)
    hf.close()
    

if __name__ == "__main__":
    
    print("Formatting Data ...")
    save_dir = format_dir_string(M_s, M_A)
    MHD_data = compile_data(save_dir)
    
    print("Measuring 3PCF Coefficients ...")
    Coefficients_Real, Coefficients_Imag = compute_NPCF_Coefficients(nbins, m_max)
    
    print("Saving Coefficients")
    save_coefficients(save_dir, Coefficients_Real, Coefficients_Imag)
    


