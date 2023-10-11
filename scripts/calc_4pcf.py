import numpy as np
import astropy.io.fits as pyf
import os
import sarabande
import time
import argparse

########################################
#       COMMAND LINE ARGUMENTS
########################################

#Command Line Arguments
parser = argparse.ArgumentParser(description='4pcf script')
parser.add_argument('-ell_max', type=int,
                    help="""This is the maximum ell value.""")
parser.add_argument('-n_bins', type=int,
                    help="""This is the number of bins.""")
parser.add_argument('-Ma', type=float,
                    help="""This is the Alfvenic number.""")
parser.add_argument('-Ms', type=float,
                    help="""This is the Mach number.""")
parser.add_argument('-t_slice', type=int,
                    help="""This is the time slice.""")

results = parser.parse_args()

print(f"We will be using ell_max {results.ell_max} in our data.")
print(f"We will be using {results.n_bins} bins in our data.")
print(f"We will be using {results.t_slice} time in our data.")
print(f"We will be using {results.Ma} Alfvenic number in our data.")
print(f"We will be using {results.Ms} Mach number in our data.")

########################################
#               FUNCTIONS
########################################

def open_file(Ms, Ma, t_slice):
    """
    Arguments:
    Ms [float]: The sonic mach number of the density field we’re using
    Ma [float]: The Alfvenic mach number of the density field we’re using
    t_slice [int]: the time slice of the simulation a value that can be [500, 550, … 900]

    returns: data [numpy array] of shape (256, 256, 256) representing the density field of a simulation after it has been reformatted to a log-normal transformation

    """

    if Ms == 0.7 and Ma == 0.7:
        Simtype = 'Ms_0p7__Ma_0p7/'
    elif Ms == 0.7 and Ma == 2.0:
        Simtype = 'Ms_0p7__Ma_2p0/'
    elif Ms == 1.2 and Ma == 0.7:
        Simtype = 'Ms_1p2__Ma_0p7/'
    elif Ms == 1.2 and Ma == 2.0:
        Simtype = 'Ms_1p2__Ma_2p0/'
    elif Ms == 7.0 and Ma == 0.7:
        Simtype = 'Ms_7p0__Ma_0p7/'
    elif Ms == 7.0 and Ma == 2.0:
        Simtype = 'Ms_7p0__Ma_2p0/'
    else:
        Simtype = None #throw error


    Folder = "/blue/zslepian/williamson.v/data/" + Simtype + f"dens_t{t_slice}.fits.gz"


    hdulist = pyf.open(Folder)
    data = hdulist[0].data.astype(np.float64) #here’s our data

    #reformatted the data to a log-normal transformation
    data = np.log(data)
    data = (data-np.mean(data))/(np.std(data))

    return data


def calc_4pcf(data, Ms, Ma, t_slice, nbins, ell_max):
    """
    Arguments:
    data [np.array] : The density field loaded in, 
    should be a log normal density field
    Ms [float] : The sonic mach number of the density 
    field we’re using
    Ma [float] : The Alfvenic mach number of the density 
    field we’re using
    t_slice [int] : the time slice of the simulation a 
    value that can be [500, 550, … 900]
    nbins [int] : number of radial bins we’re using
    ell_max [int] : the maximum ell value we’re using
	
    returns → returns nothing, but it saves the full and disconnected 4pcf measurements to numpy files

    """

    #string to directory to save temporary calculations into
    save_dir = '/blue/zslepian/williamson.v/output/'

    #make a unique identifier string
    if Ms == 0.7 and Ma == 0.7:
        Simtype = 'Ms_0p7__Ma_0p7'
    elif Ms == 0.7 and Ma == 2.0:
        Simtype = 'Ms_0p7__Ma_2p0'
    elif Ms == 1.2 and Ma == 0.7:
        Simtype = 'Ms_1p2__Ma_0p7'
    elif Ms == 1.2 and Ma == 2.0:
        Simtype = 'Ms_1p2__Ma_2p0'
    elif Ms == 7.0 and Ma == 0.7:
        Simtype = 'Ms_7p0__Ma_0p7'
    elif Ms == 7.0 and Ma == 2.0:
        Simtype = 'Ms_7p0__Ma_2p0'
    else:
        Simtype = None #throw error

    save_name = f"_{Simtype}_{t_slice}_"
    
    print(f"\nThe save name is {save_name}.\n")

    #create measure_obj
    _4PCF = sarabande.measure(nPCF=4, projected=False,
            density_field_data = data, save_dir=save_dir,
            save_name=save_name, nbins=nbins, ell_max=ell_max,
            physical_boxsize=256, rmin=1e-14, rmax=128, normalize=True)

    print("\nObject has been made, running function now.\n")

    sarabande.calc_zeta(_4PCF, verbose_flag=True, parallelized=True,
    calc_bin_overlaps=False, calc_odd_modes=True, calc_disconnected=True)

    #save files
    np.save(f"zetas/4pcf_full_{save_name}.npy", _4PCF.zeta); np.save(f"zetas/4pcf_disconnected_{save_name}.npy", _4PCF.zeta_disconnected)

    return None


##############################
#     WE RUN CODE BELOW
##############################

if __name__ == "__main__": #for scripts only

    #Loading in the Data
    data = open_file(results.Ms, results.Ma, results.t_slice)

    print(f"\nThe data has been loaded in.\n")

    start = time.time()

    # calculate the 4PCF
    calc_4pcf(data, results.Ms, results.Ma, results.t_slice, results.n_bins, results.ell_max)

    end = time.time()

    #Results
    print(f"\nIt took {end - start} seconds to run and it successfully saved!\n")


