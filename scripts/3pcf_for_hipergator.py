import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1 as axgrid
from mpl_toolkits.axes_grid1 import make_axes_locatable
import astropy.io.fits as pyf
import os
import sarabande
from google.colab import drive
drive.mount('/content/drive')
save_name = 'demo'
eps = 1e-15
ell_max = 2
nbins = 4
#don't set it to zero, because spherical harmonics are ill-defined at the origin
bin_min = 1-1e-5
#radius from the origin.
bin_max = 128+1e-5

switch = {
	'LIN' : np.linspace(bin_min, bin_max, nbins+1),
	'INV' : 1./np.linspace(1./bin_min, 1./bin_max, nbins+1),
	'LOG' : np.exp(np.linspace(np.log(bin_min), np.log(bin_max), nbins+1))
}
bin_edges = switch['LIN']
ld_one_d = 256
boxsize = ld_one_d
N_mesh = ld_one_d
x = np.linspace(-ld_one_d/2, ld_one_d/2-1 , ld_one_d)
xsq = x*x
m_center = np.where(x==0)[0][0]
X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
X = -X
Y = -Y
Z = -Z
Xsq, Ysq, Zsq = np.meshgrid(xsq, xsq, xsq,indexing='ij')
#PRECOMPUTE POWERS OF ARRAYS AND COMBINATIONS (e.g. x - iy).
Rsq = Xsq+Ysq+Zsq
R = np.sqrt(Rsq)
del Rsq
zero_ind = np.where(R==0)
R[zero_ind] = eps
X[zero_ind] = eps
Y[zero_ind] = eps
Z[zero_ind] = eps
boundsandnumber = np.zeros((2, nbins+1))
boundsandnumber[0,:] = bin_edges
for i in range(nbins):
  boundsandnumber[1,i] = np.sum(np.logical_and(R >= bin_edges[i],
                                               R < bin_edges[i+1]))
np.save('bin_bounds_and_pixel_number_'+save_name+'.npy',boundsandnumber)
def ylm_save(ylm, ell, m):
    np.save('ylm_'+save_name+'_'+str(ell)+'_'+str(m)+'.npy',ylm)#zs: to work on my laptop.
    del ylm

def ylm_transform_save(ylm_on_shell, ell, m, i):
    FT = np.fft.fftn(np.fft.fftshift(ylm_on_shell))
    np.save('YLMtilde_'+save_name+'_'+str(ell)+'_'+str(m)+'_bin_'+str(i)+'.npy',FT)
    del FT
    #COMPUTE YLMS SEQUENTIALLY AND SAVE.
#ell, m = 0,0
y00 =.5*(1./np.pi)**.5*np.ones((ld_one_d,ld_one_d,ld_one_d))
ylm_save(y00, 0, 0)
del y00

#ell, m = 1, -1
xdivr = X/R
del X
ydivr = Y/R#we'll need these as individuals later anyway.
del Y
xmiydivr = xdivr - 1j*ydivr
y1m1 = .5*np.sqrt(3./(2.*np.pi))*xmiydivr
ylm_save(y1m1, 1, 1)
del y1m1

#ell, m = 1, 0
zdivr = Z/R
del Z
y10 = .5*np.sqrt(3./np.pi)*zdivr
ylm_save(y10, 1, 0)
del y10

#ell, m = 2, -2
xmiydivrsq = xmiydivr*xmiydivr
y2m2 = .25*np.sqrt(15./(2.*np.pi))*xmiydivrsq
ylm_save(y2m2, 2, 2)
del y2m2

#ell, m = 2, -1
y2m1 = .5*np.sqrt(15./(2.*np.pi))*xmiydivr*zdivr
ylm_save(y2m1, 2, 1)
del y2m1

#ell, m = 2, 0
xdivrsq = xdivr*xdivr
ydivrsq = ydivr*ydivr
zdivrsq = zdivr*zdivr
y20 = .25*np.sqrt(5./np.pi)*(2.*zdivrsq-xdivrsq-ydivrsq)
ylm_save(y20, 2, 0)
del y20

#ell, m = 3, -3
xmiydivrcu = xmiydivr*xmiydivrsq
y3m3 = .125*np.sqrt(35./np.pi)*xmiydivrcu
ylm_save(y3m3, 3, 3)
del y3m3

#ell, m = 3, -2
y3m2 = .25*np.sqrt(105./(2.*np.pi))*xmiydivrsq*zdivr
ylm_save(y3m2, 3, 2)
del y3m2

#ell, m = 3, -1
y3m1 = .125*np.sqrt(21./np.pi)*(xmiydivr*(4.*zdivrsq-xdivrsq-ydivrsq))
ylm_save(y3m1, 3, 1)
del y3m1

#ell, m = 3, 0
y30 = .25*np.sqrt(7./np.pi)*(zdivr*(2.*zdivrsq-3.*xdivrsq-3.*ydivrsq))
ylm_save(y30, 3, 0)
del y30

#ell, m = 4, -4
xmiydivrft = xmiydivr*xmiydivrcu
y4m4 = .1875*np.sqrt(35./(2.*np.pi))*xmiydivrft
ylm_save(y4m4, 4, 4)
del y4m4

#ell, m = 4, -3
y4m3 = .375*np.sqrt(35./np.pi)*xmiydivrcu*zdivr
ylm_save(y4m3, 4, 3)
del y4m3

#ell, m = 4, -2
y4m2 = .375*np.sqrt(5./(2.*np.pi))*xmiydivrsq*(7.*zdivrsq-1)
ylm_save(y4m2, 4, 2)
del y4m2

#ell, m = 4, -1
y4m1 = .375*np.sqrt(5./np.pi)*xmiydivr*zdivr*(7.*zdivrsq-3.)
ylm_save(y4m1, 4, 1)
del y4m1

#ell, m = 4, 0
zdivrft = zdivrsq*zdivrsq
y40 = .1875*np.sqrt(1./np.pi)*(35.*zdivrft-30.*zdivrsq+3.)
ylm_save(y40, 4, 0)
del y40

#ell, m = 5, -5
xmiydivrfi = xmiydivr*xmiydivrft
y5m5 = (3./32.)*np.sqrt(77./np.pi)*xmiydivrfi
ylm_save(y5m5, 5, 5)
del y5m5

#ell, m = 5, -4
y5m4 = (3./16.)*np.sqrt(385./(2.*np.pi))*xmiydivrft*zdivr
ylm_save(y5m4, 5, 4)
del y5m4

#ell, m = 5, -3
y5m3 = (1./32.)*np.sqrt(385./np.pi)*xmiydivrcu*(9.*zdivrsq-1.)
ylm_save(y5m3, 5, 3)
del y5m3

#ell, m = 5, -2
zdivrcu = zdivr*zdivrsq
y5m2 = (1./8.)*np.sqrt(1155./(2.*np.pi))*xmiydivrsq*(3.*zdivrcu-zdivr)
ylm_save(y5m2, 5, 2)
del y5m2

#ell, m = 5, -1
y5m1 = (1./16.)*np.sqrt(165./(2.*np.pi))*xmiydivr*(21.*zdivrft-14.*zdivrsq+1.)
ylm_save(y5m1, 5, 1)
del y5m1

#ell, m = 5, 0
zdivrfi = zdivr*zdivrft
y50 = (1./16.)*np.sqrt(11./np.pi)*(63.*zdivrfi-70.*zdivrcu+15.*zdivr)
ylm_save(y50, 5, 0)
del y50
for ell in range(0, ell_max+1):
  for m in range(0, ell+1):
    print("ell, m = ", ell, m)
    #do one ylm at a time to save lots of accessing memory
    ylm = np.load('ylm_'+save_name+'_'+str(ell)+'_'+str(m)+'.npy') 
    for i in range(nbins):
      print("bin i = ", i)
      #where is radius in bin?
      rib = np.where((R >= bin_edges[i]) & (R < bin_edges[i+1]))
      ylm_on_shell = np.zeros((ld_one_d, ld_one_d, ld_one_d)) + 0j
      ylm_on_shell[rib] = ylm[rib]
      del rib
      ylm_transform_save(ylm_on_shell, ell, m, i)
      del ylm_on_shell
    if 'ylm' in globals():
      del ylm
    file_to_rm = 'ylm_'+save_name+'_'+str(ell)+'_'+str(m)+'.npy'
    call(["rm", file_to_rm])
    save_name = 'dens_t800'
kernel_name = 'demo'
file = save_name+'.fits.gz'

bin_bounds = np.load('bin_bounds_and_pixel_number_'+kernel_name+'.npy')
nbins = bin_bounds.shape[1] - 1

!wget https://users.flatironinstitute.org/~bburkhart/data/CATS/MHD/256/b.1p.32/t_800/dens_t800.fits.gz
!ls
hdulist = pyf.open(file)
data = hdulist[0].data.astype(np.float64) #pull data
data -= np.mean(data) #subtract out the mean

#FT OF SHIFTED DATA
normalized = False # do log normalization
if normalized:
  data = (np.log(data) - np.mean(np.log(data)))/np.std(np.log(data))
ft_data = np.fft.fftn(data)
#CONVOLUTION OF DATA AND SPH_KERNEL AT A GIVEN BIN
for l in range(0, ell_max+1, 1):
  for m in range(0,l+1, 1):
    for bin in range(0, nbins, 1):
      print("l, m, bin =", l, m, bin)
      #load ft of bsph_kernels
      bsph_kernel = np.load('YLMtilde_'+kernel_name+'_'+str(l)+'_'+str(m)+
                            '_bin_'+str(bin)+'.npy')
      conv = np.fft.ifftn(ft_data*bsph_kernel)
      del bsph_kernel
      np.save(save_name+'conv_data_kernel_'+kernel_name+'_'+str(l)+'_'+str(m)+
              '_bin_'+str(bin)+'.npy', conv)
      del conv
      import time
      #load alms (results of convolution at fixed ell and for all m, to combine).
zeta = np.zeros((ell_max + 1, nbins, nbins)) + 0j

start = time.time()
for l in range(0, ell_max+1, 1):
    for bin1 in range(0, nbins, 1):
        for bin2 in range(0, bin1, 1):
            for m in range(0,l+1, 1):
                #load b1 and b2 at this m

                ylm_b1 = np.load(save_name +'conv_data_kernel_'+kernel_name+'_'+str(l)+
                                  '_'+str(m)+'_bin_'+str(bin1)+'.npy').astype(np.complex128)
                ylm_b2 = np.load(save_name+'conv_data_kernel_'+kernel_name+'_'+str(l)+
                                  '_'+str(m)+'_bin_'+str(bin2)+'.npy').astype(np.complex128)
                #form half of sum at that m.
                ylm_b1 *= ylm_b2.conjugate()
                del ylm_b2

                #average spatially now.
                ylm_avg = np.sum(data*ylm_b1)
                del ylm_b1
                # add in negative m: so basically we are computing m's in sum term by
                # term but with matching between + and - m.
                if m > 0:
                    ylm_avg += ylm_avg.conjugate() 
                zeta[l, bin1, bin2] += ylm_avg #now accumulate so we get sum over m.

#----------------
# Symmetrization
#----------------
for l in range(0, ell_max+1, 1):
    zeta[l, :, :] *= ((-1)**l / np.sqrt(2 * l + 1))
    for bin1 in range(0, nbins, 1):
        for bin2 in range(bin1, nbins, 1): 
            zeta[l, bin1, bin2] = zeta[l, bin2, bin1]

stop = time.time()

print(f"Full 3PCF Calculation Took: {stop - start:0.2f} seconds")