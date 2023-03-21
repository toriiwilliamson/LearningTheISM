import os
from tqdm import tqdm

for i in tqdm(range(500, 901, 50)):
    os.system(f'wget https://users.flatironinstitute.org/~bburkhart/data/CATS/MHD/256/b1p1/t_{i}/dens_t{i}.fits.gz')
    
os.system(f'mv dens_*.fits.gz Ms_0p7__Ma_0p7')


print("Finished First Round of Data Downloading")

for i in tqdm(range(500, 901, 50)):
    os.system(f'wget https://users.flatironinstitute.org/~bburkhart/data/CATS/MHD/256/b1p.1/t_{i}/dens_t{i}.fits.gz')
    
os.system(f'mv dens_*.fits.gz Ms_2p0__Ma_0p7')

print("Finished Second Round of Data Downloading")

for i in tqdm(range(500, 901, 50)):
    os.system(f'wget https://users.flatironinstitute.org/~bburkhart/data/CATS/MHD/256/b.1p1/t_{i}/dens_t{i}.fits.gz')
    
os.system(f'mv dens_*.fits.gz Ms_0p7__Ma_2p0')

print("Finished Third Round of Data Downloading")

for i in tqdm(range(500, 901, 50)):
    os.system(f'wget https://users.flatironinstitute.org/~bburkhart/data/CATS/MHD/256/b.1p.1/t_{i}/dens_t{i}.fits.gz')
    
os.system(f'mv dens_*.fits.gz Ms_2p0__Ma_2p0')


