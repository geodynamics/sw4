import h5py
import numpy as np

# This file will be overwritten
fname = '/global/cscratch1/sd/houhun/h012_usgs_v21-1-alpha2.iter10.h5'

gmg_file = h5py.File(fname, 'r+')
block_names = ['vres25m', 'vres50m', 'vres125m', 'vres250m']
vs_loc = int(2)

niter = 10
kappa = 0.5
for dname in block_names:
    data_all = gmg_file['blocks'][dname][()]
    
    nskip = 0
    for iter in range(0, niter):
        data_smooth = np.copy(data_all)
        for i in range(1, data_all.shape[0]-1):
            for j in range(1, data_all.shape[1]-1):
                for k in range(1, data_all.shape[2]-1):
                    # skip invalid points
                    if data_all[i+1,j,k,vs_loc] < 0 or data_all[i-1,j,k,vs_loc] < 0 or \
                       data_all[i,j+1,k,vs_loc] < 0 or data_all[i,j-1,k,vs_loc] < 0 or \
                       data_all[i,j,k+1,vs_loc] < 0 or data_all[i,j,k-1,vs_loc] < 0 :
                        nskip += 1
                        continue
                    else:
                        # 0 rho, 1 vp, 2 vs
                        for c in range(0,3):
                            data_smooth[i,j,k,c] = (1-kappa)*data_all[i,j,k,c] + \
                                                    kappa / 6 * (data_all[i+1,j,k,c] + data_all[i-1,j,k,c] + \
                                                                 data_all[i,j+1,k,c] + data_all[i,j-1,k,c] + \
                                                                 data_all[i,j,k+1,c] + data_all[i,j,k-1,c])
        data_all = np.copy(data_smooth)
    
    gmg_file['blocks'][dname][:,:,:,:] = data_all[:,:,:,:]

    print('skipped', nskip, 'points', flush=True)
    
gmg_file.close()
