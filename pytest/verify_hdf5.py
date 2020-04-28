import numpy as np
import h5py
import struct
import os

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

def read_sac_usgs(fname):
    x = []
    y = []
    z = []
    t = []
    sta = open(fname,'r')
    lines = sta.readlines()
    skip = 0
    flag = 0
    count = 0

    for line in lines:
        if line[0] == '#':
            skip += 1
            continue
        elif flag == 0:
            n = len(lines) - skip
            flag = 1
            x = np.zeros(n)
            y = np.zeros(n)
            z = np.zeros(n)
            t = np.zeros(n)
            
        t[count] = np.float32(line.split()[0])
        x[count] = np.float32(line.split()[1])
        y[count] = np.float32(line.split()[2])
        z[count] = np.float32(line.split()[3])
        count += 1
    sta.close()
    return (np.array(t), np.array(x), np.array(y), np.array(z))

def read_sac_hdf5(fname, staname):
    fid = h5py.File(fname, 'r')
    grp = fid[staname]
    npts = grp['NPTS'][0]
    nsew = grp['ISNSEW'][0]
    delta = fid['DELTA'][0]
    
    if nsew == 1:
        data_names = ['EW', 'NS', 'UP']
        has_nsew = 1
    else:
        data_names = ['X', 'Y', 'Z']
        
    x=np.array(grp[data_names[0]][0:npts])
    y=np.array(grp[data_names[1]][0:npts])
    z=np.array(grp[data_names[2]][0:npts])
    
    if npts <= 0 or delta <=0:
        print('Error with station data, npts = ', npts, ', delta = ', delta)
    origintime = fid['ORIGINTIME'][0]
    t = np.linspace(0, delta*(npts-1), npts)
        
    fid.close()
    return (t, x, y, z)

def read_sw4img(fname):
    img = open(fname,'rb')
    prec = struct.unpack('i', img.read(4))[0]
    npatch = struct.unpack('i', img.read(4))[0]
    t = struct.unpack('d', img.read(8))[0]
    plane = struct.unpack('i', img.read(4))[0]
    xi = struct.unpack('d', img.read(8))[0]
    mode = struct.unpack('i', img.read(4))[0]
    grid_info = struct.unpack('i', img.read(4))[0]
    # Skip create time
    img.read(25)
    
    #print("prec=%d, npatch=%d, t=%f, plane=%d, xi=%f, mode=%d, grid_info=%d" % \
    #      (prec, npatch, t, plane, xi, mode, grid_info))
    
    # Grid size info
    h = np.zeros(npatch, dtype=np.float64)
    zmin = np.zeros(npatch, dtype=np.float64)
    i = np.zeros(npatch, dtype=np.int)
    ni = np.zeros(npatch, dtype=np.int)
    j = np.zeros(npatch, dtype=np.int)
    nj = np.zeros(npatch, dtype=np.int)
    nelem = 0
    for u in range (0, npatch):
        h[u] = struct.unpack('d', img.read(8))[0]
        zmin[u] = struct.unpack('d', img.read(8))[0]
        i[u] = struct.unpack('i', img.read(4))[0]
        ni[u] = struct.unpack('i', img.read(4))[0]
        j[u] = struct.unpack('i', img.read(4))[0]
        nj[u] = struct.unpack('i', img.read(4))[0]
        nelem += ni[u] * nj[u]
        
        #print("patch %d: h=%f, zmin=%f, i=%d, ni=%d, j=%d, nj=%d" % \
        #      (u, h[u], zmin[u], i[u], ni[u], j[u], nj[u]))

    # Read all patch data
    if prec == 4:
        pdata = struct.unpack(str(nelem)+'f', img.read(prec*nelem))
    elif prec == 8:
        pdata = struct.unpack(str(nelem)+'d', img.read(prec*nelem))
    #print("patch %d: min=%e, max=%e" % (u, np.min(pdata), np.max(pdata)))
    img.close()
    return np.array(pdata)

def read_sw4img_hdf5(fname):
        
    image_h5      = h5py.File(fname, 'r')
    npatch        = image_h5['npatch'][0]
    time          = image_h5['time'][0]
    plane         = image_h5['plane'][0]
    coordinate    = image_h5['coordinate'][0]
    grid_info     = image_h5['gridinfo'][0]
    mode          = image_h5['mode'][0]
    creation_time = str(image_h5.attrs['creationtime'], 'utf-8')

    grid_sizes    = image_h5['grid_size']
    zmins         = image_h5['zmin']
    nis           = image_h5['ni']
    njs           = image_h5['nj']
    
    readz = False
    has_grid = False
    if 'grid' in image_h5.keys():
        has_grid = True
    
    pdata = image_h5['patches']

    return np.array(pdata)


def verify(pytest_dir, tolerance):
    ref_dir = pytest_dir + '/hdf5/loh1-h100-mr-1/'
    hdf5_dir = os.getcwd() + '/loh1-h100-mr-1-hdf5/'
    verify = True
    nsta = 0
    for i in range(1,11):
        sta_name = 'sta%02d' % i
        #print(sta_name)
        usgs_fname = ref_dir + sta_name + '.txt'
        hdf5_fname = hdf5_dir + 'sta.h5'
        usgs_t, usgs_x, usgs_y, usgs_z = read_sac_usgs(usgs_fname)
        hdf5_t, hdf5_x, hdf5_y, hdf5_z = read_sac_hdf5(hdf5_fname, sta_name)
        if np.max(hdf5_t-usgs_t) > tolerance or np.min(hdf5_t-usgs_t) < -tolerance:
            verify = False
            print ("Station [%s] time data not match!" % sta_name)
        if np.max(hdf5_x-usgs_x) > tolerance or np.min(hdf5_x-usgs_x) < -tolerance:
            verify = False
            print ("Station [%s] x data not match!" % sta_name)
            print(hdf5_x-usgs_x)
        if np.max(hdf5_y-usgs_y) > tolerance or np.min(hdf5_y-usgs_y) < -tolerance:
            verify = False
            print ("Station [%s] y data not match!" % sta_name)
        if np.max(hdf5_z-usgs_z) > tolerance or np.min(hdf5_z-usgs_z) < -tolerance:
            verify = False
            print ("Station [%s] z data not match!" % sta_name)  
        nsta += 1

    # if verify == 1:
    #     print ('All %d stations data match!' % nsta)

    nimg = 0
    for filename in os.listdir(ref_dir):
        if filename.endswith(".sw4img"):
            nimg += 1
            #print(filename)
            sw4img_file = ref_dir + filename
            h5img_file  = hdf5_dir + filename + '.h5'
            sw4_pdata = read_sw4img(sw4img_file)
            h5_pdata = read_sw4img_hdf5(h5img_file)
            if len(h5_pdata) != len(h5_pdata):
                print("Image sizes are diferent! %d/%d" % (len(sw4_pdata), len(h5_pdata)))
                verify = False
            else:
                if np.max(h5_pdata-sw4_pdata) > tolerance or np.min(h5_pdata-sw4_pdata) < -tolerance:
                    print ("Image data [%s] does not match!" % sw4img_file)        
                    verify = False
            if verify == False:
                break
                    
    # if verify == 1:
    #     print ('All %d images data match!' % nimg)
    return verify
