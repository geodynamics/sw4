import numpy as np
import h5py
import sys

if len(sys.argv) != 3:
    print('Need both input and output file name')
    sys.exit()

infile  = sys.argv[1]
outfile = sys.argv[2]

print('Input SRF file:', infile)
print('Output SRF-HDF5 file:', outfile)

f = open(infile, 'r')
lines = f.readlines()
f.close()

idx = 0
# Header
version = np.float32(lines[idx])
idx += 1
print('Version:', version)

fields = lines[idx].split()
if fields[0] != 'PLANE':
    print('Expecting PLANE on second line, exit on error...')
    sys.exit()

nplane = np.int32(fields[1])
print('Plane:', nplane)
idx += 1

seg_elon=[]
seg_elat=[]
seg_nstk=[]
seg_ndip=[]
seg_len=[]
seg_wid=[]
seg_stk=[]
seg_dip=[]
seg_dtop=[]
seg_shyp=[]
seg_dhyp=[]
for i in range(0, nplane):
    # first line
    fields = lines[idx].split()
    seg_elon.append(np.float32(fields[0]))
    seg_elat.append(np.float32(fields[1]))
    seg_nstk.append(np.int32(fields[2]))
    seg_ndip.append(np.int32(fields[3]))
    seg_len.append(np.float32(fields[4]))
    seg_wid.append(np.float32(fields[5]))
    idx += 1

    # second line
    fields = lines[idx].split()
    seg_stk.append(np.float32(fields[0]))
    seg_dip.append(np.float32(fields[1]))
    seg_dtop.append(np.float32(fields[2]))
    seg_shyp.append(np.float32(fields[3]))
    seg_dhyp.append(np.float32(fields[4]))
    idx += 1

for i in range(0, nplane):
    print('Seg #', i, 'elon =', seg_elon[i], 'elat =', seg_elat[i], 'nstk =', seg_nstk[i], 'ndip =', seg_ndip[i], 'len =', seg_len[i], 'wid =', seg_wid[i])
    print('      ', 'stk =', seg_stk[i], 'dip =', seg_dip[i], 'dtop =', seg_dtop[i], 'shyp =', seg_shyp[i], 'dhyp =', seg_dhyp[i])

seg_type = np.dtype([('ELON', 'f4'), ('ELAT', 'f4'), ('NSTK', 'i4'), ('NDIP', 'i4'), ('LEN', 'f4'), ('WID', 'f4'), ('STK', 'f4'), ('DIP', 'f4'), ('DTOP', 'f4'), ('SHYP', 'f4'), ('DHYP', 'f4'), ])

points_type = np.dtype([('LON', 'f4'), ('LAT', 'f4'), ('DEP', 'f4'), ('STK', 'f4'), ('DIP', 'f4'), ('AREA', 'f4'), ('TINIT', 'f4'), ('DT', 'f4'), ('VS', 'f4'), ('DEN', 'f4'), ('RAKE', 'f4'), ('SLIP1', 'f4'), ('NT1', 'i4'), ('SLIP2', 'f4'), ('NT2', 'i4'), ('SLIP3', 'f4'), ('NT3', 'i4'), ])

h5file = h5py.File(outfile,'w')

h5file.attrs.create('VERSION', version)
# h5file.attrs.create('NPLANE', nplane)

attr_data = []

for i in range(0, nplane):
    attr_data.append( np.array([(seg_elon[i],seg_elat[i],seg_nstk[i],seg_ndip[i], seg_len[i],seg_wid[i], seg_stk[i],seg_dip[i],seg_dtop[i],seg_shyp[i],seg_dhyp[i])], dtype = seg_type))

h5file.attrs.create('PLANE', attr_data, (nplane,), seg_type)

# POINTS
fields = lines[idx].split()
if fields[0] != 'POINTS':
    print('Expecting POINTS, exit on error...')
    sys.exit()

npoints= np.int32(fields[1])
print('Points:', npoints)
idx += 1

points_data = []
sr1_data = []
nsr1 = 0
points_dset = h5file.create_dataset('POINTS', (npoints,), points_type)
for i in range(0, npoints):
    fields = lines[idx].split()
    lon    = np.float32(fields[0])
    lat    = np.float32(fields[1])
    dep    = np.float32(fields[2])
    stk    = np.float32(fields[3])
    dip    = np.float32(fields[4])
    area   = np.float32(fields[5])
    tinit  = np.float32(fields[6])
    dt     = np.float32(fields[7])
    vs     = np.float32(fields[8])
    den    = np.float32(fields[9])
    idx += 1

    fields = lines[idx].split()
    rake   = np.float32(fields[0])
    slip1  = np.float32(fields[1])
    nt1    = np.int32(fields[2])
    slip2  = np.float32(fields[3])
    nt2    = np.int32(fields[4])
    slip3  = np.float32(fields[5])
    nt3    = np.int32(fields[6])
    idx += 1

    points_dset[i, ...] = np.array([(lon,lat,dep,stk,dip,area,tinit,dt,vs,den,rake,slip1,nt1,slip2,nt2,slip3,nt3,)], dtype = points_type);

    if nt1 > 0:
        nread = 0
        while (nread < nt1):
            fields = lines[idx].split()
            for j in range (0, len(fields)):
                sr1_data.append(np.float32(fields[j]))

            nread += len(fields)
            idx += 1
        nsr1 += nread

    # Skip sr1 and sr2
    if nt2 > 0:
        nread = 0
        while (nread < nt2):
            fields = lines[idx].split()
            nread += len(fields)
            idx += 1

    if nt3 > 0:
        nread = 0
        while (nread < nt3):
            fields = lines[idx].split()
            nread += len(fields)
            idx += 1

sr1_dset = h5file.create_dataset('SR1', (nsr1,), 'f4')
sr1_dset[...] = sr1_data

    
# SR1


h5file.close()
