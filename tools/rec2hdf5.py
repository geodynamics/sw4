import sys
import h5py
import numpy as np
from datetime import datetime

# filename='HF_M7.0_3DTOPO_VSMIN250_H6.25MR.sw4input'
# filename='hayward-att-h100-cori.in'

filename = sys.argv[1]
outname  = filename+'.h5'

with open(filename) as f:
    lines = f.readlines()

outfile = h5py.File(outname,'w')

xyz=np.array([0.0,0.0,-1.0])

i = 0
is_nsew = 0
for line in lines:
    # line = line.lower()
    if line.startswith("rec") or line.startswith("sac") :
        has_sta = 0
        # Remove "sac/rec" and "\n"
        line = line[4:-1]
        llist = line.split(' ')
        # print(llist)
        nsew=0
        # Use sta as the station name
        for pairs in llist:
            kv = pairs.split('=')
            if kv[0] == "sta":
                # print(kv[1])
                grp = outfile.create_group(kv[1])
                has_sta = 1
            if kv[0] == "nsew" and kv[1] == "1":
                nsew = 1
                is_nsew = 1

        # If sta is not found, use file as station name
        if has_sta == 0:
            for pairs in llist:
                kv = pairs.split('=')
                if kv[0] == "file":
                    # print(kv[1])
                    grp = outfile.create_group(kv[1])
                    has_sta = 1

        # No sta or file is given
        if has_sta == 0:
            grp_name = "Noname " + str(i)
            i += 1
            grp = outfile.create_group(grp_name)

        downsample = 1
        for pairs in llist:
            kv = pairs.split('=')
            if kv[0] == "depth" or kv[0] == "topodepth" or kv[0] == "z":
                xyz[2] = float(kv[1])
            elif kv[0] == "x" :
                xyz[0] = float(kv[1])
            elif kv[0] == "y" :
                xyz[1] = float(kv[1])
            elif kv[0] == "lat" :
                xyz[0] = float(kv[1])
            elif kv[0] == "lon" :
                xyz[1] = float(kv[1])
            elif kv[0] == "downsample" :
                downsample = int(kv[1])
            # elif kv[0] == "  " or kv[0] == "hdf5file" or kv[0] == "hdf5format" or kv[0] == "usgsformat" or kv[0] == "sacformat" or kv[0] == "sta" or kv[0] == "file" or kv[0] == "variables" or kv[0] == "nsew":
            #     continue
            # else:
            #     print("Ignored cmd:[", kv[0], "]")

        dset = grp.create_dataset('ISNSEW', (1,), dtype='i4')
        dset[0] = nsew
        if nsew == 0:
            dset = grp.create_dataset('STX,STY,STZ', (3,), dtype='f8')
        else:
            dset = grp.create_dataset('STLA,STLO,STDP', (3,), dtype='f8')

        dset[:] = xyz

dset = outfile.create_dataset('DOWNSAMPLE', (1,), dtype='i4')
dset[0] = downsample

# outfile.attrs.create('UNIT', 'm')
outfile.attrs["UNIT"] = np.string_("m")

now = datetime.now()
dt_string = now.strftime("%Y-%m-%dT%H:%M:%S.0")
# outfile.attrs.create('DATETIME', dt_string)
outfile.attrs["DATETIME"] = np.string_(dt_string)

outfile.close()
