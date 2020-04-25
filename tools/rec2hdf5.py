#!/usr/bin/python
import sys
import h5py
import numpy as np

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
            elif kv[0] == "writeEvery" :
                grp.attrs.create("write every", int(kv[1]))
            elif kv[0] == "variables" :
                grp.attrs.create("variables", int(kv[1]))
            elif kv[0] == "usgsformat" or kv[0] == "sacformat" or kv[0] == "sta" or kv[0] == "file" or kv[0] == "nsew":
                continue
            else:
                print("Unsupported cmd:", kv[0], kv[1])

        if nsew == 0:
            grp.attrs.create('x,y,z', xyz, shape=xyz.shape)
        else:
            grp.attrs.create('lat,lon,depth', xyz, shape=xyz.shape)

outfile.attrs.create('is_nsew', is_nsew)
outfile.close()
