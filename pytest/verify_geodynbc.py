import numpy as np
import struct
import os

def whiwriter(pytest_dir):
    whidir = pytest_dir + "/geodynbc/loh1-h100-mr-whi-1/"
    fout = open(pytest_dir + "/geodynbc/loh1-h100-mr.whi", "w")
    n = 11
    h = 50
    sx = 15000
    sy = 15000
    sz = 250
    steps = 256
    dt = 0.0132681

    # Assemble master file in memory
    text = []
    # Sides 1 & 2
    for s in range(1, 3):  # side
        for k in range(n):  # z
            for j in range(n):  # y
                with open(whidir + "side%d_z%d_y%d.txt" % (s, k + 1, j + 1)) as f:
                    text.extend(f.readlines()[13 : (13 + steps)])
    # Sides 3 & 4
    for s in range(3, 5):
        for k in range(n):  # z
            for i in range(n):  # x
                with open(whidir + "side%d_z%d_x%d.txt" % (s, k + 1, i + 1)) as f:
                    text.extend(f.readlines()[13 : (13 + steps)])
    # Sides 5 & 6
    for s in range(5, 7):
        for j in range(n):  # y
            for i in range(n):  # x
                with open(whidir + "side%d_y%d_x%d.txt" % (s, j + 1, i + 1)) as f:
                    text.extend(f.readlines()[13 : (13 + steps)])

    # Write WHI file
    fout.write(
        "grid faces=6 stepsize=%d nx=%d ny=%d nz=%d x0=%d y0=%d z0=%d adjust=1\n"
        % (h, n, n, n, sx - (h * (n-1) / 2), sy - (h * (n-1) / 2), sz - (h * (n-1) / 2))
    )
    fout.write("time timestep=%f nsteps=%d\n" % (dt, steps))
    fout.write("begindata\n")
    # Step
    for t in range(steps):
        #
        for l in range(6 * n * n):
            # Just write out the x, y, z (no time in space-delimited column 1)
            fout.write(text[t + l * steps].partition(" ")[2])

    fout.close()

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

def verify(pytest_dir, tolerance):
    reference_dir = pytest_dir + '/reference/geodynbc/loh1-h100-mr-geodynbc-1/'
    geodynbc_dir = pytest_dir + '/geodynbc/loh1-h100-mr-geodynbc-1/'
    verify = True
    nsta = 0
    for i in range(1,11):
        sta_name = 'sta%02d' % i
        #print(sta_name)
        usgs_fname = reference_dir + sta_name + '.txt'
        geodynbc_fname = geodynbc_dir + sta_name + '.txt'
        usgs_t, usgs_x, usgs_y, usgs_z = read_sac_usgs(usgs_fname)
        geodynbc_t, geodynbc_x, geodynbc_y, geodynbc_z = read_sac_usgs(geodynbc_fname)
        if np.max(geodynbc_t-usgs_t) > tolerance or np.min(geodynbc_t-usgs_t) < -tolerance:
            verify = False
            print ("Station [%s] time data not match!" % sta_name)
            return False
        if np.max(geodynbc_x-usgs_x) > tolerance or np.min(geodynbc_x-usgs_x) < -tolerance:
            verify = False
            print ("Station [%s] x data not match!" % sta_name)
            print(geodynbc_x-usgs_x)
            return False
        if np.max(geodynbc_y-usgs_y) > tolerance or np.min(geodynbc_y-usgs_y) < -tolerance:
            verify = False
            print ("Station [%s] y data not match!" % sta_name)
            return False
        if np.max(geodynbc_z-usgs_z) > tolerance or np.min(geodynbc_z-usgs_z) < -tolerance:
            verify = False
            print ("Station [%s] z data not match!" % sta_name)  
            return False
        nsta += 1

    # if verify == 1:
    #     print ('All %d stations data match!' % nsta)

    return verify
