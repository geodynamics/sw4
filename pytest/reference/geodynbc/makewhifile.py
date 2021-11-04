#!/usr/bin/env python3

# Make whi file for geodynbc input
# From sw4 run to create whi records
# Number of time steps = 1000 dt: 0.0132681


def main():
    fout = open("loh1-h100-mr.whi", "w")
    n = 4
    h = 50
    sx = 15000
    sy = 15000
    sz = 100
    whidir = "loh1-h100-mr-whi/"
    steps = 256
    dt = 0.0132681

    # Assemble master file in memory
    text = []
    # Sides 1 & 2
    for s in range(1, 3):  # side
        for k in range(n + 1):  # z
            for j in range(n + 1):  # y
                with open(whidir + "side%d_z%d_y%d.txt" % (s, k + 1, j + 1)) as f:
                    text.extend(f.readlines()[13 : (13 + steps)])
    # Sides 3 & 4
    for s in range(3, 5):
        for k in range(n + 1):  # z
            for i in range(n + 1):  # x
                with open(whidir + "side%d_z%d_x%d.txt" % (s, k + 1, i + 1)) as f:
                    text.extend(f.readlines()[13 : (13 + steps)])
    # Sides 5 & 6
    for s in range(5, 7):
        for j in range(n + 1):  # y
            for i in range(n + 1):  # x
                with open(whidir + "side%d_y%d_x%d.txt" % (s, j + 1, i + 1)) as f:
                    text.extend(f.readlines()[13 : (13 + steps)])

    # Write WHI file
    fout.write(
        "grid faces=6 stepsize=%d nx=%d ny=%d nz=%d x0=%d y0=%d z0=%d adjust=1\n"
        % (h, n + 1, n + 1, n + 1, sx - (h * n / 2), sy - (h * n / 2), sz - (h * n / 2))
    )
    fout.write("time timestep=%f nsteps=%d\n" % (dt, steps))
    fout.write("begindata\n")
    # Step
    for t in range(steps):
        #
        for l in range(6 * (n + 1) * (n + 1)):
            # Just write out the x, y, z (no time in space-delimited column 1)
            fout.write(text[t + l * steps].partition(" ")[2])

    fout.close()


if __name__ == "__main__":
    main()
