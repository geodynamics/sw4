#!/usr/bin/env python3

# Make rec lines to create whi file input
# e.g., rec x=15600 y=15800 z=0 file=sta01 usgsformat=1 sacformat=0
# source x=15000 y=15000 z=2000 mxy=1e18 t0=0.36 freq=16.6667 type=Gaussian
# n = number of nodes in box
# h, sx, sy, sz are from loh1-h100-mr-whi.in file

def main():
    f = open("whireclines.txt", "w")
    n = 11
    h = 50
    sx = 15000
    sy = 15000
    sz = 250

    # Sides 1 & 2
    for s in range(2):  # side
        for k in range(n):  # z
            for j in range(n):  # y
                f.write(
                    "rec x=%1d y=%1d z=%1d file=side%1d_z%1d_y%1d usgsformat=1 sacformat=0\n"
                    % (
                        sx - (h * (n-1) / 2) + 2 * s * (h * (n-1) / 2),  # trick for sides
                        sy - (h * (n-1) / 2) + (j * h),
                        sz - (h * (n-1) / 2) + (k * h),
                        s + 1,
                        k + 1,
                        j + 1,
                    )
                )

    # Sides 3 & 4
    for s in range(2):  # side
        for k in range(n):  # z
            for i in range(n):  # x
                f.write(
                    "rec x=%1d y=%1d z=%1d file=side%1d_z%1d_x%1d usgsformat=1 sacformat=0\n"
                    % (
                        sx - (h * (n-1) / 2) + (i * h),
                        sy - (h * (n-1) / 2) + 2 * s * (h * (n-1) / 2),  # trick for sides
                        sz - (h * (n-1) / 2) + (k * h),
                        s + 3,
                        k + 1,
                        i + 1,
                    )
                )

    # Sides 5 & 6
    for s in range(2):  # side
        for j in range(n):  # y
            for i in range(n):  # x
                f.write(
                    "rec x=%1d y=%1d z=%1d file=side%1d_y%1d_x%1d usgsformat=1 sacformat=0\n"
                    % (
                        sx - (h * (n-1) / 2) + (i * h),
                        sy - (h * (n-1) / 2) + (j * h),
                        sz - (h * (n-1) / 2) + 2 * s * (h * (n-1) / 2),  # trick for sides,
                        s + 5,
                        j + 1,
                        i + 1,
                    )
                )

    f.close()


if __name__ == "__main__":
    main()
