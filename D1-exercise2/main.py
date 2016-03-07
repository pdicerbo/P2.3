import numpy as np
import matplotlib.pyplot as plt
from morton_vars import Vars
from tools import *

def make_grid(xgrid, ygrid, vr):

    if vr.levgrid + 1 >= vr.levscan:
        vr.x_pos = np.copy(xgrid)
        vr.y_pos = np.copy(ygrid)
        # print("len vr.x_pos = ", len(vr.x_pos))
        return
    else:
        vr.levgrid += 1
        n = len(xgrid)
        print("calculation for n = ", n)
        xsub = np.zeros(4 * n)
        ysub = np.zeros(4 * n)

        xswap = xgrid / 2
        yswap = ygrid / 2

        print(xgrid)
        print(ygrid)

        for i in range(0, 4):
            iad = i * n

            for j in range(0, n):
                xsub[j + iad] = xswap[j] + vr.corner[0, i] * vr.side * 0.5
                ysub[j + iad] = yswap[j] + vr.corner[1, i] * vr.side * 0.5

        make_grid(xsub, ysub, vr)

        return

def make_key(vr):
    leavmorton = 8
    cj = 2**leavmorton / vr.side

    levkey = 0
    
    for i in range(0, len(vr.x_pos)):
        ix = (int) (vr.x_pos[i] * cj)
        iy = (int) (vr.y_pos[i] * cj)
        tmp = 0
        while levkey < leavmorton:
            bitx = bit_test(ix, levkey)
            bity = bit_test(iy, levkey)

            if bitx:
                tmp = bit_set(tmp, 2*levkey, 1)

            if bity:
                tmp = bit_set(tmp, 1 + 2*levkey, 1) 
            
            levkey += 1

        vr.key_list = np.copy(np.append(vr.key_list, tmp))
        levkey = 0
    

def main():
    print("main section of the program")
    # print(vr.quad)

    vr.quad[0,0] = 0.5
    vr.quad[1,0] = 0.5
    
    vr.quad[0,1] = 1.5
    vr.quad[1,1] = 0.5
    
    vr.quad[0,2] = 0.5
    vr.quad[1,2] = 1.5
    
    vr.quad[0,3] = 1.5
    vr.quad[1,3] = 1.5

    print(vr.quad)

    xgrid = vr.quad[0,:]*0.5 * vr.side
    ygrid = vr.quad[1,:]*0.5 * vr.side

    vr.corner[0,0] = 0.
    vr.corner[1,0] = 0.

    vr.corner[0,1] = 1.
    vr.corner[1,1] = 0.

    vr.corner[0,2] = 0.
    vr.corner[1,2] = 1.
    
    vr.corner[0,3] = 1.
    vr.corner[1,3] = 1.

    print("starting grid calculation")

    make_grid(xgrid, ygrid, vr)
    make_key(vr)

    plt.figure()
    plt.plot(vr.x_pos, vr.y_pos)
    plt.show()
    
    return
    

if __name__ == "__main__":
    vr = Vars()
    print("Kick")
    main()
    print("Done")
