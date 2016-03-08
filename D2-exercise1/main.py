import numpy as np
import matplotlib.pyplot as plt
from hilbert_vars import Vars
from tools import *

def make_hilbert_grid(xgrid, ygrid, vr):

    if vr.levgrid + 1 >= vr.levscan:

        vr.x_pos = np.copy(xgrid)
        vr.y_pos = np.copy(ygrid)

        return
    else:
        print("vr.levgrid = ", vr.levgrid)
        vr.levgrid += 1
        
        n = len(xgrid)
        ntmp = (int) (n * 0.5 - 1)

        xsub = np.zeros(4 * n)
        ysub = np.zeros(4 * n)
        
        for i in range(0, 4):

            iad = i * n
            xswap = xgrid / 2
            yswap = ygrid / 2

            if i == 0:
                xtmp = np.copy(xswap)
                xswap = np.copy(yswap)
                yswap = vr.side * 0.5 - xtmp

                # reversing arrays
                xswap = xswap[::-1]
                yswap = yswap[::-1]
                
            if i == 3:
                ytmp = np.copy(yswap)
                yswap = np.copy(xswap)
                xswap = vr.side * 0.5 - ytmp

                # reversing arrays
                xswap = xswap[::-1]
                yswap = yswap[::-1]
            
            for j in range(0, n):
                xsub[j + iad] = xswap[j] + vr.corner[0, i] * vr.side * 0.5
                ysub[j + iad] = yswap[j] + vr.corner[1, i] * vr.side * 0.5

        make_hilbert_grid(xsub, ysub, vr)

        return


def make_morton_grid(xgrid, ygrid, vr):

    if vr.levgrid + 1 >= vr.levscan:

        vr.x_pos = np.copy(xgrid)
        vr.y_pos = np.copy(ygrid)

        return
    else:
        vr.levgrid += 1
        n = len(xgrid)

        xsub = np.zeros(4 * n)
        ysub = np.zeros(4 * n)

        xswap = xgrid / 2
        yswap = ygrid / 2

        for i in range(0, 4):
            iad = i * n

            for j in range(0, n):
                xsub[j + iad] = xswap[j] + vr.corner[0, i] * vr.side * 0.5
                ysub[j + iad] = yswap[j] + vr.corner[1, i] * vr.side * 0.5

        make_morton_grid(xsub, ysub, vr)

        return

def make_morton_key(vr):
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

    # WARNING: levscan cannot be greater than 6
    vr.levscan = 6
    vr.side = 16
    
    vr.quad_pos[0,0] = 0.5
    vr.quad_pos[1,0] = 0.5
    
    vr.quad_pos[0,1] = 0.5
    vr.quad_pos[1,1] = 1.5
    
    vr.quad_pos[0,2] = 1.5
    vr.quad_pos[1,2] = 1.5
    
    vr.quad_pos[0,3] = 1.5
    vr.quad_pos[1,3] = 0.5

    print(vr.quad_pos)

    xgrid = vr.quad_pos[0,:]*0.5 * vr.side
    ygrid = vr.quad_pos[1,:]*0.5 * vr.side
    print("start with xgrid = ", xgrid)
    vr.corner[0,0] = 0.
    vr.corner[1,0] = 0.

    vr.corner[0,1] = 0.
    vr.corner[1,1] = 1.

    vr.corner[0,2] = 1.
    vr.corner[1,2] = 1.
    
    vr.corner[0,3] = 1.
    vr.corner[1,3] = 0.
    
    print("starting grid calculation")

    make_hilbert_grid(xgrid, ygrid, vr)
    # make_morton_key(vr)

    print("plotting")
    
    plt.figure()
    plt.plot(vr.x_pos, vr.y_pos)
    plt.show()
    
    return
    

if __name__ == "__main__":
    vr = Vars()
    print("kick")
    main()
    print("done")
