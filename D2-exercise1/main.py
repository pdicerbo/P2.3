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
        # print("vr.levgrid = ", vr.levgrid)
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

def make_hilbert_key(vr):

    vr.quad_table[0,0,0] = 0
    vr.quad_table[0,1,0] = 3

    vr.quad_table[0,0,1] = 1
    vr.quad_table[0,1,1] = 2

    vr.quad_table[1,0,0] = 1
    vr.quad_table[1,1,0] = 0

    vr.quad_table[1,0,1] = 2
    vr.quad_table[1,1,1] = 3
    
    vr.quad_table[2,0,0] = 2
    vr.quad_table[2,1,0] = 1

    vr.quad_table[2,0,1] = 3
    vr.quad_table[2,1,1] = 0

    vr.quad_table[3,0,0] = 3
    vr.quad_table[3,1,0] = 2

    vr.quad_table[3,0,1] = 0
    vr.quad_table[3,1,1] = 1

    for i in range(0, len(vr.x_pos)):

        x = vr.x_pos[i]
        y = vr.y_pos[i]
        
        rotation = 0
        sense = 1
        
        k = (int)(vr.side * 0.5)
        num = 0

        while k > 0:
            xbit = (int)(x / k)
            ybit = (int)(y / k)
            x -= xbit * k
            y -= ybit * k
            quadh = vr.quad_table[rotation, xbit, ybit]

            if sense == -1:
                num += k*k*(3-quadh)
            else:
                num += k*k*quadh

            rotation += vr.rotation_table[quadh]
            if rotation >= 4:
                rotation -= 4

            sense *= vr.sense_table[quadh]
            k = (int)(k * 0.5)

        vr.key_list = np.copy(np.append(vr.key_list, num))
        print(i, vr.x_pos[i], vr.y_pos[i], num)

def test2D(vr):

    for j in range(0, len(vr.x_pos)):

        sidelocal = 2**vr.levscan
        k = (int) (sidelocal * 0.5)

        ix  = (int)(vr.x_pos[j] * sidelocal / vr.side)
        iy  = (int)(vr.y_pos[j] * sidelocal / vr.side)

        num = (int)(0)
        levs = vr.levscan - 1

        while k > 0:
            rx = (int)(bit_test(ix, levs))
            ry = (int)(bit_test(iy, levs))
            quad = (3 * rx)^ry
            num += k * k * quad

            if ry == 0:
                if rx == 1:
                    ix = sidelocal - 1 - ix
                    iy = sidelocal - 1 - iy
                    
                swap = ix
                ix = iy
                iy = swap

            k = (int) (k * 0.5)
            levs -= 1

        vr.key_list = np.copy(np.append(vr.key_list, num))
        print(j, vr.x_pos[j], vr.y_pos[j], num)

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
    vr = Vars(4)
    # vr.levscan = 4
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
    # make_hilbert_key(vr)
    test2D(vr)
    # print(vr.key_list)

    print("plotting")
    
    plt.figure()
    plt.plot(vr.x_pos, vr.y_pos)
    plt.show()
    
    return
    

if __name__ == "__main__":
    print("kick")
    main()
    print("done")
