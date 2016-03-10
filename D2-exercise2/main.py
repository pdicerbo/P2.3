import numpy as np
import matplotlib.pyplot as plt
from vars import *
from tools import *

def tree_sort(g, nblist, listbodies, oldcell):

    HOC = np.zeros(4)
    LLJ = np.zeros(nblist)
    sublist = np.zeros(nblist)
    g.itercell_tree += 1

    if g.itercell_tree > g.levorder:
        return

    levscan = g.levorder - g.itercell_tree

    for j in range(0, nblist):
        i = listbodies[j]
        key_of_point = np.int64(g.key_list[j])
        lpos = g.ndim * levscan

        #ju = int(bin(key_of_point)[4:2:-1],2) #extract_bits(key_of_point, lpos, 2) # subcell index!

        # MyExtractBits
        ju = int(bin(key_of_point)[2*g.levorder-lpos+2:2*g.levorder-lpos:-1],2)
        print("key_of_point = ", key_of_point, " lpos = ", lpos, " ju = ", ju)
        kprev = HOC[ju]
        LLJ[j] = kprev
        HOC[ju] = j

    # to complete...
    # for j in range(0, nsubcell):
    #     k = HOC[jsub]

    #     if k == 0:
    #         continue

    #     nsubc = 0

    #     while k > 0:
    #         nsubc += 1
    #         i = listbodies[k]
    #         sublist[nsubc] = i
    #         k = LLJ[k]
        
    

def main():
    nbodsmax = 1000
    ncells = 2*nbodsmax
    nbodcell = nbodsmax + ncells
    ndim = 2
    nsubcell = 2**ndim
    # key_list = np.empty((0,))
    rmin = np.zeros(ndim)
    
    # others global: bottom, cellsize, pointers_of_tree,
    #                levbit, root, incells, intercell_tree,
    #                rsize, rmin

    load_data = np.loadtxt("tree.dat")
    L = (int)(load_data[0,0])
    rsize = L
    nbodies = (int)(load_data[0,1])
    print(load_data.shape, load_data[0,0], load_data[0,1])
    load_data = load_data[1:]
    print(load_data.shape, load_data[0,0], load_data[0,1])

    pos = np.zeros((2, nbodsmax + ncells))
    bottom = np.zeros((2, nbodcell))
    subindex = np.zeros(nbodies)
    bodlist = np.zeros(nbodies)
    cellsize = np.zeros(nbodcell)
    iback = np.zeros((2,nbodcell))
    pm1 = np.zeros((nsubcell,ndim))

    pos[0,:-(ncells + (nbodsmax-nbodies))] = load_data[:,0]
    pos[1,:-(ncells + (nbodsmax-nbodies))] = load_data[:,1]
    
    print("last particles:")
    print(load_data.shape, load_data[-1,0], load_data[-1,1])

    # subquadrant position
    pm1[0,0] = -1
    pm1[0,1] = -1
    
    pm1[1,0] = +1
    pm1[1,1] = -1

    pm1[2,0] = -1
    pm1[2,1] = +1

    pm1[3,0] = +1
    pm1[3,1] = +1


    # MORTON KEYS calculation
    levmorton = 10

    g.levorder = levmorton
    g.ndim = ndim
    cj = 2**levmorton / rsize

    levkey = 0

    for i in range(0, load_data.shape[0]):

        # ix = (int) (load_data[i,0] * cj)
        # iy = (int) (load_data[i,1] * cj)

        ix = (int) (pos[0,i] * cj)
        iy = (int) (pos[1,i] * cj)
        tmp = 0
        while levkey < levmorton:
            bitx = bit_test(ix, levkey)
            bity = bit_test(iy, levkey)

            if bitx:
                tmp = bit_set(tmp, 2*levkey, 1)

            if bity:
                tmp = bit_set(tmp, 1 + 2*levkey, 1) 
            
            levkey += 1

        g.key_list = np.copy(np.append(g.key_list, tmp))
        levkey = 0

    print("key_list = ", g.key_list)

    # initialization of some other variables
    incells = 1
    root = nbodies

    # initialization to 0, meaning there are no cells within the tree
    pointers_of_tree = np.zeros((nsubcell, ncells))

    for j in range(0, ndim):
        pos[j, root] = rmin[j] + 0.5*rsize
        bottom[j, root] = rmin[j]

    tree_sort(g, nbodies, bodlist, root)

        
def test_init():
    v = LinkedList(5.1)
    a = LinkedList(3.2)
    b = LinkedList(2.3)
    
    print("v = ", v)
    print("v.next = ", v.next)
    print("FIRST ADD")
    v.add(a)
    print("v = ", v)
    print("v.next = ", v.next)
    print("v.next.next = ", v.next.next)
    print("SECOND ADD")
    v.add(b)
    print("v = ", v)
    print("v.next = ", v.next)
    print("v.next.next = ", v.next.next)

if __name__ == "__main__":
    g = GlobVars(0);
    main()
    print("Done!")
