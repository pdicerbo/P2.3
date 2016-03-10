import numpy as np
import matplotlib.pyplot as plt
from vars import *
from tools import *

def tree_sort(g, nblist, listbodies, oldcell):

    # HOC = np.zeros(4)
    HOC = -1 * np.ones(4)
    LLJ = np.zeros(nblist)
    sublist = np.zeros(nblist)
    g.itercell_tree += 1

    if g.itercell_tree > g.levorder:
        return

    levscan = g.levorder - g.itercell_tree

    for j in range(0, nblist):
        i = listbodies[j]
        key_of_point = (int)(g.key_list[i])
        lpos = 2 * ( g.itercell_tree - 1) # g.ndim * levscan

        # MyExtractBits
        ju = int(format(key_of_point, 'b').zfill(20)[lpos:lpos+2],2)

        # print("key_of_point = ", key_of_point, " lpos = ", lpos, " ju = ", ju, "i = ", i)

        kprev = HOC[ju]
        LLJ[j] = kprev
        HOC[ju] = j

    # print("\n\n=====================================================\n\n")
    
    for jsub in range(0, g.nsubcell):
        k = HOC[jsub]

        # if k == 0:
        if k == -1:
            continue

        nsubc = 0

        # while k > 0:
        while k >= 0:
            i = listbodies[k]
            sublist[nsubc] = i
            k = LLJ[k]
            nsubc += 1
        
        if nsubc > 1:
            g.incells += 1
            if g.incells > g.ncells:
                print("error")
                return

            newcell = g.incells + g.root - 1
            g.pointers_of_tree[jsub, oldcell] = newcell
            g.cellsize[newcell] = g.cellsize[oldcell] * 0.5

            for m in range(0, g.ndim):
                g.pos[m, newcell] = g.pos[m,oldcell] + g.pm1[jsub,m]*0.5*g.cellsize[newcell]
                g.bottom[m, newcell] = g.pos[m,newcell] - 0.5 * g.cellsize[newcell]

            g.iback[0,newcell] = jsub+1
            g.iback[1,newcell] = oldcell

            tree_sort(g, nsubc, sublist, newcell)
            
        else:

            pbody = sublist[nsubc]
            g.pointers_of_tree[jsub, oldcell] = pbody
            g.iback[0, pbody] = jsub + 1
            g.iback[1, pbody] = oldcell

    g.itercell_tree -= 1

def print_cells(g):
    out = open("treecells.dat", "w")
    print("\n\tcells:\n")
    for j in range(g.incells):
        ind = j + g.root
        tosave = "%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n" % (g.bottom[0,ind], g.bottom[1,ind],
                  g.bottom[0,ind] + g.cellsize[ind], g.bottom[0,ind],
                  g.bottom[0,ind] + g.cellsize[ind], g.bottom[1,ind] + g.cellsize[ind],
                  g.bottom[0,ind], g.bottom[1,ind] + g.cellsize[ind])
        print(tosave[:-2])
        out.write(tosave)
    print("")
    out.close()

def sorti(g):
    tmp = np.sort(g.key_list)
    j = 0
    i = 0
    for el in tmp:
        while el != g.key_list[j]:
            j += 1
        g.subindex[i] = j
        i += 1
        j = 0

    print("g.subindex = ", g.subindex)
                
def main():
    nbodsmax = 1000
    g.ncells = 2*nbodsmax
    nbodcell = nbodsmax + g.ncells
    ndim = 2
    g.nsubcell = 2**ndim
    # key_list = np.empty((0,))
    rmin = np.zeros(ndim)
    
    load_data = np.loadtxt("tree.dat")
    L = (int)(load_data[0,0])
    rsize = L
    nbodies = (int)(load_data[0,1])
    print(load_data.shape, load_data[0,0], load_data[0,1])
    load_data = load_data[1:]
    print(load_data.shape, load_data[0,0], load_data[0,1])

    g.pos = np.zeros((2, nbodsmax + g.ncells))
    g.bottom = np.zeros((2, nbodcell))
    g.subindex = np.zeros(nbodies)
    bodlist = np.zeros(nbodies)
    g.cellsize = np.zeros(nbodcell)
    g.iback = np.zeros((2,nbodcell))
    g.pm1 = np.zeros((g.nsubcell,ndim))

    g.pos[0,:-(g.ncells + (nbodsmax-nbodies))] = load_data[:,0]
    g.pos[1,:-(g.ncells + (nbodsmax-nbodies))] = load_data[:,1]
    
    print("last particles:")
    print(load_data.shape, load_data[-1,0], load_data[-1,1])

    # subquadrant position
    g.pm1[0,0] = -1
    g.pm1[0,1] = -1
    
    g.pm1[1,0] = +1
    g.pm1[1,1] = -1

    g.pm1[2,0] = -1
    g.pm1[2,1] = +1

    g.pm1[3,0] = +1
    g.pm1[3,1] = +1

    for k in range(0, nbodies):
        bodlist[k] = k

    # MORTON KEYS calculation
    levmorton = 10

    g.levorder = levmorton
    g.ndim = ndim
    cj = 2**levmorton / rsize

    levkey = 0

    for i in range(0, load_data.shape[0]):

        ix = (int) (g.pos[0,i] * cj)
        iy = (int) (g.pos[1,i] * cj)
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
    g.incells = 1
    g.root = nbodies
    g.cellsize[g.root] = rsize
    
    # initialization to 0, meaning there are no cells within the tree
    g.pointers_of_tree = np.zeros((g.nsubcell, g.ncells))

    for j in range(0, ndim):
        g.pos[j, g.root] = rmin[j] + 0.5*rsize
        g.bottom[j, g.root] = rmin[j]

    tree_sort(g, nbodies, bodlist, g.root)
    print_cells(g)
    sorti(g)
        
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
