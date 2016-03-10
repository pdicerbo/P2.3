# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 09:14:57 2015
&&
Created on Mon Feb 16 15:12:14 2015

@author: nicola
"""

import numpy as np

class MortonVars:
    def __init__(self):
        self.x_pos = np.zeros((0,),dtype=float)
        self.y_pos = np.zeros(self.x_pos.shape)
        self.quad = np.zeros((2,4),dtype=float)
        self.corner = np.zeros((2,4),dtype=float)
        self.levgrid = 0
        self.levscan = int(4)
        self.side = int(16)
        self.key_list = np.empty((0,)) 
        return


class HilbertVars:
    def __init__(self, lev):
        self.n_points = int(4**lev) #int(4096)
        self.rotation_table = np.array([3, 0, 0, 1], dtype=np.int32)
        self.sense_table = np.array([-1, 1, 1, -1], dtype=np.int32)
        self.quad_table = np.zeros((4,2,2), dtype=np.int32) #(0:3,0:1,0:1)
        self.quad_pos = np.zeros((2,4), dtype = np.float)
        self.corner = np.zeros((2,4), dtype = np.float)
        self.x_pos = np.zeros((self.n_points,),dtype=np.float)
        self.y_pos = np.zeros(self.x_pos.shape,dtype=np.float)
        self.levgrid = np.int32(0)
        self.levscan = np.int32(lev) #np.int32(0)
        self.side = np.int32(0)
        self.key_list = np.empty((0,))

class GlobVars:
    def __init__(self, i):
        self.itercell_tree = (int)(i)
        self.levorder = 0
        self.key_list = np.empty((0,))
        self.ndim = 0
        self.nsubcell = 0
        self.incells = 0
        self.ncells = 0
        self.root = 0
        self.pos = np.empty((0,))
        self.pm1 = np.empty((0,))
        self.iback = np.empty((0,))
        self.bottom = np.empty((0,))
        self.cellsize = np.empty((0,))
        self.pointers_of_tree = np.empty((0,))

        
class LinkedList:
    def __init__(self, head = None, next = None):
        self.head = head
        self.next = next

    def add(self, a):
        if self.next == None:
            self.next = a
        else:
            self.next.add(a)
            
    def __str__(self):
        return str(self.head)
