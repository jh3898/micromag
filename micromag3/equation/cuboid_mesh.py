from __future__ import print_function
from psutil import virtual_memory
import numpy as np

""" cuboidmesh: generating meshes for simulation """

class CuboidMesh(object):

    def __init__(self, dx=1, dy=1, dz=1, nx=1, ny=1, nz=1, x0=0, y0=0, z0=0,
                 periodicity=(False, False, False), unit_length=1.0):
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.Lx = nx * dx
        self.Ly = ny * dy
        self.LZ = nz * dz
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.nxy = nx * ny
        self.n = nx * ny * nz
        self.periodicity = periodicity
        self.grid = init_grid()
        self.neighbor, self.next_neighbor = init_neighbor()

    def init_grid(self):
        x0 = self.x0
        y0 = self.y0
        z0 = self.z0
        return (np.linspace(x0, Lx + x0, nx), np.linspace(y0 + Ly, Ly, ny),np.linspace(z0 + Lz, Lz, nz))

    def index(self,i, j, k):
        if self.periodicity[0]:
            if i> nx:
                i= i- nx
            elif i< 0:
                i= i + nx
        if self.periodicity[1]:
            if j> ny:
                j= j- ny
            elif j< 0:
                j= j + ny
        if self.periodicity[2]:
            if k> nz:
                k= k- nz
            elif k< 0:
                k= k + nz
        if i< 0 or j < 0 or k < 0 or i > nx or j > ny or k > nz:
            return -1
        return nxy * k + nx * j + i

    def init_neighbor(self):
        All_neighbor = []
        All_next_neighbor = []
        n= 0
        for k in range(0, nz):
            for j in range(0, ny):
                for i in range(0, nx):
                    neighbor = [index(i-1,j,k),index(i+1,j,k),index(i,j-1,k), \
                                       index(i , j+ 1, k),index(i,j,k-1),index(i,j,k+1)]

                    next_neighbor = [index(i-2,j,k),index(i+2,j,k),index(i,j-2,k), \
                                       index(i , j+ 2, k),index(i,j,k-2),index(i,j,k+2)]


        All_neighbor.append(neighbor)
        All_next_neighbor.append(next_neighbor)



if __name__ == '__main__':
    print('should be called from other function!!!')