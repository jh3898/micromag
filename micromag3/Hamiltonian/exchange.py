#import micromag3
from common.constant import mu_0, Ms_inv

class Exchange(object):
    """compute the exchange field"""

    def __init__(self,A,name= 'Exchange'):
        self.A= A
        self.name= name
        self.jac= True

    def compute_field(self,t):

        exch()


    def exch(self,m, field, energy, Ms_inv, A, dx, dy, dz, n, ngbs):
        ax = 1 / dx ** 2
        ay = 1 / dy ** 2
        az = 1 / dz ** 2
        for i in range(0,n):
            num= ngbs(i*6)
            ix= 3*i
            iy= 3*i+1
            iz= 3*i+2
            avr_x = ax*(m[ix+ngbs[i*6+1]] -2*m[ix] +m[ix + ngbs[i*6+2]])+ \
                    ay*(m[ix+ngbs[i*6+3]] -2*m[ix] +m[ix + ngbs[i*6+4]])+ \
                    az *(m[ix + ngbs[i * 6 +5]] - 2 * m[ix] + m[ix + ngbs[i * 6 +6]])
            avr_y = ax*(m[i*3+1+ngbs[i*6+1]] -2*m[3*i+1] +m[3*i +1+ ngbs[i*6+2]])+ \
                    ay*(m[i*3+1+ngbs[i*6+3]] -2*m[3*i+1] +m[3*i + ngbs[i*6+4]])+ \
                    az *(m[i * 3 + ngbs[i * 6 +5]] - 2 * m[3 * i] + m[3 * i + ngbs[i * 6 +6]])

























