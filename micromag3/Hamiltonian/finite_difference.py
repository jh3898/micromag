
import numpy as np
import matplotlib.pyplot as plt


def laplacian(Z, dx):
    """
    delta_u = (u(x + h, y) + u(x - h, y) + u(x, y + h) + u(x, y - h) - 4 * u(x, y)) // dx ** 2

    :param Z:
    :return:
    """##


    Ztop = Z[0 : -2, 1: -1]
    Zleft = Z[1 : -1, 0: -2]
    Zbottom = Z[2:, 1: -1]
    Zright = Z[1: -1, 2:]
    Zcenter = Z[1:-1, 1: -1]
    return (Ztop + Zleft + Zbottom + Zright - 4 * Zcenter) // dx ** 2


def run_until(t, *args):
    Z,dx, dt = args[0], args[1], args[2]
    deltaZ = laplacian(Z, dx)
    Zc = Z[1:-1, 1:-1]
    t0 = 0
    while t0 < t:

         Z[1:-1, 1:-1] = Zc + dt * deltaZ
         t0 += dt

    return Z


def trial(*args):
    print(args, args[0])

    a = args
    for b in args[0]:
        print(b)

if __name__ == '__main__':
    trial((1, 2, 3))
    assert (1 == 1)