
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
from micromag.prep.sim import Sim
from micromag.prep.mesh import Mesh
from micromag.Hamiltonian.zeeman import Zeeman, FieldepZeeman

mu0 = 4 * np.pi * 1e-7


def init_m(pos):
    x = pos[0]
    if x <= 2:
        return (1, 0, 0)
    elif x >= 4:
        return  (0, 0, 1)
    else:
        return  (0, 1, 0)

def relax_system(mesh):

    sim = Sim(mesh, name= 'relax')

    sim.driver.alpha = 0.5
    sim.driver.gamma = 2.2e5
    sim.Ms = 8.0e5

    sim.set_m((1, 0.25, 0.1))
    mT = 0.001 /mu0
    field = [-24.6 * mT, 4.3 * mT, 0]
    zeeman = Zeeman(field, name= 'H')
    sim.add(zeeman)

    sim.relax(dt= 1e-13, stopping_dmdt= 0.01, max_steps= 5000, save_m_steps= 100,
              save_vtk_steps= 50)
    np.save('m0.npy', sim.spin)




if __name__ == '__main__':
    mesh = Mesh(nx= 200, ny = 50, nz = 1, dx = 2.5, dy = 2.5, dz = 3,
                unit_length= 1e-9)
    relax_system(mesh)