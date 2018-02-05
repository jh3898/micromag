import micro_jinting2.micro.micro_clib
import numpy as np
from micro_jinting2.common.constant import mu_0
from micro_jinting2.micro.energy import Energy
import micro_jinting2.common.helper


class UniaxialAnisotropy(Energy):

    """
        compute the anisotropy field with the energy density E = K[1- (m.u)^2]
    """

    def __init__(self, Ku, axis=(1, 0, 0), name='Anisotropy'):
        self.Ku = Ku
        self.name = name
        self.jac = True
        self.axis = axis

    def setup(self, mesh, spin, Ms):
        super(UniaxialAnisotropy, self).setup(mesh, spin, Ms)

        self._Ku = helper.init_scalar(self.Ku, self.mesh)
        self._axis = helper.init_vector(self.axis, self.mesh, True)

    def compute_field(self, t=0, spin=None):
        if spin is not None:
            m = spin
        else:
            m = self.spin
        micro_clib.compute_anisotropy_micro(m,
                                            self.field,
                                            self.energy,
                                            self.Ms_inv,
                                            self._Ku,
                                            self._axis,
                                            self.nx,
                                            self.ny,
                                            self.nz)
        return self.field
