
import numpy as np
import os
import zipfile

class MinimizerBase():

    def __init__(self, mesh, spin, magnetization, magnetization_inv, field, pins,
                 interactions, name, data_saver):

        self.mesh = mesh
        self.spin = spin

        self._magnetization = magnetization
        self._magnetization_inv = magnetization_inv
        self.field = field

        self._pins = pins
        self.interactions = interactions
        self.name = name
        self.data_saver = data_saver

        ############# variables


    def save_m(self, ZIP = False):

        if not os.path.exists('%s_npys' % self.name):
            os.makedirs('%s_npys' % self.name)
        name = '%s_npys/m_%g.npy' % (self.name, self.step)
        np.save(name, self.spin)
        if ZIP:
            with zipfile.ZipFile('%s_m.zip' % self.name, 'a') as myzip:
                myzip.write(name)
            try:
                os.remove(name)
            except OSError:
                pass



