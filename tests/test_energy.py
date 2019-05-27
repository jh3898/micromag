
import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np


def test_energy(do_plot = False):

    mesh = Mesh(nx = 30)
    print('Mesh done')

    energy = data['E_total']

    for i in range(len(energy) - 1):
        assert energy[i] > energy[i + 1]


if __name__ == '__main__':
    test_energy(do_plot= False)