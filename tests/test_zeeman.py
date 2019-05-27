
from micromag.Hamiltonian.zeeman import Zeeman

def test_zeeman():
    zeeman = Zeeman((1, 0, 0))
    field = zeeman.compute_field()
    assert field[6] == 1.2 * (2 + 0.5)

def test_zeeman_energy():

    zeeman = Zeeman((1, 0, 0))
    field = zeeman.compute_field()

    #exp_energy = 8 * (- mu0 * H * Ms * mesh.dx * mesh.dy * mesh.dz)

    # assert np.abs(zeeman.compute_energy() - exp_energy)) < 1e-10


if __name__ == '__main__':
    test_zeeman_energy()