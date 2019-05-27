import numpy as np

from minimizerbase import MinimizerBase

class GradientDescentMinimizer(MinimizerBase):
    """
    minimization: evolution of the system
    h_i = (m_i + H) / ||(m_i + H)||
    m_i + 1 = m_i + alpha * (h_i - m_i)
    """


    def __init__(self, mesh, spin, magnetization, magnetization_inv, field, pins,
                 interactions, name, data_saver, use_jac = False, integrator = None):

        super().__init__(mesh, spin, magnetization, magnetization_inv, field, pins,
                 interactions, name, data_saver)



        self.t = 1e-4
        self._alpha = 0.1
    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, value):
        self._alpha = value
        self._alpha_field = value * np.ones_like(self.spin)


    def run_step(self):

        self.spin_last[:] = self.spin[:]
        self.update_effective_field()

        self._new_spin[self._material] = (self.spin + self.field)[self._material]
        self._new_spin[self._material] = (self.spin_last + self._alpha_field * (self._new_spin - self.spin_last))[self._material]
        self.spin[:] = self._new_spin[:]


    def minimize(self, stopping_dm = 1e-2, max_steps = 2000, save_data_steps = 10,  save_m_steps = None):

        self.step = 0

        while self.step < max_steps:
            self.run_step()

            max_dm = (self.spin - self.spin_last).reshape(-1, 3) ** 2
            max_dm = np.max(np.sqrt(np.sum(max_dm, axis = 1)))

            if max_dm < stopping_dm and self.step >0:
                print('Finished at: max_tau = {:<8.3g} max_dm = {:<10.3g} counter = {}'.format(
                    np.max(np.abs(self.tau)), max_dm, self.step) )

                self.compute_effective_field()
                break

            if self.step % save_data_steps == 0:
                self.compute_effective_field()

            if save_m_steps is not None and self.step % save_data_steps == 0:
                self.save_m()

            self.step += 1



if __name__ == '__main__':
    args = [i for i in range(9)]
    solution = GradientDescentMinimizer(*args)
    res = solution.alpha
    print(res)
    solution.alpha = 0.3
    print(solution.alpha, solution._alpha_field)