from __future__ import division
from __future__ import print_function

class LLG(object):
    """
    solve Landau Liftshitz Gilbert equation

    dm - gamma
    ---- = --------  (m X H_eff  + a * m X ( m X H_eff) )
    dt
    2
    (1 + a)
    """
    def __init__(self,mesh, spin, Ms, field, pins, interations, name, data_saver, integrator = '',
                 use_jac= False):


    def LLG_rhs(self, t, y, ydot):
        self.t =t
        self.compute_effective_field(t)
        compute_llg_rhs(ydot, self.spin, self.field, self._alpha, self._pins, self.gamma,
                        self.mesh.n, self.do_precession, self.default_c)
        return 0

    def compute_llg_rhs(self, dmdt, m, h, alpha, pins, gamma, n, do_precession, default_c):
        n= self.mesh.n

        H_eff = self.field
        m = self.spin
        #mx= m[0]


        for ia in range(0,n):
            i = 3 * ia
            j = 3 * ia + 1
            k = 3 * ia + 2
            coff = -self.gamma / (1 + (self._alpha[i]) ** 2)
            mx = m[i]
            my = m[j]
            mz = m[k]
            H_effx = H_eff[i]
            H_effy = H_eff[j]
            H_effz = H_eff[k]
            ### dot m * H_eff
            mH= mx * H_effx + my * H_effy + mz * H_effz
            mm= mx * mx + my * my + mz * mz

            #### cross product
            mXH = np.cross((mx, my, mz), (H_effx, H_effy, H_effz))

            dmdt[i] = coff * (mXH[0] + self._alpha[ia] * (mH * mx - mm * H_effx))
            dmdt[j] = coff * (mXH[1] + self._alpha[ia] * (mH * my - mm * H_effy))
            dmdt[k] = coff * (mXH[2] + self._alpha[ia] * (mH * mz - mm * H_effz))

            if default_c < 0:
                c= 6 * sqrt(dmdt[i]**2 + dmdt[j]**2 + dmdt[k]**2 )
            else:
                c= default_c

            dmdt[i] += c * (1 - mm) * mx
            dmdt[j] += c * (1 - mm) * my
            dmdt[k] += c * (1 - mm) * mz








