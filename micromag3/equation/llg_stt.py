from __future__ import division
from __future__ import print_function
import numpy as np
import common.constant as const
class LLG_STT(object):
    """
    solve Landau Liftshitz Gilbert equation

    dm - gamma
    ---- = --------  (m X H_eff  + a * m X ( m X H_eff) )
    dt
    2
    (1 + a)


    Tau = (js * del ) m
    """
    def __init__(self,mesh, spin, Ms, field, pins, interations, name, data_saver, integrator = '',
                 use_jac= False):


    def LLG_stt_rhs(self, t, y, ydot):
        self.t =t
        self.compute_effective_field(t)
        compute_stt_field(self.spin,
                          self.field_stt,
                          self._jx * self.update_j_fun(t),
                          self._jy * self.update_j_fun(t),
                          self._jz * self.update_j_fun(t),
                          self.mesh.dx * self.mesh.unit_length,
                          self.mesh.dy * self.mesh.unit_length,
                          self.mesh.dz * self.mesh.unit_length,
                          self.mesh.neighbours,
                          self.n)

        compute_llg_rhs(ydot, self.spin, self.field, self._alpha, self._pins, self.gamma,
                        self.mesh.n, self.do_precession, self.default_c)
        return 0

    def compute_stt_field(self.spin, self.field_stt,
                                   self._jx * self.update_j_fun(t),
                                   self._jy * self.update_j_fun(t),
                                   self._jz * self.update_j_fun(t),
                                   self.mesh.dx * self.mesh.unit_length,
                                   self.mesh.dy * self.mesh.unit_length,
                                   self.mesh.dz * self.mesh.unit_length,
                                   self.mesh.neighbours,
                                   self.n):

        field = np.zeros(3 *n, 1)
        for i in range(0,n):
            nn_i = 6 * i


            """in the x direction, nearest neighbors (f(x+1) -f(x-1))/2 if both
            neighbor exist, if left one exists, use f(x) -f(x-1)
            """
            if ngbs[nn_i] is not -1 and ngbs[nn_i +1 ] is not -1:
                factor_x = 2
                nn_x1 = ngbs[nn_i]
                nn_x2 = ngbs[nn_i + 1]

            elif ngbs[nn_i] == -1 and ngbs[nn_i +1 ] is not -1:
                factor_x = 1
                nn_x1 = i
                nn_x2 = ngbs[nn_i + 1]

            elif ngbs[nn_i] is not -1 and ngbs[nn_i +1 ] == -1:
                factor_x = 1
                nn_x1 = ngbs[nn_i]
                nn_x2 = i

            else:
                factor_x = 0

            if factor_x:
                field[3 * i ] += jx[i] * (spin[3* nn_x2] - spin[3* nn_x1] )/ (factor_x *dx)
                field[3 * i + 1 ] += jx[i] * (spin[3 * nn_x2 +1] - spin[3 * nn_x1 +1]) / (factor_x * dx)
                field[3 * i + 2] += jx[i] * (spin[3 * nn_x2 + 2] - spin[3 * nn_x1 + 1]) / (factor_x * dx)
            ###### in the y direction
            if ngbs[nn_i+2] is not -1 and ngbs[nn_i +3 ] is not -1:
                factor_y = 2
                nn_y1 = ngbs[nn_i+2]
                nn_y2 = ngbs[nn_i + 3]

            elif ngbs[nn_i+2] == -1 and ngbs[nn_i +3 ] is not -1:
                factor_y = 1
                nn_y1 = i
                nn_y2 = ngbs[nn_i + 3]

            elif ngbs[nn_i+2] is not -1 and ngbs[nn_i +3 ] == -1:
                factor_y = 1
                nn_y1 = ngbs[nn_i+2]
                nn_y2 = i

            else:
                factor_y = 0

            if factor_y:
                field[3 * i ] += jy[i] * (spin[3* nn_x2] - spin[3* nn_x1] )/ (factor_x *dx)
                field[3 * i + 1 ] += jy[i] * (spin[3 * nn_x2 +1] - spin[3 * nn_x1 +1]) / (factor_y * dy)
                field[3 * i + 2] += jz[i] * (spin[3 * nn_x2 + 2] - spin[3 * nn_x1 + 2]) / (factor_y * dy)
            ###### in the z direction

            if ngbs[nn_i+4] is not -1 and ngbs[nn_i +5 ] is not -1:
                factor_z = 2
                nn_z1 = ngbs[nn_i+4]
                nn_z2 = ngbs[nn_i + 5]

            elif ngbs[nn_i+4] == -1 and ngbs[nn_i +5 ] is not -1:
                factor_z = 1
                nn_z1 =i
                nn_z2 = ngbs[nn_i + 5]

            elif ngbs[nn_i+4] is not -1 and ngbs[nn_i +5 ] == -1:
                factor_z = 1
                nn_z1 = ngbs[nn_i+4]
                nn_z2 = i

            else:
                factor_z = 0

            if factor_z:
                field[3 * i ] += jz[i] * (spin[3* nn_z2] - spin[3* nn_z1] )/ (factor_z *dz)
                field[3 * i + 1 ] += jz[i] * (spin[3 * nn_z2 +1] - spin[3 * nn_z1 +1]) / (factor_z * dz)
                field[3 * i + 2] += jz[i] * (spin[3 * nn_z2 + 2] - spin[3 * nn_z1 + 2]) / (factor_z * dz)
            ###### in the z direction


    def compute_llg_stt_rhs(self, dmdt, m, h, alpha, pins, gamma, n, do_precession, default_c):
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


            ###### adding current to the standard LLG equation above
            coff_stt = u0 / (1 + (self._alpha[ia]) ** 2)
            mht = mx * h_stt[i] + my * h_stt[j] + mz * h_stt[k]
            Hpx = mm * h_stt[i] - mht * mx
            Hpy = mm * h_stt[j] - mht * my
            Hpz = mm * h_stt[k] - mht * mz

            mtH = np.cross((mx, my, mz), (Hpx, Hpy, Hpz))
            dmdt[i] += coff_stt * (1 + self._alpha[ia]*beta)* Hpx - (beta- alpha[ia]* mtH[0])
            dmdt[j] += coff_stt * (1 + self._alpha[ia]*beta)* Hpx - (beta- alpha[ia]* mtH[1])
            dmdt[k] += coff_stt * (1 + self._alpha[ia]*beta)* Hpx - (beta- alpha[ia]* mtH[2])

                c= 6 * sqrt(dmdt[i]**2 + dmdt[j]**2 + dmdt[k]**2 )

            dmdt[i] += c * (1 - mm) * mx
            dmdt[j] += c * (1 - mm) * my
            dmdt[k] += c * (1 - mm) * mz








