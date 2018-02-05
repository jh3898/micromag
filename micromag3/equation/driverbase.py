from __future__ import division
from __future__ import print_function

import os
import numpy as np
import zipfile
import micro_jinting2.common.helper as helper
from micro_jinting2.common.integrators import StepIntegrator,ScipyIntegrator

class DriverBase(object):
    """
    common methods for the micromagnetic driver class
    """
    def __init__(self):
        pass

    def initiate_variables(self,n_spins):
        """ variables for micro"""
        self._alpha = np.zeros(n_spins, dtype= np.float)
        self.t= 0
        self.spin_last= np.ones(3*n_spins, dtype= np.float)
        self.dm_dt= np.zeros(3*n_spins, dtype = np.float)
        self.integrator_tolerances_set= False
        self.step =0

    def get_alpha(self):
        """ return the array with spatially dependent gilbert damping
         per mesh/lattice site"""
        return self._alpha

    def set_alpha(self,value):
        """ """
        self._alpha[:]= helper.init_scalar(value,self.mesh)

    alpha= property(get_alpha, set_alpha)

    def set_integrator(self, integrator, use_jac):

        self.integrator= ScipyIntegrator(self.spin, self.step_rhs)

    def stat(self):
        return self.integrator.stat()

    def set_default_options(self):
        pass

    def reset_integrator(self,t=0):
        #self.integrator.reset()
        self.t= t
        self.step =0

    def compute_effective_field(self,t):
        """compute the effective field from the simulation interaction,
        calling Energy class
        """
        self.field[:] =0
        for obj in self.interactions:
            self.field += obj.compute_field(t)

    def compute_effective_field_jac(self,t,spin):
        self.field[:]=0
        for obj in self.interations:
            if obj.jac:
                self.field += obj.compute_field(t,spin)

    def compute_dmdt(self,dt):
        m0= self.spin_last
        m1= self.spin
        dm= (m1-m0).reshape((3,-1))
        max_dm = np.max(np.sqrt(np.sum(dm**2, axis =0)))
        max_dmdt= max_dm /dt
        return max_dmdt

    def run_until(self,t):
        if t <= self.t:
            if t == self.t and self.t ==0.0:
                self.compute_effective_field(t)
                self.data_saver.save()
                return
            else:
                raise ValueError("t must be >= sim.t")
        ode= self.integrator

        self.spin_last[:]= self.spin[:]

        flag = ode.run_until(t)

        if flag <0:
            raise Exception('Run cython run until failed')

        self.spin[:]= ode.y[:]

        self.t= t
        self.step += t

        self.compute_effective_field(t)
        self.data_saver.save()

    def relax(self,dt= 10e-12,stopping_dmdt= 0.01, max_steps =1000,
              save_m_steps= 100, save_vtk_stps=100):
        """evolve until dmdt < stopping_dmdt"""
        while self.step <max_steps:
            _dt = dt
            t= self.t + _dt
            self.run_until(t)

            if (save_m_steps is not None) and (self.step % save_m_steps ==0):
                self.save_m()

            dmdt= self.compute_dmdt(_dt)
            print("# {:< 4} t= {:8.3g} dt ={:.3g} max_dmdt ={:.3g}".format(
                self.step,
                self.t,
                _dt,
                dmdt / self._dmdt_factor))
            if dmdt < stopping_dmdt * self._dmdt_factor:
                break

        if save_m_steps is not None:
                self.save_m()

    #### define save functions:
    def save_m(self,ZIP= False):
        if not os.path.exists('%s_pys' % self.name):
            os.makedirs('%s_pys' % self.name,int(self.step))
        name= '%s_pys/m_%g.py' % (self.name, int(self.step))
        np.save(name,self.spin)















