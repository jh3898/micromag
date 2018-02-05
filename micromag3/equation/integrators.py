from scipy.integrate import ode
import warnings

epsilon = 1e-16

class BaseIntegrator(object):
    def __init__(self,spins,rhs_fun):
        self.y = spins
        self.t = 0
        self.rhs = rhs_fun
        self.rhs_evals_nb = 0

    def reset(self, spin, t):
        self.y= spin
        self.t= t

    def run_until(self,t):
        pass

    def set_initial_value(self, spins, t):
        pass

    def rhs_evals(self):
        return self.rhs_evals_nb

class StepIntegrator(BaseIntegrator):
    def __init__(self, spins, rhs_fun, step= "euler", stepsize= 1e-15):
        super(StepIntegrator,self).__init__(spins, rhs_fun)

    def run_until(self,t):
        while abs(self.t -t )> epsilon:
            self.t, self.y, evals= self._step(self.t, self.y, self.stepsize, self.rhs)
            self.rhs_evals_nb += evals
            if self.t >t:
                break
        return 0

    def set_options(self, rtol= 1e-8, atol= 1e-8, integrator= "dopri5"):
        self.rtol = rtol
        self.atol = atol
        self.integrator = integrator

    def set_step(self, step):
        step_choices= {'euler': euler_step, 'rk4':runge_kutta_step}
        if step not in step_choices:
            raise NotImplemented("step must be euler or rk4")
        self._step = step_choices[step]


class ScipyIntegrator(BaseIntegrator):
    def __init__(self, spins, rhs_fun):
        super(ScipyIntegrator,self).__init__(spins,rhs_fun)
        self.integrator_created = False
        self.internal_timesteps= [0]

        def rhs_wrap(y,t):
            self.rhs_evals_nb += 1
            return rhs_fun(y,t)
        self.rhs = rhs_wrap  ### overwriting rhs to count evals

    def solout(self,t,y):
        self.internal_timesteps.append(t)
        return 0

    def set_options(self, rtol =1e-8,atol =1e-8, integrator= "dopri5"):
        self.rtol= rtol
        self.atol= atol
        self.integrator = integrator

    def _create_integrator(self):
        self.ode = ode(self.rhs).set_integrator(self.integrator,rtol= self.rtol,
                                                atol= self.atol)
        self.ode.set_solout(self.solout)
        self.ode.set_initial_value(self.y, self.t)
        self.integrator_created= True

    def run_until(self,t):
        if not self.integrator_created:
            self._create_integrator()

        r = self.ode.integrate(t)
        if not self.ode.successful():
            raise RuntimeError("Integration with ode not sucessful")
        self.y[:] = r
        self.t = t
        return 0

def runge_kutta_step(t,y,h,f):
    k1= f(t,y)
    k2= f(t, h/2.0, y+h *k1/2.0)
    k3= f(t, h/2.0 ,y+h *k2/2.0)
    k4= f(t, h/2.0, y+h *k3)

    tp = t + h
    yp= y + h *(k1 + 2*k2 + 2*k3 + k4)/ 6.0
    evals= 4
    return tp, yp, evals













