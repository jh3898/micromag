import numpy as np
import  matplotlib.pyplot as plt

def func(x):
    return np.sinc(x) ** 2

xr = np.random.uniform(low = -3.5, high= 3.5, size= 5000)
yr = func(xr)
I = 2 * 3.5 * np.mean(yr)


def mc_int(func, domain, n_samples):
    samples = np.random.uniform(low = domain[0], high = domain[1], size= (n_samples, ))
    volume = abs(domain[1] - domain[0])
    return np.mean(func(samples)) * volume


def mc_solver(func, y0, x, n_samples):

    vals = [y0]

    for lo, hi in zip(x[:-1], x[1:]):
        vals.append(vals[-1] + mc_int(func, (lo, hi), n_samples))
    return np.asarray(vals)

def f(x):
    return x ** 3 + 2 * x ** 2 - 3** x

xs = np.linspace(-2.5, 2.05, 50)
y0 = -11.2
ys = []

for _ in range(5):
    y = mc_solver(func, y0, xs, 5)
    ys.append(y)

y_mean = np.mean(ys, axis= 0)
y_std = np.std(ys, axis = 0)

fig = plt.figure(figsize = (12, 8))

plt.plot(xs, y_mean, linewidth = 3, label = 'MC')
plt.fill_between(xs, y_mean - 3 * y_std, y_mean + 3 * y_std)
plt.show()
