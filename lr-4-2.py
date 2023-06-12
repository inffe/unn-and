import numpy as np
import random as rn
import matplotlib.pyplot as plt

def function(x, r):
    return r*x*(1-x)

def lyapunov(x0, r, n):
    x = x0
    sum = 0
    for i in range(n):
        zxc = r - 2*r*x
        x = function(x, r)
        sum += np.log(abs(zxc))
    return sum / n

if __name__ == "__main__":

    x = 0.1
    n = 1000

    result = []
    rarray = []

    for r in np.arange(0, 4, 0.0001):
        result.append(lyapunov(x, r, n))
        rarray.append(r)

plt.plot(rarray, result)
plt.grid()
plt.title('Ляпуновский показатель')
plt.xlabel('r')
plt.ylabel(r'$\lambda$')
plt.show()

