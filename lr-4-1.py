import numpy as np
import random as rn
import matplotlib.pyplot as plt


if __name__ == "__main__":

    ac = 0.001
    n = 1000
    numtoplot = 200

    result = np.zeros(n)
    result[0] = rn.random()

    for r in np.arange(0, 4, ac):
        for i in range(n - 1):
            result[i+1] = r * result[i] * (1 - result[i])
        plt.plot([r] * numtoplot, result[n - numtoplot:], 'b.', markersize=.03)


plt.title('Бифуркационная диаграмма')
plt.xlabel('r')
plt.ylabel('x')
plt.show()
