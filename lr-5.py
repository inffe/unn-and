import numpy as np
import random as rn
import matplotlib.pyplot as plt

def ref1(protein_count, c, maxT, h, x, t):
    cur_t = 0.0
    while (cur_t < maxT):
        x.append(protein_count * np.exp(-c * cur_t))
        cur_t += h
        t.append(cur_t)

def model1(protein_count, c, maxT, t, x):
    cur_t = 0.0
    cur_pc = protein_count

    while(cur_t < maxT):
        if (cur_pc == 0):
            break
        a = c*cur_pc
        rand_value = rn.random()
        tau = (1.0 / a) * np.log(1.0/rand_value)

        cur_pc -= 1
        x.append(cur_pc)

        cur_t += tau
        t.append(cur_t)

def functionOne(beta, gamma, m_prev, x_prev):
    return beta * m_prev - gamma * x_prev

def functionTwo(K, n, m_prev, x_prev):
    return K / (1.0 + pow(x_prev, n)) - m_prev

def ref2(K, beta, gamma, n, maxT, h, x, m, t):
    cur_t = 0
    while (cur_t < maxT):
        cur_t += h

        x_k1 = functionOne(beta, gamma, m[-1], x[-1])
        m_k1 = functionTwo(K, n, m[-1], x[-1])

        x_k2 = functionOne(beta, gamma, m[-1] + h / 2 * m_k1, x[-1] + h / 2 * x_k1)
        m_k2 = functionTwo(K, n, m[-1] + h / 2 * m_k1, x[-1] + h / 2 * x_k1)

        x_k3 = functionOne(beta, gamma, m[-1] + h / 2 * m_k2, x[-1] + h / 2 * x_k2)
        m_k3 = functionTwo(K, n, m[-1] + h / 2 * m_k2, x[-1] + h / 2 * x_k2)

        x_k4 = functionOne(beta, gamma, m[-1] + h / 2 * m_k3, x[-1] + h / 2 * x_k3)
        m_k4 = functionTwo(K, n, m[-1] + h / 2 * m_k3, x[-1] + h / 2 * x_k3)

        x_value = x[-1] + (x_k1 + 2*x_k2 + 2*x_k3 + x_k4) * h/6
        m_value = m[-1] + (m_k1 + 2*m_k2 + 2*m_k3 + m_k4) * h/6

        x.append(x_value)
        m.append(m_value)
        t.append(cur_t)

def model2(protein_count, mrnk_count, K, beta, gamma, n, maxT, t, x, m):
    cur_t = 0
    t.append(cur_t)

    cpt = protein_count
    x.append(cpt)

    cmk = mrnk_count
    m.append(cmk
             )
    while(cur_t < maxT):
        a1 = K / (1 + pow(x[-1], n))
        a2 = m[-1]
        a3 = beta * m[-1]
        a4 = gamma * x[-1]
        A = a1 + a2 + a3 + a4

        p1 = a1/A
        p2 = a2/A
        p3 = a3/A
        p4 = a4/A

        tau = 1 / A * (np.log(1 / rn.random()))

        if (rn.random() < p1):
            cmk += 1
        elif (rn.random() < p1 + p2):
            cmk -= 1
        elif (rn.random() < p1 + p2 + p3):
            cpt += 1
        elif (rn.random() < p1 + p2 + p3 + p4):
            cpt -= 1

        t.append(cur_t)
        x.append(cpt)
        m.append(cmk)
        cur_t += tau


if __name__ == "__main__":

    protein_count = 25
    c = 3.0
    maxT = 5
    h = 0.1

    mrnk_count = 40;
    K = 250
    beta = 1.0
    gamma = 2
    n = 6

    xOne = []

    xTwo = []
    xTwo.append(protein_count)
    tOne = []
    tOne.append(0)

    tTwo = []
    tTwo.append(0)

    xThree = []
    xThree.append(protein_count)
    mThree = []
    mThree.append(mrnk_count)
    tThree = []
    tThree.append(0)

    xFour = []
    mFour = []
    tFour = []

ref1(protein_count, c, maxT, h, xOne, tOne)
xOne.append(0)
model1(protein_count, c, maxT, tTwo, xTwo)
tTwo.append(5)
xTwo.append(0)

ref2(K, beta, gamma, n, maxT, h, xThree, mThree, tThree)
model2(protein_count, mrnk_count, K, beta, gamma, n, maxT, tFour, xFour, mFour)

#plt.plot(tOne, xOne, label="ref")
#plt.plot(tTwo, xTwo, label="model")

plt.plot(tFour, mFour, label="model")
plt.plot(tThree, mThree, label="ref")

plt.legend(loc="upper right")
plt.title('Сравнение двух моделей авторепрессора для молекул мРНК')
plt.xlabel('t')
plt.ylabel('x')
plt.show()

