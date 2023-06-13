import numpy as np
import matplotlib.pyplot as plt

def filler(yFix, y1):
    for i in range(y1[0].size):
        yFix = np.append(yFix, None)
        yFix[-1] = (y1[0][i], y1[1][i])
    return yFix

def error(yFix, y, yError):
    for i in range(len(yFix)):
        yError = np.append(yError, None)
        yError[-1] = ()

def result (y, yFix, error, errorE):
    print(len(yFix))
    for i in range(len(yFix)):
        error = np.append(error, np.sqrt((y[i][0] - yFix[i][0])*(y[i][0] - yFix[i][0]) + (y[i][1] - yFix[i][1])*(y[i][1] - yFix[i][1])))
        errorE = np.append(errorE, np.sqrt(
            (yE[i][0] - yFix[i][0]) * (yE[i][0] - yFix[i][0]) + (yE[i][1] - yFix[i][1]) * (yE[i][1] - yFix[i][1])))
    return error, errorE

def function(t, y): #функция
    return (y[1], -y[0])

def anal_function(t, A, fi0):
    t = t.astype(float)
    return (A * np.cos((t + fi0*np.ones_like(t))*1.0), -A * np.sin((t + fi0*np.ones_like(t))*1.0))

def anal_function2(t, A, fi0):
    return (A * np.cos(t + fi0), -A * np.sin(t + fi0))

def euler(h, t, y): #метод эйлера
    while t[-1] < 12:
        y = np.append(y, None)
        y[-1] = (tuple(map(sum, zip(y[-2], tuple(h*x for x in function(t[-1], y[-2]))))))
        t = np.append(t, None)
        t[-1] = (t[-2] + h)
    return t, y
def rk4(h, t, y):

    while t[-1] < 12:
        k1 = function(t[-1], y[-1])
        k2 = function(t[-1] + h/2, tuple(map(sum, zip(y[-1], tuple(h/2*x for x in k1)))))
        k3 = function(t[-1] + h/2, tuple(map(sum, zip(y[-1], tuple(h/2*x for x in k2)))))
        k4 = function(t[-1] + h, tuple(map(sum, zip(y[-1], tuple(h*x for x in k3)))))
        y = np.append(y, None)
        y[-1] = (tuple(map(sum, zip(y[-2], tuple(h/6*x for x in k1), tuple(2*h/6*x for x in k2), tuple(2*h/6*x for x in k3), tuple(h/6*x for x in k4)))))
        t = np.append(t, None)
        t[-1] = (t[-2] + h)
    return t, y

if __name__ == "__main__":

    A = 1.0
    fi0 = 0.0

    t0 = 0.0
    y0 = anal_function2(t0, A, fi0)

    h = 0.1
    h2 = 0.01
    h3 = 0.001

    y = np.array([])
    t = np.array([])

    y2 = np.array([])
    t2 = np.array([])

    y3 = np.array([])
    t3 = np.array([])

    yE = np.array([])
    tE = np.array([])

    y = np.append(y, None)
    y[-1] = y0
    t = np.append(t, None)
    t[-1] = t0

    y2 = np.append(y2, None)
    y2[-1] = y0
    t2 = np.append(t2, None)
    t2[-1] = t0

    y3 = np.append(y3, None)
    y3[-1] = y0
    t3 = np.append(t3, None)
    t3[-1] = t0

    yE = np.append(yE, None)
    yE[-1] = y0
    tE = np.append(tE, None)
    tE[-1] = t0

    t, y = rk4(h, t, y)
    t2, y2 = rk4(h2, t2, y2)
    t3, y3 = rk4(h3, t3, y3)

    tE, yE = euler(h3, tE, yE)

    y1 = anal_function(t, A, fi0)
    y12 = anal_function(t2, A, fi0)
    y13 = anal_function(t3, A, fi0)

    yFix = np.array([])
    yFix2 = np.array([])
    yFix3 = np.array([])

    yFix = filler(yFix, y1)
    yFix2 = filler(yFix2, y12)
    yFix3 = filler(yFix3, y13)

    yError = np.array([])
    yError2 = np.array([])
    yError3 = np.array([])

    error = np.array([])
    error2 = np.array([])
    error3 = np.array([])

    errorE = np.array([])
    errorE2 = np.array([])
    errorE3 = np.array([])

    error, errorE = result(y, yFix, error, errorE)
    error2, errorE2 = result(y2, yFix2, error2, errorE2)
    error3, errorE3 = result(y3, yFix3, error3, errorE3)

    #for i in range(len(y1)):
        #error[i] = np.sqrt((y[0][i] - yFix[0][i])*(y[0][i] - yFix[0][i]) + (y[1][i] - yFix[1][i])*(y[1][i] - yFix[1][i]))

    plt.plot(t, error, label="h = 0.1")
    plt.plot(t2, error2, label="h = 0.01")
    plt.plot(t3, error3, label="h = 0.001")

    plt.semilogy()


    #plt.plot(tE, errorE, label="ME")
    plt.legend(loc="upper right")

    plt.title('Ошибка численного решения при различых шагах интегрирования')
    plt.xlabel('t')
    plt.ylabel(r'$\delta$x(t)')
    plt.show()
