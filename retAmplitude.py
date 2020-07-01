from expressions import *
import matplotlib.pyplot as plt


def retInt(t, k):
    return np.log(np.abs(G(t, k)) ** 2)


# Simpson's rule maximum grid index
GN = 200


def simpsonGInt(t):
    step = np.pi * 2 / GN
    oddStepInd = range(1, GN, 2)
    evenStepInd = range(2, GN, 2)
    k00 = 0
    oddStep = [retInt(t, k00 + j * step) for j in oddStepInd]
    evenStep = [retInt(t, k00 + j * step) for j in evenStepInd]

    rst = step / 3 * (retInt(t, k00) + 4 * sum(oddStep) + 2 * sum(evenStep) + retInt(t, k00 + step * GN))
    return rst


def gMinus(t):
    return -1 / (2 * np.pi) * simpsonGInt(t)


allTime = [t0 + dt * j for j in timeIndexV1]
gm = [gMinus(tj) for tj in allTime]
plt.figure()
plt.plot(allTime, gm)
plt.show()
