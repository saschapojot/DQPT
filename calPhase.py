from consts import *
import pandas as pd
import matplotlib.pyplot as plt

cutOff = 0.8


def jumpDesision(x):
    '''
       because the angle change value is set to (-pi,pi), we need to
       transform it to the true value
       :param x: input phase
       :return:
       '''
    x /= (np.pi)
    # the angle decreases, crosses -pi
    if x >= cutOff:
        x -= 1
    # the angle increases, crosses +pi
    if x <= -cutOff:
        x += 1
    return x


V1Tab = pd.read_csv('V1.csv', index_col=0)
U2Tab = pd.read_csv('U2.csv', index_col=0)
SU2Tab = U2Tab.cumsum()
# create table V2
V2Tab = pd.DataFrame(columns=colNamesV1, index=rowNamesV1)
# create tabel beta
betaTab = pd.DataFrame(columns=colNamesV1, index=rowNamesV1)
for n in range(0, NN):
    V2Tab.iloc[0, n] = 0

for m in range(1, M + 1):
    for n in range(0, NN):
        V2Tab.iloc[m, n] = SU2Tab.iloc[m - 1, n + 1] - SU2Tab.iloc[m - 1, n]

for m in range(0, M + 1):
    for n in range(0, NN):
        # print(V1Tab.iloc[m,n])
        betaTab.iloc[m, n] = jumpDesision(
            V1Tab.iloc[m, n] + V2Tab.iloc[m, n])

windingNum = betaTab.sum(axis=1)
Ts = [t0 + dt * tInd for tInd in timeIndexV1]
plt.plot(Ts, windingNum)
plt.xlabel('Time')
plt.ylabel('Winding Number')
plt.xticks(np.arange(Ts[0], Ts[-1], 2))
plt.title('$\mu = $' + str(mu) + ', $r = $' + str(r) + ', $\gamma = $' + str(gam))
plt.savefig(outFigName)
