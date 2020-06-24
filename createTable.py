from consts import *
from expressions import *
import pandas as pd
from multiprocessing import Pool
from datetime import datetime

# this script produces table V1 and U2
# np has been included
# starting time is 0

# create dataframe for V1
'''colNamesV1 = ['k' + str(n) for n in range(0, NN)]
rowNamesV1 = ['t' + str(m) for m in range(0, M + 1)]'''
V1Tab = pd.DataFrame(columns=colNamesV1, index=rowNamesV1)
# create dataframe for U2
colNamesU2 = [str(n) for n in range(0, NN + 1)]
rowNamesU2 = [str(j) for j in range(0, M)]
U2Tab = pd.DataFrame(columns=colNamesU2, index=rowNamesU2)


# calculate entries of V1
def entryV1(tInd, kInd):
    '''
   calculate V1Tab[tInd,kInd]
    :param tInd:
    :param kInd:
    :param simpN:
    :return:
    '''
    tTmp = t0 + tInd * dt
    kTmp = k0 + kInd * dk
    rst = V1(tTmp, kTmp, kTmp + dk)
    return rst


def fillV1(tkInd):
    tInd = tkInd[0]
    kInd = tkInd[1]
    rst = entryV1(tInd, kInd)
    return [tInd, kInd, rst]


beforeCalV1 = datetime.now()

# fill in V1 with multithreading
pool1 = Pool(10)
V1Vals = pool1.map(fillV1, tkIndexV1)
pool1.close()
pool1.join()
for item in V1Vals:
    tInd = item[0]
    kInd = item[1]
    val = item[2]
    V1Tab.iloc[tInd, kInd] = np.real(val)

# write V1 to csv
afterCalV1 = datetime.now()
print('Computation time for V1: ', afterCalV1 - beforeCalV1)
V1Tab.to_csv('V1.csv')


# calculate entries for U2
def entryU2(tInd, kInd, simpN):
    '''
    calculate U2[tInd, kInd]
    :param tInd:
    :param kInd:
    :param simpN:
    :return:
    '''
    tTmp = t0 + tInd * dt
    kTmp = k0 + kInd * dk
    rst = simpsonT2(tTmp, dt, simpN, kTmp)
    return rst


def fillU2(tkInd):
    tInd = tkInd[0]
    kInd = tkInd[1]
    rst = entryU2(tInd, kInd, gridSimpson)
    return [tInd, kInd, rst]


beforeCalU2 = datetime.now()
# fill in U2 with multithreading
pool2 = Pool(10)
U2Vals = pool2.map(fillU2, tkIndexU2)
pool2.close()
pool2.join()
for item in U2Vals:
    tInd = item[0]
    kInd = item[1]
    val = item[2]
    U2Tab.iloc[tInd, kInd] = np.real(val)

# write U2Tab to csv
afterCalU2 = datetime.now()
print('Computation time for U2: ', afterCalU2 - beforeCalU2)
U2Tab.to_csv('U2.csv')
