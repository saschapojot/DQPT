import numpy as np

mu = 0.3
r = 0.6
gam = 1
# starting time is 0
t0 = 0
# time grids, t0,...,tM
M = 150
# timeIndex = range(0, M + 1)
# time difference between adjacent grid
totalTime = 30
dt = totalTime / M
# momentum grids, k0,...,kN
NN = 500
# kIndex = range(0, NN + 1)
k0 = 0
# momentum difference between adjacent grid
dk = 2 * np.pi / NN
# (timeIndex, kIndex)


# indices for V1
timeIndexV1 = range(0, M + 1)
kIndexV1 = range(0, NN)
tkIndexV1 = [(tInd, kInd) for tInd in timeIndexV1 for kInd in kIndexV1]
colNamesV1 = ['k' + str(n) for n in range(0, NN)]
rowNamesV1 = ['t' + str(m) for m in range(0, M + 1)]
# indices for U2
timeIndexU2 = range(0, M)
kIndexU2 = range(0, NN + 1)
tkIndexU2 = [(tInd, kInd) for tInd in timeIndexU2 for kInd in kIndexU2]
# maximum grid index in Simpson integration, start from 0
gridSimpson = 20
# create a data frame, row :t, col:k.
outFigName = 'jump1.png'
