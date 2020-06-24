from sympy import *

from consts import *

# this script calculates the expressions
# a, b, t1, t2, z, z1, z2, z3 = symbols('a,b,t1,t2,z,z1,z2,z3', cls=Symbol)
t1, t2 = symbols('t1,t2', cls=Symbol, real=True)
k = symbols('k', cls=Symbol, real=True)

# consts
z1 = mu + r * cos(k)
z2 = 0
z3 = r * sin(k) + I * gam / 2
w1 = mu + r * cos(k)
w2 = 0
w3 = r * sin(k)
z = sqrt(z1 ** 2 + z2 ** 2 + z3 ** 2)
a = I * w2 - w1
b = w3 + sqrt(w1 ** 2 + w2 ** 2 + w3 ** 2)
# x component of vector psi, at time t1 and t2
psiX1 = a * cos(t1 * z) + (-z2 * b - I * z1 * b - I * z3 * a) / z * sin(t1 * z)
psiX2 = a * cos(t2 * z) + (-z2 * b - I * z1 * b - I * z3 * a) / z * sin(t2 * z)
# y component of psi
psiY1 = b * cos(t1 * z) + 1 / z * (I * z3 * b - I * z1 * a + z2 * a) * sin(t1 * z)
psiY2 = b * cos(t2 * z) + 1 / z * (I * z3 * b - I * z1 * a + z2 * a) * sin(t2 * z)
# x component of d psi/dt
dPsiXDt = diff(psiX1, (t1, 1))
dPsiYDt = diff(psiY1, (t1, 1))

# assemble H^{f}
# Pauli matrices
s1 = Matrix([[0, 1], [1, 0]])
s2 = Matrix([[0, -I], [I, 0]])
s3 = Matrix([[1, 0], [0, -1]])
Hf = z1 * s1 + z2 * s2 + z3 * s3
# Term1, -ilog(<t1|t2>/|<t1|t2>|), assemble the numerator and lambdify it
# it is a function of t1, t2 and k
numerator1 = conjugate(psiX1) * psiX2 + conjugate(psiY1) * psiY2
# lambdify of numerator1
Num1Func = lambdify([t1, t2, k], numerator1, 'numpy')


def T1(t1, t2, k):
    numTmp = Num1Func(t1, t2, k)
    return -1j * np.log(numTmp / np.abs(numTmp))


# Term3, assemble the numerator <t2|t2>
numerator3 = conjugate(psiX2) * psiX2 + conjugate(psiY2) * psiY2
# Term3, assemble the denominator<t1|t1>
denominator3 = conjugate(psiX1) * psiX1 + conjugate(psiY1) * psiY1
term3Tmp = -I * log(numerator3 / denominator3) / 2
T3 = lambdify([t1, t2, k], term3Tmp, 'numpy')


# V11 is the first part of V1, V11=T1(0,mt,n1k)-T1(0,mt,nk)
def V11(mt, nk, n1k):
    return T1(0, mt, n1k) - T1(0, mt, nk)


def V12(mt, nk, n1k):
    return T3(0, mt, n1k) - T3(0, mt, nk)


# V1=V11+V12
def V1(mt, nk, n1k):
    return T1(0, mt, n1k) - T1(0, mt, nk) + T3(0, mt, n1k) - T3(0, mt, nk)


# Term2, assemble the numerator
# <t1|Hf|t1>
leftPsiTmp = Matrix(1, 2, [conjugate(psiX1), conjugate(psiY1)])
rightPsiTmp = Matrix(2, 1, [psiX1, psiY1])
numerator2 = (leftPsiTmp * Hf * rightPsiTmp)[0]
# Term2, assemble the denominator<t1|t1>
denominator2 = conjugate(psiX1) * psiX1 + conjugate(psiY1) * psiY1
term2Tmp = numerator2 / denominator2
T2 = lambdify([t1, k], term2Tmp)


def simpsonT2(jt, dt, N, nk):
    '''
    Simpson integration of T2. SimpsonT2 is U2
    :param jt:
    :param dt:
    :param N:
    :param nk:
    :return:
    '''
    t0 = jt
    tn = jt + dt
    tStep = dt / N
    evenStep = range(2, N, 2)
    oddStep = range(1, N, 2)
    evenStepVal = [T2(t0 + tStep * j, nk) for j in evenStep]
    oddStepVal = [T2(t0 + tStep * j, nk) for j in oddStep]
    rst = tStep / 3 * (T2(t0, nk) + 4 * sum(oddStepVal) + 2 * sum(evenStepVal) + T2(tn, nk))
    return rst
