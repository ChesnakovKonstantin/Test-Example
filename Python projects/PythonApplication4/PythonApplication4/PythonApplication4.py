import matplotlib . pyplot as plt
import numpy as np
from scipy.integrate import odeint
from sympy import Symbol , solve , lambdify , Matrix

def GetPyFunc(y, eq, valK1m, valK3m, valK2, valK3):
    return lambdify(y, eq.subs(k2, valK2).subs(km1, valK1m).subs(km3, valK3m).subs(k3, valK3))

def PlotsFor1Param(k1, x, y, yVal, eq1, eq2, eqDet, eqTr, valK1m, valK3m, valK2 = 2, valK3 = 0.0032):
    func1 = GetPyFunc(y, eq1, valK1m, valK3m, valK2, valK3)
    plt.plot(func1(yVal), yVal, color = 'b', label ="$y_{k_2}$")
    xVal = func1(yVal)
    func2 = GetPyFunc(y, eq2, valK1m, valK3m, valK2, valK3)
    plt.plot(func1(yVal), func2(yVal), color = 'g', label =" $x_{k_2}$")
    print(eqDet.subs(x, eq2))
    funcDet = GetPyFunc(y, eqDet.subs(x, eq2).subs(k1, eq1), valK1m, valK3m, valK2, valK3)
    detA = funcDet(yVal)
    funcTr = GetPyFunc(y, eqTr.subs(x, eq2).subs(k1, eq1), valK1m, valK3m, valK2, valK3)
    trA = funcTr(yVal)
    for i in range(1, len(yVal)):
        if (detA[i-1]*detA[i] <= 0):
            print(i)
            print(func1(i))
            plt.plot(func1(yVal[i]), yVal[i], color ='b', linestyle = ' ', marker ='o')
            plt.plot(func1(yVal[i]), func2(yVal[i]), color ='g', linestyle = ' ', marker ='o')
        if (trA[i-1]*trA[i] <= 0):
            print(i)
            print(func1(i))
            plt.plot(func1(yVal[i]), yVal[i], color ='b', linestyle = ' ', marker ='s')
            plt.plot(func1(yVal[i]), func2(yVal[i]), color ='g', linestyle = ' ', marker ='s')
    plt.xlim([-2, 2])
    plt.ylim([0,1])
    plt.xlabel('$k_2$')
    plt.ylabel('x,y')
    plt.show()

def PlotFor2Param(k1, k2, x, y, yVal, eq1, eq2, eqDet, eqTr, valK1m, valK3m, valK2 = 2, valK3 = 0.0032):
    eqk1_1Det = solve(eqDet.subs(x, eq2), k1)
    eqForK1Det = eq1 - eqk1_1Det[0]
    eqK2Det = solve(eqForK1Det.subs(km1, valK1m).subs(km3, valK3m).subs(k3, valK3), k2)
    funcK2Det = lambdify(y, eqK2Det[0])
    funcK1Det = lambdify(y, eqk1_1Det[0].subs(x, eq2).subs(k2, eqK2Det[0]).subs(km1, valK1m).subs(km3, valK3m).subs(k3, valK3))
    hopf, = plt.plot(funcK1Det(yVal), funcK2Det(yVal), linestyle = '--', color = 'r', label = 'hopf line')

    eqk1_1Tr = solve(eqTr.subs(x, eq2), k1)
    eqForK1Tr = eq1 - eqk1_1Tr[0]
    eqK2Tr = solve(eqForK1Tr.subs(km1, valK1m).subs(km3, valK3m).subs(k3, valK3), k2)
    funcK2Tr = lambdify(y, eqK2Tr[0])
    funcK1Tr = lambdify(y, eqk1_1Tr[0].subs(x, eq2).subs(k2, eqK2Tr[0]).subs(km1, valK1m).subs(km3, valK3m).subs(k3, valK3))
    sadle, = plt.plot(funcK1Tr(yVal), funcK2Tr(yVal), color = 'g', label = 'sadle-nodle line')

    plt.xlim([-1, 2])
    plt.ylim([-1,5])
    plt.legend(handles=[hopf, sadle])
    plt.show()

k1 = Symbol("k1")
k2 = Symbol("k2")
k3 = Symbol("k3")
km1 = Symbol("km1")
km3 = Symbol("km3")
x = Symbol("x")
y = Symbol("y")

#eq1 = k1*(1 - x - y) - km1*x - k3*x*y
#eq2 = k2*(1 - x - y)**2 - km2*y**2 - k3*x*y

eq1 = k1*(1 - x - y) - km1*x - k3*x + km3*y - k2*((1 - x - y)**2)*x
eq2 = k3*x - km3*y

#eq1 = k1*(1 - x - y) - km1*x - k2*((1 - x - y)**2)*x
#eq2 = k3*(1 - x - y)**2 - km3*y*y

res = solve([eq1, eq2], k1, x)
print(res)
k1Y = res[0][0]
xY= res[0][1]
A = Matrix([eq1, eq2])
var_vector = Matrix([x, y])
jacA = A.jacobian(var_vector)
detA = jacA.det()
traceA = jacA.trace()
DA = traceA**2 - 4* detA
yVal = np.linspace(0, 1, 100000)

PlotsFor1Param(k1, x, y, yVal, k1Y, xY, detA, traceA, valK1m = 0.01, valK3m = 0.002)
PlotFor2Param(k1, k2, x, y, yVal, k1Y, xY, detA, traceA, valK1m = 0.01, valK3m = 0.002)
#plt.plot(GetPyFunc(k2Y, 0.04, 0.02)(y), y, color = 'b', label ="$y_{k_2}$")
#plt.plot(funck2(y), funcxy(y), color = 'g', label =" $x_{k_2}$")