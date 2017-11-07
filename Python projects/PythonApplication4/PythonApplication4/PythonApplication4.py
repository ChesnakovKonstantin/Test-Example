import matplotlib . pyplot as plt
import numpy as np
from scipy.integrate import odeint
from sympy import Symbol , solve , lambdify , Matrix

def NonNeg(mylist):
    return [elem for elem in mylist if elem > 0]

def vectorfield(w, t, p):
    x, y = w
    k1, k2, k3, km1, km3 = p

    # Create f = (x1',y1',x2',y2'):
    #f = [k1*(1 - x - y) - km1*x - k3*x + km3*y - k2*((1 - x - y)**2)*x,
    #     k3*x - km3*y]
    f = [k1*(1 - x - 2*y) - km1*x - k3*x*(1 - x - 2*y) + km3*y - k2*((1 - x - 2*y)**2)*x,
        k3*x*(1 - x - 2*y) - km3*y]

    return f

def GetPyFunc(y, eq, valK1m, valK3m, valK2, valK3):
    return lambdify(y, eq.subs(k2, valK2).subs(km1, valK1m).subs(km3, valK3m).subs(k3, valK3))

def PlotsFor1Param(k1, x, y, yVal, eq1, eq2, eqDet, eqTr, valK1m, valK3m, valK2 = 0.95, valK3 = 0.0032):
    func1 = GetPyFunc(y, eq1, valK1m, valK3m, valK2, valK3)
    yGraph, = plt.plot(func1(yVal), yVal, color = 'b', label ="y")
    xVal = func1(yVal)
    func2 = GetPyFunc(y, eq2, valK1m, valK3m, valK2, valK3)
    xGraph, = plt.plot(func1(yVal), func2(yVal), color = 'g', label ="x")
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
    plt.xlim([0, 2])
    plt.ylim([0,1])
    plt.xlabel('$k_1$')
    plt.ylabel('x,y')
    plt.legend(handles=[xGraph, yGraph])
    plt.show()



def PlotFor2Param(k1, k2, x, y, yVal, eq1, eq2, eqDet, eqTr, valK1m, valK3m, valK2 = 0.95, valK3 = 0.0032):
    eqk1_1Det = solve(eqDet.subs(x, eq2), k1)
    eqForK1Det = eq1 - eqk1_1Det[0]
    eqK2Det = solve(eqForK1Det.subs(km1, valK1m).subs(km3, valK3m).subs(k3, valK3), k2)
    funcK2Det = lambdify(y, eqK2Det[0])
    funcK1Det = lambdify(y, eqk1_1Det[0].subs(x, eq2).subs(k2, eqK2Det[0]).subs(km1, valK1m).subs(km3, valK3m).subs(k3, valK3))
    k1DetAll = funcK1Det(yVal)
    k2DetAll = funcK2Det(yVal)
    k1DetPos = [k1DetAll[i] for i in range(len(k1DetAll)) if k1DetAll[i] > 0 and k2DetAll[i] > 0]
    k2DetPos = [k2DetAll[i] for i in range(len(k2DetAll)) if k1DetAll[i] > 0 and k2DetAll[i] > 0]
    hopf, = plt.plot(k1DetPos, k2DetPos, linestyle = '--', color = 'r', label = 'hopf line')

    eqk1_1Tr = solve(eqTr.subs(x, eq2), k1)
    eqForK1Tr = eq1 - eqk1_1Tr[0]
    eqK2Tr = solve(eqForK1Tr.subs(km1, valK1m).subs(km3, valK3m).subs(k3, valK3), k2)
    funcK2Tr = lambdify(y, eqK2Tr[0])
    funcK1Tr = lambdify(y, eq1.subs(x, eq2).subs(k2, eqK2Tr[0]).subs(km1, valK1m).subs(km3, valK3m).subs(k3, valK3))
    k1AllTr = funcK1Tr(yVal)
    k2AllTr = funcK2Tr(yVal)
    k1PosTr = [k1AllTr[i] for i in range(len(k1AllTr)) if k1AllTr[i] > 0 and k2AllTr[i] > 0]
    print(len(k1PosTr))
    k2PosTr = [k2AllTr[i] for i in range(len(k2AllTr)) if k1AllTr[i] > 0 and k2AllTr[i] > 0]
    print(len(k2PosTr))
    sadle, = plt.plot(k1PosTr, k2PosTr, color = 'g', label = 'sadle-nodle line')

    plt.xlim([0, 2])
    plt.ylim([0,5])
    plt.legend(handles=[hopf, sadle])
    plt.xlabel('k1')
    plt.ylabel('k2')
    plt.show()

def PhasePortrait(p):
    Y, X = np.mgrid[0.1:0.3:5000j, 0:0.72:5000j]
    K1, K2, K3, KM1, KM3 = p
    U, V = [K1*(1 - X - 2*Y) - KM1*X - K3*X*(1 - X - 2*Y) + KM3*Y - K2*((1 - X - 2*Y)**2)*X,
            K3*X*(1 - X - 2*Y) - KM3*Y]
    #U, V = [K1*(1 - X - Y) - KM1*X - K3*X + KM3*Y - K2*((1 - X - Y)**2)*X,
    #        K3*X - KM3*Y]
    plt.streamplot(X, Y, U, V, density = [2, 2])
    plt.xlim([0, 0.72])
    plt.ylim([0.1, 0.3])
    plt.xlabel('x')
    plt.ylabel('y')
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

#eq2 = k3*x - km3*y
#eq1 = k1*(1 - x - y) - km1*x - k3*x + km3*y - k2*((1 - x - y)**2)*x

#eq2 = k3*y - km3*x
#eq1 = k1*(1 - x - y) - km1*y - k3*y + km3*x - k2*((1 - x - y)**2)*y

#y and x and changed

eq1 = k1*(1 - y - 2*x) - km1*y - k3*y*(1 - y - 2*x) + km3*x - k2*((1 - y - 2*x)**2)*y
eq2 = k3*y*(1 - y - 2*x) - km3*x

res = solve([eq1, eq2], k1, x)
print(res)
k1Y = res[0][0]
xY= res[0][1]
A = Matrix([eq1, eq2])
var_vector = Matrix([y, x])
jacA = A.jacobian(var_vector)
detA = jacA.det()
traceA = jacA.trace()
DA = traceA**2 - 4* detA
yVal = np.linspace(0, 1, 100000)

PlotsFor1Param(k1, x, y, yVal, k1Y, xY, detA, traceA, valK1m = 0.01, valK3m = 0.002)
PlotFor2Param(k1, k2, x, y, yVal, k1Y, xY, detA, traceA, valK1m = 0.01, valK3m = 0.002)
p = [0.12, 0.98, 0.0032, 0.01, 0.002]
w0 = [0.25, 0.25]
t = [float(i) for i in range(10000)]
wsol = odeint(vectorfield, w0, t, args=(p,))
print(wsol[:, 0])
x_t, = plt.plot(t, wsol[:, 0], color = 'xkcd:sky blue', label = 'X(t)')
y_t, = plt.plot(t, wsol[:, 1], color = 'xkcd:peach', label = 'Y(t)')
plt.legend(handles=[x_t, y_t])
plt.show()
PhasePortrait(p)
#Y, X = np.mgrid[0.1:0.3:5000j, 0:0.72:5000j]
#K1, K2, K3, KM1, KM3 = p
#U, V = [K1*(1 - X - 2*Y) - KM1*X - K3*X*(1 - X - 2*Y) + KM3*Y - K2*((1 - X - 2*Y)**2)*X,
#        K3*X*(1 - X - 2*Y) - KM3*Y]
#plt.streamplot(X, Y, U, V, density = [2, 2])
#plt.xlim([0, 0.72])
#plt.ylim([0.1, 0.3])
#plt.show()