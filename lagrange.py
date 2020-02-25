import numpy as np
import matplotlib.pyplot as plt

def basis(i,m,xi,x):
    if m==len(xi):
        return 1
    elif i==m:
        return 1*basis(i,m+1,xi,x)
    else:
        return (x-xi[m])/(xi[i]-xi[m])*basis(i,m+1,xi,x)

def dbasis(i,m,j,xi,x,value):
    if j==len(xi):
        return 0
    elif j==i:
        return dbasis(i,len(xi),j+1,xi,x,value)
    elif m==-1:
        return 1
    else:
        if m==i or m==j:
            return dbasis(i,m-1,j,xi,x,value)
        elif m==len(xi):
            value = 1/(xi[i]-xi[j])*dbasis(i,m-1,j,xi,x,value)
            value += dbasis(i,len(xi),j+1,xi,x,value)
            return value
        else:
            return (x-xi[m])/(xi[i]-xi[m])*dbasis(i,m-1,j,xi,x,value)

def lagrange(xi,x):
    N = len(xi)
    phi = []
    dphi = []
    for i in range(N):
        phi.append(basis(i,0,xi,x))
        dphi.append(dbasis(i,len(xi),0,xi,x,0))
    return phi, dphi


GLL = 10
N = 8
xi,w  = np.polynomial.legendre.leggauss(GLL)
x = np.linspace(-1, 1, N)
f = np.sin(xi**2)
der = np.cos(x**2)*2*x
phi, dphi = lagrange(xi,x)
dphi = np.array(dphi)
der_ap = f@dphi

plt.figure()
for i in range(GLL):
    plt.plot(x,phi[i], label=r'$\phi_{%i}$' %(i+1))
plt.legend()
plt.show()

plt.figure()
plt.plot(x,der,label='Analytical Derivative')
plt.plot(x,der_ap,'--',label='Interpolation')
plt.legend()
plt.show()
