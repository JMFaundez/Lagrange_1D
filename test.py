import numpy as np
import matplotlib.pyplot as plt
from lagrange import lagrange1d, lagrange2d

def test_1d():
    GL = 8
    #xi,w  = np.polynomial.legendre.leggauss(GL)
    xi = np.linspace(20,70,GL)

    fun = lambda x: 8*x**5 - np.sqrt(x) + 3*x**2 -10
    dfun = lambda x: 40*x**4 - 0.5*1/np.sqrt(x) + 6*x

    N = 100
    x = np.linspace(20, 70, N)

    phi, dphi = lagrange1d(xi,x)
    dphi = np.array(dphi)

    fi = fun(xi)

    der = dfun(x)
    f_an = fun(x)

    der_ap = fi@dphi
    f_ap = fi@phi
    print(len(phi))
    plt.figure(figsize=(15,8))
    plt.subplot(2,2,1)
    for i in range(GL):
        plt.plot(x,phi[i], label=r'$\phi_{%i}$' %(i+1))
        plt.legend()
        plt.title('Basis Functions')

    plt.subplot(2,2,3)
    for i in range(GL):
        plt.plot(x,dphi[i], label=r"$\phi_{%i}'$" %(i+1))
        plt.legend()
        plt.title('Basis Functions Derivatives')


    plt.subplot(2,2,2)
    plt.plot(x,f_an,label='Analytical Function')
    plt.plot(x,f_ap,'--',label='Interpolation')
    plt.plot(xi, fi, 'x', label = r'$f_i$')
    plt.title(r'$f(x)$')
    plt.legend()

    plt.subplot(2,2,4)
    plt.plot(x,der,label='Analytical Derivative')
    plt.plot(x,der_ap,'--',label='Interpolation')
    plt.title(r"$f'(x)$")
    plt.legend()

    plt.show()

def test_2d():
    GL = 12
    #xi,w  = np.polynomial.legendre.leggauss(GL)
    xi = np.linspace(-2*np.pi,2*np.pi,GL)
    eta = np.linspace(0, 4*np.pi, GL)
    XI,ETA = np.meshgrid(xi,eta)
    N = 100
    x = np.linspace(-2*np.pi, 2*np.pi, N)
    y = np.linspace(0, 4*np.pi, N)
    X,Y = np.meshgrid(x,y)

    fun = lambda x,y: np.sin(x) + np.cos(y) + x
    dfundx = lambda x,y: np.cos(x) + 1
    dfundy = lambda x,y: -np.sin(y)

    phi, dphi_xi, dphi_eta = lagrange2d(xi,eta,X,Y)
    phi = np.asarray(phi)
    dphi_xi = np.asarray(dphi_xi)
    dphi_eta = np.asarray(dphi_eta)

    fi = fun(XI,ETA).reshape(GL**2)

    f_an = fun(X,Y)
    dfx_an = dfundx(X,Y)
    dfy_an = dfundy(X,Y)

    f_ap = np.zeros((N,N))
    dfx_ap = np.zeros((N,N))
    dfy_ap = np.zeros((N,N))
    for i in range(GL**2):
        f_ap += phi[i,:,:]*fi[i]
        dfx_ap += dphi_xi[i,:,:]*fi[i]
        dfy_ap += dphi_eta[i,:,:]*fi[i]

    plt.figure()
    plt.subplot(3,2,1)
    plt.contourf(X,Y,f_an)
    plt.title(r'Analytical')
    plt.ylabel(r'$f(x,y)$')

    plt.subplot(3,2,2)
    plt.contourf(X,Y,f_ap)
    plt.title(r'Approximation')

    plt.subplot(3,2,3)
    plt.contourf(X,Y,dfx_an)
    plt.ylabel(r'$\partial f(x,y)/ \partial x$')

    plt.subplot(3,2,4)
    plt.contourf(X,Y,dfx_ap)

    plt.subplot(3,2,5)
    plt.contourf(X,Y,dfy_an)
    plt.ylabel(r'$\partial f(x,y)/ \partial y$')

    plt.subplot(3,2,6)
    plt.contourf(X,Y,dfy_ap)
    plt.show()

#test_1d()
test_2d()
