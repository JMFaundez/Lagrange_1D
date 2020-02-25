import numpy as np
import matplotlib.pyplot as plt
from lagrange import lagrange1d


GLL = 8
xi,w  = np.polynomial.legendre.leggauss(GLL)

f = np.sin(xi**2)

N = 50
x = np.linspace(-1, 1, N)

phi, dphi = lagrange1d(xi,x)
dphi = np.array(dphi)

der = np.cos(x**2)*2*x
f_an = np.sin(x**2)
der_ap = f@dphi
f_ap = f@phi

plt.figure(figsize=(15,8))
plt.subplot(2,2,1)
for i in range(GLL):
    plt.plot(x,phi[i], label=r'$\phi_{%i}$' %(i+1))
plt.legend()
plt.title('Basis Functions')

plt.subplot(2,2,2)
plt.plot(x,f_an,label='Analytical Function')
plt.plot(x,f_ap,'--',label='Interpolation')
plt.title(r'$f(x)=sin(x^2)$')
plt.legend()

plt.subplot(2,2,4)
plt.plot(x,der,label='Analytical Derivative')
plt.plot(x,der_ap,'--',label='Interpolation')
plt.title(r"$f'(x)=2xcos(x^2)$")
plt.legend()

plt.show()
