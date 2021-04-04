import numpy as np
import matplotlib.pyplot as plt
from math import * 
from scipy.integrate import *
from pylab import * 
from scipy.integrate import quad


MHD = np.zeros((80, 90, 5), dtype=float)
BGI = np.zeros((80, 90, 5), dtype=float)
Fp = np.zeros((80), dtype=float) 
AngMHD = np.zeros((90,2), dtype=float)
AngBGI = np.zeros((90,2), dtype=float)  
B0 = [0.5, 1.5, 3, 5, 10]
V = [0.3, 0.3, 0.2, 0.1, 0.1]


def PMHD(p, chi, b):
    return b**2/p*(1 +(sin(chi))**2)

def xMHD(p, chi, b):
    return -b**2/p**2*sin(chi)*cos(chi)

def PBGI(p, chi, b):
    Q = 0.7*p/b**0.57/sqrt(cos(chi))
    if Q > 1:
        A = 1
    else:
        A = Q
    return b**2/p*(A*(cos(chi))**2 + 0.01/sqrt(p))

def xBGI(p, chi, b):
    Q = 0.7*p/b**0.57/sqrt(cos(chi))
    if Q > 1:
        A = 1
    else:
        A = Q
    return A*b**2/p**2*sin(chi)*cos(chi)

P0 = 0.3
Pend = 1
B12 = 4

dx = 0.0001


for i in range(450):
   xi0 = i/5 + 0.1
   x0 = pi/180*xi0
   P = P0
   x = x0
   while 0.7*P/B12**0.57/sqrt(cos(x)) < 2:
      P = P + PMHD(P, x, B12)*dx
      x = x + xMHD(P, x, B12)*dx
      gx = 180/pi*x
      iP = int(P/0.1)
      ix = int(gx)
      if iP < 80:
        MHD[iP, ix, 0] = MHD[iP, ix, 0] + 1


for i in range(450):
   xi0 = i/5 + 0.1
   x0 = pi/180*xi0
   P = P0
   x = x0
   while 0.7*P/B12**0.57/sqrt(cos(x)) < 2:
      P = P + PBGI(P, x, B12)*dx
      x = x + xBGI(P, x, B12)*dx
      gx = 180/pi*x
      iP = int(P/0.1)
      ix = int(gx)
      if iP < 80:
        BGI[iP, ix, 0] = BGI[iP, ix, 0] + 1

#for j in range(80):
#   for i in range(90):
#      Fp[j] = Fp[j] + PxiB[j, i, 0] 
#   print(j/10, Fp[j])      


for i in range(90):
    j = int(10*Pend)
    AngMHD[i,0] = i
    AngBGI[i,0] = i
    AngMHD[i,1] = MHD[j, i, 0]
    AngBGI[i,1] = BGI[j, i, 0]
#    print(i, PxiB[10, i, 0])


ymax = np.max(AngBGI)

fig, ax = plt.subplots()
x = np.linspace(0, 90)
plt.xlim(1, 90)
plt.ylim(0, 1.2*ymax)
data1 = np.array(AngMHD)
data2 = np.array(AngBGI)
X1,Y1 = data1.T
X2,Y2 = data2.T
plt.scatter(X1,Y1, color = 'blue', s=15, label="MHD")
plt.scatter(X2,Y2, color = 'red', s=15, label="BGI")
plt.title('$P_0$ = '+str(P0)+', P = '+str(Pend)+', $B_{12}$ = '+str(B12)+'')
plt.grid(True,which="both", ls="-")
plt.grid(True,which="both", ls="-")
plt.xlabel('$\chi$')
#plt.ylabel('$\lambda g(x_{0})$')
plt.legend()
plt.show()    


#fig, ax = plt.subplots()
#x = np.linspace(0, 1)
#plt.xlim(0.0001, 1.0)
#plt.ylim(0, 0.1)
#plt.plot(x, x**2*(cos(ch)*(1 - x**2) + 1/2*sin(ch)*(x - x**3))**3, label="fitting")
#plt.title(''+str(PSR)+', $n_{\pm}$ (P = '+str(P)+', $B_{12}$ = '+str(B12)+', $\chi$ = '+str(chi)+'$^{\circ}$), $\lambda = 92$')
#plt.grid(True,which="both", ls="-")
#plt.grid(True,which="both", ls="-")
##ax.vlines(xcr, 0, 8, color = 'black', linewidth = 1.5, linestyle = '--')
#plt.xlabel('$r_{0}/R_0$')
#plt.ylabel('$n_{\pm}$')
#plt.legend()
#plt.show() 