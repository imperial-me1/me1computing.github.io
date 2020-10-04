#
# Golden ratio
# Vincent Legat - 2018
# Ecole Polytechnique de Louvain
#

import matplotlib 
from matplotlib import pyplot as plt
matplotlib.rcParams['toolbar'] = 'None'
plt.rcParams['figure.facecolor'] = 'silver'

#
# -0- Definition du nombre d'or
#

print("\n -0- Golden ratio number \n")
from math import sqrt

phi = (1 + sqrt(5))/2
print(phi)

# 
# -1- Polynomes avec python/numpy
#

print("\n -1- Golden ratio number as roots of a polynomial... \n")
from numpy import *

p = [1,-1,-1]
r = roots(p)
print(r)

plt.figure("Golden ratio : roots of a polynomial")
x = linspace(-2,2,100)
plt.plot(x,polyval(p,x))

#
# -2- Calcul symbolique avec python/sympy
#

print("\n -2- Golden ratio from a symbolic calculation  :-) \n")
from sympy import symbols,solve,evalf

x = symbols('x')
r = solve(1/x-x+1,x)
print(r)

phi = r[0]
print(phi)
print(phi.evalf())
print(phi.evalf(50))

#
# -3- Calcul numerique avec python/scipy.optimize
#

print("\n -3- Golden ratio from a numerical calculation (:) \n")
from scipy.optimize import fsolve 

f = lambda x : 1.0/x - x + 1.0
phi = fsolve(f,0.8)
print(phi[0])

plt.figure("Golden ratio : numerical computation")
x = linspace(0.1,4,100)
plt.xlim(( 0,4))
plt.ylim((-3,7))
plt.plot(x,f(x),'-b')
plt.plot(phi,0,'or',markersize=12)
plt.plot(0.8,0,'og',markersize=12)

#
# -4- Fractions continues
#

print("\n -4- Golden ratio from a continued fraction (./.) \n")
n = 6 
p = '1' 
for k in range(n):
  p = '1+1/(' + p + ')'
print(p)

phi = eval(p)
print(phi)
err = (1 + sqrt(5))/2 - phi
print(err)

#
# -5- Et encore un autre algorithme....
#

print("\n -5- And the last algorithm (../..) \n")
n = 6
p = 1
q = 1
for k in range(n):
    s = p
    p = p + q
    q = s
p = phi 

phi = eval('%d/%d' % (p,q))
print(phi)
err = (1 + sqrt(5))/2 - phi
print(err)

plt.show()