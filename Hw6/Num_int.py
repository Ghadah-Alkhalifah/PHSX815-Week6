import numpy as np
import matplotlib.pyplot as plt
import math as mt
from scipy import integrate
import sys


# Define the function to be integrated
def fun(x):
    return (mt.sqrt(2*x))*np.exp(-x)

# Trapezoidal rule
def trap(a, b, n, f):
    h=(b-a)/n
    y=f(a)+f(b)
    for i in range (1, n):
        x=a+i*h
        y+=2*f(x)
        Tot=0.5*h*y
    return Tot

# Gaussian quadrature with Legendre 
def gaus_quad(a, b, n, f):
    x, w = np.polynomial.legendre.leggauss(n)
    Integral=0
    for i in range (n):
        x_val=f((x[i]*(b-a)/2) + (b+a)/2)
        Integral+=x_val*0.5*(b-a)*w[i]
    return Integral

# Analytic solution        
def exac_f(a,b):
    Analytic = integrate.quad(fun,a,b)
    return Analytic

# Integral limits
a=0
b=5
n=100

# if the user includes the flag -h or --help print the options
if '-h' in sys.argv or '--help' in sys.argv:
    print ("Usage: %s [-a] lower limit [-b] upperlimit [-n] partition" % sys.argv[0])
    print
    sys.exit(1)
# read the user-provided seed from the command line (if there)
if '-a' in sys.argv:
    p = sys.argv.index('-a')
    a = int(sys.argv[p+1])
    
if '-b' in sys.argv:
    p = sys.argv.index('-b')
    b = int(sys.argv[p+1])
    
if '-n' in sys.argv:
    p = sys.argv.index('-n')
    n = int(sys.argv[p+1])

# Compute the numerical integrals and the error
trapzoid=trap(a,b,n,fun)
gaussian=gaus_quad(a, b, n, fun)
exact=exac_f(a,b)
erro_trap = abs(exact[0]- trapzoid)
erro_gaus= abs(exact[0]- gaussian)

# Print the results
print()
print('Analytic Integral  Trapezoidal Integral  Gaussian Integral')

print("{:.8f}           {:.8f}           {:.8f} ".format(exact[0], trapzoid, gaussian))
print('===========================================================')
print('Trapezoidal Error  Gaussian Error')
print("{:.8f}            {:.8f}".format(erro_trap, erro_gaus))
print('===========================================================')
