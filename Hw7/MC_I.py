import numpy as np
import matplotlib.pyplot as plt
import math as mt
#from scipy import integrate
import sys
from sympy import *




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

# MC integral
def MCI(a, b, n, f, seed):
    #seed=33333 defult seed
    np.random.seed(seed)
    sum1=0
    for i in range(n):
        xrandom=np.random.uniform(a,b)
        sum1+=f(xrandom)
    mean_sum=sum1/int(n)
    integralm=(b-a)*mean_sum
    return integralm
    


# Analytical solution
def exac_2(x,a,b):
    
    return integrate((sqrt(2*x)*exp(-x)),(x, a, b)).evalf(10)
    

# Integral limits
a=0
b=5
n=100
seed=33333
x=symbols('x')

# if the user includes the flag -h or --help print the options
if '-h' in sys.argv or '--help' in sys.argv:
    print ("Usage: %s [-a] lower limit [-b] upperlimit [-n] partition  [-seed] seed number" % sys.argv[0])
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

if '-seed' in sys.argv:
        p = sys.argv.index('-seed')
        seed =int(sys.argv[p+1])

# Compute the numerical integrals and the error
trapzoid=trap(a,b,n,fun)
gaussian=gaus_quad(a, b, n, fun)
mc=MCI(a,b,n,fun, seed)
exact2=exac_2(x,a,b)

# compute errors
erro_trap = abs(exact2- trapzoid)
erro_gaus= abs(exact2- gaussian)
erro_mc= abs(exact2- mc)

# Print the results
print('==========================================================================')
print('Analytic Integral  Trapezoidal Integral  Gaussian Integral   Monte carlo')

print("{:.9f}           {:.9f}           {:.9f}       {:.9f} ".format(exact2, trapzoid, gaussian, mc))
print('==========================================================================')
print('Trapezoidal Error  Gaussian Error      Monte carlo Error')
print("{:.9f}            {:.9f}         {:0.9f}".format(erro_trap, erro_gaus, erro_mc ))
print('==========================================================================')

