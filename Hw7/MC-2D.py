import numpy as np
import matplotlib.pyplot as plt
import math as mt
import sys
from sympy import *

# Define the function to be integrated
def fun(x,y):
    return (mt.sqrt(x))*np.exp(-y)

# MC integral
def MCI(a, b, c, d, n, f, seed):
    np.random.seed(seed)
    sum1=0
    for i in range(n):
        xrandom=np.random.uniform(a,b)
        yrandom=np.random.uniform(c,d)
        sum1+=f(xrandom,yrandom)
    mean_sum=sum1/int(n)
    integralm=(b-a)*(d-c)*mean_sum
    return integralm
    
# Analytical solution
def exac_2(x,a,b,y,c,d):
    return integrate((sqrt(x)*exp(-y)),(x, a, b),(y, c, d) ).evalf(10)

# Integral limits
a=0
b=5
c=1
d=4
n=100
seed=33333
x=symbols('x')  
y=symbols('y')

exact=exac_2(x,a,b,y,c,d)
print('The Analytic solution=', exact)
print('===========================================================')

# List of sample sizes
n_list = [100, 1000, 10000, 1000000]

# Compute Monte Carlo approximation and error for each sample size
for n in n_list:
    mc = MCI(a, b, c, d, n, fun, seed)
    error = abs(mc - exact)
    print(f"n = {n}: Monte Carlo approximation = {mc:.9f}, error = {error:.9f}")


