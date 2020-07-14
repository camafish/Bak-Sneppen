import numpy as np
from math import exp
from math import log

def printout(output, name):
    print(name +  "\n   minimizer = (" + str(float(output[0])) +")\n   minimum = " + str(float(output[1]))+ "\n   number of iterations = " + str(len(output[2])-1))

def NewtonMethod(f, df, startpt, precision):
    x = startpt                             # the first x is set to the starting point
    history = [float(f(x))]                 # create list to record all f(x) values
    while abs(f(x)) > precision:         # stopping criterion
        #print(x, f(x), df(x))
        x = x - f(x)/df(x)           
        history.append(float(f(x)))         

    return (x, f(x), history) 

def make_4th_poly(a,b,c,d,e):
    def f(x):
        return a*(x)**4 + b*(x)**3 + c*(x)**2 + d*(x) + e
    return f

def make_4th_poly_prime(a,b,c,d,e):
    def df(x):
        return (4*a*(x)**3 + 3*b*(x)**2 + 2*c*(x) + d)
    return df

def make_4th_poly_prime_prime(a,b,c,d,e):
    def ddf(x):
        return 12*a*(x)**2 + 6*b*(x) + 2*c
    return ddf

a = .035
b = -.393
c = 2.01
d = -4.721
e = 7

p = make_4th_poly(a,b,c,d,e)
pp = make_4th_poly_prime(a,b,c,d,e)
ppp = make_4th_poly_prime_prime(a,b,c,d,e)

f = lambda x : exp(p(log(x)))
fp = lambda x : (1/x)*f(x)*pp(log(x))
fpp = lambda x : (1/x)*fp(x)*((ppp(log(x))/pp(log(x))) + pp(log(x)) - 1)

startpt = 16
precision = 10**(-6)

NM = NewtonMethod(fp, fpp, startpt, precision)
printout(NM, "Newton's Method")

# now NM[0] contains the zero of the derivative. if we plug that into f, we should get the minimum of f

print("via NM f is minimized at x = " + str(NM[0]) +" with value f(x) = "+str(f(NM[0])))
