#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import scipy
import scipy.misc
import math

########################################################################
def factorial(n):
    #return scipy.special.gamma(n+1.0)
    return scipy.misc.factorial(n, exact=0)

########################################################################
def Laguerre(x, l_int):
    # http://en.wikipedia.org/wiki/Generalized_Laguerre_polynomials#Recursive_definition
    l = np.float64(l_int)
    L = 0.0
    if   (l_int == 0):
        L = 1.0
    elif (l_int == 1):
        L = 1.0 - x
    else:
        L = ( (2.0 * l - 1.0 - x) * Laguerre(x, l-1) - (l-1.0) * Laguerre(x, l-2) ) / l

    return L
# def Laguerre(x, l):

########################################################################
def GLaguerre(x, nu_int, alpha_int):
    # As defined by Friedrich, "Theoretical Atomic Physics", Section A.2
    nu = np.float64(nu_int)
    alpha = np.float64(alpha_int)
    GL = 0.0
    assert(alpha_int >= 0)

    # http://mathworld.wolfram.com/LaguerrePolynomial.html
    if   (alpha_int == 0):
        GL = Laguerre(x, nu_int)
    elif (nu_int == 0):
        GL = 1.0
    elif (nu_int == 1):
        #GL = 1.0 - x
        GL = 1.0 - x + alpha
    else:
        ## Friedrich, equation (A.15), page 441
        #GL = (
        #    (2.0*nu + alpha - 1.0 - x) * GLaguerre(x, nu_int-1, alpha_int)
        #    + (nu - 1.0 + alpha) * GLaguerre(x, nu_int-2, alpha_int)
        #) / nu

        # http://en.wikipedia.org/wiki/Generalized_Laguerre_polynomials#Recurrence_relations (Second one)
        #for i in range(nu_int+1):
        #    GL += GLaguerre(x, i, alpha_int-1)
        GL = GLaguerre(x, nu_int, alpha_int-1) + GLaguerre(x, nu_int-1, alpha_int)

    return GL
# def GLaguerre(x, l, alpha):

########################################################################
def phi_nl(r, n_int, l_int):
    assert(n_int > l_int)

    n = np.float64(n_int)
    l = np.float64(l_int)
    a = 1.0
    two_r_over_na = (2.0 * r) / (n * a)

    #print "n = ", n_int
    #print "l = ", l_int
    #print "n-l-1 = ", n_int-l_int-1
    #print "2*l+1 = ", 2*l+1
    #laguerre = special.genlaguerre(n_int-l_int-1, 2*l_int+1)

    # Friedrich, equation 1.138, p. 25
    phi_nl = \
          (1.0 / n) \
        * np.sqrt( factorial(n-l-1.0) / (a * factorial(n+l)) ) \
        * two_r_over_na**(l+1.0) \
        * GLaguerre(two_r_over_na, n_int-l_int-1, 2*l_int+1) \
        * np.exp(-two_r_over_na / 2.0)
        #* laguerre(two_r_over_na) \

    return phi_nl
# def phi_nl(r, n, l):

#x = np.linspace(0, 3*2**math.pi, 1000)
#y = np.sin(x)
#plt.plot(x,y)

r = np.linspace(0, 100.0, 1000)
n = 5
l = 1
phi = phi_nl(r, n, l)


matplotlib.rc('text', usetex=True)


fig = plt.figure(figsize=(6, 3))
ax  = fig.add_subplot(1, 1, 1)
plt.plot(r, phi, '-r', label = r'$|n,l>=|' + str(n) + ',' + str(l) + '>$', lw=3)
plt.plot(26.3582, 0.0, 'xm', ms=20, mew=3)
plt.legend(loc='best', shadow=True, fancybox=True)
plt.xlabel(r"$r$")
plt.ylabel(r"$\phi$")

plt.show()


