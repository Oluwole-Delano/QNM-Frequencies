"""
19th August 2025
KCL Lionel London Project

This script aims to calculate the Quasi-Normal mode frequencies of Kerr 
(and Scwharzchild as a special case of Kerr) black hole pertubations 
using the equations outlined in Leaver 1985. 

"""

import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from scipy import optimize 
from mpmath import findroot 
import math 

# Scipi/numpy root finder

# 1. Definitions and Inputs

    # Define quantities
s = -2 # Field Spin-Weight for gravitaional fields 
l = 2
m = 0
a = 0 # Angular momentum per unit mass, 0 <= a <= 1/2
Alm = l*(l+1) - s*(s+1) # Value of Alm for a = 0 (Schwarzchild)
b = (1 - 4*a**2)**0.5 # Auxiliary rotation parameter 
j = complex(0,1)
N = 10

# Check if parameters have reasonable values
def check_parameters(a,l,m,s):
    if (a > 0.5) or (a < 0): return False
    if (abs(m) > l): return False
    if (abs(s) != 0.5) and (abs(s) != 1) and (abs(s) != 2): return False
    else: return True

    # Define equations

# Define as quantities

# Leaver 20
k_1 = 0.5 * abs(m-s) 
k_2 = 0.5 * abs(m+s) 

def alpha_theta_n(n):
    return -2*(n+1)*(n + 2*k_1 + 1)
def beta_theta_n(n,w):
    return n*(n-1) * 2*n*(k_1 + k_2 + 1 - 2*a*w) - [2*a*w*(2*k_1 + s + 1) \
             - (k_1 + k_2)*(k_1 + k_2 + 1) - a**2 * w**2 + s*(s+1) + Alm]
def gamma_theta_n(n,w):
    return 2*a*w*(n + k_1 +k_2 + s)

# Leaver 26
def c_0(w):
    return 1 - s - w*j - (2*j/b) * (w/2 - a*m) 
def c_1(w):
    return -4 + 2*j*w*(2 + b) + (4*j/b) * (w/2 - a*m)
def c_2(w):
    return s + 3 - 3*j*w - (2*j/b)*(w/2 - a*m)
def c_3(w):
    return w**2 * (4 + 2*b - a**2) - 2*a*m*w - s - 1 + (2+b)*j*w - Alm + ((4*w + 2*j)/b)*(w/2 - a*m)
def c_4(w):
    return s + 1 -2*w**2 - (2*s+3)*j*w - ((4*w + 2*j)/b)*(w/2 - a*m)

# Leaver 25
def alpha_r_n(w,n):
    return n**2 +(c_0(w) + 1)*n + c_0(w) 
def beta_r_n(w,n):
    return -2*n**2 + (c_1(w) + 2)*n + c_3(w) 
def gamma_r_n(w,n):
    return n**2 + (c_2(w) - 3)*n + c_4(w) - c_2(w) + 2 

# 2. Implement method to calculate roots of continued fraction equation

    # Define continued fraction equations
def continued_fraction_Alm():
    pass

# Leaver 13
def continued_fraction(w,n):
    if (n == 0):
        return ( alpha_r_n(w,0) * gamma_r_n(w,1) ) / beta_r_n(w,0)
    
    return alpha_r_n(w,n-1) * ( gamma_r_n(w,n) / ( beta_r_n(w,n-1) - continued_fraction(w,n-1) ) )

def w_equation(w,n=N):
    return beta_r_n(w,0) - continued_fraction(w,n) 

def complex_w_equation(w,n=N):
    return [np.real(w_equation(w,n)), np.imag(w_equation(w,n))] 

# Implement leaver 14

    # Countour plot in imaginary plane of abs(function)
I = 500
lim = 5
x, y = np.meshgrid(np.linspace(-1,1,I),np.linspace(-5,0,I))
w = x + 1j*y
values = w_equation(w,N)

fig, ax = plt.subplots()
contour_plot = ax.contourf(x,y,np.log(abs(values)),50)
fig.colorbar(contour_plot)

fig2, ax2 = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax2.plot_surface(w.real, w.imag, np.log(abs(values)), cmap=cm.coolwarm, 
                       linewidth=0, antialiased=False)
# Customize the z axis.
ax2.set_zlim(-2, 2)
ax2.zaxis.set_major_locator(LinearLocator(10))
ax2.zaxis.set_major_formatter('{x:.02f}')
# Add a color bar which maps values to colors.
fig2.colorbar(surf, shrink=0.5, aspect=5)

    # Root finding for w
w_initial_guess = -1.28 - 1j*0.8
w_roots = findroot(w_equation, w_initial_guess, solver='muller')
print("Did this finally frickin work: ", w_roots)
print(w_equation(w_roots))
#w_roots = optimize.root(complex_w_equation, w_initial_guess, args=(N), method='lm')
#print("The roots of the equation are: ", w_roots)

    # Plot w_root on the two graphs
#ax2.scatter(w_roots.real,w_roots.imag,0,color='red',marker='x')
w_reference = 0.7473433688360835-0.17792463137787093*j
ax.scatter(w_reference.real,w_reference.imag,color='red',marker='x')
#ax.scatter(w_roots.real,w_roots.imag,color='red',marker='x') 

plt.show()

# .real/imag on a func
# root

# 3. Output the results of the calculation

    # Store real and Imaginary parts of the frequency in a csv file

    # Plot an argand diagram for the calculated frequencies, may be able to find some cool 
    # visuals for this 