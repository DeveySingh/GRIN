# Main file

# from geometry import *
import numpy as np
import matplotlib.pyplot as plt
uM = 1e-6 #micro metre
nM = 1e-9 #nano metre
epsilon = 1*uM 
n0 = 1.3
sqrt_A = np.sqrt(184615384.61538467)
A = sqrt_A**2
""" 
Temporary implementation of the algorithm
"""
# Given/known quantities
def n(x, y, n0=1.3, sqrt_a=sqrt_A):
    """ 
    Gradient index profile
    modified Selfoc SML W10
    refractive index changes from 1 to 1.607 in the radius of 39.5 uM
    """
    r = np.sqrt(x**2 +  y**2)
    return n0 * (1 - 0.5*(sqrt_a*r)**2)

def D(x, y):
    """ 
    Gradient of the the refractive index
    """
    return n0*A*np.array([-x, -y])

def phi(x, y)

# print(n(0.5))
x = np.arange(-50, 50.1, 0.1)*uM
y = np.arange(-50, 50.1, 0.1)*uM
X, Y = np.meshgrid(x, y)
# r = np.sqrt(X**2 + Y**2)

n_profile = n(X, Y)
n_profile = np.clip(n_profile, 1, 1.607)
# print(X, Y)
# print(X[0, 0], Y[1, 0])
plt.title("Index Profile")
plt.imshow(n_profile, extent=[-50, 50, -50, 50])
plt.colorbar()
plt.xlabel("x(uM)")
plt.ylabel("y(uM)")
plt.show()


# n_profile[500, 0]
# 2*(1-1/1.607)/(((50*uM)**2))


if __name__ == "__main__":
    pass
    # c = Cuboid(10,4,10)
    # print(c.GetSurfaceArea())
    # print(c.GetVolume())