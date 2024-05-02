import numpy as np
import matplotlib.pyplot as plt


x = np.linspace(-0.2, 0.2, 10000)
f = 5**3/(5**2 + x**2)
g = 5**3/(5**2 + (x - 0.025*0.2)**2)

plt.figure()
# plt.plot(x, f)
# plt.plot(x, g)
plt.plot(x, g-f)
plt.grid()
plt.show()