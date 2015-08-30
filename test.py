from numpy import *
from matplotlib.pyplot import *
p = logspace(-6, -0.5, 101)
y = (1-exp(-p))/p
semilogx(p, y)
show()