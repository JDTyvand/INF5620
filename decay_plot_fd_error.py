from numpy import *
from matplotlib.pyplot import *
p = logspace(-6, -0.5, 101)
FE = (1-exp(-p))/p
BE = (exp(p)-1)/p
CN = (exp(p/2) - exp(-p/2))/p
semilogx(p, FE, p, BE, p, CN)
show()

from sympy import *

s = Symbol('s')
FE = (1-exp(-s))/s
BE = (exp(s)-1)/s
CN = (exp(s/2) - exp(-s/2))/s

print FE.series(s, 0, 6)
print BE.series(s, 0, 6)
print CN.series(s, 0, 6)