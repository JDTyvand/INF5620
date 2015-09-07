from numpy import *
from matplotlib.pyplot import *


def mesh_function(f, t):
	tmp = zeros(len(t))
	for i in range(len(t)):
		tmp[i] = f(t[i])
	return tmp

def f(t):
	if t >= 0 and t <= 3:
		return exp(-t)
	if t > 3 and t <= 4:
		return exp(-3*t)
	else:
		return 0

def main(T, dt):
	dt = float(dt)
	Nt = int(round(T/dt))
	T = Nt*dt
	t = linspace(0, T, Nt+1)

	m_a = mesh_function(f, t)

	plot(t, m_a)
	show()

main(4, 0.1)