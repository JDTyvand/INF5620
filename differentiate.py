from numpy import *

def differentiate(u, dt):
	d = zeros(len(u))
	for i in range(1, len(u)-1):
		d[i] = (u[i+1]-u[i-1])/(2*dt)
	d[0] = (u[1]-u[0])/dt
	d[-1] = (u[-1]-u[-2])/dt
	return d

def test_differentiate():
	t = linspace(0,6,13)
	u = 10*t + 5
	dt = t[1]-t[0]	
	d_c = differentiate(u, dt)
	d_e = 10
	error = abs(d_c - d_e).max()
	Tol=1E-15
	assert error<Tol,"Differentiation formula is incorrect"
	print "Test passed with error %15f" %error

def main():
	test_differentiate()

main()
