from numpy import *
from matplotlib.pyplot import *

def solver(I, a , T, dt, theta):
	dt = float(dt)
	Nt = int(round(T/dt))
	T = Nt*dt
	u = zeros(Nt+1)
	t = linspace(0, T, Nt+1)

	u[0] = I
	for n in range(0, Nt):
		u[n+1] = (1 - (1-theta)*a*dt)/(1 + theta*dt*a)*u[n]
	return u, t

def exact_solution(t_e, I, a):
	return I*exp(-a*t_e)

def main(I, a, T):
	theta_values = [0, 0.5, 1]
	for theta in theta_values:
		print theta
		dt_values = [1.25]
		for dt in dt_values:
			u, t = solver(I, a, T, dt, theta)
			t_e = linspace(0, T, 1001)
			u_e = exact_solution(t_e, I, a)

			plot(t_e, u_e, 'b-', t, u, 'r--o')
			legend(['exact', 'numerical'])
    		xlabel('t')
    		ylabel('u')
    		title('theta=%g, dt=%g' % (theta, dt))
    		show()
main(1, -1, 8)
