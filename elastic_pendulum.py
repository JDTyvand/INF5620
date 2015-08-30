from numpy import *
import matplotlib.pyplot as plt

def simulate(
    beta=0.9,                 # dimensionless parameter
    Theta=30,                 # initial angle in degrees
    epsilon=0,                # initial stretch of wire
    num_periods=6,            # simulate for num_periods
    time_steps_per_period=60, # time step resolution
    plot=True,                # make plots or not
    ):

    P = 2*pi
    T = P*num_periods
    dt = P/time_steps_per_period
    Nt = T/dt
    x = zeros(Nt+1)
    y = zeros(Nt+1)
    t = linspace(0, T, Nt+1)
    theta = zeros(Nt+1)
    x[0] = (1 + epsilon)*sin(Theta)
    y[0] = 1 - (1 + epsilon)*cos(Theta)
    theta[0] = Theta
    L = sqrt(x[0]**2 + (y[0] - 1)**2)
    x[1] = x[0] - (dt**2/2)*(beta/(1 - beta) * (1 - beta/L))*x[0]
    y[1] = y[0] - (dt**2/2)*((beta/(1 - beta) * (1 - beta/L))*(y[0] - 1)) - (dt**2/2)*beta

    for i in range(len(x)-1):
   		L = sqrt(x[i]**2 + (y[i] - 1)**2)
   		x[i+1] = 2*x[i] - x[i-1] - dt**2*(beta/(1 - beta) * (1 - beta/L))*x[i]
   		y[i+1] = 2*y[i] - y[i-1] - dt**2*(beta/(1 - beta) * (1 - beta/L))*(y[i] - 1) - (dt**2*beta)
   		theta[i+1] = 1/tan(x[i+1]/(1-y[i+1]))

    if plot == True:
   		plt.figure()
   		plt.plot(x, y)
   		plt.gca().set_aspect('equal')
   		plt.figure()
   		plt.plot(t, theta)
   		plt.show()

    return x, y, theta, t

def zero_test(Theta):
	x, y , theta, t = simulate(Theta=Theta)
	print x.max()
	print y.max()
	if abs(x.max()) > 0 or abs(y.max()) > 0:
		print "Error in the numerical scheme!"
	else:
		print "Theta = 0 and epsilon = 0 gives x = y = 0 for all times, as intended."
def main():
	simulate()
	zero_test(0)
main()
