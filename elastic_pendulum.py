# coding: utf-8
import odespy
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
    
    Theta = Theta*pi/180
    """
    P = 2*pi
    T = P*num_periods
    dt = P/time_steps_per_period
    Nt = T/dt
    x = zeros(Nt+1)
    y = zeros(Nt+1)
    t = linspace(0, T, Nt+1)
    x[0] = (1 + epsilon)*sin(Theta)
    y[0] = 1 - (1 + epsilon)*cos(Theta)
    L = sqrt(x[0]**2 + (y[0] - 1)**2)
    x[1] = x[0] - 0.5*(dt**2)*(beta/(1 - beta))*(1 - (beta/L))*x[0]
    y[1] = y[0] - 0.5*(dt**2)*(beta/(1 - beta))*(1 - (beta/L))*(y[0] - 1) - (dt**2/2)*beta

    for i in range(len(x)-1):
   		L = sqrt(x[i]**2 + (y[i] - 1)**2)
   		x[i+1] = 2*x[i] - x[i-1] - dt**2*(beta/(1 - beta))*(1 - beta/L)*x[i]
   		y[i+1] = 2*y[i] - y[i-1] - dt**2*(beta/(1 - beta))*(1 - beta/L)*(y[i] - 1) - (dt**2*beta)
   	"""

    ic = [0,                          # x’=vx
         (1 + epsilon)*sin(Theta),       # x
         0,                              # y’=vy
         1 - (1 + epsilon)*cos(Theta),   # y
         ]

    def f(u, t, beta):
        vx, x, vy, y = u
        L = np.sqrt(x**2 + (y-1)**2)
        h = beta/(1-beta)*(1 - beta/L)  # help factor
        return [-h*x, vx, -h*(y-1) - beta, vy]

    # Non-elastic pendulum (scaled similarly in the limit beta=1)
    # solution Theta*cos(t)
    P = 2*pi
    dt = P/time_steps_per_period
    T = num_periods*P
    omega = 2*pi/P
    time_points = np.linspace(
        0, T, num_periods*time_steps_per_period+1)
    
    solver = odespy.EulerCromer(f, f_args=(beta,))
    solver.set_initial_condition(ic)
    u, t = solver.solve(time_points)
    x = u[:,1]
    y = u[:,3]

    theta = arctan(x/(1-y))

    if plot == True:
   		plt.figure()
   		plt.plot(x, y)
   		plt.gca().set_aspect('equal')
   		plt.figure()
   		plt.plot(t, theta*180/pi)
   		plt.show()
    if abs(Theta) < 10*pi/180:
      theta_exact = Theta*cos(t)
      plt.figure()
      plt.plot(t, theta, t, theta_exact)

    return x, y, theta, t

def pure_vertical():
    beta = 0.9
    omega = sqrt(beta/(1-beta))
    P = 2*pi/omega
    N = 5
    num_periods = N/omega
    n = 600
    time_steps_per_period = omega*n

    def y_exact(t):
      return -0.1*cos(omega*t)

    x, y, theta, t = simulate(Theta=0, epsilon=0.1, num_periods=num_periods, time_steps_per_period=time_steps_per_period, plot=False)
    tol = 0.0005
    assert abs(x.max()) < tol
    y_e = y_exact(t)
    err = abs(y_e - y).max()
    print err
    assert err < tol

def zero_test():
    x, y , theta, t = simulate(Theta=0)
    print x.max()
    print y.max()
    if abs(x.max()) > 0 or abs(y.max()) > 0:
		  print "Error in the numerical scheme!"
    else:
		  print "Theta = 0 and epsilon = 0 gives x = y = 0 for all times, as intended."

def demo(beta=0.999, Theta=40, num_periods=3):
    x, y, theta, t = simulate(beta=beta, Theta=Theta, num_periods=num_periods, time_steps_per_period=600, plot=True)
def main():
    #simulate()
    #zero_test()
    #pure_vertical()
    demo()
main()
