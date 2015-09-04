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
    
    Theta_s = Theta
    Theta = Theta*pi/180
    
    omega = 1.
    P = 2*pi/omega
    dt = P/time_steps_per_period
    T = num_periods*P
    Nt = num_periods*time_steps_per_period
    x = zeros(Nt+1)
    y = zeros(Nt+1)
    t = linspace(0, T, Nt+1)
    x[0] = (1 + epsilon)*sin(Theta)
    y[0] = 1 - (1 + epsilon)*cos(Theta)
    L = sqrt(x[0]**2 + (y[0] - 1)**2)
    C = beta/(1-beta)*(1 - beta/L)
    x[1] = (1 - 0.5*C*dt**2)*x[0]
    y[1] = (1 - 0.5*C*dt**2)*y[0] + 0.5*(C - beta)*dt**2

    for i in range(1,len(x)-1):
      L = sqrt(x[i]**2 + (y[i] - 1)**2)
      C = beta/(1-beta)*(1 - beta/L)
      x[i+1] = (2 - C*dt**2)*x[i] - x[i-1]
      y[i+1] = (2 - C*dt**2)*y[i] - y[i-1] + (C - beta)*dt**2
    
    theta = arctan(x/(1-y))

    if plot == True:
      plt.figure()
      plt.plot(x, y, 'b-')
      plt.title('Motion of pendulum')
      plt.xlabel('x')
      plt.ylabel('y')
      plt.savefig('Motion_%d.png' % Theta_s)
      #plt.gca().set_aspect('equal')
      plt.figure()
      plt.plot(t, theta*180/pi, 'b-')
      plt.title('Angle in degrees')
      plt.xlabel('t')
      plt.ylabel('theta')
      plt.savefig('Angle_%d.png' % Theta_s)
    if abs(Theta) < 10*pi/180:
      theta_ne = Theta*cos(t)
      plt.figure()
      plt.plot(t, theta, t, theta_ne)
      plt.title('Angle elastic vs. non-elastic')
      plt.xlabel('t')
      plt.ylabel('theta')
      plt.savefig('Elastic_vs_non_%d.png' % Theta_s)

    return x, y, theta, t

def pure_vertical():
    beta = 0.9
    omega = sqrt(beta/(1-beta))
    P = 2*pi/omega
    N = 6
    num_periods = N/omega
    n = 600
    time_steps_per_period = omega*n
    epsilon = 0.1

    def y_exact(t):
      """
      Returns the exact solution for the position of y,
      with a minus sign due to a positive epsilon stretching
      the pendulum in the negative y direction.
      """
      return -epsilon*cos(omega*t)

    x, y, theta, t = simulate(Theta=0, epsilon=epsilon, num_periods=num_periods, time_steps_per_period=time_steps_per_period, plot=True)
    tol = 0.0001
    y_e = y_exact(t)
    assert abs(x.max()) < tol
    err = abs(y_e - y).max()
    assert err < tol

def zero_test():
    x, y , theta, t = simulate(Theta=0)
    if abs(x.max()) > 0 or abs(y.max()) > 0:
		  print "Error in the numerical scheme!"
    else:
		  print "Theta = 0 and epsilon = 0 gives x = y = 0 for all times, as intended."

def demo(beta=0.9, Theta=60, num_periods=3):
    x, y, theta, t = simulate(beta=beta, Theta=Theta, num_periods=num_periods, time_steps_per_period=600, plot=True)
def main():
    simulate()
    simulate(Theta=8)
    zero_test()
    pure_vertical()
    demo()
main()
