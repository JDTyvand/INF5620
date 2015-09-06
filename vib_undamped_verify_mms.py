import sympy as sym
from numpy import *
from matplotlib.pyplot import *
b, c, V, t, I, w, dt = sym.symbols('b c V t I w dt')  # global symbols
f = None  # global variable for the source term in the ODE

def ode_source_term(u):
    """Return the terms in the ODE that the source term
    must balance, here u'' + w**2*u.
    u is symbolic Python function of t."""
    return sym.diff(u(t), t, t) + w**2*u(t)

def residual_discrete_eq(u):
    """Return the residual of the discrete eq. with u inserted."""
    num = DtDt(u, dt) + w**2*u(t)
    R = ode_source_term(u) - num
    return sym.simplify(R)

def residual_discrete_eq_step1(u):
    """Return the residual of the discrete eq. at the first
    step with u inserted."""
    f_0 = ode_source_term(u).subs(t,0)
    exact = (ode_source_term(u).subs(t, dt) - sym.diff(u(t), t, t).subs(t, dt))/w**2
    R = exact - ((1 - (w**2*dt**2)/2)*I + dt*V + dt**2*f_0/2)
    return sym.simplify(R)

def DtDt(u, dt):
    """Return 2nd-order finite difference for u_tt.
    u is a symbolic Python function of t.
    """
    return (u(t+dt) - 2*u(t) + u(t-dt))/(dt**2)

def solver(I, W, V, T, step):
    """
    Solver function for computing the numerical solution,
    given initial conditions u(0)=I, u'(0)=V and frequency W. 
    """
    step = float(step)              # avoid integer division
    Nt = int(round(T/step))         # no of time intervals
    T = Nt*step                     # adjust T to fit time step dt
    u_values = zeros(Nt+1)          # array of u[n] values
    times = linspace(0, T, Nt+1)

    def f(time):
        return sin(time)

    u_values[0] = I
    u_values[1] = ((1 - (W**2*step**2)/2)*I + step*V + step**2*f(0)/2)
    for i in range(1,len(times)-1):
        u_values[i+1] = (2 - W*step**2)*u_values[i] - u_values[i-1] + step**2*f(i*step)

    plot(times, u_values)
    show()

def main(u):
    """
    Given some chosen solution u (as a function of t, implemented
    as a Python function), use the method of manufactured solutions
    to compute the source term f, and check if u also solves
    the discrete equations.
    """
    print '=== Testing exact solution: %s ===' % u
    print "Initial conditions u(0)=%s, u'(0)=%s:" % (u(t).subs(t, 0), sym.diff(u(t), t).subs(t, 0))

    # Method of manufactured solution requires fitting f
    global f  # source term in the ODE
    f = sym.simplify(ode_source_term(u))

    # Residual in discrete equations (should be 0)
    print 'residual step1:', residual_discrete_eq_step1(u)
    print 'residual:', residual_discrete_eq(u)

def linear():
    main(lambda t: V*t + I)
def quadratic():
    main(lambda t:b*t**2 + V*t + I)
def cubic():
    main(lambda t:c*t**3 + b*t**2 + V*t + I)

if __name__ == '__main__':
    #linear()
    #quadratic()
    #cubic()
    solver(I=1, W=2, V=0, T=4, step=0.1)
