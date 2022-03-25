from numpy import *
def solver(I, a, T, dt, theta):
    """Soluciona u'(t) = -a*u(t) con la condicion inicial u(0)=I, para t en el intervalo de tiempo (0, T]. El tiempo de divide en intervalos de longitud igual a dt. theta = 1 para el metodo hacia atrás de Euler, theta = 0 para el metodo hacia adelante de Euler y thena = 0.5 para el metodo de Crank-Nicolson."""
    dt = float(dt)
    Nt = int(round(T/dt))
    T = Nt*dt
    u = zeros(Nt+1)
    t = linspace(0, T, Nt+1)
    u[0]=I
    for n in range(0, Nt):
        u[n+1] = (1 - (1-theta)*a*dt)/(1 + theta*a*dt)*u[n]
    for i in range(len(t)):
        print ('t=%3.2f u=%1.9f' % (t[i], u[i]))
    return u, t
def u_exact(t, I, a):
    return I*exp(-a*t)
from matplotlib.pyplot import *
def plot_numerical_and_exact(theta, I, a, T, dt):
    """Compara las soluciones numérica y exacta en un grafico."""
    u, t = solver(I=I, a=a, T=T, dt=dt, theta=theta)
    t_e = linspace(0, T, 1001)
    u_e = u_exact(t_e, I, a)
    plot(t, u, 'r--o', t_e, u_e, 'b-')
    legend(['numerical', 'exact'])
    xlabel('t')
    ylabel('u')
    title('theta=%g, dt=%g' % (theta, dt))
    savefig('plot_%s_%g.png' % (theta, dt))
plot_numerical_and_exact(I=1, a=2, T=8, dt=0.8, theta=0)
show()