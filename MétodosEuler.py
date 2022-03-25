from pandas import *
from matplotlib import *
from numpy import *


def solver(I, a, T, dt, theta):
    """Soluciona u'(t) = -a*u(t) con la condicion inicial u(0)=I, 
    para t en el intervalo de tiempo (0, T]. El tiempo de divide 
    en intervalos de longitud igual a dt. theta = 1 para el metodo 
    hacia atrás de Euler, theta = 0 para el metodo hacia adelante 
    de Euler y theta = 0.5 para el metodo de Crank-Nicolson."""
    dt = float(dt)
    Nt = int(round(T/dt))
    T = Nt*dt
    u = zeros(Nt+1)
    t = linspace(0, T, Nt+1)
    u[0]=I
    for n in range(0, Nt):
        u[n+1] = (1 - (1-theta)*a*dt)/(1 + theta*a*dt)*u[n]
    print ('...............................')
    print ('tiempo   incremento')    
    print ('...............................')
    for i in range(len(t)):
        print ('t = %4.3f u = %12.6f' % (t[i], u[i]))
    return u, t
def u_exact(t, I, a):
    return I*exp(-a*t)
from matplotlib.pyplot import *
def explore(I, a, T, dt, theta, makeplot=True):
    """ejecuta, determina la medida del error y dibuja las soluciones numérica y exacta (si makeplot = True)."""
    u, t = solver(I, a, T, dt, theta)
    u_e = u_exact(t, I, a)
    e = u_e - u
    E = sqrt(dt*sum(e**2))
    if makeplot:
        figure()
        t_e = linspace(0, T, 1001)
        u_e = u_exact(t_e, I, a)
        plot(t, u, 'r--o', t_e, u_e, 'b-')
        legend(['numerica', 'exacta'])
        xlabel('t')
        ylabel('u')
        theta2name = {0: 'Forward_Euler', 1: 'Backward_Euler', 0.5: 'Crank-Nicolson'}
        title('Grafico del metodo %s para theta=%g, dt=%g' % (theta2name[theta],theta, dt))
        savefig('%s_%g.png' % (theta2name[theta], dt))
        show()
    return E
def main(I, a, T, dt_values, theta_values=(0, 1, 0.5)):                                 
    theta2name = {0: 'Forward Euler', 1: 'Backward Euler', 0.5: 'Crank-Nicolson'}
    for theta in theta_values:
        print ('Calculos para el metodo %s' % (theta2name[theta]))                           
        for dt in dt_values:
            E = explore(I, a, T, dt, theta, makeplot=True)
            print ('---------------------------------')
            print ('theta   dt    error')
            print ('%4.3f  %6.3f  %12.5E' % (theta, dt, E))
            print ('---------------------------------')            
        print ('Fin de los calculos para el metodo %s' % (theta2name[theta])) 
        print (' ')           
main(I=1, a=2, T=8, dt_values=[0.2, 0.4, 0.8]) #Entre corchetes ponga los valores de theta separados por coma