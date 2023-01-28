import numpy as np
import matplotlib.pyplot as plt

N=10**5
tb=0
te=2
phi0=0
omega0=0
dt=(tb-te)/N
A=-11067
a=3.05*10**-2
b=4.366
c=4.449
#d=0.03037
d=0.0001

x0=0
y0=0
vx0=0
vy0=0

t=np.linspace(tb,te,N)
omega=np.zeros(N)
omega[0]=omega0
phi=np.zeros(N)
phi[0]=phi0

x=np.zeros(N)
x[0]=x0
y=np.zeros(N)
y[0]=y0
vx=np.zeros(N)
vx[0]=vx0
vy=np.zeros(N)
vy[0]=vy0

def f(u):
    return A*np.cos(u)

def fx(u,du):
    return (a*du**2-b-c*np.sin(u))*np.cos(u)+d*f(u)*np.sin(u)

def fy(u,du):
    return (a*du**2-b-c*np.sin(u))*np.sin(u)-d*f(u)*np.cos(u)

for i in range(1,N):
    omega_dumb=omega[i-1]+f(phi[i-1])*dt
    phi_dumb=phi[i-1]+omega_dumb*dt
    omega[i]=omega[i-1]+f((phi_dumb+phi[i-1])/2)*dt
    phi[i]=phi[i-1]+(omega[i]+omega[i-1])*dt/2
    
    _phi=(phi[i]+phi[i-1])/2
    _omega=(omega[i]+omega[i-1])/2
    vx[i]=vx[i-1]+fx(_phi,_omega)*dt
    x[i]=x[i-1]+(vx[i]+vx[i-1])*dt/2
    vy[i]=vy[i-1]+fy(_phi,_omega)*dt
    y[i]=y[i-1]+(vy[i]+vy[i-1])*dt/2

plt.plot(t,x)
plt.show()
