import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

N=10**5
tb=0
te=2
dt=(tb-te)/N
#A=-11067*np.pi/180
A=-200 
a=4.05*10**-2
b=4.366
c=4.449+2
#d=0.03037
d=0.0001
l2=0.016
t=np.linspace(tb,te,N)

phi0=0
omega0=0
x0=0
y0=0
vx0=0
vy0=0

def _autism(t,u):
    phi,phidot=u
    phi2dot=A*np.cos(phi)
    return [phidot,phi2dot]

def autism(t,u):
    phi,phidot,x,xdot,y,ydot=u
    phi2dot=A*np.cos(phi)
    x2dot=(a*phidot**2-b-c*np.sin(phi))*np.cos(phi)+d*phi2dot*np.sin(phi)
    y2dot=(a*phidot**2-b-c*np.sin(phi))*np.sin(phi)-d*phi2dot*np.sin(phi)
    return [phidot,phi2dot,xdot,x2dot,ydot,y2dot]

u0=[phi0,omega0,x0,vx0,y0,vy0]
res=solve_ivp(autism,(tb,te),u0,rtol=1e-10,method="RK45",t_eval=t)

print(len(res.t))

x2=res.y[2,:]+l2*np.cos(res.y[0,:])
y2=res.y[4,:]+l2*np.sin(res.y[0,:])
#plt.plot(res.y[2,:],res.y[4,:])
plt.title("c={}".format(c))
plt.plot(x2,y2)
plt.show()
