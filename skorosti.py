import numpy as np
import pandas
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data=pandas.read_csv('danni_magnitomahalo.csv')

def parabola(x,a,b,c):
    return a*x**2+b*x+c

def V(T,X,n):
    res=np.zeros(len(X))
    for i in range(n,len(X)-n+1):
        popt,pcov=curve_fit(parabola,T[i-n:i+n+1],X[i-n:i+n+1])
        res[i]=2*popt[0]*T[i]+popt[1]
    return res

plt.plot(data.t,V(data.t,data.xA,10))
plt.show()
