import numpy as np
from scipy.integrate import tplquad,IntegrationWarning

mu_0=4*np.pi*10**-7

def B(r0,h,R,H,M):
    '''(r0,h) - координати на точката
        R - радиус на магнита
        H - височина на магнита
        M=Br/mu_0 - магнетизация на магнита, където Br е остатъчната индукция, 
        която може да се види на сайта на прамиз
        Всичко трябва да бъде в основни мерни единици!'''
    def rho(z,r,phi): # разстоянието до точката
        return np.sqrt((r0-r*np.cos(phi))**2+(r*np.sin(phi))**2+(H+h-z)**2)

    def dBx_r(z,r,phi): #диференциал на х-компонентата
        return (3*(H+h-z)*(r0-r*np.cos(phi))/rho(z,r,phi)**5)*r

    def dBy_r(z,r,phi): #диференциал на у-компонентата
        return (-3*(H+h-z)*r*np.sin(phi)/rho(z,r,phi)**5)*r

    def dBz_r(z,r,phi): #диференциал на z-компонентата
        return (3*(H+h-z)**2/rho(z,r,phi)**5-1/rho(z,r,phi)**3)*r
    try:
        Bx=mu_0*M/(4*np.pi)*tplquad(dBx_r,0,2*np.pi,lambda phi: 0,lambda phi: R,lambda phi,r: 0,lambda phi,r: H)[0]
    except IntegrationWarning: #тук не помня защо ми се е налагало да оправям тази грешка
        Bx=0
    try:
        By=mu_0*M/(4*np.pi)*tplquad(dBy_r,0,2*np.pi,lambda phi: 0,lambda phi: R,lambda phi,r: 0,lambda phi,r: H)[0]
    except IntegrationWarning:
        By=0
    Bz=mu_0*M/(4*np.pi)*tplquad(dBz_r,0,2*np.pi,lambda phi: 0,lambda phi: R,lambda phi,r: 0,lambda phi,r: H)[0]
    return np.array((Bx,By,Bz))

if __name__=="__main__":
    print(B(0.001,0.003,0.2,0.1,1.26/mu_0)) #тук стойностите са леко произволни