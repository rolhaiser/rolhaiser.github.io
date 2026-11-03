#definimos el metodo de runge kutta de cuarto orden y un ejemplo
import numpy as np
import matplotlib.pyplot as plt
print("Inicio del programa :) \n")
def EDO(x,y):
    ydv = np.cos(x)
    return ydv
def sol_exa(x):
    return np.sin(x)
#condicion inicial
xo= 0
yo= 0
h=-0.1
#si queremos predecir valores anteriores a la condicion inicial
#podemos definir h negativo, y cuidar de definir X=np.arange(inicio, final<inicial, h)
def rk4(x,y,f=EDO,h=h):
    k1 = f(x,y)
    k2= f(x+h/2,y+h/2*k1)
    k3 = f(x+h/2,y+h/2*k2)
    k4 = f(x+h/2,y+h/2*k3)
    y_sig = y + (k1+2*k2+2*k3+k4)*h/6
    return y_sig

"""
X = np.arange(xo,10*np.sign(h),h) #la fn signo es por si h es negativo :>
Y= np.zeros(len(X))
E= sol_exa(X)

error_prom =0
for i in range(len(Y)):
    error_prom += abs(Y[i]-E[i])
error_prom = error_prom/len(Y)

plt.plot(X,E+1,c="b") #+1 para que no se tapen
Y[0] = yo
for i in range(len(Y)-1):
    Y[i+1]=rk4(X[i],Y[i])
plt.plot(X,Y,c="g")
plt.title(f"error promedio = {error_prom}")
plt.show()
"""

def EDOS(z,t,y):
    return np.array([[y],[2*np.e**t-2*y-z]])
t0 = 0
z0=0
zdv0 = 1
#y0 = np.array([[0],[1]]) #z(0), z'(0)
#V[0] = y0[1][0]
#print(rk4(0,y0))


T = np.arange(t0,10*np.sign(h),h) #la fn signo es por si h es negativo :>
P = np.zeros(len(T))
V = np.zeros(len(T))
P[0] = z0
V[0] = zdv0

for i in range(len(T)-1):
    P[i+1]=rk4(T[i],P[i],EDOS)
    V[i+1]=rk4(T[i],V[i],EDOS)
plt.plot(T,P)
plt.plot(T,V)
plt.show()

print("programa ejecutado :)")

#ejemplos para tratar
#   ydv = -2y,    xo = 0, yo=1,  y_exa= e**-2x
#   ydv = x  ,    xo = 0, yo=0,  y_exa=x**2/2
#   ydv = cos x,  xo = 0, yo=0,  y_exa= sen x
