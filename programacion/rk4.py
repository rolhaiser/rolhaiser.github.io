#!/usr/bin/env python3
#definimos el metodo Runge-Kutta de 4to orden con ejemplos
import numpy as np #numerical python
import matplotlib.pyplot as plt #graficar
import datetime #medir tiempo de ejecucion
#pagina con todos los tipos de graficos con codigos :D
#    https://matplotlib.org/stable/gallery/index.html

inicio = datetime.datetime.now()
print(f"\nInicio del programa a las {inicio.strftime("%H:%M")} :) \n")
#-----
def rk4(x,y,f,h=0.01): #actualiza el valor una vez
    k1 = f(x,y)
    k2= f(x+h/2,y+h/2*k1)
    k3 = f(x+h/2,y+h/2*k2)
    k4 = f(x+h/2,y+h*k3)
    y_sig = y + (k1+2*k2+2*k3+k4)*h/6
    return y_sig

# Ejemplo 1): y'=f(x,y), con y(x0)=y0
#"""
#   ejemplos para tratar
#       ydv = -2y,    xo = 0, yo=1,  y_exa= e**-2x
#       ydv = x  ,    xo = 0, yo=0,  y_exa=x**2/2
#       ydv = cos x,  xo = 0, yo=0,  y_exa= sen x

xo = 0; yo=0
h = 0.01
def ydv1(x,y): return np.cos(x)
def yexa1(x): return np.sin(x)
#   vamos a graficar la solucion en el siguiente intervalo
xmin = -4; xmax = 4
X_sig = np.arange(xo,xmax+h,h)
Y_sig = np.zeros(len(X_sig)); Y_sig[0] = yo
for i in range(len(Y_sig)-1):
    Y_sig[i+1]=rk4(X_sig[i],Y_sig[i],ydv1)

if xmin < xo:
    X_ant = np.arange(xo,xmin-h,-h)
    Y_ant = np.zeros(len(X_ant)); Y_ant[0] = yo
    for i in range(len(Y_ant)-1):
        Y_ant[i+1] = rk4(X_ant[i], Y_ant[i], ydv1, h=-h)
    X_sig = np.concatenate((X_ant[::-1], X_sig[1:])) #invertimos el orden y unimos
    Y_sig = np.concatenate((Y_ant[::-1], Y_sig[1:]))

plt.plot(X_sig,yexa1(X_sig),label="exacta")
plt.plot(X_sig,Y_sig,label="rk4", color="green")
plt.scatter(xo,yo,color="red") #punto condicion inicial

#   configuracion extra de los graficos
plt.title("Ejemplo 1 RK4: y'=cos(x)")
plt.xlabel("x"); plt.ylabel("y")
plt.xlim(xmin, xmax)
plt.axhline(0, color='black', linewidth=1, linestyle='-')  # línea en y=0
plt.axvline(0, color='black', linewidth=1, linestyle='-')  # línea en y=0
plt.grid(True, which='both', linestyle='--', linewidth=0.7, color='gray', alpha=0.7)
plt.legend(loc="upper right"); plt.show()
#"""

# Ejemplo 2): y''= f(x,y,y')
#"""
#   ejemplos para tratar
#       y''+2y'+y = 2e^x, con y(0)=0, y'(0) = 1
#   como es de orden 2, definimos 2-1 variables auxiliares z = y'
#   z'=2e^x-2z-y,   vecY=[y,z],    vecY'=[y',z']= f(x,vecY)
#No hace falta que el sistema sea lineal, basta con que la
#derivada de mayor orden se pueda despejar
xo = 0; yo = 0; dvyo = 1
vecYo= np.array([yo,dvyo]) #vector condicion inicial
h=0.01
xmin = -4; xmax= 4

def ydv2(x,vecY):
    y = vecY[0]; z= vecY[1]
    vecYdv = np.array([z,2*np.e**x-2*z-y]) #vector derivada
    return vecYdv
def yexa2(x): return np.sinh(x)

X_sig = np.arange(xo,xmax+h,h)
vecY_sig = np.zeros([2,len(X_sig)]); vecY_sig[:,0] = vecYo
for i in range(len(X_sig)-1):
    vecY_sig[:,i+1] = rk4(X_sig[i], vecY_sig[:,i], ydv2, h=h)

if xmin < xo:
    X_ant = np.arange(xo,xmin-h,-h)
    vecY_ant = np.zeros([2,len(X_ant)]); vecY_ant[:,0] = vecYo
    for i in range(len(X_ant)-1):
        vecY_ant[:,i+1] = rk4(X_ant[i], vecY_ant[:,i], ydv2, h=-h)
    X_sig = np.concatenate((X_ant[::-1], X_sig[1:])) #invertimos el orden y unimos
    vecY_sig = np.hstack((vecY_ant[:, ::-1], vecY_sig[:, 1:])) #une columnas


posicion= vecY_sig[0]; velocidad= vecY_sig[1]
plt.plot(X_sig,yexa2(X_sig),label="exacto")
plt.plot(X_sig,posicion,label="rk4",color="green")
plt.scatter(xo,vecYo[0],color="red") #punto condicion inicial

plt.title("Ejemplo 2 RK4: $y''+2y'+y = 2e^x$")
plt.xlabel("x"); plt.ylabel("y")
plt.xlim(xmin, xmax)
plt.axhline(0, color='black', linewidth=1, linestyle='-')  # línea en y=0
plt.axvline(0, color='black', linewidth=1, linestyle='-')  # línea en y=0
plt.grid(True, which='both', linestyle='--', linewidth=0.7, color='gray', alpha=0.7)
plt.legend(loc="lower right"); plt.show()
#"""

# Ejemplo 3): pendulo
#"""
#   y''= -g/Lsen(y)    con y(0)=45°;    y'(0)=0
xo = 0; yo = (179)*np.pi/180; dvyo = 0
g = 9.8; L = 1
h = 0.01
vecYo= np.array([yo,dvyo]) #vector condicion inicial
xmin = -10; xmax= 10
def ydv3(x,vecY):
    y = vecY[0]; z= vecY[1]
    vecYdv = np.array([z,-g/L*np.sin(y)]) #vector derivada
    return vecYdv

X_sig = np.arange(xo,xmax+h,h)
vecY_sig = np.zeros([2,len(X_sig)]); vecY_sig[:,0] = vecYo
for i in range(len(X_sig)-1):
    vecY_sig[:,i+1] = rk4(X_sig[i], vecY_sig[:,i], ydv3, h=h)

if xmin < xo:
    X_ant = np.arange(xo,xmin-h,-h)
    vecY_ant = np.zeros([2,len(X_ant)]); vecY_ant[:,0] = vecYo
    for i in range(len(X_ant)-1):
        vecY_ant[:,i+1] = rk4(X_ant[i], vecY_ant[:,i], ydv3, h=-h)
    X_sig = np.concatenate((X_ant[::-1], X_sig[1:])) #invertimos el orden y unimos
    vecY_sig = np.hstack((vecY_ant[:, ::-1], vecY_sig[:, 1:])) #une columnas

posicion= vecY_sig[0]; velocidad= vecY_sig[1]
plt.plot(X_sig,posicion,label="rk4",color="green")
plt.scatter(xo,vecYo[0],color="red") #punto condicion inicial

plt.title(rf"Ejemplo 3 RK4: péndulo $\theta_0 = {int(yo*180/np.pi)}°$")
plt.xlabel("t [s]"); plt.ylabel(r"$\theta(t)\:\:[rad]$")
plt.xlim(xmin, xmax)
plt.axhline(0, color='black', linewidth=1, linestyle='-')  # línea en y=0
plt.axvline(0, color='black', linewidth=1, linestyle='-')  # línea en y=0
plt.grid(True, which='both', linestyle='--', linewidth=0.7, color='gray', alpha=0.7)
plt.legend(loc="lower right"); plt.show()
#"""
#-----
fin = datetime.datetime.now()
print(f"programa ejecutado en {fin-inicio} \n")