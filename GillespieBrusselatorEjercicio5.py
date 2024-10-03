# Se importan las librerias necesarias
import numpy as np # Para operaciones con vectores
import math # Para funciones matemáticas
import random as rnd # Generación de números aleatorios
import matplotlib.pyplot as plt # Para graficar

# Condiciones iniciales
t0 = 0.0
X0=0.0
Y0=0.0

# Parámetros
alpha = 2.
beta = 1.
K_R = 1
n = 2
Omega_rnd = [rnd.randint(0, 1000) for _ in range(5)] # Lista de valores Omega
Steps = 5000

# Se define la Matriz estequiométrica vista en clase
S =[[1.,0.,-1.,0.],[0.,1.,0.,-1.]]

# Función que genera número aleatorio de distribución exponencial
def dist_exp(a):
    r = rnd.random()  # Se genera un número aleatorio de distribución uniforme
    return -(1./a) * math.log(r)  # Se transforma en distribución exponencial

# Elige aleatoriamente la reacción definida en el vector de propensiones ni y una tasa de reacción global A
def dist_reaction(ni,A):
	r = rnd.random() # Número aleatorio de distribucion uniforme
	if (r < ni[0]/A): # Se elige reaccion 0
		return 0
	elif (r < (ni[0] + ni[1])/A): #Se elige reacción 1
		return 1
	elif (r < (ni[0] + ni[1] + ni[2])/A): #Se elige reacción 2
		return 2
	elif (r <= (ni[0] + ni[1] + ni[2] + ni[3])/A): #Se elige reacción 3
		return 3

# Se inicializan los vectores solución y el tiempo para cada valor de X0_rnd
for Omega in Omega_rnd:
    # Se inicializan las condiciones iniciales
    Y = np.zeros([2, Steps + 1])
    t = np.zeros(Steps + 1)
    Y[0][0] = X0
    Y[1][0] = Y0
    

    # Ciclo del algoritmo
    for i in range(Steps):
        # Se actualiza vector de propensiones
        ni = np.array([alpha*(1/(1+(Y[1][i]/(K_R*Omega))**n)), alpha*(1/(1+(Y[0][i]/(K_R*Omega))**n)), beta*(Y[0][i]/Omega), beta*(Y[1][i]/Omega)])
        # Se actualiza tasa de reacción global
        a = sum(ni)
        # Se elige aleatoriamente tiempo de reacción
        tau = dist_exp(a)
        # Se elige la reacción que se llevará a cabo
        mu = dist_reaction(ni, a)
        # Se actualiza el estado del sistema y se avanza el tiempo
        Y[0][i + 1] = Y[0][i] + S[0][mu]
        t[i + 1] = t[i] + tau
        Y[1][i + 1] = Y[1][i] + S[1][mu]

    plt.plot(t, Y[1], label=f'r_2, Omega = {Omega}')
    plt.plot(t, Y[0], label=f'r_1, Omega = {Omega}')
	
# Se muestra la trafica
plt.xlabel('Tiempo')
plt.ylabel('Moléculas de la especie')
plt.title('Simulación estocástica para distintos Omega')
plt.legend()
plt.show()