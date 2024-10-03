# Se importan las librerias necesarias
import numpy as np # Para operaciones con vectores
import math # Para funciones matemáticas
import random as rnd # Generación de números aleatorios
import matplotlib.pyplot as plt # Para graficar

# Condiciones iniciales
t0 = 0.0
X0_rnd = [rnd.randint(10, 1000) for _ in range(10)] # Lista de valores iniciales
print(X0_rnd)

# Parámetros
alpha = 20.
beta = .5
Steps = 1050
s_s_v = alpha / beta  # Steady-state value

# Se define la Matriz estequiométrica vista en clase
S = [[1., -1.]]  # Solo una reacción reversible

# Función que genera número aleatorio de distribución exponencial
def dist_exp(a):
    r = rnd.random()  # Se genera un número aleatorio de distribución uniforme
    return -(1./a) * math.log(r)  # Se transforma en distribución exponencial

# Elige aleatoriamente la reacción definida en el vector de propensiones ni y una tasa de reacción global A
def dist_reaction(ni, A):
    r = rnd.random()  # Número aleatorio de distribución uniforme
    if r < ni[0] / A:  # Se elige reacción 0
        return 0
    elif r < (ni[0] + ni[1]) / A:  # Se elige reacción 1
        return 1

# Se inicializan los vectores solución y el tiempo para cada valor de X0_rnd
for X0 in X0_rnd:
    # Se inicializan las condiciones iniciales
    Y = np.zeros([1, Steps + 1])
    t = np.zeros(Steps + 1)
    F = np.zeros([1, Steps + 1])
    Y[0][0] = X0

    # Ciclo del algoritmo
    for i in range(Steps):
        # Se actualiza vector de propensiones
        ni = np.array([alpha, beta * Y[0][i]])
        # Se actualiza tasa de reacción global
        a = sum(ni)
        # Se elige aleatoriamente tiempo de reacción
        tau = dist_exp(a)
        # Se elige la reacción que se llevará a cabo
        mu = dist_reaction(ni, a)
        # Se actualiza el estado del sistema y se avanza el tiempo
        Y[0][i + 1] = Y[0][i] + S[0][mu]
        t[i + 1] = t[i] + tau
        F[0][i + 1] = Y[0][i + 1] - s_s_v

    # Se grafica la evolución en el tiempo
    plt.plot(t, Y[0] / s_s_v, label=f'X0 = {X0}')  # Normalizado y etiquetado por X0
    plt.xlabel('Tiempo')
    plt.ylabel('Número de moléculas normalizado')
    plt.title('Simulación estocástica para distintos números de moléculas X0')

# Se muestra la trafica
plt.legend()
plt.figure()
plt.hist(F[0], bins=100) #Para el último X0
plt.show()