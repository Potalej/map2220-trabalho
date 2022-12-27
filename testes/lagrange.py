def L(X, x, k):
    prod = 1
    for i, xi in enumerate(X):
        if i != k:
            prod *= (x-xi)/(X[k]-xi)
    return prod

def polinomioLagrange (f, X, x):
    soma = 0
    for i, xi in enumerate(X):
        soma += f(xi)*L(X, x, i)
    return soma

from math import sin, pi

f = lambda x: sin(x)

X = [0.5*pi, pi, 1.5*pi, 2*pi]
Y = [f(x) for x in X]

Psx = [1.5*pi*i/50 + pi/2 for i in range(51)]
Psy = [polinomioLagrange(f, X, x) for x in Psx]
Fy = [f(x) for x in Psx]

import matplotlib.pyplot as plt

plt.title("Interpolação de Lagrange para $f(x)=\sin{x}$ com \n quatro pontos de base")
plt.plot(Psx, Fy, c="black", label="$f(x)=sin(x)$")
plt.plot(Psx, Psy, c="black", label="$P_4(x)$", linestyle="dashed")
plt.scatter(X, Y, c="black")


plt.legend()
plt.show()