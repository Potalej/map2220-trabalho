R = 10
R5 = 0.193
R6 = 4.10622*10**(-4)
R7 = 5.45177*10**(-4)
R8 = 4.4975*10**(-7)
R9 = 3.40735*10**(-5)
R10 = 9.615*10**(-7)

f1 = lambda x1, x2, x3, x4, x5: x1*x2+x1-3*x5
f2 = lambda x1, x2, x3, x4, x5: 2*x1*x2+x1+3*R10*x2**2+x2*x3**2+R7*x2*x3+R9*x2*x4+R8*x2-R*x5
f3 = lambda x1, x2, x3, x4, x5: 2*x2*x3**2+R7*x2*x3+2*R5*x3**2+R6*x3-8*x5
f4 = lambda x1, x2, x3, x4, x5: R9*x2*x4+2*x4**2-4*R*x5
f5 = lambda x1, x2, x3, x4, x5: x1*x2+x1+R10*x2**2+x2*x3**2+R7*x2*x3+R9*x2*x4+R8*x2+R5*x3**2+R6*x3+x4**2-1

Jac = [
    [
        lambda x1, x2, x3, x4, x5: x2+1,
        lambda x1, x2, x3, x4, x5: x1,
        lambda x1, x2, x3, x4, x5: 0,
        lambda x1, x2, x3, x4, x5: 0,
        lambda x1, x2, x3, x4, x5: -3,
    ],
    [
        lambda x1, x2, x3, x4, x5: 2*x2+1,
        lambda x1, x2, x3, x4, x5: 2*x1+6*R10*x2+x3**2+R7*x3+R9*x4+R8,
        lambda x1, x2, x3, x4, x5: 2*x2*x3+R7*x2,
        lambda x1, x2, x3, x4, x5: R9*x2,
        lambda x1, x2, x3, x4, x5: -R,
    ],
    [
        lambda x1, x2, x3, x4, x5: 0,
        lambda x1, x2, x3, x4, x5: 2*x3**2+R7*x3,
        lambda x1, x2, x3, x4, x5: 4*x2*x3+R7*x2+4*R5*x3+R6,
        lambda x1, x2, x3, x4, x5: 0,
        lambda x1, x2, x3, x4, x5: -8,
    ],
    [
        lambda x1, x2, x3, x4, x5: 0,
        lambda x1, x2, x3, x4, x5: R9*x4,
        lambda x1, x2, x3, x4, x5: 0,
        lambda x1, x2, x3, x4, x5: R9*x2+4*x4,
        lambda x1, x2, x3, x4, x5: -4*R,
    ],
    [
        lambda x1, x2, x3, x4, x5: x2+1,
        lambda x1, x2, x3, x4, x5: x1+2*R10*x2+R7*x3+R9*x4+R8,
        lambda x1, x2, x3, x4, x5: 2*x2*x3+R7*x2+2*R5*x3+R6,
        lambda x1, x2, x3, x4, x5: R9*x2+2*x4,
        lambda x1, x2, x3, x4, x5: 0,
    ],

]

F = [f1, f2, f3, f4, f5]
p0 = [10, 10, 10, 10, 10]

print('n \t erro \t k(J)')
i = 0
J_conds = []


from metodo_newton import Newton
from diferencasFinitas import DiferencasFinitas
from numpy.linalg import solve

dif = DiferencasFinitas(h=10e-5).derivada_definicao
newton = Newton(metodo_sistema_linear=solve, metodo_diferenciacao_parcial=dif, h=10e-5)

while True:
    p0, y, J_cond = newton.passo(F,p0)
    erro = newton.norma_infinito(y)
    print(i, '\t', erro)
    if erro < 10**(-8):
        break
    i += 1
print()
print("Solução: ", p0)
print()

sol = [0.003431, 31.325636, 0.068352, 0.859530, 0.036963]
difs = [abs(p0[i] - sol[i]) for i in range(len(p0))]
print('Diferenças: ', difs)
