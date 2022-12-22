from broyden import Broyden
import numpy as np

R = np.float64(10)
R5 = np.float64(0.193)
R6 = np.float64(4.10622*10**(-4))
R7 = np.float64(5.45177*10**(-4))
R8 = np.float64(4.4975*10**(-7))
R9 = np.float64(3.40735*10**(-5))
R10 = np.float64(9.615*10**(-7))

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

import numpy as np

p0 = [np.float64(10), np.float64(10), np.float64(10), np.float64(10), np.float64(10)]

print('n \t erro \t k(J)')
i = 0
J_conds = []

broyden = Broyden()

# while True:
erro_aceitavel = 10e-2
p0, erro = broyden.passo_broyden(F,p0, erro=erro_aceitavel)
print()
print("Solução: ", p0)
print()

sol = [0.003431, 31.325636, 0.068352, 0.859530, 0.036963]
difs = [abs(p0[i] - sol[i]) for i in range(len(p0))]
print('Diferenças: ', difs)

print()
print('Erro:', erro)