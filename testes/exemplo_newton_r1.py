"""
    Exemplo do método de Newton para funções univariadas

    Aqui consta um exemplo de aplicação do Método de Newton para
    encontrar a raiz de uma função f:R->R não linear, que no caso
    é f(x) = cos(x)-x³. Uma interpretação equivalente é resolver
    a equação não linear cos(x)=x³.

    O ponto inicial considerado será p0=1/2.
"""
from tabulate import tabulate
from numpy import float64, cos, sin

f = lambda x: cos(x)-x**3       #  f (x)
f_ = lambda x: -sin(x)-3*x**2   #  f'(x)
f__ = lambda x: -cos(x)-6*x     # f''(x)

p0 = float64(3/4) # ponto inicial

p = float64(0.86547403310161444662) # solução exata

# método de Newton
def newton (p0, f, f_p0): return p0 - f(p0)/f_p0

"""
    É importante também medirmos o erro real e o erro estimado para
    verificar a teoria de convergência quadrática.

    Primeiramente, podemos ver que a condição M |e0| < 1. Para isso,
    considere que estamos pensando no intervalo (.5, 1) e que então,
    já que tomamos p0 = 3/4, no pior dos casos e_0 = 1/4. Assim:
"""
e_0 = float64(1/4)

M = abs(f__(float64(1)))/(2*abs(f_(float64(.5))))

print(f'\nGarantindo a 3ª condição:\nM|e_0| = {M*abs(e_0)} < 1\n')

"""
    Agora podemos ir aos erros
"""
e_n = e_0
erro_absoluto = []
erro_estimado = []

for i in range(8):
    der = f_(p0)
    p0 = newton(p0, f, der)
    
    en1 = abs(p-p0)
    
    er_est = M*(e_n**2)

    erro_absoluto.append(en1)
    erro_estimado.append(er_est)

    e_n = er_est

tabela = [
    [i for i in range(len(erro_absoluto))],
    erro_absoluto,
    erro_estimado
]
tabela = list(zip(*tabela))
print(tabulate(tabela, headers=['n', 'absoluto', 'estimativa']))