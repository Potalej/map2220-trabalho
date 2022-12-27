"""
    Testes de Diferenciação Numérica

    Aqui testaremos a diferenciação numérica via diferenças progressivas
    e via diferenças centradas. A ideia é diferenciar a função f(x)=sen(x),
    cuja derivada é conhecida e é a função f'(x)=cos(x), e conferir a margem
    de erro da derivação conforme o tamanho do passo.
"""
from tabulate import tabulate
from ..diferencasFinitas import DiferencasFinitas
from math import sin, cos

f = lambda x: sin(x)
f_ = lambda x: cos(x)
f__ = lambda x: -sin(x)
f___ = lambda x: -cos(x)

x0 = 1/2
valor_real = f_(x0)

epsilon_estrela = 7e-17

"""
    Primeiramente, vamos encontrar a derivada usando valores de h cada vez
    menores só para mostrar que funciona.
"""
hs = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]

difs_progressivas = []
difs_centradas = []

for h in hs:
    dif = DiferencasFinitas(h)

    dif_prog = abs(valor_real - dif.progressivas(f, [x0], 0))
    dif_cent = abs(valor_real - dif.centradas(f, [x0], 0))

    difs_progressivas.append(dif_prog)
    difs_centradas.append(dif_cent)

tabela = [hs, difs_progressivas, difs_centradas]
tabela = list(zip(*tabela))

print("Diferença entre valor real e calculado para alguns valores de h")
print(tabulate(tabela, headers=["h", "progr.", "centr."], tablefmt="presto"))

"""
    Também é preciso estimar o erro. Nesse caso, para diferenças progressivas
    o erro será -hf''(C)/2, onde C habita entre x0 e x0+h, e para as diferenças
    centradas o erro será -h²f'''(D)/6, onde D habita entre x0 e x0+2h.

    Os valores C e D devem ser encontrados de modo a maximizar as funções |f''| e 
    |f'''| no intervalo. Observe que f''(x)=-sin(x) e f'''(x)=-cos(x). Podemos
    ver via gráfico que enquanto |f''| aumenta conforme x aumenta no intervalo 
    [x0, x0+h], |f'''| diminui. Então devemos escolher x0 para maximizar |f'''|
    e x0+h para maximizar |f''|.
"""
csi_maximo_prog = lambda h: x0+h
erro_estimado_prog = lambda h: 0.5*h*abs(f__(csi_maximo_prog(h))) + 2*epsilon_estrela/h * abs(f(x0 + h))

csi_maximo_cent = lambda h: x0
erro_estimado_cent = lambda h: (1/6)*h*h*abs(f___(csi_maximo_cent(h))) + epsilon_estrela/h * abs(f(x0 + h))

estimado_difs_progressivas = []
estimado_difs_centradas = []

for h in hs:
    dif = DiferencasFinitas(h)

    dif_prog = erro_estimado_prog(h)
    dif_cent = erro_estimado_cent(h)

    estimado_difs_progressivas.append(dif_prog)
    estimado_difs_centradas.append(dif_cent)

tabela = [hs, difs_progressivas, estimado_difs_progressivas, difs_centradas, estimado_difs_centradas]
tabela = list(zip(*tabela))

print("\nErro máximo estimado para alguns valores de h")
print(tabulate(tabela, headers=["h", "prog real", "prog estim", "cent real", "cent estim"], tablefmt="presto"))

"""
    Agora, é importante observar a questão do tamanho ótimo do passo. Conseguimos
    estimar para as progressivas e para as centradas.
"""

h_prog = dif.h_otimo_progressivas_1(epsilon_estrela, f, f__, x0)
h_cent = dif.h_otimo_centradas_1(epsilon_estrela, f, f___, x0)

erro_estimado_h_otimo_prog = erro_estimado_prog(h_prog)
erro_estimado_h_otimo_cent = erro_estimado_cent(h_cent)

dif.h = h_prog
dif_progr = abs(valor_real - dif.progressivas(f, [x0], 0))
dif.h = h_cent
dif_centr = abs(valor_real - dif.centradas(f, [x0], 0))

print("\nErro estimado para um valor ótimo de h")
print(f"h* = {h_prog} | estimado = {erro_estimado_h_otimo_prog} | erro = {dif_progr}")
print(f"h* = {h_cent} | estimado = {erro_estimado_h_otimo_cent} | erro = {dif_centr}")