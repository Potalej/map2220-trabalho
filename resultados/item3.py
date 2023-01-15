from math import sqrt
from numpy.random import ranf
from tabulate import tabulate

kappa = sqrt(2)/2

f1 = lambda x1, x2, x3, x4, x5, x6: x1**2 + x3**2 - 1
f2 = lambda x1, x2, x3, x4, x5, x6: x2**2 + x4**2 - 1
f3 = lambda x1, x2, x3, x4, x5, x6: x5*x3**3 + x6*x4**3
f4 = lambda x1, x2, x3, x4, x5, x6: x5*x1**3 + x6*x2**3
f5 = lambda x1, x2, x3, x4, x5, x6: x5*x1*x3**2 + x6*x2*x4**2
f6 = lambda x1, x2, x3, x4, x5, x6: x5*x1**2*x3+x6*x2**2*x4
F = [f1,f2,f3,f4,f5,f6]

def testar(ponto):
  residuo = []
  for f in F:
    residuo.append(round(f(*ponto), 10))
  return max(abs(i) for i in residuo)

def primeiroCaso (alpha:float):
  # para cada alpha há 4 combinações possíveis
  lista = [
    [kappa, kappa, kappa, kappa, -alpha, alpha],
    [kappa, kappa, -kappa, -kappa, -alpha, alpha],
    [-kappa, -kappa, kappa, kappa, -alpha, alpha],
    [-kappa, -kappa, -kappa, -kappa, -alpha, alpha]
  ]
  return lista

def segundoCaso (alpha:float):
  beta = sqrt(1 - alpha**2)
  lista = [
    [kappa, alpha, kappa, beta, 0, 0],
    [-kappa, alpha, kappa, beta, 0, 0],
    [-kappa, alpha, -kappa, beta, 0, 0],
    [kappa, alpha, -kappa, beta, 0, 0]
  ]
  return lista

def terceiroCaso (alpha:float):
  beta = sqrt(1 - alpha**2)
  lista = [
    [alpha, kappa, beta, kappa, 0, 0],
    [alpha, -kappa, beta, kappa, 0, 0],
    [alpha, -kappa, beta, -kappa, 0, 0],
    [alpha, kappa, beta, -kappa, 0, 0]
  ]
  return lista

def primeiroSubCaso (alpha:float):
  return [
    [1,1,0,0,alpha,-alpha],
    [-1,-1,0,0,alpha,-alpha],
    [1,-1,0,0,alpha,alpha],
    [-1,1,0,0,alpha,alpha],
  ]

def segundoSubCaso (alpha:float):
  return [
    [0,0,1,1,alpha,-alpha],
    [0,0,-1,-1,alpha,-alpha],
    [0,0,1,-1,alpha,alpha],
    [0,0,-1,1,alpha,alpha],
  ]

def terceiroSubCaso ():
  return [
    [1,0,0,1,0,0],
    [-1,0,0,-1,0,0],
    [1,0,0,-1,0,0],
    [-1,0,0,1,0,0],
  ]

def quartoSubCaso ():
  return [
    [0,1,1,0,0,0],
    [0,-1,-1,0,0,0],
    [0,1,-1,0,0,0],
    [0,-1,1,0,0,0],
  ]
  
def quintoSubCaso1 (L):
  delta = -4*(-1+L**4)*(1-L**2)
  x3 = sqrt(delta)/(2*(L**4 - 1))
  x4 = L * x3
  x1 = -sqrt(1-(L**2-1)/(L**4-1))
  x2 = x1/L
  x6 = 0
  x5 = 0
  return [[x1,x2,x3,x4,x5,x6]]

def quintoSubCaso2 (L):
  delta = -4*(1/L**2 - 1/L**4 - L**2 + 1)
  x2 = sqrt(delta)/(2*L**2 + 2/L**2)
  x3 = x2 / L
  x1 = sqrt(1 - delta/(L**2 * 4*(1/L**2 + L**2)**2))
  x4 = sqrt(1 - x2**2)
  x5 = 0
  x6 = 0
  return [[x1,x2,x3,x4,x5,x6]]

def quintoSubCaso3(L):
  delta = -4*(L**2-1/L**2)*(1-L**2)
  x3 = sqrt(delta)/(2*L**2 - 2/L**2)
  x2 = x3 / L
  x1 = sqrt(1 - delta/(2*L**2 - 2/L**2)**2)
  x4 = sqrt(1 - x2**2)
  return [[x1,x2,x3,x4,0,0]]

def quintoSubCaso4(L):
  delta = -4*(1-L**4)*(L**2-1)
  x4 = sqrt(delta)/(2-2*L**4)
  x3 = L * x4
  x2 = sqrt(1-x4**2)
  x1 = sqrt(1-x3**2)
  return [[x1,x2,x3,x4,0,0]]

sols1 = primeiroCaso # R
sols2 = segundoCaso # [-1, 1]
sols3 = terceiroCaso # [-1, 1]
sols41 = primeiroSubCaso # R
sols42 = segundoSubCaso # R
sols43 = terceiroSubCaso # R 
sols44 = quartoSubCaso # R

sols451 = quintoSubCaso1 # nor -1, 1
sols452 = quintoSubCaso1 # nor -1, 0, 1
sols453 = quintoSubCaso1 # nor -1, 0, 1
sols454 = quintoSubCaso1 # nor -1, 1

funcoes = [
  ['Primeiro caso',sols1, 0],
  ['Segundo caso', sols2, 'lim1'],
  ['Terceiro caso', sols3, 'lim1'],
  ['Primeiro sub-caso', sols41, 0],
  ['Segundo sub-caso', sols42, 0],
  ['Terceiro sub-caso', sols43, 'nao'],
  ['Quarto sub-caso', sols44, 'nao'],
  ['Quinto Sub-caso 1', sols451, [-1,1]],
  ['Quinto Sub-caso 2', sols451, [-1,0,1]],
  ['Quinto Sub-caso 3', sols451, [-1,0,1]],
  ['Quinto Sub-caso 4', sols451, [-1,1]]
]

def numerosAleatorios (a, b, qntd):
  lista = ranf(qntd)
  for i, x in enumerate(lista):
    lista[i] = x * (b - a) + a
  return lista

def gerarNumeros (qntd:int, intervalo:list=[], restricoes:list=[]):
  if len(intervalo) > 0:
    lista = numerosAleatorios(intervalo[0], intervalo[1],qntd)
  if len(restricoes) > 0:
    if len(intervalo) > 0:
      parar = False
      for rest in restricoes:
        if rest in lista:
          parar = True
          break
      if parar:
        lista = gerarNumeros (qntd, intervalo, restricoes)
    else:
      lista = numerosAleatorios(-100, 100,qntd)
      parar = False
      for rest in restricoes:
        if rest in lista:
          parar = True
          break
      if parar:
        lista = gerarNumeros (qntd, intervalo, restricoes)

  return lista

qntdPorFuncao = 15

megastring = "EXEMPLOS DE SOLUÇÕES\n\n"

for nome, funcao, limites in funcoes:
  intervalo = [-100, 100]
  restricoes = []
  parametros = True

  if limites == 'lim1':
    intervalo = [-1,1]

  elif limites == 'nao':
    parametros = False

  elif limites != 0:
    restricoes = limites

  if parametros:
    # gera numeros aleatorios para parâmetros
    lista = gerarNumeros(qntdPorFuncao, intervalo, restricoes)
    # aplica
    solucoes = []
    for param in lista:
      sol = funcao(param)
      for s in sol:
        solucoes.append([str(x) for x in s] + [str(testar(s))])
        # solucoes += "\n" + '\t'.join([str('%.3f' % x) for x in s]) + '\t' + str(testar(s))
    solucoes = tabulate(solucoes, headers=['x1','x2','x3','x4','x5','x6','res.'], tablefmt='presto')
  else:
    sol = funcao()
    solucoes = []
    for s in sol:
      solucoes.append([str(x) for x in s] + [str(testar(s))])
    solucoes = tabulate(solucoes, headers=['x1','x2','x3','x4','x5','x6','res.'], tablefmt='presto')
  
  megastring += F"# MÉTODO: {nome}\n\n"
  megastring += solucoes + '\n\n'

with open('item3_exemplos_solucoes.txt', 'w', encoding='utf8') as arq:
  arq.write(megastring)
