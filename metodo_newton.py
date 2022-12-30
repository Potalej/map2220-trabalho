"""
  Método de Newton Multivariado

  O Método de Newton para funções de uma variável real pode ser usado para encontrar raízes
da função. Já para aplicações de maiores dimensões, havendo tantas quanto variáveis, é possível
aplicar o mesmo princípio do método para encontrar as soluções do sistema não-linear que 
descreve a aplicação.
"""
from baseMetodoNumerico import BaseMetodoNumerico
from diferenciacao import Diferenciacao
from time import time
from numpy import matrix

class Newton (BaseMetodoNumerico):
  """
    Método de Newton Multivariado.

    Parâmetros
    ----------
    metodo_sistema_linear : function
      Uma função que resolva sistemas lineares, informada a matriz de
      coeficientes e a de vetores.
    metodo_diferenciacao_parcial : function
      Uma função que será utilizada para calcular as derivadas parciais. Deve receber os
      parâmetros `f` (function), `p` (list | np.array) e `i` (int).
    h : float
      Tamanho do passo a ser considerado na diferenciação.
  """
  def __init__ (self, metodo_sistema_linear, metodo_diferenciacao_parcial=False, h=0):
    super().__init__(metodo_sistema_linear=metodo_sistema_linear)
    if metodo_diferenciacao_parcial:
      self.Dif = Diferenciacao(metodo_diferenciacao_parcial=metodo_diferenciacao_parcial, h=h)
      self.jacobiana = self.Dif.jacobiana
    pass

  def passo_tempo (self, F, p0, Jac=[]):
    """
      Aplica um passo do Método de Newton Multivariado da exata mesma forma que a função
      `passo`. Porém, retorna também o tempo demandado em cada etapa. Pode ser mais lento 
      que a função `passo`.

      Parâmetros
      ----------
      F : list
        Aplicação a ser desenvolvida no método.
      p0 : np.array
        Posição inicial.
      Jac : list = []
        Matriz jacobiana. Se não for passada nenhuma, será 
        calculada numericamente.

      Retorna
      -------
      (np.array, np.array, np.array, dict)
        Respectivamente, novo ponto `p`, vetor `y` (resultado 
        de `Jy - F`), `J` e informação de tempo.
    """
    # dicionário para armazenar os tempos
    tempo = dict()

    t0 = time()
    # matriz jacobiana
    if len(Jac) == 0: J = matrix(self.jacobiana(F, p0))
    else: J = matrix([[ self.funcao(Jij, p0) for Jij in Ji] for Ji in Jac])
    tempo["jacobiana"] = time() - t0

    # - F(x0)
    t0 = time()
    Fx0 = - matrix([ [self.funcao(f, p0)] for f in F ])
    tempo["Fx0"] = time() - t0

    # J(x0)*y = -F(x0)
    t0 = time()
    y = self.resolver_sistema(J, Fx0)
    tempo["y"] = time() - t0

    # aplica o passo
    t0 = time()
    p = p0 + y
    tempo["p"] = time() - t0

    # o tempo total será a soma de todos os tempos
    tempo["total"] = sum(tempo[i] for i in tempo)
    
    return p, y, J, tempo

  def passo (self, F, p0, Jac=[]):
    """
      Aplica um passo do Método de Newton Multivariado.

      Parâmetros
      ----------
      F : list
        Aplicação a ser desenvolvida no método.
      p0 : np.array
        Posição inicial.
      Jac : list = []
        Matriz jacobiana. Se não for passada nenhuma, será 
        calculada numericamente.

      Retorna
      -------
      (np.array, np.array, np.array)
        Respectivamente, novo ponto `p`, vetor `y` (resultado 
        de `Jy - F`) e `J`.
    """
    # matriz jacobiana
    if len(Jac) == 0: J = matrix(self.jacobiana(F, p0))
    else: J = matrix([[ self.funcao(Jij, p0) for Jij in Ji] for Ji in Jac])

    # - F(x0)
    Fx0 = - matrix([ [self.funcao(f, p0)] for f in F ])

    # J(x0)*y = -F(x0)
    y = self.resolver_sistema(J, Fx0)

    # aplica o passo
    p = p0 + y
    
    return p, y, J, ()

  def aplicar (self, F, p0, Jac=[], erro_admitido:float=10e-5, qntd_maxima_passos:int=10e2, solucao_exata=[], qntd_exata_passos:int=-1, medir_tempo:bool=False)->tuple:
    """
      Para facilitar a aplicação, pode-se utilizar esta função.

      Parâmetros
      ----------
      F : list
        Aplicação a ser desenvolvida no método.
      p0 : np.array
        Posição inicial.
      Jac : list | np.array | np.matrix = []
        Matriz jacobiana. Se não for passada nenhuma, será 
        calculada numericamente.
      erro_admitido : float
        É o erro máximo admitido. Assim que o erro se tornar estritamente 
        menor a função encerrará.
      qntd_maxima_passos : int
        Caso o método não venha a convergir para a solução, ao atingir esta
        quantidade de passos a função encerrará.
      solucao_exata : list = []
        Caso se possua a solução exata, é possível informar que seja calculado
        o erro real do método.
      qntd_exata_passos : int = -1
        Se passado valor inteiro positivo, a função encerrá somente ao atingir 
        esta quantidade exata de passos.
      medir_tempo : bool = False
        Caso se deseje medir o tempo necessário para calcular cada passo, deve-se
        passar `True` e as informações de tempo serão retornadas no dicionário de 
        informações.
    """
    # se o ponto inicial for do tipo list, precisa converter
    if type(p0) == list: p0 = matrix(p0)

    # dicionário de informações
    info = {
      "erros": [], 
      "erro_real": [],
      "passo": 0
    }
    
    # caso se deseje armazenar informações de tempo
    if medir_tempo: 
      info["tempo"] = dict()
      metodo = self.passo_tempo
    else:
      metodo = self.passo


    while True:
      # aplica o método
      p0, y, J_cond, tempo = metodo(F, p0, Jac)
      if medir_tempo:
        for tipo in tempo: 
          try:
            info["tempo"][tipo] += [tempo[tipo]]
          except:
            info["tempo"][tipo] = [tempo[tipo]]
      
      # medição do erro
      erro = self.norma_infinito(y)
      info["erros"].append(erro)

      # caso possua solução exata informada, a calcula
      if len(solucao_exata) > 0: 
        info["erro_real"].append(self.norma_infinito(solucao_exata - p0))

      # adiciona à quantidade de passos
      info["passo"] += 1

      # verifica se há quantidade exata de passos para parar 
      if qntd_exata_passos > 0:
        # se tiver, verifica se já bateu
        if info["passo"] == qntd_exata_passos:
          break 
      else:
        # verifica se o passo ultrapassou o limite ou o erro ficou abaixo do admitido
        if info["passo"] >= qntd_maxima_passos or erro < erro_admitido:
          break

    return p0, info