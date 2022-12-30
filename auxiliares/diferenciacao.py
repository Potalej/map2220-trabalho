"""
  MÉTODOS DE DIFERENCIAÇÃO NUMÉRICA

  Aqui consta uma classe para diferenciação numérica a partir de métodos informados.
  A classe conta com o cálculo da matriz jacobiana para aplicações multivariadas e 
  também o cálculo do gradiente.
"""
from auxiliares.vetores import Vetores

class Diferenciacao (Vetores):
  """
    Diferenciação numérica.

    Parâmetros
    ----------
    metodo_diferenciacao : function
      Uma função que será utilizada para calcular as derivadas parciais. Deve receber
      os parâmetros `f` (function), `p` (list || numpy.matrix) e `i` (int).
    h : float = 1e-5
      Tamanho do passo na diferenciação.
  """
  def __init__ (self, metodo_diferenciacao, h : float = 1e-5):
    super().__init__()
    self.metodo_diferenciacao = metodo_diferenciacao
    self.h = h 
  
  def derivada_parcial (self, f, p, i:int)->float:
    """
      Calcula a derivada parciala de uma função `f` no ponto `p` com respeito à
      coordenada `i`, assumindo que esta exista.
    """
    return self.metodo_diferenciacao(f, p, i)

  def gradiente_local (self, f, p):
    """
      Calcula o gradiente de uma função `f` em um dado ponto `p`.
    """
    return self.matriz([
      [self.derivada_parcial(f, p, i)] for i in range(len(p))
    ])

  def jacobiana_local (self, F:list, p):
    """
      Calcula a matriz Jacobiana de uma aplicação `F` no ponto `p`.
    """
    return self.matriz([
      [ self.derivada_parcial(f, p, i) for i in range(len(p))]
      for f in F
    ])

  def Jacobiana (self, J:list):
    """
      Converte uma matriz Jacobiana passada na forma de matriz de funções
      para uma função de matriz.

      Parâmetros
      ----------
      J : list
        Matriz de funções da matriz Jacobiana.

      Retorna
      -------
      function 
        Função matriz Jacobiana.
    """
    return lambda x: self.matriz([
      [ self.funcao(J[i][j], x) for j in range(len(J)) ] 
      for i in range(len(J))
    ])