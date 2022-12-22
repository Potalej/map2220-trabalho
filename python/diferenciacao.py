"""
  Métodos numéricos para diferenciação numérica

  Aqui consta uma classe para diferenciação numérica a partir de métodos informados.
A classe também conta com o cálculo da matriz jacobiana para aplicações multivaridas.
"""

class Diferenciacao:
  """
    Diferenciação numérica via diferenças finitas.

    Parâmetros
    ----------
    metodo_diferenciacao_parcial : function
      Uma função que será utilizada para calcular as derivadas parciais. Deve receber os
      parâmetros `f` (function), `p` (list | np.array) e `i` (int).
    h : float = 10e-5
      Tamanho do passo na diferenciação.
  """
  def __init__ (self, metodo_diferenciacao_parcial, h : float = 10e-5):
    self.metodo_diferenciacao_parcial = metodo_diferenciacao_parcial
    self.h = h

  def gradiente (self, f, p):
    """
      Calcula o gradiente de uma função `f` em um dado ponto `p`.

      Parâmetros
      ----------
      f : function
        Função a ser derivada.
      p : list | np.array 
        Vetor de posição na qual a função será derivada.

      Retorna
      -------
      list
        Vetor gradiente calculado em `p`.
    """
    grad = [
      self.derivada_parcial(f, p, i)
      for i in range(len(p))
    ]
    return grad

  def derivada_parcial (self, f, p, i:int)->float:
    """
      Calcula a derivada parcial de uma função `f` no ponto `p` com respeito à coordenada `i`.

      Parâmetros
      ----------
      f : function
        Função a ser derivada.
      p : list | np.array
        Vetor de posição na qual a função será derivada.
      i : int
        Índice da coordenada a ser derivada.
    """
    return self.metodo_diferenciacao_parcial(f, p, i)

  def jacobiana (self, F:list, p)->list:
    """
      Calcula a matriz jacobiana de uma aplicação `F` no ponto `p`.

      Parâmetros
      ----------
      F : list
        Aplicação a ser diferenciada, isto é, uma lista de funções reais.
      p : list | np.array
        Vetor de posição na qual a função será derivada.

      Retorna
      -------
      list
        Matriz jacobiana de `F` em `p`.
    """
    Jac = [
      self.gradiente(f, p)
      for f in F
    ]
    return Jac
