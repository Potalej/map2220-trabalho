"""
  CLASSE PARA APLICAÇÃO DE MÉTODOS NUMÉRICOS

  Esta classe contém facilidades numéricas básicas para servir de base
  na criação de classes de métodos numéricos para a resolução de sistemas
  de equações não-lineares.
"""
from auxiliares.diferenciacao import Diferenciacao

class BaseMetodoNumerico(Diferenciacao):
  """
    Classe para aplicação de métodos de numéricos.

    Parâmetros
    ----------
    metodo_sistema_linear : bool || function = False
      Uma função que resolva sistemas lineares, informada a matriz de 
      coeficientes e a de vetores. Como nem sempre é necessária, por
      padrão está definida como False.
    metodo_diferenciacao : bool || function = False
      Uma função que será utilizada para calcular as derivadas parciais
      numericamente caso necessário.
    h : float = 0
      Tamanho do passo a ser considerado na diferenciação numérica, se for
      o caso.
  """
  def __init__ (self, metodo_sistema_linear=False, metodo_diferenciacao=False, h:float=0):
    if metodo_diferenciacao:
      super().__init__(metodo_diferenciacao, h)
    self.metodo_sistema_linear = metodo_sistema_linear

  def resolver_sistema_linear (self, A, b):
    """
      Para resolver sistemas lineares, independente de uqal método tenha
      se informado no instanciamento.

      Parâmetros
      ----------
      A : numpy.matrix | list
        Matriz de coeficientes do sistema linear.
      b : numpy.matrix | list
        Vetor de resultado.

      Retorna
      -------
      numpy.matrix
        Vetor de solução.
    """
    return self.metodo_sistema_linear(A, b)  