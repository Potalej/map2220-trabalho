"""
  Método de Newton Multivariado

  O Método de Newton para funções de uma variável real pode ser usado para encontrar raízes
da função. Já para aplicações de maiores dimensões, havendo tantas quanto variáveis, é possível
aplicar o mesmo princípio do método para encontrar as soluções do sistema não-linear que 
descreve a aplicação.
"""
from python.baseMetodoNumerico import BaseMetodoNumerico
from python.diferenciacao import Diferenciacao
from numpy import array

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
  def __init__ (self, metodo_sistema_linear, metodo_diferenciacao_parcial, h):
    super().__init__(metodo_sistema_linear=metodo_sistema_linear)
    self.Dif = Diferenciacao(metodo_diferenciacao_parcial=metodo_diferenciacao_parcial, h=h)
    self.jacobiana = self.Dif.jacobiana
    pass

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
    if len(Jac) == 0: J = array(self.jacobiana(F, p0))
    else: J = array([[ Jij(*p0) for Jij in Ji] for Ji in Jac])

    # - F(x0)
    Fx0 = - array([ f(*p0) for f in F ])

    # J(x0)*y = -F(x0)
    y = self.resolver_sistema(J, Fx0)

    # aplica o passo
    p = p0 + y

    return p, y, J