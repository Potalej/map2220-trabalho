"""
    CLASSE BÁSICA PARA APLICAÇÃO DE MÉTODOS NUMÉRICOS

    Esta classe contém métodos numéricos básicos de diferenciação
  que possibilita a fácil aplicação de métodos numéricos elaborados
  para a resolução de sistemas não lineares, como o método de Newton
  multivariado e os métodos Quasi-Newton, como o de Broyden

"""

class BaseMetodoNumerico:
  """
    Classe básica para aplicação de métodos numéricos.

    Parâmetros
    ----------
    metodo_sistema_linear : function
      Uma função que resolva sistemas lineares, informada a matriz de 
      coeficientes e a de vetores.
  """
  def __init__ (self, metodo_sistema_linear=False):
    if metodo_sistema_linear:
      self.metodo_sistema_linear = metodo_sistema_linear

  def resolver_sistema (self, A, b):
    """
      Para resolver sistemas lineares, independente de qual método
      se decidir utilizar no instanciamento.

      Parâmetros
      ----------
      A : numpy.array | numpy.matrix
        Matriz de coeficientes do sistema linear.
      b : numpy.array | numpy.matrix
        Vetor de resultados.
      
      Retorna
      -------
      numpy.array
        Vetor de solução.
    """
    return self.metodo_sistema_linear(A, b)

  def numero_condicao (self, M, norma)->float:
    """
      Calcula o número de condição de uma determinada matriz a partir
      de uma norma matricial informada.

      Parâmetros
      ----------
      M : list | np.array | np.matrix
        Matriz que se deseja calcular o número de condição.
      norma : function
        Norma matricial a ser utilizada.
    """
    return norma(M)

  def funcao (self, f, x):
    """
      Aplica uma função `f` sobre um vetor `x`.

      Parâmetros
      ----------
      f : function
        Função a ser aplicada.
      x : np.matrix
        Vetor no qual a função será aplicada.
    """
    return f(*[x[i,0] for i in range(len(x))])

  def norma_infinito (self, y):
    return max(abs(y[i,0]) for i in range(len(y)))

  def norma_2 (self, y):
    return sum(y[i,0]**2 for i in range(len(y)))**.5