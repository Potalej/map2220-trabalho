"""
  CLASSE BÁSICA PARA APLICAÇÃO DE FUNÇÕES VETORIAIS

  Criar funções vetoriais de maneira generalizada pode ser complicado
  no Python, então esta classe ajuda a aplicar funções em vetores 
  matriciais, ou seja, que possuem dois índices. 
"""
from numpy import matrix

class Vetores:
  """
    Classe básica para aplicação de funções vetoriais
  """
  def matriz (self, M:list):
    """
      Recebe uma matriz na forma de lista e retorna na forma de 
      numpy.matrix.
    """
    return matrix(M)
  
  def funcao (self, f, x):
    """
      Aplica uma função `f` sobre um vetor `x` matricial.

      Parâmetros
      ----------
      f : function
        Função a ser aplicada.
      x : np.matrix | list (n x 2) | dict ((i,2), i=1,..,n)
        Vetor matricial.

      Retorna
      -------
      function
        Função que agora pode receber um vetor matricial.
    """
    return f(*[x[i,0] for i in range(len(x))])

  def funcoes_f (self, F:list):
    """
      A partir de uma lista de funções multivariadas, cria uma lista de funções
      que recebam facilmente vetores matriciais.
    """
    lista = []
    for f in F:
      lista.append(lambda x, fi=f: self.funcao(fi,x))
    return lista

  def funcao_F (self, F:list):
    """
      A partir de uma lista de funções que recebem vetores matriciais monta um
      campo vetorial.
    """
    return lambda x: self.matriz([
      [ (lambda x, fi=f: fi(x))(x) ] for f in F
    ])

  def norma_infinito (self, y):
    return max(abs(y[i,0]) for i in range(len(y)))

  def norma_2 (self, y):
    return sum(y[i,0]**2 for i in range(len(y)))**.5