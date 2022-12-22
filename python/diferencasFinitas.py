"""
  Métodos de diferenças finitas

  Para fazer a diferenciação numérica, é frequente o uso do método de
diferenças finitas. Este pode ser calculado de várias formas, então a
classe aqui presente contém alguns deles.
"""

class DiferencasFinitas:
  """
    Métodos de diferenças finitas para diferenciação
    numérica.

    Parâmetros
    ----------
    h : float
      Tamanho do passo a ser considerado na diferenciação.
  """
  def __init__ (self, h:float):
    self.h = h

  def derivada_definicao (self, f, p, i:int)->float:
    """
      Calcula a derivada a partir da definição clássica, que com limites
      é exata mas que aqui assume a forma:
      
      `f'(x) = (f(x+h)-f(x))/h`

      Parâmetros
      ----------
      f : function
        Função a ser derivada.
      p : list | np.array
        Ponto em que a derivada será calculada.
      i : int
        Coordenada que será derivada.
      
      Retorna
      -------
      float
        Derivada parcial de `f` em `p` em relação à coordenada `i`.
    """
    ph = [*p] # faz uma cópia do ponto
    ph[i] += self.h # adiciona o passo h
    return (f(*ph)-f(*p))/self.h