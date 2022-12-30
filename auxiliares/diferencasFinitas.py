"""
  Métodos de diferenças finitas

  Para fazer a diferenciação numérica, é frequente o uso do método de
diferenças finitas. Este pode ser calculado de várias formas, então a
classe aqui presente contém alguns deles.
"""

from math import sqrt

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

  def progressivas (self, f, p, i:int)->float:
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
    ph = p.copy() # faz uma cópia do ponto
    ph[i,0] += self.h # adiciona o passo h
    return (f(ph)-f(p))/self.h
  
  def centradas (self, f, p, i:int)->float:
    """
      Calcula a derivada a partir da expansão com três pontos considerando
      o ponto do meio.

      `f'(x) = [f(x+h)-f(x-h)]/2h`

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
    ph_pos = p.copy() # faz uma cópia do ponto
    ph_neg = p.copy() # faz uma cópia do ponto
    ph_pos[i,0] += self.h # adiciona o passo h
    ph_neg[i,0] -= self.h # adiciona o passo h
    return (f(ph_pos)-f(ph_neg))/(2*self.h)

  def h_otimo_progressivas_1 (self, epsilon_estrela, f, f__, x0):
    """
      Calcula o valor ótimo para o passo `h`.

    """
    return 2*sqrt(epsilon_estrela*abs(f(x0)))/sqrt(abs(f__(x0)))

  def h_otimo_centradas_1 (self, epsilon_estrela, f, f___, x0):
    cbrt = lambda x: x**(1/3)
    return cbrt(3*epsilon_estrela*abs(f(x0)))/cbrt(abs(f___(x0)))