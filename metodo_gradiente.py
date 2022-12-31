"""
  Método do Gradiente

  O Método do Gradiente para aproximar soluções de sistemas de equações 
  não lineares. Por sua pequena taxa de convergência, o utilizamos para
  conseguir melhores pontos iniciais para aplicar o método de Newton ou 
  de Broyden.
"""
from auxiliares.baseMetodoNumerico import BaseMetodoNumerico

class Gradiente (BaseMetodoNumerico):
  """
    Método do Gradiente.
  
    Parâmetros
    ----------
    F : list
      Uma lista de funções das equações do sistema.
    metodo_diferenciacao : bool || function = False
      Uma função que será utilizada para calcular as derivadas parciais
      numericamente caso necessário.
    h : float = 0
      Tamanho do passo a ser considerado na diferenciação numérica, se
      for o caso.
  """
  def __init__ (self, F, metodo_diferenciacao=False, h=0):
    super().__init__(metodo_diferenciacao=metodo_diferenciacao, h=h)

    # adapta as funções
    self.F_lista = self.funcoes_f(F)
    self.F = self.funcao_F(self.F_lista)

    # cria a função g
    self.g = self.funcao_g()

  def funcao_g (self):
    """
      Cria a função `g(x) = f1(x)² + f1(x)² + ... + fn(x)²`.
    """
    return lambda x: sum( f(x)**2 for f in self.F_lista )

  def grad_g (self, Jac, F):
    """
      Calcula o gradiente analítico da função `g`.
    """
    return lambda x: 2 * Jac(x).T * F(x)

  def valores_h (self, lista_g: list, lista_alpha: list)->list:
    """
      Coeficientes do polinômio interpolador a partir das diferenças
      divididas de Newton.
    """
    h_1 = (lista_g[1] - lista_g[0])/(lista_alpha[1] - lista_alpha[0])
    h_2 = (lista_g[2] - lista_g[1])/(lista_alpha[2] - lista_alpha[1])
    h_3 = (h_2 - h_1)/(lista_alpha[2] - lista_alpha[0])

    return [h_1, h_2, h_3]

  def raiz_derivada_polinomio (self, h_1:float, h_3:float, alpha_2:float)->float:
    """
      Calcula a raiz da derivada do polinômio quadrático.
    """
    return (h_3 * alpha_2 - h_1)/(2 * h_3)

  def passo (self, p0, gradiente_g=False):
    """
      Aplica um passo do Método do Gradiente.
    """
    # aplica g em p0
    g_p0 = self.g(p0)

    # calcula ou aplica o gradiente
    if not gradiente_g:
      gradg_p0 = self.gradiente_local(self.g, p0)
    else:
      gradg_p0 = gradiente_g(p0)
    
    # calcula o versor relativo ao gradiente de g em p0
    z0 = self.norma_2(gradg_p0)
    z = gradg_p0/z0

    # o alpha1, por padrão, é 0
    alpha_1 = 0
    g_1 = g_p0

    # o alpha3 ainda preciso ver
    alpha_3 = 1
    g_3 = self.g(p0 - alpha_3 * z)

    g_3_anterior = g_3

    if g_3 >= g_1:
      alpha_3 = alpha_3/2
      g_3 = self.g(p0 - alpha_3 * z)

      coeficiente = 0.5 if g_3 < g_3_anterior else g_3 > g_3_anterior

      while g_3 >= g_1:
        alpha_3 *= coeficiente
        g_3 = self.g(p0 - alpha_3 * z)
    
    # alpha2 será a metade de alpha3
    alpha_2 = alpha_3/2
    g_2 = self.g(p0 - alpha_2*z)

    # cálculo dos coeficientes
    lista_h = self.valores_h([g_1,g_2,g_3],[alpha_1,alpha_2,alpha_3])

    # caso h_3 seja positiva, então o mínimo global será o ponto crítico
    if lista_h[2] > 0:
      # cálculo da raiz da derivada do polinômio
      alpha_0 = self.raiz_derivada_polinomio(lista_h[0], lista_h[2], alpha_2)

      # se a raiz estiver no intervalo [alpha_1, alpha_3], a aceitamos
      if alpha_1 <= alpha_0 <= alpha_3:
        alpha_hat = alpha_0

      # caso contrário, ficamos com alpha_3
      else: alpha_hat = alpha_3

    # caso h_3 < 0, não existe mínimo global e devemos ficar com alpha_3
    else:
      alpha_hat = alpha_3

    # calcula o novo ponto
    p1 = p0 - alpha_hat * z
    
    return p1

  def aplicar (self, p0, Jac=[], qntd_exata_passos:int=2):
    """
      Para facilitar a aplicação, pode-se utilizar esta função.

      Parâmetros
      ----------
      p0 : list || np.matrix
        Ponto inicial.
      Jac : list = []      
        Matriz Jacobiana. Se não for passada nenhuma, será calculada
        numericamente.
      qntd_exata_passos : int = -1
        Se passado um valor inteiro positivo, a função encerrá somente ao 
        atingir esta quantidade de passos.
    """

    # se o ponto inicial for do tipo lista, precisa converter
    if type(p0) == list: 
      if type(p0[0]) == list: p0 = self.matriz(p0)
      else: p0 = self.matriz([[p] for p in p0])

    # se passado algum gradiente, é possível montar a função grad(g) diretamente
    gradiente_g = False
    if len(Jac) > 0: 
      Jac = self.Jacobiana(Jac)
      gradiente_g = self.grad_g(Jac, self.F)

    # para salvar os pontos
    pontos = []    
    # aplica o método
    passo = 0
    while True:
      # aplica um passo do método do Gradiente
      p0 = self.passo(p0, gradiente_g)
      # salva o novo ponto
      pontos.append(p0)
      # adiciona ao passo
      passo += 1
      # se atingiu o número de passos, encerra
      if passo == qntd_exata_passos: break
    return pontos