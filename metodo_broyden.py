""""
  Método de Broyden

  O método de Broyden resolve sistemas de equações não lineares, sendo uma generalização 
  do método da Secante. Sua vantagem sobre o Método de Newton é não precisar calcular as
  derivadas em cada passo, mas somente uma vez. Sua desvantagem é sua convergência apenas
  superlinear.
"""
from auxiliares.baseMetodoNumerico import BaseMetodoNumerico
from numpy.linalg import inv
from time import perf_counter as time

class Broyden (BaseMetodoNumerico):
  """
    Método de Broyden.

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
  def __init__ (self, F:list, metodo_diferenciacao=False, h=0):
    if metodo_diferenciacao: super().__init__(metodo_diferenciacao=metodo_diferenciacao, h=h)
    else: super().__init__()
    
    # adapta as funções
    self.F_lista = self.funcoes_f(F)
    self.F = self.funcao_F(self.F_lista)

  def novoA_inv (self, p0, p1, Fp, A0_inv):
    """
      Calcula a nova matriz A^-1.

      Parâmetros
      ----------
      p0 : numpy.matrix
        Ponto anterior.
      p1 : numpy.matrix
        Ponto atual.
      A0_inv : numpy.matrix
        Matriz A inversa atual.
    """
    y1 = self.F(p1) - Fp
    s1 = p1 - p0
    s1t = s1.T
    A1_inv = A0_inv + (s1 - A0_inv * y1)*s1t*A0_inv / (s1t * A0_inv * y1)
    return A1_inv

  def passo_tempo (self, p, A_inv):
    """
      Aplica um passo do Método de Broyden, armazenando
      o tempo decorrido para calcular cada parte.
    """
    tempo = dict()

    # calcula o novo p
    t0 = time()
    Fp = self.F(p)
    p1 = p - A_inv * Fp
    tempo["p1"] = time() - t0

    return p1, tempo, Fp

  def passo (self, p, A_inv):
    """
      Aplica um passo do Método de Broyden.
    """
    # calcula o novo p
    p1 = p - A_inv * self.F(p)
    return p1, 0

  def aplicar (self, p0, Jac=[], erro_admitido:float=1e-5, qntd_maxima_passos=1e2, solucao_exata=[], qntd_exata_passos:int=-1, medir_tempo:bool=False, limitacao_float:bool=False, exibir_causa_fim:bool=True)->tuple:
    """
      Para facilitar a aplicação, pode-se utilizar esta função.

      Parâmetros
      ----------
      p0 : numpy.matrix
        Ponto inicial.
      Jac : list = []
        Matriz Jacobiana. Se não for passada nenhuma, será calculada
        numericamente.
      erro_admitido : float
        É o erro máximo admitido. Assim que o erro se tornar estritamente 
        menor a função encerrá.
      qntd_maxima_passos : int
        Caso o método não venha a convergir para a solução, ao atingir
        esta quantidade de passos a função encerrará.
      solucao_exata : list = []
        Caso se possua a solução exata, é possível informá-la para que
        seja calculado o erro real do método.
      qntd_exata_passos : int = -1
        Se passado um valor inteiro positivo, a função encerrá somente ao 
        atingir esta quantidade de passos.
      medir_tempo : bool = False
        Caso se deseje medir o tempo necessário para calcular cada passo,
        deve-se passar `True` e as informações de tempo serão retornadas
        no dicionário de informações.
      limitacao_float : bool = False
        Encerra o método somente quando atingir o número máximo de passos
        informado ou quando a norma infinito das diferenças de passo para
        passo atingir seu mínimo fixo.
      exibir_causa_fim : bool = True
        Exibirá a causa do encerramento do método quando `True`.
    """
    # se o ponto inicial for do tipo lista, precisa converter
    p0 = self.ponto_matriz(p0)

    # se a matriz Jacobiana for do tipo lista, então precisa converter
    if len(Jac) > 0: Jac = self.Jacobiana(Jac)(p0)
    # se ela não for informada, precisa calcular numericamente
    else: 
      Jac = self.jacobiana_local(self.F_lista, p0)

    # precisa calcular a inversa da Jacobiana
    J_inv = inv(Jac)

    # precisamos definir que A_1 = J_1 e calcular F(p0)
    A_inv = J_inv
    
    # arruma a solução exata, caso tenha
    if len(solucao_exata) > 0:
      if type(solucao_exata[0]) == list:solucao_exata = self.matriz(solucao_exata)
      else: solucao_exata = self.matriz([[p] for p in solucao_exata])

    # dicinário de informações
    info = { "erro": [], "erro real": [], "passo": 0, "x": [] , "jacs": [], "residuo": []}

    # caso se deseje armazenar informações de tempo
    if medir_tempo:
      info["tempo"] = dict()
      metodo = self.passo_tempo
    else: metodo = self.passo

    # aplicando
    while True:
      # aplica o método
      p1, tempo, Fp = metodo(p0, A_inv)

      # salva o valor obtido
      info["x"].append(p1)

      # salva o resíduo
      info["residuo"].append(self.norma_infinito(self.F(p1)))

      # erro
      erro = self.norma_infinito(p1 - p0)
      info["erro"].append(erro)

      # erro real
      if len(solucao_exata) > 0:
        info["erro real"].append(self.norma_2(solucao_exata - p0))

      # adicionada à quantidade de passos
      info["passo"] += 1
      
      # verifica se há quantidade exata de passos para parar
      if qntd_exata_passos > 0:
        # se tiver, verifica se já bateu
        if info["passo"] == qntd_exata_passos: 
          if exibir_causa_fim: print('[ quantidade exata de passos atingida ]')
          break
      else:
        # verifica se o passo ultrapassou o limite ou o erro ficou abaixo do admitido
        if info["passo"] >= qntd_maxima_passos:
          if exibir_causa_fim: print('[ quantidade máxima de passos atingida ]')
          break
        # caso queira verificar por limitação de ponto flutuante, verifica
        elif limitacao_float:
          if len(info["x"]) >= 2:
            if info["erro"][-1] == info["erro"][-2]:
              if exibir_causa_fim: print('[ limitação de ponto flutuante ]')
              break
        elif erro <= erro_admitido:
          if exibir_causa_fim: print('[ erro inferior ao erro admitido ]')
          break
        # se por algum acaso o erro zerar, então convergiu
        if erro == 0:
          if exibir_causa_fim: print('[ a norma da diferença zerou ]')
          break

      # agora calculamos o novo A_inv
      if medir_tempo:
        t0 = time()
        A_inv = self.novoA_inv(p0, p1, Fp, A_inv)
        tempo["A_inv"] = time() - t0
        tempo["total"] = sum(tempo[i] for i in tempo)
        
        # salva o tempo
        for etapa in tempo:
          try:    info["tempo"][etapa] += [tempo[etapa]]
          except: info["tempo"][etapa] = [tempo[etapa]]

      else: A_inv = self.novoA_inv(p0, p1, A_inv)

      # armazena a matriz jacobiana
      info["jacs"].append(A_inv)

      p0 = p1
  
    return p0, info