"""
  Método de Newton Multivariado

  O Método de Newton para funções de uma variável real pode ser usado para encontrar raízes
  da função. Já para aplicações de maiores dimensões, havendo tantas quanto variáveis, é possível
  aplicar o mesmo princípio do método para encontrar as soluções do sistema não-linear que 
  descreve a aplicação.
"""
from auxiliares.baseMetodoNumerico import BaseMetodoNumerico
from time import perf_counter as time

class Newton (BaseMetodoNumerico):
  """
    Método de Newton Multivariado.

    Parâmetros
    ----------
    F : list
      Uma lista de funções das equações do sistema.
    metodo_sistema_linear : function
      Uma função que resolva sistemas lineares, informada a matriz de
      coeficientes e a de vetores.
    metodo_diferenciacao : bool || function = False
      Uma função que será utilizada para calcular as derivadas parciais
      numericamente caso necessário.
    h : float = 0
      Tamanho do passo a ser considerado na diferenciação numérica, se
      for o caso.
  """
  def __init__ (self, F:list, metodo_sistema_linear, dominio:list=[], metodo_diferenciacao=False,h=0):
    if len(dominio) > 0:
      self.dominio = dominio
    if metodo_diferenciacao:
      super().__init__(metodo_sistema_linear, metodo_diferenciacao, h)
    else:
      super().__init__(metodo_sistema_linear)

    # adapta as funções
    self.F_lista = self.funcoes_f(F)
    self.F = self.funcao_F(self.F_lista)
  
  def passo_tempo (self, p0, Jac=[]):
    """
      Aplica um passo do Método de Newton, armazenando
      o tempo decorrido para calcular cada parte.
    """
    # dicionário para armazenar os tempos
    tempo = dict()

    t0 = time()
    # matriz Jacobiana
    if type(Jac) == list: J = self.jacobiana_local(self.F_lista, p0)
    else: J = Jac(p0)
    tempo["jacobiana"] = time() - t0

    # - F(x0)
    t0 = time()
    Fx0 = - self.F(p0)
    tempo["Fx0"] = time() - t0

    # J(x0) * y = - F(x0)
    t0 = time()
    y = self.resolver_sistema_linear(J, Fx0)
    tempo["y"] = time() - t0

    # aplica o passo
    tempo["p"] = time() - t0
    p = p0 + y

    # tempo total é a soma de todos
    tempo["total"] = sum(tempo[i] for i in tempo)

    return p, y, J, tempo

  def passo (self, p0, Jac=[]):
    """
      Aplica um passo do Método de Newton.
    """
    # matriz Jacobiana
    if type(Jac) == list: J = self.jacobiana_local(self.F_lista, p0)
    else: J = Jac(p0)

    # - F(x0)
    Fx0 = - self.F(p0)

    # J(x0) * y = - F(x0)
    y = self.resolver_sistema_linear(J, Fx0)

    # aplica o passo
    p = p0 + y

    return p, y, J, 0

  def aplicar (self, p0, Jac=[], erro_admitido:float=1e-5, qntd_maxima_passos:int=1e2, solucao_exata=[], qntd_exata_passos:int=-1, residuo_admitido:float=0, medir_tempo:bool=False, limitacao_float:bool=False, exibir_causa_fim:bool=True, adaptado:bool=False)->tuple:
    """
      Para facilitar a aplicação, pode-se utilizar esta função.

      Parâmetros
      ----------
      p0 : list || np.matrix
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
      residuo_admitido : float = 0
        Caso queira limitar via resíduo.
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
      adaptado : bool = False
        Parece delimitar o método sobre um domínio. Necessita do parâmetro 
        `dominio` na inicialização.
    """
    # se o ponto inicial for do tipo lista, precisa converter
    p0 = self.ponto_matriz(p0)

    # se a matriz Jacobiana for do tipo lista, então precisa converter
    if len(Jac) > 0: Jac = self.Jacobiana(Jac)

    # dicinário de informações
    info = { "erro": [], "erro real": [], "passo": 0, "x": [], "jacs": [], "residuo": [] }

    # caso se deseje armazenar informações de tempo
    if medir_tempo:
      info["tempo"] = dict()
      metodo = self.passo_tempo
    else:
      metodo = self.passo
    
    # arruma a solução exata, caso tenha
    if len(solucao_exata) > 0:
      if type(solucao_exata[0]) == list:solucao_exata = self.matriz(solucao_exata)
      else: solucao_exata = self.matriz([[p] for p in solucao_exata])

    # armazena as informações de erro
    self.qntd_exata_passos = qntd_exata_passos
    self.qntd_maxima_passos = qntd_maxima_passos
    self.limitacao_float = limitacao_float
    self.erro_admitido = erro_admitido
    self.residuo_admitido = residuo_admitido

    # começa o método
    while True:
      # aplica o método
      p1, y, J_cond, tempo = metodo(p0, Jac)

      # caso esteja delimitado, verifica se está dentro do domínio
      if adaptado:
        for i in range(len(p1)):
          if p1[i,0] < self.dominio[i][0]:
            # p1[i,0] = self.dominio[i][0]
            p1[i,0] = p0[i,0]
          elif p1[i,0] > self.dominio[i][1]:
            # p1[i,0] = self.dominio[i][1]
            p1[i,0] = p0[i,0]
        # nesse caso, calcula o novo y
        y = p1 - p0

      # salva o ponto
      p0 = p1

      # salva o valor obtido
      info["x"].append(p0)

      # armazena a matriz jacobiana
      info["jacs"].append(J_cond)

      # salva o resíduo
      info["residuo"].append(self.norma_infinito(self.F(p0)))
      
      # medição de tempo
      if medir_tempo:
        for etapa in tempo:
          try:    info["tempo"][etapa] += [tempo[etapa]]
          except: info["tempo"][etapa]  = [tempo[etapa]]

      # medição de erro
      erro = self.norma_infinito(y)
      info["erro"].append(erro)

      # caso possua solução exata informada, a calcula
      if len(solucao_exata) > 0:
        info["erro real"].append(self.norma_infinito(solucao_exata - p0))
      
      # adiciona à quantidade de passos
      info["passo"] += 1

      # verifica se precisa parar
      msg = self.parada(info)
      if msg:
        if exibir_causa_fim: print(f'Causa da parada: {msg}')
        break
    
    return p0, info