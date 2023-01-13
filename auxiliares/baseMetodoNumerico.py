# -*- coding: UTF-8 -*-
"""
  CLASSE PARA APLICAÇÃO DE MÉTODOS NUMÉRICOS

  Esta classe contém facilidades numéricas básicas para servir de base
  na criação de classes de métodos numéricos para a resolução de sistemas
  de equações não-lineares.
"""
from auxiliares.diferenciacao import Diferenciacao
from auxiliares.exibir import grafico_tempo, percentual

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
      Para resolver sistemas lineares, independente de qual método tenha
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

  def bateria_testes_tempo (self, qntd_vezes:int, params:list, etapas:list=[], exibirGrafico:bool=False):
    """
      Para facilitar a medição do tempo consumido pelo método instanciado na
      aplicação.
    """
    tempo_medio = {}
    qntd_total_passos = 0
    for vez in range(qntd_vezes):
      # aplica o método
      x_vez, info_vez = self.aplicar(*params)
      qntd_total_passos += info_vez["passo"]
      for indice in info_vez['tempo']:
        try:    tempo_medio[indice] += sum(info_vez['tempo'][indice])
        except: tempo_medio[indice]  = sum(info_vez['tempo'][indice])
    
    percentuais = {}
    for indice in tempo_medio:
      if indice != 'total': percentuais[indice] = tempo_medio[indice] / tempo_medio['total']

    if etapas == []: etapas = [*percentuais.keys()]

    print(f"> TEMPO TOTAL: {round(tempo_medio['total'], 4)}s")
    print(f"> TEMPO TOTAL MÉDIO: {round(tempo_medio['total']/qntd_total_passos, 6)}s")
    print(f"> QNTD MÉDIA DE PASSOS: {round(qntd_total_passos/qntd_vezes, 2)}")
    print("> PERCENTUAIS")
    for indice in percentuais:
      print(f"* {indice}: {percentual(percentuais[indice])}%")

    if exibirGrafico: grafico_tempo(etapas, percentuais.values())
  
  def parada (self, info:dict)->bool:
    msg = False
    # verifica se há quantidade exata de passos para parar
    if self.qntd_exata_passos > 0:
      # se tiver, verifica se já bateu
      if info["passo"] == self.qntd_exata_passos:
        msg = 'quantidade exata de passos atingida'
    else:
      # verifica se o passo ultrapassou o limite ou o erro ficou abaixo do admitido
      if info["passo"] >= self.qntd_maxima_passos:
        msg = "quantidade máxima de passos atingida"
      # caso queira verificar por limitação de ponto flutuante, verifica
      elif self.limitacao_float:
        msg = self.criterio_float(info)
      elif info["erro"][-1] < self.erro_admitido:
        msg = 'erro inferior ao erro admitido'
      # se por algum acaso o erro zerar, então convergiu
      if info["erro"][-1] == 0:
        msg = 'a norma da diferença zerou'
    return msg     

  def criterio_float (self, info:dict)->bool:
    """
      Aplica o critério de parada por ponto flutuante. Levanta em conta o resíduo e
      a diferença relativa.
    """
    if len(info["x"]) >= 3:
      if info["erro"][-1] == info["erro"][-3]:
        return 'limitação de ponto flutuante (dif.)'
      elif info["residuo"][-1] == info["residuo"][-3]:
        return 'limitação de ponto flutuante (res.)'
    else:
      return False