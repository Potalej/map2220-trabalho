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

  def bateria_testes_tempo (self, qntd_vezes:int, params:list, etapas:list=[]):
    """
      Para facilitar a medição do tempo consumido pelo método instanciado na
      aplicação.
    """
    tempo_medio = {}
    for vez in range(qntd_vezes):
      # aplica o método
      x_vez, info_vez = self.aplicar(*params)
      for indice in info_vez['tempo']:
        try:    tempo_medio[indice] += sum(info_vez['tempo'][indice])
        except: tempo_medio[indice]  = sum(info_vez['tempo'][indice])
    
    percentuais = {}
    for indice in tempo_medio:
      if indice != 'total': percentuais[indice] = tempo_medio[indice] / tempo_medio['total']

    if etapas == []: etapas = [*percentuais.keys()]

    print(f"> TEMPO TOTAL: {round(tempo_medio['total'], 4)}s")
    print(f"> TEMPO TOTAL MÉDIO: {round(tempo_medio['total']/qntd_vezes, 6)}s")
    print("> PERCENTUAIS")
    for indice in percentuais:
      print(f"* {indice}: {percentual(percentuais[indice])}%")

    grafico_tempo(etapas, percentuais.values())
  