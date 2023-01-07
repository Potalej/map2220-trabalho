from tabulate import tabulate
import matplotlib.pyplot as plt

def exibir (cabecalho:list, colunas:list, formato:str="presto")->str:
  """ Exibe uma tabela formatada com as colunas e cabeçalho passados. """
  tabela = list(zip(*colunas))
  tabela = tabulate(
    tabular_data=tabela,
    headers=cabecalho,
    tablefmt=formato
  )
  print(tabela)
  return tabela

def exibir_vetores (info:dict, formato:str="presto")->str:
  
  # tamanho do vetor para facilitar
  n = len(info["x"][0])

  tabela = []
  for i in range(info["passo"]):
    tabela.append([i])
    tabela[-1] += [ info["x"][i][j,0] for j in range(n)]
    tabela[-1] += [ info["erro"][i] ]
    tabela[-1] += [ info["residuo"][i] ]

  tabela = list(zip(*tabela))

  cabecalho = ["i"] + [f"x{i+1}" for i in range(n)] + ["erro", "resíduo"]
  exibir(cabecalho, tabela, formato)

def grafico_tempo (etapas:list, valores:list):
  plt.bar(etapas, valores, width=.6, color='black')
  plt.xlabel('Etapa', fontsize=14)
  plt.ylabel('Razão com o tempo total', fontsize=14)
  plt.show()

def percentual (x):
  return round(x*100, 2)