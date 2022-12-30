from tabulate import tabulate

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