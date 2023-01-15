# Trabalho computacional final de MAP2220 - Fundamentos de Análise Numérica

Na natureza, é muito comum encontrar fenômenos cuja descrição matemática se dá na forma de equações não-lineares, seja em dinâmica dos fluidos, mecânica celeste ou mesmo em bioinformática. Porém, a resolução de sistemas de equações não-lineares é uma tarefa árdua, quando não impossível, de ser feita analiticamente, ou seja, com ferramentais matemáticos exatos. Dessa forma, surge a necessidade de ter métodos numéricos, que se aproximam das soluções sempre com alguma margem de erro. Neste trabalho são estudados brevemente dois destes método: o Método de Newton e o Método de Broyden. Após a teoria também são resolvidas os três itens solicitados no enunciado do trabalho.

## Regras e observações

- Os módulos do python relativos a solução de problemas não-lineares e otimização **não podem** ser usados (*scipy.optimize*);
- Os módulos relativos a soluções de sistemas linares **podem** ser utilziados (*scipy.linalg* ou *numpy.linalg*);
- Pesquise formas de medir o tempo computacional e explique em seu relatório o que foi utilziado. A medida do tmepo computacional pode requerer várias execuções do código para que a medida de tempo seja confiável. Esteja atento a isso.
- Dependendo das condições iniciais pode ser necessário modificar as variáveis para que elas permaneçam dentro do domínio da solução. No caso de raízes de números negativos utilize o valor absoluto com prevenção.
- Inclua comentários no seu código esclarecendo suas escolhas e decisões de implementação.
- A entrega é um relatório, jupyter notebook pode ser usado, mas a entrega é um relatório em pdf que pode conter recortes, trechos, do jupyter notebook.
- A avaliação do trabalho não se limita apenas a obtenção dos resultados numéricos. Os resultados devem ser usados para suportar as análises que permitem comprovar as expectativas teóricas.