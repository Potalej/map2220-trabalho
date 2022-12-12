# map2220-trabalho
Trabalho computacional final de MAP2220 - Fundamentos de Análise Numérica.

## Coisas para fazer
- [ ] (provavelmente) Cálculo de derivada de funções;
- [ ] Método de Newton multivariado para resolver equações não lineares;
  - [ ] Norma do resíduo;
  - [ ] Critério de convergência;
- [ ] Método de Broyden;
- [ ] Aproximação numérica da Jacobiana via Diferenças Finitas por expressões de 1ª e 2ª ordem;
  - [ ] Definir o passo h;
  - [ ] Usar essa Jacobiana para aplicar o Método de Newton;
  - [ ] Repetir a ideia usando o Método de Broyden (nesse caso apenas 1ª Jacobiana é aproximada);
- [ ] Resolver aleatoriamente um sistema doidão;

## Regras e observações

- Os módulos do python relativos a solução de problemas não-lineares e otimização **não podem** ser usados (*scipy.optimize*);
- Os módulos relativos a soluções de sistemas linares **podem** ser utilziados (*scipy.linalg* ou *numpy.linalg*);
- Pesquise formas de medir o tempo computacional e explique em seu relatório o que foi utilziado. A medida do tmepo computacional pode requerer várias execuções do código para que a medida de tempo seja confiável. Esteja atento a isso.
- Dependendo das condições iniciais pode ser necessário modificar as variáveis para que elas permaneçam dentro do domínio da solução. No caso de raízes de números negativos utilize o valor absoluto com prevenção.
- Inclua comentários no seu código esclarecendo suas escolhas e decisões de implementação.
- A entrega é um relatório, jupyter notebook pode ser usado, mas a entrega é um relatório em pdf que pode conter recortes, trechos, do jupyter notebook.
- A avaliação do trabalho não se limita apenas a obtenção dos resultados numéricos. Os resultados devem ser usados para suportar as análises que permitem comprovar as expectativas teóricas.
