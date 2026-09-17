# Não desistir valeu: a estatística estava certa, a população é que estava errada

**2026-09-08.** Três experimentos, ~18 s de CPU, e um resultado que muda o que o papel mede.

---

## Onde estávamos

`r_val` tem 97,7 % da massa em `{1, ∞}` — quase binário. Das quatro direções propostas, duas
morreram por lema (L2 tirou "identificar as afirmações falsas" da tabela de largura) e duas
pareciam mortas por herança (um quantil de variável quase binária continua quase binário).

A aposta de largura passa com um erro (0,710) e falha com dois (0,903), e o casco seletivo é
idêntico ao blanket em 67,2 % dos problemas. Não era uma semana boa.

---

## O que eu tentei

**A tradução literal de Rosenbaum.** Ele nunca perguntou "o conjunto de ajuste ainda é válido?".
Ele pergunta **quanto viés seria preciso para virar a conclusão**. A tradução direta:

> **Raio de nulidade** $r_0$: o menor orçamento de retração em que o casco dos efeitos possíveis
> passa a conter zero. Em palavras: *quantas afirmações do especialista teriam de cair antes de eu
> deixar de poder afirmar que existe um efeito?*

Primeira medição, sobre todos os problemas: **84 % "nunca"**. Parecia pior que `r_val`.

---

## 🔑 O que estava errado não era a estatística

Só **21,2 %** dos problemas sorteados são informativos, no sentido de a classe de equivalência
deixar alguma incerteza sobre o efeito. Nos outros 79 % o casco é um ponto e **nenhuma estatística
de robustez tem o que medir**. Toda medição anterior estava diluída por eles.

Condicionando no estrato informativo, as duas estatísticas **trocam de lugar**:

| | todos os problemas | **informativos** |
|---|---|---|
| `r_val` = 1 | 43,9 % | **89,5 %** |
| `r_val`, valores efetivos | 2,28 | **1,51** |
| raio de nulidade, valores efetivos | 2,66 | **3,46** |
| raio de nulidade, níveis com ≥5 % da massa | 3 | **4** |

`n = 1800` problemas, 381 informativos. Valores efetivos $=2^H$, com $H$ a entropia da
distribuição.

![](estrato.pdf)

**Condicionar piora `r_val` e melhora o raio de nulidade.** No estrato em que a pergunta faz
sentido, `r_val` colapsa para quase constante (89,5 % num único valor, 1,51 valores efetivos) e o
raio de nulidade abre (14,2 / 22,8 / 8,4 / 2,1 / 52,5 sobre `{0,1,2,3,nunca}`, 3,46 valores
efetivos, **três vezes mais informativo em bits**: 0,60 → 1,79).

---

## Por que isso é a coisa certa, e não um conserto cosmético

1. **É o objeto de Rosenbaum, não uma variante do nosso.** Γ pergunta o que derruba a conclusão;
   $r_0$ pergunta o que derruba a conclusão. A linhagem que o papel invoca desde julho passa a ser
   a linhagem que ele de fato instancia.
2. **Sobrevive à crítica que matou o resto.** A tela de um bit ("alguma afirmação toca a
   consulta?") prediz quebra de *validade*. Ela não diz nada sobre em que orçamento o efeito deixa
   de ser distinguível de zero.
3. **Explica o achado da direção em vez de ser atropelado por ele.** Dois terços da largura vêm da
   classe discordar se $X$ causa $Y$. Um raio que mede exatamente *quando a existência do efeito
   deixa de estar garantida* está medindo aquilo de propósito, com o nome certo.
4. **É acionável.** "Sua conclusão de que há efeito sobrevive à retirada de 2 das 4 afirmações"
   é uma frase que um clínico contesta. "$r_{val} = 1$" não é.

---

## O que muda no papel

- **A quantidade principal passa a ser $r_0$**, com `r_val` reportado ao lado como o objeto de
  validade que ele é. Não é abandonar o raio; é medir a ruptura da *conclusão* em vez da ruptura
  de uma propriedade técnica do conjunto de ajuste.
- **O estrato informativo vira condição de escopo declarada**, não uma nota. Todo número condiciona
  nele, e a fração (21,2 %) é ela própria um resultado: diz quando conhecimento prévio importa.
- **`r_val` quase constante no estrato informativo é um achado**, não um defeito a esconder: dentro
  dos problemas com ignorância real, a validade do conjunto ótimo quebra na primeira retração
  quase sempre. Isso é uma afirmação estrutural sobre a classe de problemas.

## O que ainda não foi feito

- $r_0$ não foi levado à tabela de largura: falta medir a largura que ele exige para cobertura
  honesta de 95 %, que é a métrica única escolhida.
- Tudo continua sobre a classe de equivalência **verdadeira** e mecanismos lineares-gaussianos.
- O estrato informativo foi definido por largura de blanket positiva. Um referee vai perguntar se
  isso é seleção pós-tratamento; a resposta é que a largura do blanket não depende de $K$, mas
  precisa estar escrita.

## Arquivos

`signradius.py` o raio de nulidade · `stratum.py` a comparação condicionada, com entropia ·
`fig_stratum.py` → `estrato.pdf` · roda tudo em ~18 s.

---

# Segunda rodada: a lei de localidade, medida na moeda que importa

A tese empírica da linha desde julho é *o dano é local*. Ela sempre foi medida em taxa de falha
silenciosa. **Nunca foi medida em cobertura**, que é a moeda da tabela principal, e nunca teve
controle negativo. Agora tem.

## O teste

Todas as políticas retraem **no máximo uma** afirmação — orçamento casado, para que a comparação
não seja sobre quanto se retrai:

| política | retrai |
|---|---|
| `PLAIN` | nada |
| `RH1_far` | qualquer uma das que **não** tocam a consulta |
| `RH1_near` | qualquer uma das que **tocam** a consulta |
| `RH_1` | qualquer uma |
| `BLANKET` | a classe inteira |

600 problemas informativos por braço. `local_hedge.py`.

## O resultado

**Um erro:**

| política | cobre | largura/blanket | enumerações |
|---|---|---|---|
| `PLAIN` | 0,685 | 0,014 | 1,0 |
| **`RH1_far`** | **0,687** | 0,023 | 2,2 |
| `RH1_near` | **0,998** | 0,621 | **3,1** |
| `RH_1` | 1,000 | 0,630 | 4,3 |
| `BLANKET` | 1,000 | 1,000 | 1,0 |

**Dois erros:** `PLAIN` 0,363 · **`RH1_far` 0,370** · `RH1_near` 0,760 · `RH_1` 0,763 ·
`BLANKET` 1,000.

![](localidade.pdf)

## 🟢 O controle negativo é o resultado mais forte da página

**Retrair as afirmações que não tocam a consulta move a cobertura em $+0{,}2$ ponto percentual
com um erro e $+0{,}7$ com dois.** Elas custam largura e não compram nada.

Isso é a lei de localidade medida em cobertura, com controle negativo, pela primeira vez. Até
hoje a linha a media em taxa de falha silenciosa, que é um evento intermediário; agora está
medida na quantidade que o praticante reporta.

## 🔴 E ela não paga em largura

`RH1_near` recupera praticamente toda a cobertura de `RH_1` (0,998 contra 1,000; 0,760 contra
0,763) — mas por **0,621 contra 0,630 de largura**, uma economia de **1,5 %**. Restringir à
vizinhança da consulta é essencialmente **de graça**, e não é **barato**.

⚠️ E isso não é uma surpresa: **é L2 outra vez.** O lema dizia que uma tela é neutra em largura
porque as afirmações filtradas não estavam alargando o casco. Aqui isso aparece como medida, em
outra moeda. **A lei de localidade e L2 são o mesmo fato visto de dois lados**, e o papel deve
dizer isso em vez de apresentá-los como dois resultados.

## 🎯 Em que moeda a localidade paga, então

| moeda | ganho |
|---|---|
| largura | **1,5 %** — desprezível |
| **custo** | **28 %** menos enumerações (3,1 contra 4,3) |
| **elicitação** | não peça ao especialista as afirmações longe da consulta: elas movem a cobertura em $\le 0{,}7$ pp |

**A segunda linha é um resultado de sistemas; a terceira é um resultado de protocolo, e é a que
um clínico usa.** Nenhuma das duas é uma linha na tabela de largura, e insistir em pôr localidade
ali é o erro que este experimento fecha.

## E a hipótese que eu mesmo tinha, morta

Eu supus que um hedge restrito à vizinhança seria **mais estreito** ao mesmo nível de cobertura.
A primeira versão do teste (`LOCAL`, retraindo qualquer subconjunto das próximas) deu
**0,822 de largura contra 0,630 de `RH_1`** — mais larga e sem cobrir melhor. A hipótese estava
errada, e o erro era de desenho: `LOCAL` retrai mais, não melhor. A versão de orçamento casado é
a comparação correta e ela dá o resultado acima.
