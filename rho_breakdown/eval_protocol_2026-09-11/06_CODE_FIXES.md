# O diagnóstico do código, e os nove passos

**Escrito 2026-09-11**, depois que a apostila expôs um erro meu e a investigação dele mudou o
diagnóstico duas vezes. Todo número aqui sai de `diagnostico_amenabilidade.py`, que roda read-only
contra `~/bkrobust` e grava `diagnostico_amenabilidade.json`.

## O que eu afirmei e tive de retirar

**Afirmei** que `g0_undirected_edges = 0` nas 831 instâncias admissíveis, e que portanto a bola
estava centrada na verdade no corpus inteiro. **Errado.** O campo só é atribuído no ramo
`O_INTRACTABLE` de `measure.py` (linha 240); em todos os outros ele fica no valor padrão do
`dataclass`. Reconstruindo `G0` instância a instância: **543** totalmente orientadas (cobertura 1,0,
logo iguais ao DAG verdadeiro) e **288** de cobertura parcial, **todas** com aresta livre.

**Afirmei** que Lever 0 devolveria a população apagada. **Errado, e medido.**

## O diagnóstico das 518 rejeições

| das 518 `o_g0_not_identified` | |
|---|---|
| não amenáveis: o efeito não é identificável por ajuste | **518** |
| amenáveis, com conjunto válido recuperável pela forma fechada | **0** |

Não há nada a certificar nelas. Trocar o critério de ótimo pela forma fechada recupera **zero**.
A deleção está substantivamente certa; o **rótulo** é que está errado, porque diz
`o_g0_not_identified` quando o fato é não-amenabilidade.

## Onde o defeito está de verdade

Na regra de desbaste do `select_knowledge`: ela ordena o conjunto recuperador **alfabeticamente** e
retém um passo constante. Testando, para as rejeições das oito redes menores, todas as escolhas
possíveis de `|K|` afirmações do mesmo tamanho (1.238 subconjuntos):

| entre 40 casos não degenerados | |
|---|---|
| existe outro subconjunto do mesmo tamanho que **preserva** a amenabilidade | **29** (72%) |
| nenhum subconjunto daquele tamanho funciona | **11** (28%) |
| degenerados à parte (`round(1 × 0,5) = 0`, o analista fica sem conhecimento) | 61 |

Na maioria dos casos testáveis a instância não morre porque o analista sabe **pouco**, e sim porque
sabe as coisas **erradas** — e quais são as erradas é decidido por `sorted()`. ⚠️ Escopo: 40 casos
em 8 redes pequenas, não as 518.

## O reenquadramento

As 518 não são dados faltando: são o estrato em que o raio já vale **zero**. Apagá-las é remover a
cauda esquerda da distribuição e reportar a média do que sobrou. Isso sugere medir o **raio de
identificabilidade** — quantas afirmações podem ser retiradas antes de o efeito deixar de ser
estimável por ajuste — que é monótono (retrair só acrescenta arestas livres, e amenabilidade só
piora) e está definido exatamente onde o raio atual morre.

## Os nove passos

### Camada 1 — horas de trabalho, muda o denominador do estudo

1. **Trocar o desbaste determinístico por amostragem com semente**, e reportar a distribuição sobre
   subconjuntos em vez de um ponto. `measure.py#select_knowledge`. Maior efeito por linha escrita.
2. **Separar o rótulo de rejeição**: `not_amenable` contra `optimal_not_identified`. Hoje as 518
   saem com o rótulo errado. `measure.py#evaluate`.
3. **Atribuir `g0_undirected_edges` em todos os ramos**, não só no de enumeração intratável.
   `measure.py:240`.
4. **Veredito de três valores para `o`**: `if not o` confunde `None` com conjunto vazio, então as
   518 não podem ser decompostas nos dois estados que o protocolo declara.
5. **Desfazer a confusão estimador/resposta**: subir `DEFAULT_SEARCH_BUDGET` (hoje 3, em
   `hybrid.py`) acima do maior `|K_G0|`, ou imprimir a perna de despacho como coluna.
   Medido: `local_up_fast` respondeu 811 instâncias, todas em raio {1,2,3}; `e1_ladder` 20, todas
   em {4..14}.
6. **Commitar o script que escreve `results/axisa2/random_control.jsonl`**. Nenhum script na árvore
   o produz, então a regra de seleção de `K` daquele braço não é verificável.

### Camada 2 — a semana

7. **Reportar amenabilidade como estatística de escopo de primeira classe**, por cobertura e por
   rede. É o número que decide se o arcabouço tem público, e hoje está escondido dentro de uma
   contagem de rejeição.
8. **Definir e medir o raio de identificabilidade** nas 518, para ver se a cauda esquerda tem
   estrutura ou é toda zero.

### Camada 3 — não é código, e nenhum patch substitui

9. **`select_knowledge` recebe o DAG verdadeiro.** Consertar o desbaste melhora o realismo do
   *tamanho* de `K`, não da sua *origem*. A origem só muda com conhecimento elicitado: o
   questionário e os quatro fornecedores de `02_PROTOCOL.md`.

## Reproduzir

```bash
python3 diagnostico_amenabilidade.py     # ~5 s, read-only contra ~/bkrobust
```
