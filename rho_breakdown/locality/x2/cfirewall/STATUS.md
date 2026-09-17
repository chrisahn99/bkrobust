# C-FIREWALL — estado em 2026-08-23

## Veredito: **não refutado**, e agora parcialmente provado

**0 violações em 70 375 325 ensaios / 103 105 398 verificações (ensaio × D).**

| p | modo | CPDAGs | ensaios | reports | moveram `O*` | violações |
|---|---|---|---|---|---|---|
| 3 | **censo** | 11 | 78 | 72 | 18 | **0** |
| 4 | **censo** | 185 | 10 996 | 10 276 | 2 232 | **0** |
| 5 | **censo** | 8 782 | 2 161 055 | 2 053 715 | 402 840 | **0** |
| 6 | amostra | 400 000 sorteios | 22 851 273 | 22 537 373 | 4 400 761 | **0** |
| 7 | amostra | 120 000 | 16 435 189 | 16 252 091 | 2 870 115 | **0** |
| 8 | amostra | 40 000 | 10 839 188 | 10 735 797 | 1 774 275 | **0** |

Limite da regra de três sobre a taxa de violação: **< 4,3 × 10⁻⁸**.
Poder do detetor (braço de controlo, dentro do censo): **3 005 488 / 13 791 476 = 21,8 % inválidos**.

⚠️ **p=6 como censo é mensuravelmente intratável**: 5,40 s/CPDAG numa sonda uniforme ⇒ ≈160
core-hours para os seus 1 067 825 CPDAGs — e o custo é quase todo desperdiçado, porque K₆ sozinho
alcança 130 023 MPDAGs em 198 s e **produz zero ensaios** (um esqueleto completo não tem par não
adjacente). Daí amostra a 6/7/8, rotulada como tal.

## 🔴 Três achados estruturais que mudam a redação

### 1. Os três testes de coerência colapsam num só — **verificado por mim, de forma independente**

Censo escrito de raiz, p=3/4/5, **791 484 avaliações**, só duas células alguma vez povoadas:

| p | avaliações | `conflict` dispara | `cycle` ⟺ `not extendable` discorda |
|---|---|---|---|
| 3 | 84 | **0** | **0** |
| 4 | 5 940 | **0** | **0** |
| 5 | 785 460 | **0** | **0** |

⇒ **A hipótese é ACICLICIDADE. Dor–Tarsi não compra nada neste operador.** Isto simplifica o
teorema, **e torna o argumento de custo mais forte**: um teste de ciclo é O(V+E), muito mais barato
do que extensibilidade — e é ele que separa 0 % de ~35 %.

### 2. O caminho de prova oferecido é **falso**, nas duas metades
`cn_H` é subconjunto **estrito** de `cn_{G0}` em 45 720 ensaios e incomparável em 6 360; `forb`
também encolhe (160); `pa(cn)` encolhe (4 360). E **282 480 de 513 560 (55,0 %)** dos pais que a
aresta acrescenta **escapam** a `forb_H` e entram em `O*_H`. **`cn` não é monótona no conjunto de
arestas** — o fecho de Meek dentro de `H` pode orientar `v_{i+1} → v_i` numa aresta não dirigida em
`G0`, **destruindo** um caminho possivelmente causal.

✅ **Substituto, 0 exceções em 8,8 M ensaios:** `O*(X,Y,H) ∩ forb(X,Y,G₀) = ∅`. Como
`forb_D ⊆ forb_{G₀}` para `D ∈ [G₀]`, isto **liquida a cláusula (1)** do critério de ajuste para a
classe inteira de uma vez. Falta a cláusula (2).

### 3. A forma forte «todo `D ∈ [G₀]`» **vale**, e é empiricamente equivalente à forma de um só `D`
28,1 % dos reports vieram de `|[G₀]| ≥ 2` (máx. 156), 2 954 238 dos quais moveram — logo a forma
forte é exercitada, não vazia. Decisivo: no braço de controlo, onde 3 005 488 conjuntos **são**
inválidos, `ctrl_violating_someD = 0` — **a validade nunca se parte dentro de `[G₀]`**.
⇒ O sorteio único de `move2.py` não estava a subestimar nada. Agora está medido.

## Estado da prova
- **Cláusula (1)** — provada, sob uma hipótese extra **(M)**: `poss_de_{G₀}(X) ⊆ poss_de_H(X)`,
  que vale em **100 %** dos ensaios do censo p≤5 em que `O*` se moveu.
- **Cláusula (2)** — condição suficiente fecha **98,0 %** das verificações (ensaio × D).
- 🔴 **Os ~2 % restantes da cláusula (2) estão genuinamente em aberto**, e mostrou-se que
  **nenhum argumento de contenção lá chega** (`O*_H ⊊ O*_{G0}` em 123 000 ensaios, incomparável em
  24 480).

## Bónus: sobrevive a k ≥ 2
18 077 546 ensaios, 6 372 313 moveram, **0 violações**, com censos completos em (p=4,k=2),
(p=4,k=3), (p=5,k=2). **Não segue por indução** — depois de uma asserção, `H` já não é MPDAG de
nenhum CPDAG que contenha `D`.

## Higiene
`cfirewall/proof/` tem **87 ficheiros .py de vários agentes concorrentes**, todos de hoje. Um
agente sinalizou que `cat >` pode ter sobrescrito ficheiros homónimos. Os resultados vivem em
`results/` com nomes por agente e não foram afetados, mas **a pasta `proof/` deve ser tratada como
rascunho, não como código canónico**.
