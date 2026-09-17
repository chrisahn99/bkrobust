# O braço Ĉ estimado — a eficácia do firewall é uma função monótona da qualidade da descoberta

Corrido em **betelgeuse**, 2026-08-23. 400 SCM por célula, p ∈ {6,7,8,9,10}, ER grau {1.5,2.0,2.5},
iSCM, α = 0,01. Código: `x2/chat/code/{pc.py,chat_arm.py}` — **PC escrito de raiz**, sem
`causal-learn` nem `gcastle`.

## Portão, corrido antes de cada célula
`test_population_recovers_cpdag`: com o teste de IC do oráculo, o PC tem de devolver **exatamente**
`dag_to_cpdag(D)`. **150/150** e **60/60 PASS** em todas as execuções. Se falhar, nada neste
ficheiro significa coisa alguma.

## O condicionamento, que na primeira tentativa estava errado
As taxas S1/S1rej só são comparáveis ao baseline se forem medidas **no mesmo conjunto**. Tudo
abaixo está condicionado a **o relatório-base já ser válido em `D`** — ou seja, mede-se
*newly invalid*, o dano que o enunciado falso acrescenta, e não o que a descoberta já tinha feito.

## A curva

| `n` | SHD(Ĉ,C) | baseline inválido | **S1** (passa Dor–Tarsi) | **S1rej** (carimbado sem teste) | separação |
|---|---|---|---|---|---|
| 200 | 3,16 | 35,9 % | **10,87 %** [9,10 – 12,78] | 17,44 % [10,61 – 24,92] | 1,6× |
| 500 | 2,35 | 28,8 % | **6,36 %** [4,92 – 8,33] | 18,53 % [13,22 – 23,99] | 2,9× |
| 2 000 | 1,97 | 28,6 % | **4,34 %** [3,37 – 5,49] | 25,18 % [18,65 – 32,33] | 5,8× |
| 10 000 | 1,20 | 22,9 % | **3,13 %** [2,30 – 4,04] | 24,70 % [19,56 – 30,55] | 7,9× |
| 50 000 | 0,78 | 17,1 % | **2,44 %** [1,71 – 3,21] | 30,20 % [25,09 – 35,31] | 12,4× |
| **∞** (IC oráculo) | **0,00** | **0,0 %** | **0,00 %** [0,0000 – 0,0000] | 34,77 % [29,74 – 39,87] | **∞** |

IC95 por **bootstrap agrupado ao nível do SCM** (os ensaios são agrupados por grafo).

## Três leituras

🎯 **1. O firewall replica de forma independente.** O braço de IC-oráculo é o mesmo objeto que
`move2.py` mediu (0 / 93 373, 40,4 %), mas por **outro código, outro sorteio de SCM e outra gama de
p**: **0 / 22 282 e 34,8 %**. A conjetura C-FIREWALL já não assenta num único script escrito por um
agente numa tarde.

🎯 **2. A eficácia do firewall é monótona na qualidade da descoberta.** 1,6× a `n`=200 → 12,4× a
`n`=50 000 → perfeita no limite do oráculo. Não é «funciona» ou «não funciona»: é uma curva, e a
curva é a figura.

⚠️ **3. E o erro de descoberta domina.** O baseline inválido — sem enunciado falso nenhum — vai de
**35,9 %** a **17,1 %** e só chega a zero com o CPDAG do oráculo. Onde o baseline já está inválido,
os enunciados falsos deixam-no inválido em ~85–93 % dos casos. **A ordem de grandeza do erro de
BK que o firewall deixa passar (2–11 %) é menor do que a do erro de descoberta (17–36 %).**
Isso tem de estar no papel: um relator aplicado vai fazer exatamente esta pergunta, e a resposta
honesta favorece o teste sem exagerar o seu alcance.

## Uma anomalia a não vender como resultado
`S1 TRUE stmt` tem taxa **superior** a `S1 FALSE stmt` em todo `n` (p.ex. 16,6 % vs 3,1 % a
`n`=10 000). Não é paradoxo: um enunciado «verdadeiro» aqui é um par não adjacente em Ĉ **mas
adjacente em `D`** — isto é, um **falso negativo do PC**. É uma subamostra selecionada para
dificuldade, logo a comparação é confundida. Reportar como diagnóstico, nunca como contraste.

## Limitações declaradas
- Só ER; sem scale-free nem cadeias de cliques (a família rica em caminhos sem cordas continua por testar).
- p ≤ 10; α = 0,01 fixo, sem varrimento.
- `bk_assert` só, um enunciado de cada vez (ρ = 1).
- X e Y de nó único.
