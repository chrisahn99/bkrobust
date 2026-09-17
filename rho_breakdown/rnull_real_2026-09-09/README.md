# `r₀` sobre estrutura real, a partir de `Σ̂`

Responde à pergunta de 09/09: *podemos calcular o raio de nulidade nos dados reais?*
**Sim.** Nenhuma função no caminho do `r₀` recebe o DAG verdadeiro; ele só aparece
para gerar a amostra, como a natureza faria, e no fim para auditar.

## Substrato

Quatro redes bayesianas **gaussianas** do corpus de 39 redes do `~/bkrobust`
(`results/axisa3/networks/example_models/`). Os coeficientes não são simulados: são
ajustes a dados reais. `ecoli70` vem de microarranjos de *E. coli*.

| rede | nós | arcos | não orientadas |
|---|---|---|---|
| ecoli70 | 46 | 70 | 25 (35,7%) |
| arth150 | 107 | 150 | 45 (30,0%) |
| magic-niab | 44 | 66 | 10 (15,2%) |
| magic-irri | 64 | 102 | 12 (11,8%) |

## Arquivos

| | |
|---|---|
| `load_real.py` | carrega uma rede gaussiana como `LinearSEM` com os coeficientes ajustados |
| `rnull_real.py` | o estimando e o erro-padrão a partir de `Σ` (populacional **ou** amostral); o casco; os raios |
| `sweep.py` / `sweep_fast.py` | a varredura, 2 coberturas × 4 tamanhos × 5 réplicas |
| `direcao.py` / `direcao_fast.py` | a direção do erro em amostra finita — a tabela que decide o protocolo |
| `probe.py`, `diag.py`, `diag2.py` | os três diagnósticos que acharam os erros abaixo |

```bash
python3 -u sweep_fast.py > sweep_fast.log   # ~90 s
python3 -u direcao_fast.py
```

## Resultados

`sweep_fast.log`, `direcao_fast.py`. 226 instâncias, 5 réplicas por tamanho.

- `r₀ = 1` em **34%** das consultas a cobertura plena; `r₀ = ∞` em 66%.
- Entre os `r₀` finitos, o zero é **estrutural** em **82%**: alguma extensão dá
  efeito exatamente nulo, isto é, nela `X` nem é causa de `Y`.
- **O casco pontual erra para o lado perigoso 24,7% (n=200) a 18,1% (n=20.000)**;
  a banda de 95% erra ≤ 1,7% em todo `n`. O protocolo tem de ser a banda.
- **O critério de validade em MPDAG custa 30,7 min por chamada** numa rede de 46
  nós, contra 0,02–0,04 s para um `r₀` inteiro. O gargalo é o portão.

## Três erros meus, e o que os pegou

1. **`has_zero` num fio de navalha.** O casco do oráculo termina em `-0.0000`
   (≈ −1e-17), então `lo <= 0 <= hi` dava `False` num zero que é estrutural.
   Pegou: `diag2.py`, comparando instância a instância oráculo contra `n=20.000` e
   achando 11 divergências em 60 com todos os extremos exatamente em zero.
   Corrigido com tolerância `1e-9` e com a coluna `structural_zero`, que virou
   resultado em vez de defeito.
2. **Seleção induzida por mim.** Filtrar as instâncias por `Kloc` não vazio mantinha
   só as componentes que o conhecimento orientava por completo, o que fixava
   `|Θ(G₀)| = 1`. Isso está *correto* como premissa (com conhecimento pleno o efeito
   é identificado num ponto) mas tem de ser dito, não descoberto pelo leitor.
3. **`r_val` foi tirado do experimento, não esquecido.**
   `is_valid_adjustment_set_mpdag` levou **1.839,78 s (30,7 min) numa única
   chamada** numa rede de 46 nós (`probe.log`) — a primeira estimativa, de
   ">120 s", era só o meu limite de espera. Uma chamada por grafo da bola torna o
   experimento inviável. O Chris já mede `r_val` em estrutura real; o que faltava
   era `r₀`.

## Censurado

`arth150` a cobertura 0,5 excedeu 20 minutos e **não** conta como medição: a
enumeração de extensões explode em componentes sub-determinadas de 107 nós. As
tabelas do `sweep_fast` são das três redes que terminaram.
