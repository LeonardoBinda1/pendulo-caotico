# Novas figuras 3 e 5 do pêndulo

O manuscrito não foi alterado. Este pacote contém somente as novas famílias de figuras, o código completo, um notebook independente e os resultados dos controles numéricos.

## Abrir e executar

Abra `pendulo_figuras_jupyter.ipynb` no Jupyter. O código completo está dentro do notebook, que também pode ser usado sem o arquivo `.py`. As saídas gráficas já estão incluídas para consulta. Execute as células na ordem. Quando não há dados numéricos na pasta, o programa os calcula e os grava; a execução completa pode demorar vários minutos e produzir 70–100 MB de arquivos locais. As amostras grandes não foram incluídas no ZIP para reduzir seu tamanho; são reprodutíveis pelo programa.

Dependências: `numpy`, `scipy`, `matplotlib`, `numba`. No Jupyter, se necessário: `%pip install numpy scipy matplotlib numba`.

O arquivo `pendulo_figuras.py` contém os mesmos métodos e também pode ser executado diretamente: `python pendulo_figuras.py`. A opção `--recompute` ignora os dados já gerados. No notebook, a variável `RECALCULAR` desempenha esse papel. Para uma experiência rápida, pode reduzir `N_FRACTAL`; para a comparação quantitativa, mantenha `2**21`.

## Figura 3

O arquivo principal `figura_3_seis_regimes` compara:

| Painel | q | f | Resultado observado |
|---|---:|---:|---|
| a | 0.5 | 0.9 | Período 1 |
| b | 0.5 | 1.07 | Período 2 |
| c | 0.5 | 1.2 | Caos filamentar |
| d | 0.2 | 0.8 | Atrator caótico mais preenchido |
| e | 0.1 | 1.2 | Atrator caótico mais preenchido |
| f | 0.05 | 0.6 | Caos com amortecimento fraco |

Em todos: `Omega=2/3`, fase zero, estado inicial `(theta,omega)=(0.2,0)`, 8000 períodos descartados, 65536 pontos retidos e 480 passos por período. Os quatro casos classificados como caóticos têm expoentes positivos em todos os oito blocos da amostra. Os três casos com menor amortecimento também foram conferidos com 960 passos por período; os valores estão em `gallery_refined.json`.

A comparação varia dois parâmetros; não constitui uma única cascata de bifurcações. Para manter `q=0.5` fixo, use `figura_3_amplitudes_q_fixo`, também incluída, ou `figure3(fixed_damping=True)` no programa. Essa alternativa mostra períodos 1, 2 e 3 e janelas caóticas, sem afirmar uma cascata 1→2→4.

O artigo de referência [Rafat, Wheatland e Bedding, arXiv:0812.0393](https://arxiv.org/abs/0812.0393) usa um sistema conservativo e várias condições iniciais para preencher suas seções. Aqui não foram misturados transientes, ruído ou fases diferentes para produzir curvas mais densas. A redução de `q` muda fisicamente o atrator e permite a apresentação mais preenchida.

## Figura 5: os dois regimes estão identificados

Nos painéis a–b, as variedades são calculadas para `q=0.05, f=0.12`. Nos painéis c–f, a geometria fractal continua sendo calculada para `q=0.5, f=1.2`, o atrator do manuscrito anterior. Não se deve atribuir a dimensão deste último aos cruzamentos do primeiro.

Os cruzamentos deixaram de ser apenas inferências visuais de polilinhas. Foram resolvidos por shooting com equações variacionais, depois conferidos aumentando de 2 para 3 e 4 períodos as pernas de aproximação à sela. Uma formulação em deslocamentos em relação à órbita periódica evita perda de precisão nas sementes minúsculas. O integrador DOP853 também foi repetido com tolerância relativa dez vezes menor. Os resultados são:

| Cruzamento | theta | omega | Ângulo no plano adimensional |
|---|---:|---:|---:|
| C1 | 1.9163977629 | 1.0322926314 | 10.9025 graus |
| C2 | 2.0184937116 | 0.9109710777 | 8.3003 graus |

Nos controles em deslocamentos, os resíduos de encontro são inferiores a `2e-13` e a variação de posição entre pernas é inferior a `1e-11`. Esses números documentam convergência interna do cálculo; não equivalem a uma prova matemática com aritmética intervalar. Os ângulos não nulos dão evidência numérica clara de transversalidade. O arquivo `crossing_checks.json` contém todos os testes.

## Contagem de caixas: o que foi melhorado

Foram usados 2097152 pontos, 16 origens de grade e uma convenção periódica também para os índices angulares. Esta última evita que a costura `theta=-pi` introduza células artificiais adicionais. As coordenadas espaciais são fixas: `u=(theta_bar+pi)/(2*pi)`, `v=omega/4`. A janela principal contém cinco escalas: `delta=2**(-j)`, `j=6,...,10`. A função também calcula todas as janelas vizinhas de quatro escalas e as exibe em `diagnostico_convergencia`.

Na janela principal, as medianas entre grades são aproximadamente `D0=1.291` e `D1=1.204`. As faixas entre origens são `1.281–1.295` e `1.200–1.215`, respectivamente. O expoente medido é aproximadamente `0.13174`, dando `D_KY=1.20854`. A comparação com passo reduzido usa o mesmo número de pontos em ambas as integrações, para separar o efeito do passo do efeito do tamanho da amostra.

As bandas e barras representam a dispersão entre grades, não intervalos estatísticos de confiança. A janela principal foi escolhida examinando o comportamento intermediário, não por um protocolo cego. Caixas muito grandes ainda medem a forma global; caixas muito pequenas ainda sofrem com amostragem finita. A dependência não desaparece matematicamente em uma nuvem finita: o avanço é sua redução e quantificação explícita.

`D0` é uma dimensão do suporte. `D1` envolve as probabilidades de visitação das caixas e é a comparação mais próxima da construção de Kaplan–Yorke, sob hipóteses adicionais. Não se deve forçar igualdade entre `D0` e `D_KY`. Os expoentes negativos usados em Kaplan–Yorke são inferidos da identidade exata de contração; não foram medidos independentemente por um cálculo completo do espectro.

## Conteúdo

- PNG de alta resolução e PDF com texto e linhas vetoriais para as duas figuras principais.
- Uma versão alternativa da figura 3 com amortecimento fixo.
- Um gráfico adicional de controle da janela e do passo.
- Código `.py`, notebook independente, legendas em inglês e resultados de validação em JSON.

O código `.py` foi executado com sucesso. Todas as células do notebook foram executadas em sequência no mesmo processo Python e suas saídas capturadas. O ambiente de preparação não permitia abrir os canais de comunicação de um kernel Jupyter; portanto, não se alega um teste do transporte de kernel neste ambiente.
