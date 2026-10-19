# Notas da revisão

Foi usado como base o `main.tex` do projeto do Overleaf enviado nesta conversa. O artigo permanece em inglês e em duas colunas; o título foi preservado. A versão completa tem 10 páginas, incluindo as cinco figuras, o apêndice analítico e as referências. O crescimento em relação às 7 páginas da versão anterior decorre principalmente da inclusão das figuras e da documentação dos cálculos.

## Escolha e verificação das figuras

| Material fornecido | Decisão | Motivo e posição no artigo |
|---|---|---|
| Espaço de fases do pêndulo amortecido | Substituído por curvas de energia conservativas calculadas exatamente | A seção II trata de `q=f=0`, com centros, selas e separatriz. Uma figura dissipativa não representa esse limite. A nova Fig. 1 distingue explicitamente os dois casos. |
| Melnikov e variedades | Mantido o conteúdo do painel analítico, recalculado; painéis de variedades não incluídos | A nova Fig. 2 usa `q=0.05`, `Omega=2/3` e três razões `f/f_M`. As interseções dos painéis fornecidos não podem ser certificadas somente pelo PNG: faltam método de obtenção da órbita periódica, resíduos, sementes, interpolação e teste de refinamento das variedades. A ausência dessa validação não demonstra que os painéis estejam errados. |
| Figura anterior de Poincaré e Lyapunov | Estrutura mantida e cálculos refeitos | A nova Fig. 3 compara `f=0.9` e `1.2`, incluindo médias de Lyapunov e 32000 períodos. Rótulos em inglês e parâmetros documentados. |
| Bifurcação e Lyapunov | Recalculada | A nova Fig. 4 usa reinicialização independente em cada amplitude. Isso torna explícita a escolha de bacia; pode produzir diferenças em relação a uma varredura por continuação. O marcador analítico está em `f_M=1.018774...` e é identificado na legenda. |
| Geometria fractal | Estrutura dos quatro painéis mantida e cálculos refeitos | A nova Fig. 5 mostra o atrator, ampliação, contagem de caixas e dependência com tamanho da amostra. A normalização espacial e as escalas ajustadas são explícitas. |

Os cinco gráficos são gerados de novo a partir da equação. Isso permite testar os cálculos e entregar versões vetoriais, em vez de atribuir aos PNGs parâmetros de integração desconhecidos. As legendas identificam os painéis e suas condições; o texto interpreta a física e os limites de cada diagnóstico.

## Principais correções científicas

- Definidos o tempo físico e o tempo adimensional, torques, coeficiente viscoso, frequência da força, velocidade angular e convenção angular no cilindro.
- Explicitados os autovalores de centros e selas no limite conservativo. Um centro é estável no sentido conservativo, mas não é um atrator.
- Separadas a órbita heteroclínica no plano desenrolado e a homoclínica no cilindro.
- Definidas variedades estável e instável, órbita periódica de sela, interseção transversal, tangência, zero simples e significado de uma ferradura.
- Introduzido o balanço de energia para justificar o integrando de Melnikov. A quantidade mostrada é `epsilon M`, uma separação normalizada por energia, não uma distância física.
- Acrescentada a derivação da transformada de `sech` pelo integral beta no apêndice. Esclarecido por que a igualdade no limiar não é um zero simples.
- Separados o regime perturbativo e o estudo de atratores com `q=0.5`. Melnikov não prevê, por si só, a bacia ou a existência de um atrator caótico.
- Definidos fase de amostragem, período mínimo e limitação da projeção no diagrama de bifurcações.
- Explicitados cálculo tangente, norma, renormalização, transiente tangente, unidades dos expoentes e diferença entre mapa e fluxo autônomo.
- Distinguido o parâmetro perturbativo `epsilon` da escala de observação `delta`.
- Explicado que a dimensão literal de uma nuvem finita é zero; as inclinações medem o atrator subjacente apenas em uma faixa intermediária de resolução.
- Derivada a contração por período. O expoente negativo é inferido de `lambda_-=-q-lambda_+`; não foi apresentado como se viesse de uma medição independente de todo o espectro.

## Resultados e limites de precisão

Com 240 passos por período, os exemplos longos dão `lambda_max=-0.250001` e `+0.131480`. Com 480 passos por período, dão `-0.250001` e `+0.131176`. Os oito blocos do caso caótico permanecem positivos, com médias entre `0.126760` e `0.136524`. A comparação independente com DOP853 por um período dá diferença máxima de aproximadamente `1.21e-7` no estado. A matriz variacional também reproduz `det Phi(T_D)=exp(-q T_D)`.

O Kaplan--Yorke da seção é aproximadamente `1.208`. Para a contagem de caixas com as coordenadas e janela especificadas, a inclinação é `1.163`, mas janelas vizinhas dão `1.093` e `1.206`, e deslocamentos da grade dão valores entre `1.075` e `1.163`. Portanto, a versão revisada não afirma que `D_0=1.25` esteja convergido. A proximidade entre dimensões diferentes não substitui um teste de convergência. Esse é o principal limite quantitativo que permanece no estudo.

A varredura em amplitudes não foi refinada em todos os pontos. Ela ilustra janelas periódicas e regiões irregulares para uma condição inicial; não fornece limiares exatos nem inventaria todas as bacias. As figuras de variedades também exigiriam validação adicional antes de serem usadas como prova numérica de interseções.

## Texto e bibliografia

A introdução mantém a sequência histórica: Galileo; debate sobre isocronismo; Huygens; Euler e Young; Poincaré; van der Pol, Cartwright e Littlewood; Lorenz; Peixoto e Smale. As transições conectam essa história à pergunta central do trabalho, evitando uma sucessão isolada de nomes.

Foram corrigidos título e páginas de De Jong, autor e periódico de Buyse, títulos dos trabalhos de Bistafa e a identificação do texto de Young de 1813. A referência a Young sobre interferência acústica não sustentava a afirmação sobre resistência em oscilações forçadas. Uma fonte histórica de Bleck-Neuhaus foi acrescentada para esse contexto.

As entradas atribuídas a Gutiérrez (2016) e Wickert (2018) não foram identificadas com segurança nas buscas pelos títulos fornecidos; foram retiradas da bibliografia ativa, sem afirmar que tais autores não tenham outros trabalhos relevantes. A entrada incompleta de Littlewood de 1945 foi substituída na citação pelo trabalho conjunto verificável. Pesin e Richman--Moorman foram retirados porque não são usados nesta versão, que não calcula entropia.

Os campos `autores` foram mantidos para preenchimento. Foram removidos os emails fictícios e normalizadas as afiliações; conferir os nomes e emails reais antes de encaminhar o artigo.

## Fontes consultadas para as correções bibliográficas

- De Jong: https://advlabs.aapt.org/document/ServeFile.cfm?DocID=4686&ID=14323
- Buyse: https://socpol.uvvg.ro/docs/2017-2/1.INTRO_Buyse.pdf
- Bistafa (2021): https://arxiv.org/abs/2105.03319
- Bistafa (2022): https://scholarlycommons.pacific.edu/euleriana/vol2/iss1/4/
- Young e contexto histórico: https://arxiv.org/pdf/1811.08353
- Cartwright--Littlewood: https://doi.org/10.1112/jlms/s1-20.3.180
- Smale: https://boletin.math.org.mx/pdf/2/5/5.html

O artigo não está apresentado como um novo método de detecção de caos. A contribuição proposta é a exposição integrada e reproduzível dos cálculos e das limitações dos diagnósticos.
