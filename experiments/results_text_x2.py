"""Commentary of section 18 of RESULTS.md (experiment X2, Gate 2 of CHECKLIST_AA.md). Numbers come from the CSVs."""

X2_TEXT = r"""---

## 18. X2 — Puerta 2: ¿cierra el método con pocas fuentes?

Es la Parte 2 de `CHECKLIST_AA.md`. Las decisiones de la Parte 0 siguen siendo del grupo; aquí se usaron los
valores provisionales anotados en el checklist.

**Modelo.** K símbolos de zona. Todas las cargas de una zona escalan por (1 + δ·ε_k), P y Q juntas. Las zonas
salen solo de la topología (bisección espectral, igual número de barras por zona, anidadas al duplicar K); no
son áreas reales. δ es el de diseño: 15 % en IEEE-33 y 5 % en PEGASE. En PEGASE se corrieron dos modelos de
balance: *slack única* (la barra slack absorbe todo el desbalance, como en el manuscrito) y *participación* (la
generación sigue al cambio de demanda en partes proporcionales al despacho).

**Referencia.** Con K ≤ 16, flujo AC exacto en los 2^K vértices de la caja de ruido. Con K = 32 y 64, en los dos
vértices extremos de cada barra y en el vértice que minimiza su modelo de segundo orden. Además, 20.000 muestras
AC por caso. Ningún vértice ni muestra quedó sin solución AC.

**Intervalo que se prueba.** [V* − r_i + lo_i, V* + r_i + hi_i], con r_i el radio afín y [lo_i, hi_i] una cota de
εᵀH_iε. Es una cota de segundo orden: el resto de orden δ³ se mide, no se acota, y no hay redondeo dirigido.

**Cotas que se comparan.** *AF2*: la forma cuadrática ensamblada, monomio a monomio; es también lo que hace
Liao et al. (2017) con intervalos. *RAA-K*: el operador del grupo, portado de `ibex-lib/raak/affine/ibex_KNomial.h`
a `cdhim/knomial.py` (reproduce los 80 casos de `knomial_golden.txt` sin ninguna diferencia), con cliques de
hasta 2, 3, 6 y K variables; con 6 o más, la búsqueda de cliques tiene un tope de 20.000 nodos por forma.
*Bloques espectrales*: los autovectores de H_i como bloques, acotados con el mismo teorema k-nomial. Ninguna
redondea hacia fuera.

### 18.1 El escenario

{{x2_setting}}

*Script: `experiments/x2_few_sources.py` · CSV: `results/x2_summary.csv`, `results/x2_bus.csv`*

- **Con pocas fuentes la caja 3σ deja de tener sentido.** En PEGASE con K = 8 es más ancha que la caja afín en
  <<X2_S3_WIDER_P8>> de <<X2_NBUS_P8>> barras, es decir, se sale del soporte, y aun así el <<X2_S3_SYS_P8>> % de
  las muestras tiene alguna barra fuera, porque está centrada en V* y el voltaje real cae por debajo.
- La caja afín queda a entre <<X2_AFF_MC_LO>> y <<X2_AFF_MC_HI>> veces el rango de Monte Carlo al 99,73 %. Con
  2.951 símbolos independientes eran 3 a 4 veces (C8).
- **El coste baja a K resoluciones.** El primer orden tarda <<X2_T1_MAX>> s como máximo en PEGASE; el modelo de
  segundo orden para todas las barras, <<X2_T2_MAX>> s, y las cotas, hasta <<X2_TB_MAX>> s (tabla de 18.6).
- El modelo de participación deja el radio del ángulo en <<X2_ANG_RATIO_LO>> a <<X2_ANG_RATIO_HI>> veces el de la
  slack única, pero no cambia el cuadro: las barras que más se apartan son barras débiles, donde el voltaje cae
  hasta <<X2_DV_MAX>> p.u. en un vértice.

### 18.2 La caja afín de primer orden

{{x2_affine}}

- **Con pocas fuentes la caja afín falla también por muestreo, no solo en los vértices.** En PEGASE con K = 8,
  el <<X2_AFF_MCS_P8>> % de las muestras (participación) y el <<X2_AFF_MCS_S8>> % (slack única) tiene alguna barra
  fuera. Con un símbolo por carga el Monte Carlo no lo veía (sección 17.2).
- En los vértices el exceso llega a <<X2_AFF_SHORT_MAX>> p.u. de voltaje, y entre <<X2_AFF_N3_LO>> y
  <<X2_AFF_N3_HI>> barras de PEGASE exceden su caja en más de 0,001 p.u.
- En IEEE-33 el exceso es <<X2_AFF_EXC_I_MED>> % del radio en mediana y <<X2_AFF_EXC_I_MAX>> % como máximo, y ninguna
  muestra sale.

### 18.3 Tamaño de cada orden

{{x2_orders}}

- En la barra mediana de PEGASE el segundo orden es <<X2_D2_MED_LO>> a <<X2_D2_MED_HI>> % del radio y el tercero
  <<X2_D3_MED_LO>> a <<X2_D3_MED_HI>> %. El tercer orden es el <<X2_32_MED_LO>> a <<X2_32_MED_HI>> % del segundo en
  mediana y el <<X2_32_P90_LO>> a <<X2_32_P90_HI>> % en el percentil 90.
- **En las peores barras la serie converge despacio.** El segundo orden supera varias veces el radio y el
  tercero llega a <<X2_D3_PU_LO>>–<<X2_D3_PU_HI>> p.u. Son barras de radio grande, no artefactos de radios
  diminutos: en la peor, el radio afín es <<X2_WB_R>> p.u., el segundo orden agrega hasta <<X2_WB_D2>> y el tercero
  hasta <<X2_WB_D3>> (participación, K = 8; `results/x2_bus.csv`, barra <<X2_WB_BUS>>).
- **El resto ya no es de un solo signo en PEGASE.** H_i es semidefinida negativa en todas las barras de IEEE-33,
  pero solo en el <<X2_NSD_LO>> a <<X2_NSD_HI>> % de las de PEGASE, y la proporción baja al subir K.

### 18.4 Criterio 1: ¿la cota de segundo orden contiene la solución AC?

{{x2_gate1}}

**No se cumple.**

- **Una cota de segundo orden ajustada falla exactamente por el tercer orden.** En IEEE-33 el lado inferior de
  AF2 y de RAA-K coincide con el mínimo exacto de la forma cuadrática, y la solución AC queda
  <<X2_MONO_SHORT_I>> p.u. por debajo en las 32 barras. Los bloques espectrales "contienen" solo porque
  sobreestiman la profundidad (18.6); eso es holgura, no garantía.
- **En PEGASE fallan todas.** Con participación y K = 8, los bloques espectrales dejan <<X2_SPEC_OUT_P8>> barras con
  algún vértice fuera, <<X2_SPEC_N4_P8>> por más de 10⁻⁴ p.u. y <<X2_SPEC_N3_P8>> por más de 10⁻³ p.u.; el déficit
  máximo es <<X2_SPEC_SHORT_P8>> p.u. Con slack única: <<X2_SPEC_OUT_S8>>, <<X2_SPEC_N4_S8>>, <<X2_SPEC_N3_S8>> y
  <<X2_SPEC_SHORT_S8>> p.u.
- **Por muestreo, en cambio, casi no falla.** Como máximo el <<X2_SPEC_MCS_MAX>> % de las muestras tiene alguna
  barra fuera de la cota de segundo orden (<<X2_SPEC_MCB_MAX>> barras a lo más), contra <<X2_AFF_MCS_P8>> % y
  <<X2_AFF_MCS_S8>> % de la caja afín. La corrección de segundo orden es lo que hace usable el intervalo; lo que
  no da es contención en el peor caso.

Barrido de δ con el menor K de cada sistema:

{{x2_sweep}}

- El déficit crece como δ³: en PEGASE con participación, <<X2_SW_P1>> p.u. con δ = 1 %, <<X2_SW_P2>> con 2 % y
  <<X2_SW_P5>> con 5 %. Con δ = 1 % ninguna barra pasa de 10⁻⁴ p.u.; con 2 % ya son <<X2_SW_N4_P2>>.

### 18.5 Criterio 2: ancho frente al rango de referencia

{{x2_width}}

**Se cumple.** El intervalo con bloques espectrales queda entre <<X2_W_SPEC_LO>> y <<X2_W_SPEC_HI>> veces el rango AC
de los vértices en mediana, y <<X2_W_SPEC_MAX>> veces en la peor barra. Con AF2 la peor barra llega a
<<X2_W_MONO_MAX>> veces. El producto afín estándar, término a término, da <<X2_W_NAT_P_LO>> a <<X2_W_NAT_P_HI>>
veces en PEGASE con participación y <<X2_W_NAT_S_LO>> a <<X2_W_NAT_S_HI>> con slack única.

### 18.6 Criterio 3: ¿cuánto mejoran los bloques a AF2?

> **Aviso (15-11-2026).** Las cifras de RAA-K de esta sección se calcularon sin normalizar las formas
> cuadráticas. El operador usa una tolerancia absoluta y con formas de este tamaño pierde bloques, así que
> estas cifras **subestiman a RAA-K** (sección 25.8, con la medición). Las cotas siguen siendo válidas; solo son
> más anchas de lo que el operador puede dar. Falta rehacer la sección con la forma normalizada.

Ancho de la cota sobre el rango exacto de εᵀH_iε, mediana; entre paréntesis, el lado positivo espurio como
porcentaje de la profundidad. Solo barras con H_i semidefinida negativa y K ≤ 16, donde el rango exacto es
[mínimo en los vértices, 0]:

{{x2_operators}}

Sobre todas las barras, ancho de AF2 dividido por el de cada alternativa (mediana barra a barra), y tiempos:

{{x2_gain}}

**Se cumple, y el operador del grupo es el mejor cuando las fuentes son pocas.**

- **Con hasta 8 zonas gana RAA-K.** En IEEE-33 con K = 4 da el rango exacto (<<X2_EX_RKN_I4>>, contra
  <<X2_EX_SPEC_I4>> de los bloques espectrales y <<X2_EX_AF2_I4>> de AF2). Con K = 8 queda a
  <<X2_EX_RKN_K8_LO>>–<<X2_EX_RKN_K8_HI>> veces el rango exacto en los tres casos; los bloques espectrales, a
  <<X2_EX_SPEC_K8_LO>>–<<X2_EX_SPEC_K8_HI>>, y AF2, a <<X2_EX_AF2_K8_LO>>–<<X2_EX_AF2_K8_HI>>.
- **Con más zonas se invierte.** Con K = 16, RAA-K queda a <<X2_EX_RKN_K16_LO>>–<<X2_EX_RKN_K16_HI>> y los bloques
  espectrales a <<X2_EX_SPEC_K16_LO>>–<<X2_EX_SPEC_K16_HI>>. Sobre todas las barras de PEGASE, la ventaja de RAA-K
  frente a AF2 baja de <<X2_GRKN_P8_LO>>–<<X2_GRKN_P8_HI>> veces con K = 8 a <<X2_GRKN_P64_LO>>–<<X2_GRKN_P64_HI>>
  con K = 64, mientras la de los bloques espectrales se mantiene en <<X2_G_LO>>–<<X2_G_HI>>. Con muchas variables
  el presupuesto diagonal se agota antes de absorber lo que está fuera de la diagonal, y la búsqueda solo une
  variables cuya diagonal tiene el mismo signo.
- **La ganancia está en el lado superior, no en el que manda.** En formas semidefinidas negativas, el lado
  inferior de AF2 ya es casi exacto (<<X2_LOW_AF2>> veces la profundidad) y RAA-K lo conserva
  (<<X2_LOW_RKN>>). Lo que RAA-K quita es el lado positivo espurio, que con AF2 vale
  <<X2_RISE_LO>> a <<X2_RISE_HI>> % de la profundidad y con RAA-K <<X2_RISE_RKN_LO>> a <<X2_RISE_RKN_HI>> %. Los
  bloques espectrales lo eliminan, pero aflojan el lado inferior (<<X2_LOW_SPEC>> veces).
  Para la cota inferior de voltaje, que es la que importa, AF2 y RAA-K dan lo mismo; RAA-K evita inflar la
  superior.
- **El tamaño del clique importa.** Con pares (el método binomial) la ganancia es la menor: <<X2_GRK2_LO>> a
  <<X2_GRK2_HI>> veces. Pero dejar crecer el clique sin límite no siempre ayuda: en <<X2_KN_WORSE>> de los
  <<X2_NCASES>> casos, k = K queda peor que k = 6, por el orden voraz y el tope de nodos.
- **Costo.** Las cotas de todas las barras tardan hasta <<X2_TB_MAX>> s (K = 64), casi todo en la búsqueda de
  cliques.
- **Lectura para el grupo.** El teorema vale para bloques cualesquiera. Una combinación (los autovectores
  dominantes primero, RAA-K sobre el residuo) podría juntar lo mejor de las dos; no está probada.

### 18.7 Veredicto de la Puerta 2

| Criterio | Resultado | Evidencia |
|---|---|---|
| 1. La cota contiene la solución AC en todos los vértices y muestras | **No** | Falla por el término de tercer orden: hasta <<X2_MONO_SHORT_I>> p.u. en IEEE-33 y <<X2_SPEC_SHORT_MAX>> p.u. en PEGASE (18.4) |
| 2. Ancho ≤ 1,5 veces el rango de referencia, en mediana | **Sí** | <<X2_W_SPEC_LO>> a <<X2_W_SPEC_HI>> (18.5) |
| 3. El operador por bloques mejora a AF2 | **Sí; cuánto depende del número de fuentes** | RAA-K: <<X2_GRKN_I_LO>>–<<X2_GRKN_I_HI>> veces en IEEE-33; en PEGASE, de <<X2_GRKN_P8_LO>>–<<X2_GRKN_P8_HI>> con K = 8 a <<X2_GRKN_P64_LO>>–<<X2_GRKN_P64_HI>> con K = 64. Bloques espectrales: <<X2_G_LO>>–<<X2_G_HI>>. La ganancia es en el lado superior (18.6) |

Qué cambia en el plan:

1. **El punto 3.4 deja de ser opcional.** Sin control del tercer orden no se puede afirmar contención, ni siquiera
   en IEEE-33. O se hace el punto fijo r = −J⁻¹q(Cε + r), o el paper declara una cota de segundo orden con el
   déficit medido y su ley en δ.
2. **Hay dos afirmaciones que sí se sostienen hoy.** La corrección de segundo orden baja la fracción de muestras
   con alguna barra fuera de <<X2_AFF_MCS_P8>>–<<X2_AFF_MCS_S8>> % a <<X2_SPEC_MCS_MAX>> % como máximo, y cuesta
   segundos en 2.869 barras.
3. **El rango de δ importa.** Con zonas coherentes, δ = 5 % ya es un caso duro para un método de segundo orden en
   PEGASE; con δ = 1 % el déficit es de <<X2_SW_P1>> p.u.
4. **El resto de un solo signo no es general.** Vale en IEEE-33; en PEGASE con zonas, no. La salida "solo se
   corrige la cota inferior" (punto 3.2) hay que plantearla con H_i indefinida.
5. **La elección de la slack pesa menos de lo que parecía** en estas medidas: el modelo de participación baja los
   ángulos y la sobreestimación del producto estándar, pero no el déficit en las barras débiles.
"""

X3_TEXT = r"""---

## 19. X3 — ¿El modelo de segundo orden encuentra el peor escenario de cada barra?

Para cada barra se proponen los vértices de la caja de ruido donde su voltaje sería extremo: con el modelo de
primer orden (los signos de c_i) y con el de segundo orden (el vértice que minimiza o maximiza c_iᵀx + xᵀH_ix).
Se resuelve el flujo AC exacto solo en esos vértices y se compara con el extremo verdadero sobre los 2^K
vértices, que X2 calculó por enumeración. No hay números aleatorios.

{{x3}}

*Script: `experiments/x3_worst_vertex.py` · CSV: `results/x3_worst_vertex.csv` · "acierta": el voltaje alcanzado
coincide con el extremo de la enumeración dentro de 10⁻⁹ p.u. · "brecha": cuánto le falta al peor candidato*

- En IEEE-33 el peor vértice es el mismo para todas las barras (todas las zonas al mismo lado) y los dos modelos
  lo encuentran.
- **En PEGASE el segundo orden acierta el voltaje mínimo en <<X3_LO2_LO>> a <<X3_LO2_HI>> % de las barras; el
  primer orden, en <<X3_LO1_LO>> a <<X3_LO1_HI>> %.** Para el máximo: <<X3_HI2_LO>> a <<X3_HI2_HI>> % contra
  <<X3_HI1_LO>> a <<X3_HI1_HI>> %.
- Bastan <<X3_SOLVES_LO>> a <<X3_SOLVES_HI>> flujos AC para todas las barras, porque muchas comparten vértice;
  la enumeración usa 256 o 65.536.
- **No es exacto.** Donde falla, al mínimo le faltan hasta <<X3_GAP_LO2_LO>>–<<X3_GAP_LO2_HI>> p.u. según el caso
  (con el primer orden, <<X3_GAP_LO1_LO>>–<<X3_GAP_LO1_HI>>). Lo que entrega es un valor alcanzado, es decir, una
  cota interior del rango: sirve para medir cuán conservadora es una cota exterior, no para reemplazarla.
"""
