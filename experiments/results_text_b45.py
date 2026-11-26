"""Sections 27 to 29 of RESULTS.md: second-order correction (B4), outputs from the affine forms (B5), claims."""

B45_TEXT = r"""
---

## 27. B4 — Corrección de segundo orden de las cotas CD-HIM

**Por qué.** En B1 todo el exceso sobre el 0,27 % nominal viene del error de linealización: las mismas cotas
probadas contra el modelo lineal no lo superan. El modelo de segundo orden de |V| (sección 17, `SecondOrder`) da
ese error de forma explícita, d|V_i| = c_iᵀε + εᵀH_iε, sin formar H_i.

**Qué se agrega** (`cdhim/correction.py`). Con símbolos independientes, de media cero y simétricos, la parte
lineal y la cuadrática no están correlacionadas, así que:

- **media:** E[q_i] = Var(ε)·tr(H_i), exacta y barata;
- **varianza y covarianza de rama** de q, por muestreo del término cuadrático solo (sin flujo de potencia);
- **rango de q**, por el mismo muestreo; para una elipse, proyectado sobre sus ejes principales, y la elipse se
  reemplaza por una elipse exterior de la suma de Minkowski de la elipse y ese rango;
- **para la caja afín determinista:** q evaluado en los dos vértices donde la parte lineal es extrema.

La columna "sustituto de 2.º orden" prueba las mismas cotas contra cᵀε + q(ε) en vez de AC: mide lo que dejan
fuera el tercer orden y los superiores.

### 27.1 Punto de diseño (mismas muestras AC que B1)

{{b4_design}}

*Script: `experiments/b4_second_order.py` · CSV: `results/b4_design_bus.csv`, `results/b4_design_branch.csv` ·
4.000 muestras del término cuadrático (semilla de B1 + 1)*

- **El intervalo CD-HIM por barra queda resuelto:** en PEGASE la peor barra pasa de <<B4_LMAX_0>> % a
  <<B4_LMAX_R>> %, con un ancho <<B4_LMAX_RW>> % mayor en mediana.
- **La elipse necesita el rango, no los momentos.** Con media y covarianza de q la peor elipse baja de
  <<B4_ELL_0>> % (<<B4_ELL_0N>> ramas sobre 0,27 %) a <<B4_ELL_M>> % (<<B4_ELL_MN>> ramas). A lo largo del eje menor de
  las elipses delgadas el término cuadrático es casi de un solo signo y de cola pesada (asimetría cercana a ±2,7 y
  curtosis cercana a 13 en las ramas peores; diagnóstico en esta sesión, sin CSV), lo que un ajuste de dos momentos no representa. Con el rango
  proyectado la peor elipse queda en <<B4_ELL_K>> %, con un área <<B4_ELL_KA>> % mayor en mediana y hasta
  <<B4_ELL_KAM>> veces en la peor rama.
- **El tercer orden no pesa en el punto de diseño:** el sustituto de segundo orden reproduce las violaciones AC.
- **La cota rigurosa de q término a término no sirve:** es válida para el modelo de segundo orden pero multiplica
  el ancho de la caja afín por más de cien en mediana en PEGASE.

### 27.2 Barrido carga × δ (las 50 celdas de C1, mismas muestras)

{{b4_sweep}}

*Script: `experiments/b4_second_order.py` · CSV: `results/b4_sweep.csv`. El ancho de la elipse no está en el
barrido.*

- Con el rango de q, ninguna de las <<B4_SW_CELLS>> celdas tiene una barra o una elipse sobre el nominal. Con solo
  media y varianza siguen fallando <<B4_SW_MOM_BAD>> de las 25 celdas de PEGASE.
- La caja 3σ corregida sigue fallando: lo que falla ahí es la normalidad, no la linealización.
- "Sobre 0,27 %" compara una frecuencia estimada con 20.000 muestras con el nivel nominal: una cota exacta al
  nivel nominal supera el umbral por azar en una parte de los objetivos. La elipse corregida queda bajo el umbral
  porque es algo más ancha.

### 27.3 Costo

{{b4_timing}}

*Script: `experiments/b4_second_order.py` (modo `timing`) · CSV: `results/b4_timing.csv` · un proceso*

En PEGASE la corrección completa suma <<B4_T_COR>> s a los <<B4_T_LIN>> s del método lineal. La sección 28.5 muestra
que bastan muchas menos muestras.

---

## 28. B5 — Salidas desde las formas afines, sin elipse

**Por qué.** Las formas afines ΔV_i = c_iᵀε conservan la dependencia entre barras a través de los símbolos
compartidos. Cualquier combinación lineal, por ejemplo la caída |V_i| − |V_j| de una rama, tiene su propia forma
afín y su propio modelo de segundo orden (`cdhim/affine_forms.py`, sin formar H). La elipse, G y γ_th dejan de
ser necesarios.

**Radios.** Para la parte lineal: la caja afín Σ|c_k| (peor caso), 3σ (supone normalidad) y el radio de
Bernstein para una suma de términos independientes, acotados y de media cero, r = ML/3 + √((ML/3)² + 2σ²L) con
L = ln(2/α), M = max|c_k|, topado por la caja afín. A todos se les suma el rango de q como en B4.

**Peor caso.** Además del vértice de primer orden, se busca el vértice extremo del modelo de segundo orden por
volteos de signo guiados por la ganancia exacta de cada volteo. En una prueba con 16 zonas en PEGASE, frente a la
enumeración de los 65.536 vértices del modelo, la búsqueda da el extremo exacto en 86 % (máximo) y 97 % (mínimo)
de las barras, contra 58 % y 78 % del vértice de primer orden (medición en esta sesión, sin CSV; la prueba
unitaria `tests/test_affine_forms.py` lo exige exacto en IEEE-33 con 8 zonas). No es una cota: entrega un punto
alcanzado.

### 28.1 Punto de diseño, Monte Carlo AC (mismas muestras que B1)

{{b5_design}}

*Script: `experiments/b5_affine_forms.py` · CSV: `results/b5_mc.csv`, `results/b5_width.csv`*

- **Barras.** Bernstein con el rango de q tiene su peor barra en <<B5_BERN_BUS>> % en PEGASE, con un ancho de
  <<B5_BERN_AFF>> veces la caja afín y <<B5_BERN_SIG>> veces 3σ en mediana. El intervalo CD-HIM corregido mide
  <<B5_CDHIM_BERN>> veces el de Bernstein en mediana, sin una justificación de cobertura propia.
- **Caídas por rama.** La forma afín de la caída con Bernstein no supera el nominal (peor rama <<B5_BERN_DROP>> %).
  Restar los intervalos de las dos barras da caídas <<B5_SUB_CDHIM>> veces más anchas en mediana con CD-HIM,
  <<B5_SUB_BERN>> con Bernstein y <<B5_SUB_DET>> con las cajas afines (PEGASE; en IEEE-33, <<B5_SUB_CDHIM_33>> veces con
  CD-HIM).

### 28.2 El peor caso en los vértices: la caja de segundo orden no encierra la solución AC

{{b5_vertex}}

*Script: `experiments/b5_affine_forms.py` · CSV: `results/b5_vertex.csv` · flujo AC exacto en el vértice propuesto
para cada barra o rama*

- **PEGASE-2869:** con todas las fuentes en un extremo, hasta <<B5_VX_BUS>> de <<B5_VX_NBUS>> barras y <<B5_VX_DROP>> de
  <<B5_VX_NDROP>> caídas quedan fuera incluso de la cota que toma el extremo de segundo orden hallado por búsqueda; el
  exceso llega a <<B5_VX_EXC>> p.u., el <<B5_VX_REL>> % del ancho de la cota. Pesan los órdenes superiores al segundo.
- **IEEE-33:** también quedan fuera, por <<B5_VX_33>> p.u. como máximo (el tercer orden de la sección 18).
- El Monte Carlo no lo detecta: las muestras no llegan a los vértices.
- **Consecuencia:** sin una cota rigurosa del resto, que la sección 22 mostró que no escala a transmisión, no se
  puede afirmar una cota determinista a esta escala. La caja afín queda como referencia del peor caso del modelo
  lineal, no como garantía.

### 28.3 Barrido carga × δ (las 50 celdas de C1, mismas muestras)

{{b5_sweep}}

*Script: `experiments/b5_affine_forms.py` (modo `sweep`) · CSV: `results/b5_sweep.csv` · sin búsqueda de vértices*

En las <<B5_SW_CELLS>> celdas, Bernstein con el rango de q no supera el nominal en ninguna barra (peor
<<B5_SW_BERN_BUS>> %) ni en ninguna caída (peor <<B5_SW_BERN_DROP>> %). La versión 3σ corregida falla en <<B5_SW_SIG_BAD>>
celdas.

### 28.4 Qué radio: lo que aporta suponer cada cosa

{{b5_radii}}

*Script: `experiments/b5_affine_forms.py` (modo `radii`) · CSV: `results/b5_radii.csv` · 1.000 muestras de q*

- **Los radios que usan solo la varianza no sirven:** Chebyshev pide <<B5_R_CHEB>>σ en PEGASE y queda topado por la
  caja afín en el <<B5_R_CHEB_CAP>> % de las barras; la desigualdad de Gauss para distribuciones unimodales, en el
  <<B5_R_UNI_CAP>> %. Equivalen al peor caso.
- **Los que usan el soporte acotado sí:** Hoeffding (<<B5_R_HOEF>>σ) y Bernstein (<<B5_R_BERN>>σ) quedan en la mitad de la
  caja afín. La varianza que agrega Bernstein gana solo <<B5_R_BH>> % sobre Hoeffding.
- Todos necesitan la corrección de segundo orden en las barras con una sola fuente dominante, donde el radio queda
  topado por la caja afín y el error de linealización lo cruza.

### 28.5 Cuántas muestras del término cuadrático hacen falta, y el costo del método

{{b5_surrogate}}

*Script: `experiments/b5_affine_forms.py` (modo `surrogate`) · CSV: `results/b5_surrogate.csv` · tres semillas por
tamaño · un proceso*

Con 250 muestras la cobertura no cambia: el margen del radio absorbe el error de estimar la cola del rango. El
método completo en PEGASE (estudio lineal, modelo de segundo orden y 250 muestras) toma <<B5_T_METHOD>> s, unas
<<B5_SPEEDUP>> veces menos que el Monte Carlo AC por lotes de C8 (<<B5_T_MC>> s). El barrido 28.3 usó 4.000 muestras;
las 250 están probadas solo en el punto de diseño. El muestreo necesita doble precisión: q es una suma de términos
grandes que casi se cancelan.

### 28.6 Fuentes correlacionadas (B6)

**Cómo se modela.** A la manera de la aritmética afín, con símbolos compartidos: cada fuente es
ε_k = w·ζ_z(k) + (1 − w)·η_k, con un factor ζ por zona (o uno global) y un ruido propio η_k, todos independientes
y uniformes en [−1, 1] (`FactorModel`). El soporte de ε_k sigue en [−1, 1]; la correlación entre dos fuentes de
la misma zona es w²/(w² + (1 − w)²), y su varianza baja (1/6 con ρ = 0,5). La forma afín reescrita sobre los
símbolos independientes admite Bernstein sin cambios, y el rango de q se muestrea con la misma ley. La cópula
gaussiana de C2 no se puede escribir así: para ella se prueba un modelo de un factor global con la misma
correlación de Pearson. Balance con la slack sola o con los generadores en partes fijas (`participation_factors`).
Punto de diseño, 20.000 muestras AC por caso, 1.000 muestras de q.

**Barras**

{{b6_bus}}

**Caídas por rama**

{{b6_drop}}

*Script: `experiments/b6_correlated.py` · CSV: `results/b6_correlated.csv` · semilla por caso en la columna `seed`*

- **Suponer independencia con fuentes correlacionadas no sirve:** en PEGASE la peor barra queda fuera entre el
  <<B6_IND_MIN>> % y el <<B6_IND_MAX>> % de las muestras.
- **Con la correlación escrita como símbolos compartidos el método cubre** (peor barra <<B6_FAC_MAX>> % en PEGASE), pero
  la concentración se pierde: el ancho pasa de <<B6_IND_W>> veces la caja afín con fuentes independientes a entre
  <<B6_FAC_WMIN>> y <<B6_FAC_WMAX>> veces. Un factor común no se promedia. Con un factor global en IEEE-33 el radio queda
  topado por la caja afín.
- **Un modelo de factores mal especificado falla:** el de un factor con la misma correlación que la cópula deja fuera
  hasta el <<B6_COP_MATCH>> % de las muestras en una barra de PEGASE. Igualar la correlación no basta: la cópula tiene otra
  dependencia en las colas. 3σ con la Σ_ε verdadera llega a <<B6_COP_S3>> %.
- **El peor caso lineal es lo único robusto a la especificación**, con una condición: con la slack sola y la cópula
  fuerte, en PEGASE incluso la caja afín con el rango de q queda fuera hasta el <<B6_COP_AFF_SLACK>> % (el error de
  linealización de C2, que el segundo orden no alcanza a cubrir); con los generadores repartiendo el desbalance,
  <<B6_COP_AFF_PART>> %.
- **Consecuencia para el paper:** la correlación hay que modelarla como estructura de símbolos compartidos y
  declararla como supuesto; el aporte probabilístico es grande con fuentes independientes o débilmente
  correlacionadas y se reduce con un factor común fuerte. Con correlación fuerte, el balance por la slack sola no es
  un supuesto razonable a escala de transmisión.

### 28.7 Cobertura simultánea y dos sistemas más (B7)

**Montaje.** Con una cota de unión (Bonferroni), cada objetivo recibe el nivel α/m, con m el número de barras (o de
caídas) con incertidumbre; "barras y caídas" usa m = barras + caídas. Se mide la fracción de muestras AC con
**algún** objetivo fuera, que es la cobertura del sistema. Radios: Bernstein, Gauss (cuantil normal), Chebyshev,
todos topados por la caja afín y con el rango de q (1.000 muestras). Se agregan IEEE-118 y PEGASE-1354, ambos con
barras PV, en una grilla de carga (0,8; 1,0; 1,2; 95 % del límite de cargabilidad) × δ (5, 15 y 30 %), con 20.000
muestras por celda.

**Punto de diseño de cada sistema** (IEEE-118 y PEGASE-1354: carga 1,0 y δ = 5 %)

{{b7_design}}

**Grilla de IEEE-118 y PEGASE-1354**

{{b7_grid}}

*Script: `experiments/b7_generality.py` · CSV: `results/b7_generality.csv` · las muestras AC sin converger se excluyen,
como en C1*

- **La cobertura simultánea con un radio sin distribución cuesta lo mismo que el peor caso.** En PEGASE-2869 el radio
  de Bernstein pasa de <<B7_BERN_PT>> veces la caja afín por barra a <<B7_BERN_BF>> con Bonferroni. Lo mismo ocurre en
  los cuatro sistemas. El crecimiento en √ln(m/α) es lento, pero el radio por barra ya estaba cerca del tope.
- **Gauss con Bonferroni cubre en la práctica** con <<B7_GAUSS_BF>> veces la caja afín (muestras con alguna barra fuera:
  <<B7_GAUSS_BF_ANY>> %), pero sin garantía: por barra, Gauss deja <<B7_GAUSS_PT_ANY>> % de las muestras con alguna
  barra fuera, y en PEGASE-1354 falla por barra en <<B7_GAUSS_BAD1354>> de las 12 celdas.
- **Bernstein por barra ya da una cobertura del sistema cercana al nominal, sin estar diseñado para eso:** muestras con
  alguna barra o caída fuera, <<B7_BERN_PT_ANY>> % en PEGASE-2869 y como máximo <<B7_BERN_PT_ANY_GRID>> % en la grilla.
  Es una observación empírica, no una garantía.
- **La ventaja de la concentración depende de cuántas fuentes comparables mueven cada barra.** Bernstein por barra mide
  <<B7_W1354>> veces la caja afín en PEGASE-1354 y <<B7_W118>> en IEEE-118, donde queda topado en el <<B7_CAP118>> % de las
  barras: allí pocas cargas mueven cada barra y no hay qué promediar.
- **Generalidad de la corrección:** en los dos sistemas nuevos Bernstein con el rango de q no supera el nominal en
  ninguna celda. Cerca del límite de cargabilidad con δ alto no converge hasta el <<B7_NC_MAX>> % de las muestras AC; esas
  celdas hay que leerlas con la advertencia de C1.

### 28.8 Frente a métodos publicados de flujo probabilístico (B8)

**Métodos.** (a) El método de cumulantes sobre el modelo lineal, con la expansión de Cornish-Fisher de cuarto orden
para los cuantiles α/2 y 1 − α/2. (b) La estimación por puntos de Hong (2m+1) sobre el flujo AC exacto: dos flujos por
fuente, en ε_k = ±0,775. Los cumulantes del aporte de cada fuente se suman, como corresponde a fuentes independientes;
la suma de momentos crudos del esquema básico pierde los términos cruzados del cuarto momento cuando hay miles de
fuentes. Ambos con Cornish-Fisher. Mismas celdas para los cuatro sistemas: diseño, carga 1,2 con δ = 30 %, y 95 % del
límite de cargabilidad con δ = 15 %. 20.000 muestras AC por celda.

{{b8}}

*Script: `experiments/b8_ppf_baselines.py` · CSV: `results/b8_ppf_baselines.csv`. Las caídas por rama están en el CSV.*

- **Los dos métodos publicados no cubren.** La peor barra queda fuera hasta el <<B8_CUM_MAX>> % (cumulantes) y el
  <<B8_PEM_MAX>> % (estimación por puntos) de las muestras; en el punto de diseño de PEGASE-2869, <<B8_CUM_N>> y <<B8_PEM_N>>
  de sus <<B8_NB>> barras superan el 0,27 %. Sus intervalos miden alrededor de <<B8_W>> veces los de Bernstein.
- La estimación por puntos usa el flujo AC exacto, pero fuente por fuente: no ve las interacciones entre fuentes,
  que son las que el término cuadrático recoge. Cuesta <<B8_PEM_SOLVES>> flujos AC (<<B8_PEM_T>> s) en PEGASE-2869.
- Bernstein con el rango de q no supera el nominal en ninguna barra de ninguna celda (peor <<B8_BERN_MAX>> %).

### 28.9 Corrientes de rama (B9)

**Por qué.** En transmisión lo que se vigila en una rama es su carga, es decir, la corriente. La corriente en cada
extremo es I = Y_f V (o Y_t V). Su forma afín sale del voltaje complejo de primer orden, d|I| = Re(conj(I₀)/|I₀| ·
Y_f E) ε. El resto no lineal se muestrea sobre el estado complejo de segundo orden ΔV₁ + ΔV₂, con
ΔV₂ = −J⁻¹q(ΔV₁) (`SecondOrder.state`), y su rango se suma como en las barras. Ambos extremos de todas las ramas en
servicio con |I₀| > 10⁻⁶ p.u.

{{b9}}

*Script: `experiments/b9_branch_currents.py` · CSV: `results/b9_branch_currents.csv`*

- **Aquí la corrección es imprescindible.** El módulo |I| es muy no lineal cerca de corriente cero y en ramas con
  cambios grandes: el error del modelo lineal llega a <<B9_ERRLIN>> p.u. Sin el rango del término no lineal, Bernstein
  deja fuera hasta el <<B9_LIN_MAX>> % de las muestras; cumulantes y estimación por puntos, hasta el <<B9_CUM_MAX>> %.
- **Con el rango, Bernstein cubre en todas las celdas** (peor extremo <<B9_BERN_MAX>> %).
- En las corrientes la concentración sí rinde frente al peor caso: en PEGASE-2869 la caja afín corregida mide
  <<B9_AFF_W>> veces el intervalo de Bernstein en mediana.

### 28.10 Distribuciones no uniformes (B10)

Las doce distribuciones de C3, en el punto de diseño de tres sistemas. Hoeffding supone solo media cero y soporte en
[−1, 1]; Bernstein además una cota de la varianza (1/3 por defecto, o la verdadera); 3σ, normalidad. El rango de q se
muestrea con la distribución verdadera.

{{b10}}

*Script: `experiments/b10_marginals.py` · CSV: `results/b10_marginals.csv` · barras; las caídas están en el CSV*

- **Con media cero y soporte en [−1, 1], los tres radios sin distribución cubren** (peor barra <<B10_IN_MAX>> %), incluida
  la distribución en forma de U, cuya varianza 1/2 supera el 1/3 que supone Bernstein por defecto (<<B10_U_MAX>> %).
  3σ falla con casi todas (hasta <<B10_S3_MAX>> %).
- **Fuera del soporte, nada cubre**, ni siquiera la caja afín (hasta <<B10_OUT_MAX>> % en PEGASE-2869). El soporte
  acotado es el supuesto que sostiene al método: δ tiene que ser el semiancho real de cada inyección, no una
  desviación típica.

### 28.11 El barrido con 250 muestras del término cuadrático

{{b5_n250}}

*Script: `experiments/b5_affine_forms.py` (modo `sweep250`) · CSV: `results/b5_sweep_n250.csv`, `results/b5_sweep.csv`*

Con 250 muestras ninguna celda supera el nominal, pero la peor barra de PEGASE sube a <<B5_N250_MAX>> %, cerca del límite.
La cola del rango queda subestimada en las celdas más exigentes. Se recomienda 1.000 muestras, que cuestan alrededor
de un segundo en PEGASE-2869 (28.5).

---

## 29. Lista de afirmaciones para el paper de formas afines

Cada afirmación con su respaldo y su límite. Sin fila de respaldo, no entra (CHECKLIST_AA.md, 6.1).

| # | Afirmación | Respaldo | Límite que hay que declarar |
|---|---|---|---|
| 1 | Las formas afines de |V| se obtienen de una factorización LU dispersa en el punto AC y escalan a 9.241 barras | C7, `results/c7_timing.csv` | El bloque de sensibilidades es denso: el método completo escala como N^1,5, no N^1,2 |
| 2 | El error de linealización es lo único que separa las cotas del 0,27 % nominal en el punto de diseño | B1 y 27.1 (columnas lineal y AC) | En la peor barra de PEGASE |
| 3 | El rango del término cuadrático corrige ese error en todo el rango carga × δ; un ajuste de media y varianza no basta | 27.2 y 28.3 | Validado por Monte Carlo (20.000 muestras por celda), en dos sistemas; el rango es muestreado, no acotado |
| 4 | Las caídas por rama desde la forma afín de la diferencia son varias veces más angostas que restar intervalos de barra, con la misma cobertura | 28.1 y 28.3 | Caídas de magnitud de voltaje |
| 4b | Las corrientes de rama se acotan con la misma construcción; ahí el rango del término no lineal es imprescindible | 28.9 | Validado en las tres celdas de B8 de cada sistema |
| 4c | Los métodos publicados de flujo probabilístico (cumulantes, estimación por puntos) no cubren a este nivel, ni en barras ni en corrientes | 28.8, 28.9 | Con Cornish-Fisher; otras expansiones (Gram-Charlier) no se probaron |
| 4d | El método cubre con cualquier distribución de media cero y soporte en [−1, 1]; fuera de ese soporte no | 28.10 | δ debe ser el semiancho real de las inyecciones |
| 5 | Un radio que usa el soporte acotado (Bernstein o Hoeffding) cubre sin suponer normalidad; los que usan solo la varianza equivalen al peor caso | 28.4 | Supone símbolos independientes con soporte [−1, 1] |
| 5b | Fuentes correlacionadas: escritas como símbolos compartidos, el mismo método cubre; suponer independencia o un modelo mal especificado no | 28.6 | La ganancia sobre el peor caso cae con un factor común fuerte; la estructura de correlación es un supuesto del modelo; con la slack sola y correlación fuerte falla la linealización en PEGASE |
| 6 | 3σ, la reformulación gaussiana habitual, no cubre por objetivo aun corregida | 27.2, 28.3, 28.7 | Con Bonferroni sí cubre en la práctica, sin garantía |
| 6b | La cobertura simultánea con un radio sin distribución cuesta lo mismo que el peor caso; la ventaja de la concentración es por objetivo y depende del número de fuentes por barra | 28.7 | En IEEE-118 la concentración no aporta |
| 7 | La cota de peor caso de segundo orden no encierra la solución AC en los vértices de una red de transmisión | 28.2 | Es un resultado negativo; en IEEE-33 el exceso es despreciable |
| 8 | El método completo es unas 15 veces más rápido que un Monte Carlo AC por lotes en 2.869 barras | 28.5 y C8 | Depende de la implementación del Monte Carlo; con 1.000 muestras de q, en vez de 250, algo menos (28.11) |

Lo que este informe no respalda todavía: una garantía simultánea más angosta que el peor caso, correlación estimada de datos
reales de carga y generación, la corrección en otros sistemas con el barrido completo de C1 (B7 usa una grilla de 12 celdas), y una revisión
bibliográfica en bases indexadas (RELATED_WORK.md, quinta pasada).
"""
