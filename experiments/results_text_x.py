"""Commentary of section 17 of RESULTS.md (exploratory experiment X1). Tables and <<numbers>> come from the CSVs."""

X_TEXT = r"""---

## 17. X1 — Exploración: qué puede aportar la aritmética afín

No lo pidieron los revisores. Responde a la pregunta de los autores de si el aporte se puede plantear desde la
aritmética afín. Es un estudio de factibilidad, con tres límites que conviene tener presentes al leerlo:

- Solo se acota el término de segundo orden del resto. Un encierro verificado necesita además el punto fijo
  r = −J⁻¹q(Cε + r) y redondeo dirigido; nada de eso está hecho.
- Las cotas se comparan contra valores *alcanzados* de la forma cuadrática (búsqueda local), no contra su rango
  exacto. Los cocientes "ancho / alcanzado" son cotas superiores de la sobreestimación.
- La cota con bloques binomiales y presupuesto diagonal es una lectura simplificada mía de la idea del grupo
  (una pasada voraz por |H_kl| decreciente), no la implementación de los autores.

**Planteamiento.** Lo que CD-HIM calcula ya es una forma afín por barra, V̂_i = V_i* + Σ_k c_ik ε_k: la matriz C
son las desviaciones parciales, la "caja afín" Σ_k |c_ik| es su rango, y la correlación entre barras vecinas es
la dependencia que llevan los símbolos compartidos. En coordenadas rectangulares las ecuaciones de flujo de
potencia son exactamente cuadráticas, así que el resto de esa forma afín no tiene serie truncada:

    d|V_i| = c_iᵀε + εᵀH_iε + O(δ³)
    εᵀH_iε = −Σ_m Re(γ_im · dV_m · conj(dI_m)) + (c⊥_iᵀε)² / (2|V_i|)      (más un término en |dV_m|² en barras PV)

dV_m y dI_m son las formas afines de la desviación de voltaje y de corriente inyectada en la barra m, y γ_im sale
de una resolución traspuesta con la misma factorización LU. El resto es una suma de N productos de dos formas
afines: una forma cuadrática en los símbolos de ruido, que es el objeto que acotan los operadores bilineales de
AF2 y RAA, la reducción por presupuesto diagonal y el teorema k-nomial.

### 17.1 El error de linealización es esa forma cuadrática

{{x1_form}}

*Script: `experiments/x1_aa_second_order.py` · CSV: `results/x1_summary.csv`, `results/x1_bus_remainder.csv` ·
N = 20.000 flujos AC por sistema · δ de diseño (15 % y 5 %)*

- Sumar εᵀH_iε a la forma afín baja el error contra el flujo AC en un factor mediano de <<X_RED_33>> (IEEE-33) y
  <<X_RED_PEG>> (PEGASE). El error de linealización medido en C1 es, casi entero, este término.
- **El resto tiene un solo signo.** H_i es semidefinida negativa en las 32 barras de IEEE-33 y lo es casi por
  completo en PEGASE (parte negativa del espectro: <<X_NEG_PEG>> % en mediana). La forma afín sobreestima el
  voltaje: es una envolvente superior, y solo la cota inferior necesita corrección.
- Pocas direcciones concentran el espectro (<<X_EIG90_33>> y <<X_EIG90_PEG>> autovalores para el 90 %, en mediana),
  con N_u = 64 y 2.951.

### 17.2 La caja afín no encierra la solución AC

Para cada barra se resuelve el flujo AC exacto en los dos vértices de la caja de ruido donde su forma afín es
extrema: ε = −sign(c_i) (todas las fuentes en el extremo que baja el voltaje de esa barra) y ε = +sign(c_i).
No hay muestreo ni semilla.

{{x1_vertex}}

*CSV: `results/x1_bus_remainder.csv` (columnas `vertex_*`)*

- **En el vértice inferior, la solución AC queda por debajo de la caja afín en <<X_VX_N_33>> de <<X_N_ACT_33>>
  barras de IEEE-33 y en <<X_VX_N_PEG>> de <<X_N_ACT_PEG>> de PEGASE.** El exceso es chico en IEEE-33
  (<<X_VX_MED_33>> % del radio en mediana, <<X_VX_MAX1_33>> % en la peor barra) y grande en PEGASE:
  <<X_VX_MED_PEG>> % en mediana, <<X_VX_P90_PEG>> % en el percentil 90 y <<X_VX_MAX_PEG>> % en la peor barra.
  <<X_VX_10_PEG>> barras exceden el radio en más de 10 % y <<X_VX_100_PEG>> en más de 100 %.
- El Monte Carlo de B1 y C1 no lo ve: su mayor error es <<X_ERR_AFF_MAX_PEG>> % del radio afín. Con 2.951 símbolos
  independientes los vértices son inalcanzables por muestreo. La caja afín "pasa" el Monte Carlo porque es entre
  tres y cuatro veces más ancha que el rango al 99,73 %, no porque sea una cota.
- La forma cuadrática predice ese exceso: en mediana deja sin explicar el <<X_VX_MISS_33>> % en IEEE-33 y el
  <<X_VX_MISS_PEG>> % en PEGASE.
- En el vértice superior la solución AC queda dentro en todas las barras de IEEE-33 y en todas menos
  <<X_VXH_N_PEG>> de PEGASE, como corresponde a un resto negativo.
- **Cuidado con la lectura física del vértice en PEGASE.** Es el escenario en que todas las cargas se mueven a la
  vez hacia el mismo lado y la barra *slack* absorbe todo el desbalance. El radio afín del ángulo es
  <<X_ANG_PEG>> rad en mediana, contra <<X_MAG_PEG>> p.u. de la magnitud. El caso es válido dentro de la caja que
  el propio modelo declara, pero es extremo; con *slack* distribuida sería bastante menor. No lo probé.

### 17.3 Cuánto cuesta acotar el resto

Sobre todas las barras, sin ensamblar H_i (costo N·N_u por barra):

{{x1_width}}

Sobre un subconjunto con H_i ensamblada (las 32 barras de IEEE-33; en PEGASE, 8 con el mayor error de
linealización relativo y 16 al azar):

{{x1_bounds}}

*CSV: `results/x1_quadratic_bounds.csv` · "ancho / alcanzado": ancho de la cota sobre el ancho del rango alcanzado
por búsqueda local · "lado positivo / profundidad": cota superior sobre el mayor valor negativo alcanzado; como
el resto es casi siempre negativo, lo que aparezca ahí es sobreestimación pura*

- **La evaluación afín estándar, producto a producto, no sirve a escala.** En IEEE-33 queda a <<X_NAT1_33>> veces
  el rango alcanzado. En PEGASE, a <<X_NAT_PEG_K>> veces en mediana: la cota del resto es <<X_NATW_PEG_K>> veces el
  radio afín. Escribir cada producto como diferencia de cuadrados de binomios escalados (dos bloques k-nomiales
  por producto, nunca peor que el producto estándar) casi no cambia eso.
- **La causa son las cancelaciones entre términos, no cada producto.** Con ángulos que se mueven
  <<X_ANG_PEG>> rad, la rotación de los fasores aparece en coordenadas rectangulares como términos de segundo
  orden grandes (cuadrados en barras PV y curvatura de |V|: <<X_SH_PV_PEG>> % y <<X_SH_CURV_PEG>> % de la cota
  estándar, en mediana) que se anulan entre sí al sumar. Acotar término a término pierde esa cancelación. En IEEE-33, sin barras PV y con ángulos
  de <<X_ANG_33>> rad, el efecto no existe.
- **Con H_i ensamblada el problema es tratable.** Monomio a monomio (AF2): <<X_MONO_33>>× y <<X_MONO_PEG>>× el
  rango alcanzado. Bloques binomiales con presupuesto diagonal: <<X_K2_33>>× y <<X_K2_PEG>>×; ganan poco porque
  cada bloque solo puede absorber lo que le alcance la diagonal y aquí domina lo que está fuera de ella.
- **El teorema k-nomial con bloques grandes es lo que mejora.** Un solo bloque espectral más el residuo:
  <<X_KN1_33>>× y <<X_KN1_PEG>>×. Ocho bloques: <<X_KN8_33>>× y <<X_KN8_PEG>>×. Y como el cuadrado de un k-nomio
  tiene rango de un solo lado, desaparece el lado positivo espurio: de <<X_RISE_MONO_PEG>> % de la profundidad
  con AF2 a <<X_RISE_SPEC_PEG>> % con todos los bloques, en PEGASE.
- **Lo que no está resuelto es el costo.** Ensamblar H_i es N_u² por barra (70 MB en PEGASE) y los bloques
  espectrales necesitan autovectores. Para todas las barras de un sistema grande hace falta elegir los bloques
  sin formar H_i, o cambiar la descomposición para que las rotaciones no aparezcan (diferencias por rama, o
  coordenadas polares con resto acotado). Ese es el problema abierto.
- **Tamaño del resto en el peor caso.** En IEEE-33 la profundidad alcanzada es <<X_DEPTH_33>> % del radio afín. En
  PEGASE, entre las 8 barras de mayor error, <<X_DEPTH_WORST_MED>> veces el radio en mediana y
  <<X_DEPTH_WORST_MAX>> en la peor; entre las 16 al azar, <<X_DEPTH_RND_MED>> veces en mediana y
  <<X_DEPTH_RND_MAX>> en la mayor.

### 17.4 Capa probabilística sobre la forma afín

{{x1_prob}}

{{x1_radii}}

*CSV: `results/x1_summary.csv` · mismas 20.000 muestras AC que 17.1*

- **Correr el centro en la media del término cuadrático, tr(H_i)/3, ayuda poco.** La peor barra de la caja 3σ pasa
  de <<X_S3_MAX_PEG>> % a <<X_S3S_MAX_PEG>> % en PEGASE y el nivel sistema casi no se mueve (<<X_S3_SYS_PEG>> % y
  <<X_S3S_SYS_PEG>> %). El segundo orden no arregla la caja 3σ: su problema a nivel sistema es que el 99,73 % es
  por barra.
- **Radio de Hoeffding**, ‖c_i‖₂·√(2 ln(2/α)): vale para cualquier distribución de símbolos independientes, de
  media cero y soporte en [−1, 1], sin teorema central del límite. Es <<X_HOEF_3S>> veces el radio 3σ y
  <<X_HOEF_33>> (IEEE-33) y <<X_HOEF_PEG>> (PEGASE) veces el radio afín. Contra AC: <<X_HOEF_SYS_PEG>> % de las
  muestras de PEGASE con alguna barra fuera, frente a <<X_S3_SYS_PEG>> % de la caja 3σ. No cubre símbolos
  correlacionados ni soportes no acotados.
- Regla M = 6 de la Ref. [10]: <<X_CLTM6_PEG>> barras de PEGASE tienen menos de 6 términos efectivos (ninguna en
  IEEE-33). En esas barras la Ref. [10] no usaría el intervalo 3σ.

### 17.5 Dependencia por rama: el equivalente afín de la ganancia G

{{x1_dep}}

*CSV: `results/x1_branch_dependency.csv`*

- El conjunto exacto de la pareja (V̂_i, V̂_j) es un zonotopo, no una caja ni una elipse. Ocupa <<X_ZON_33>> %
  (IEEE-33) y <<X_ZON_PEG>> % (PEGASE) menos área que la caja afín, en mediana. Es la versión determinista de lo
  que G quería medir, depende de todos los generadores y no solo de ρ, y no necesita Σ_ε ni independencia: la
  objeción de las fuentes correlacionadas (C2) no le aplica.
- A cambio es más grande que la elipse al 99,73 %: <<X_ZON_ELL_33>> y <<X_ZON_ELL_PEG>> veces en mediana (más chico
  en <<X_ZON_LT_PEG>> de <<X_NBR_PEG>> ramas de PEGASE).
- **Donde la dependencia rinde es en cantidades de rama.** El radio de V_i − V_j calculado con la forma afín es
  <<X_DROP_33>> % y <<X_DROP_PEG>> % menor, en mediana, que restando los dos intervalos. Los intervalos por barra
  no pueden mejorar con ningún conjunto por pares (B2); las caídas de voltaje y los flujos sí.

### 17.6 Lectura

1. Hay una orientación desde la aritmética afín que se sostiene con datos: el resto de la forma afín del flujo
   de potencia es una forma cuadrática en los símbolos de ruido, de un solo signo y dominada por pocas
   direcciones, y acotarla bien es justo el problema del teorema k-nomial.
2. La motivación es medible: la caja afín, que solo es determinista para el modelo lineal, no encierra la
   solución AC (17.2), y la aritmética afín estándar aplicada al resto sobreestima en tres órdenes de magnitud
   en un sistema de transmisión (17.3).
3. La ganancia medida del k-nomial sobre AF2, barra a barra, es un factor mediano de <<X_GAIN8_33>> (IEEE-33) y
   <<X_GAIN8_PEG>> (PEGASE) en el ancho con ocho bloques, y de <<X_GAIN1_33>> y <<X_GAIN1_PEG>> con uno, más la
   eliminación del lado espurio. No es un orden de magnitud. El salto grande está entre acotar término a término
   y acotar la forma ensamblada, y ese salto hoy cuesta N_u² por barra.
4. La elipse y G no entran en esta orientación. Lo que se conserva de CD-HIM es el motor de sensibilidades
   dispersas y la observación de que las barras vecinas comparten casi todos sus símbolos.
"""
