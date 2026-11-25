"""Commentary of the C experiments (inserted into RESULTS.md by build_report.py)."""

C_TEXT = r"""## 6. C1 — Validación contra AC no lineal: carga × incertidumbre

Carga: P_d, Q_d y la generación programada se multiplican por el factor de carga; la barra slack cubre las
pérdidas. Límite de cargabilidad por continuación: **3,62 (IEEE-33) y 1,59 (PEGASE-2869)**; el quinto nivel
es el 95 % de ese límite. N = 20.000 flujos AC por celda, las mismas muestras de ruido en todas las celdas
de un sistema. En cada celda: (a) error de linealización |V_AC − V_lineal| por barra, en p.u. y relativo al
radio del intervalo de salida; (b) miscobertura de la elipse dentro del modelo lineal; (c) violación total
de las muestras AC.

**IEEE-33**

{{c1_ieee33}}

**PEGASE-2869**

{{c1_pegase}}

**Primera celda en que cada criterio supera 0,27 %:**

{{c1_first}}

**Otras cotas contra AC (por barra), IEEE-33:**

{{c1_other_ieee33}}

**Otras cotas contra AC (por barra), PEGASE-2869:**

{{c1_other_pegase}}

*Script: `experiments/c1_ac_sweep.py` · CSV: `results/c1_ac_sweep.csv`, `results/c1_first_exceedance.csv` · N = 20.000 por celda · semilla en la columna `seed`*

Lectura:

<<C1_LECTURA>>

---

## 7. C2 — Incertidumbre correlacionada entre fuentes

El código acepta ahora una covarianza de entrada: Q = C·Σ_ε·Cᵀ (`cdhim.SourceCovariance`, argumento
`sigma_eps` de `Study.run`). El valor por defecto es Σ_ε = I/3 y da exactamente los resultados anteriores
(filas r = 0, idénticas entre sí; además lo comprueba una prueba unitaria). Las muestras se generan con una
cópula gaussiana de un factor y marginales U[−1, 1]; la correlación de Pearson resultante es
ρ = (6/π)·asin(r/2). Punto de diseño, N = 20.000, modelo lineal y AC.

{{c2}}

*Script: `experiments/c2_correlated_sources.py` · CSV: `results/c2_correlated_sources.csv` · N = 20.000 · semilla en la columna `seed`*

Lectura:

- **Con Q = ⅓CCᵀ y fuentes correlacionadas el método no cubre.** Ya con r = 0,3 la elipse falla en el 26 %
  (IEEE-33) y el 35 % (PEGASE) de las muestras; con r = 0,9, en el 57 % y el 59 %. Los intervalos por barra
  fallan en promedio 7–32 % y 12–29 %. La crítica R2.2 / R3.4 describe un riesgo real, no una formalidad.
- **Con la Σ_ε verdadera la cobertura vuelve en el modelo lineal** (elipse: 0,05–0,30 %), pero cambian tres cosas:
  1. *La elipse casi desaparece.* Un factor común empuja |ρ| sobre 0,999 y las ramas pasan a Caja: de 67,7 %
     a 3,2 % en IEEE-33 y de 79,1 % a 23,0 % en PEGASE con r = 0,9.
  2. *El intervalo probabilístico queda más ancho que la caja afín.* La "reducción" de ancho es −52 %, −108 %
     y −214 % (IEEE-33). Con fuentes coherentes σ tiende a Σ|c_k|/√3 y el radio 3,44σ supera la envolvente
     determinista. Sin `clip_to_affine_box` el método entrega algo peor que la aritmética afín.
  3. *En PEGASE la linealización deja de ser válida.* Contra AC la propia caja afín se viola en hasta el 7 %,
     22 % y 37 % de las muestras en la peor barra. Con cargas que suben juntas el desbalance total crece como N
     y no como √N, la slack absorbe una variación grande y el término de segundo orden deja de ser despreciable.
     En IEEE-33 no ocurre (máximo 0,32 %).
- Conclusión para el paper: la correlación que CD-HIM modela es la que la red induce en las salidas a partir de
  fuentes independientes. La correlación entre fuentes requiere Σ_ε y, en redes de transmisión, volver a validar
  la linealización.
- **Datos reales:** bloqueado hasta que los autores indiquen la serie. La función
  `cdhim.uncertainty.estimate_source_covariance(serie)` ya recibe una matriz (tiempo × fuentes) y devuelve
  centros, semianchos, media de ε y Σ_ε.

---

## 8. C3 — Robustez no gaussiana, rehecha

Todas las distribuciones tienen media cero (salvo la transformación original del paper, que se deja tal
cual). Momentos exactos, verificados por muestreo en `tests/test_uncertainty.py`:

{{c3_marginals}}

Para cada distribución, dos covarianzas: *true_variance* usa Q = Var(ε)·CCᵀ y *uniform_1/3* usa Q = ⅓CCᵀ
(el supuesto del paper, sea cual sea la distribución real). Violación sobre **todas** las ramas en modo
Elipse: modelo lineal con N = 200.000 y AC con N = 5.000.

**IEEE-33 (δ = 15 %, 21 ramas en modo Elipse)**

{{c3_ieee33}}

**PEGASE-2869 (δ = 5 %, 2.543 ramas en modo Elipse)**

{{c3_pegase}}

**Distribución sesgada con fuentes correlacionadas** (Beta(2,5) con soporte respetado, cópula r = 0,6, N = 20.000):

{{c3_combo}}

*Script: `experiments/c3_nongaussian.py` · CSV: `results/c3_marginals.csv`, `results/c3_nongaussian.csv`, `results/c3_skewed_correlated.csv` · semilla en la columna `seed`*

Lectura:

- **La cifra del paper no se sostiene.** Con la transformación original, ε = (u − 0,285)·4, la elipse calibrada
  con 1/3 se viola en promedio 0,84 % (IEEE-33) y 0,83 % (PEGASE), con máximo 0,90 % y 0,95 %: **todas** las
  ramas en modo Elipse superan 0,27 %. El 0,19 % del paper salía de la regularización de 10⁻⁹ y de evaluar
  tres ramas que el algoritmo manda a Caja. Además el 9,7 % de la masa cae fuera de [−1, 1], y por eso en PEGASE
  hay barras en modo Caja (caja afín, determinista solo dentro de ese soporte) que se violan en el 9,7 % de las muestras.
- **Beta(2,5) bien planteada (soporte respetado, varianza real 0,05):** la elipse queda en 0,31 % y 0,29 %, apenas
  sobre el nivel de diseño, en casi todas las ramas (21 de 21 y 2.082 de 2.543). El sesgo opuesto, Beta(5,2), da lo
  mismo. Con la varianza 1/3 del paper la elipse es 2,6 veces más grande de lo necesario y la violación es 0.
- **Lo que decide la cobertura es la varianza, y después la curtosis:**
  - Varianza real mayor que 1/3 (Beta(0,5; 0,5), en forma de U): con 1/3 la violación sube a 1,1 % (máximo
    1,5 % y 1,9 %). Con la varianza correcta baja a 0,07–0,08 %.
  - Varianza real menor que 1/3 (triangular, normal truncada, Beta sesgadas): usar 1/3 da violación ≈ 0 a costa
    de intervalos entre 1,4 y 3,5 veces más anchos de lo necesario.
  - Con la varianza correcta: triangular 0,16 %, normal truncada 0,22–0,24 %, Beta sesgadas 0,29–0,31 %.
- **Colas pesadas, fuera del modelo:** t de Student (3 g.l.) 1,4–1,5 %, lognormal 1,1 %. Sus versiones truncadas
  y reescaladas a [−1, 1] tienen soporte acotado pero siguen dando 0,85–1,1 % con la varianza correcta.
  Soporte acotado no implica cobertura.
- **El tamaño del sistema no ayuda.** Las tasas de PEGASE (2.951 fuentes) son las mismas que las de IEEE-33
  (64 fuentes) para cada distribución. Cada par de barras depende de pocas fuentes cercanas, así que el teorema
  central del límite no actúa a escala de rama. La frase del paper "CD-HIM inherently rewards system scale"
  no tiene respaldo en estos datos.
- **Sesgo y correlación a la vez:** aun con la Σ_ε correcta la elipse queda en 0,45 % en ambos sistemas (el
  factor común sesgado no se promedia). Con el modelo del paper (uniforme, independiente) la violación es
  2,6 % (IEEE-33) y 12 % (PEGASE).
- "Heavy-tailed" no describe a Beta(2,5): su soporte es acotado. Conviene retirar el término del manuscrito.

---

## 9. C4 — Regla basada en G frente a umbral directo de correlación

Se implementó la regla alternativa (`selection_rule: rho`): Elipse si ρ_th ≤ |ρ| < 1 − δ_stab.

**κ por definición de caja** (G = 1 − κ·√(1−ρ²)):

{{c4_kappa}}

**Las dos reglas, con γ_th = 15 %:**

{{c4}}

**Efecto de las ramas que cambian (solo aparece con la caja afín):**

{{c4_effect}}

*Script: `experiments/c4_rho_rule.py` · CSV: `results/c4_kappa.csv`, `results/c4_rho_rule.csv` · tiempos: mediana de 20 corridas de la etapa por ramas*

Lectura:

- **Con las dos cajas de covarianza, la regla de G y la de ρ eligen exactamente las mismas ramas** (0 cambios en
  los dos sistemas y las dos linealizaciones), con el mismo tiempo. G no aporta información que ρ no tenga.
- Con la caja del código original (κ = π/4), ρ_th = 0: la regla de ganancia acepta todo y solo actúa δ_stab.
- Con κ = 1,032 (elipse χ² frente a caja 3σ) el umbral sí actúa: ρ_th = 0,567. En IEEE-33 no cambia nada
  (|ρ| mínimo 0,73); en PEGASE la selección baja de 2.543 a 2.339 ramas.
- **Con la caja afín κ cambia de rama en rama** (0,10–0,29 en IEEE-33; 0,005–2,7 en PEGASE) y G ya no es función
  de ρ. Es la única variante con información propia: κ mide cuánto se aleja Σ|c_k| de la escala gaussiana, es
  decir, cuántas fuentes efectivas pesan en esa rama. Pero G frente a la caja afín es casi siempre mayor que
  15 % (mínimo 89 % en IEEE-33), así que el umbral vuelve a no filtrar: cambia 1 rama de 3.216 en PEGASE,
  sin efecto medible en violaciones.
- Para el paper: o se presenta la regla como umbral de correlación (más simple y equivalente), o se define G
  contra la caja afín y se elige un γ_th que realmente separe ramas. Llamarla "topológica" o "espectral" no se
  justifica con la definición actual.

---

## 10. C5 — Violación a nivel sistema y corrección de Bonferroni

Punto de diseño, N = 20.000 (las mismas muestras AC que B1). "≥ 1 barra fuera" es la fracción de muestras en
que al menos una barra PQ sale de su intervalo. Bonferroni usa el cuantil χ²₂ al nivel 1 − α/m, con
m = número de ramas en modo Elipse y α = 0,0027.

{{c5}}

*Script: `experiments/c5_system_level.py` · CSV: `results/c5_system_level.csv` · N = 20.000 · semilla en la columna `seed`*

Lectura:

- **99,73 % por rama no es 99,73 % para el sistema.** En PEGASE, casi la mitad de las muestras (48 % en el modelo
  lineal, 50 % contra AC) viola al menos una de las 2.543 elipses, con 2,5 elipses violadas por muestra en
  promedio. En IEEE-33, el 1 %.
- Los intervalos por barra con λ_max aguantan mejor porque están sobredimensionados: 3,1 % de las muestras AC
  tiene alguna barra fuera en PEGASE (9,3 % con la proyección exacta; 42 % con la caja 3σ).
- **Bonferroni devuelve la cobertura simultánea**, con costo: los radios crecen 23 % (IEEE-33) y 53 % (PEGASE).
  La reducción de ancho frente a la caja afín baja de 42,8 % a 29,6 % en IEEE-33 y **de 34,1 % a −0,5 % en
  PEGASE con λ_max**: el intervalo simultáneo queda tan ancho como la caja afín, y más ancho en 308 barras. Con la
  proyección exacta todavía se conserva una reducción de 45,9 %.
- El CRI de área frente a la caja de covarianza no cambia con Bonferroni (71,08 %) porque es un cociente libre
  de escala: otra señal de que no mide el conservadurismo de la salida. Frente a la caja afín sí baja (96,0 → 90,9 %).
- El residuo de 0,3–0,4 % contra AC con Bonferroni en PEGASE es error de linealización (la barra en modo Caja
  de B1), que ningún cuantil corrige.
- La frase "no branch is under-covered" solo es defendible por rama y dentro del modelo lineal.

---

## 11. C6 — Transferibilidad de γ_th

Barrido γ_th ∈ {5, 10, 15, 20, 25} % en cinco sistemas de pandapower (δ = 15 % en `case33bw`, 5 % en el resto),
para las tres cajas que pueden ir dentro de G. `case2869pegase` es aquí la versión de MATPOWER que trae
pandapower (punto de operación resuelto), no el archivo PGLib del resto del informe. Violación contra AC con
N = 5.000 por sistema.

{{c6}}

**Óptimo según el criterio del paper (CRI sobre tiempo normalizado) y estabilidad:**

{{c6_opt}}

*Script: `experiments/c6_gamma_transfer.py` · CSV: `results/c6_gamma_transfer.csv`, `results/c6_optimum.csv` · N = 5.000 · tiempos: mediana de 10 corridas*

Lectura:

- **Con la caja del paper no hay nada que transferir.** En los cinco sistemas, γ_th = 5, 10, 15 y 20 % dan la misma
  selección; solo 25 % cambia algo, y solo en los tres sistemas grandes. En `case33bw` y `case118` los cinco
  valores son indistinguibles.
- **El criterio "CRI sobre tiempo" no define un óptimo.** El tiempo no depende de γ_th: la etapa por ramas calcula
  los espectros de todas las ramas antes de decidir, y la variación medida en el barrido es ruido de
  cronometraje. El "óptimo" salta entre 5 y 25 % de un sistema a otro por ese ruido.
- **El criterio "solo CRI" es degenerado**: como el CRI se promedia sobre las ramas elegidas, siempre gana el
  umbral más alto, que se queda con las ramas más correlacionadas.
- Con la caja 3σ (κ = 1,032) el umbral sí actúa y la fracción de ramas en Elipse depende del sistema (64–89 %):
  ahí γ_th **es dependiente de la topología**, como sospechaba R2.7. No hay un valor común que se pueda
  justificar con estos datos.
- La violación máxima por barra contra AC no depende de γ_th, pero sí del sistema: 1,06–1,26 % en `case300`
  y 0,60 % en `case1354pegase`, ambos sobre 0,27 % a δ = 5 %. Son barras donde manda el error de linealización.
  El rango de validez δ ≤ 20 % del paper no se puede afirmar en general.

---

## 12. C7 — Escalabilidad

Un hilo, 10 repeticiones por sistema (media ± desviación estándar). Etapas: construcción del Jacobiano
(dS/dV disperso), factorización LU (SuperLU, COLAMD), resoluciones triangulares (N_u lados derechos; se
guarda solo el bloque de |V|), etapa por ramas (covarianzas, autovalores 2×2, decisión y radios). El flujo
de potencia base se reporta aparte: el método lo necesita, pero no es parte del método.

{{c7}}

**Relleno de la factorización:**

{{c7_fill}}

**Exponentes en escala log-log (y ∝ N^k, N = número de barras):**

{{c7_fit}}

*Script: `experiments/c7_scalability.py` · CSV: `results/c7_timing.csv`, `results/c7_fill.csv`, `results/c7_fit.csv`*

Lectura:

<<C7_LECTURA>>

---

## 13. C8 — Comparación con otros métodos

Todos los métodos de intervalo usan el mismo modelo lineal corregido. Punto de diseño, N = 20.000 flujos AC
(las mismas muestras que B1). "Rango AC 99,73 %" es el semiancho empírico entre los cuantiles 0,135 % y
99,865 % de las muestras AC.

{{c8}}

**Tiempos (un hilo):**

{{c8_time}}

*Script: `experiments/c8_methods.py` · CSV: `results/c8_methods.csv`, `results/c8_timing.csv` · N = 20.000*

Lectura:

<<C8_LECTURA>>

**Verificación de CLTM (Ref. [10]).** En el repositorio no hay ninguna implementación de CLTM. Lo que el
manuscrito llama "CLTM Box" es la caja 3σ del propio modelo lineal de CD-HIM (`cdhim.geometry.cltm_box_halfwidth`).
Contrastado con el texto completo de la Ref. [10], que los autores enviaron después de la primera entrega:

- *Lo que coincide.* El intervalo es [μ − 3σ, μ + 3σ] (p = 0,9973) y cada término aporta la varianza de una
  uniforme, ancho²/12. Con símbolos en [−1, 1] eso es 3·√(Σ c_k²/3), que es lo que calcula
  `cltm_box_halfwidth`. Las formas afines parten de las sensibilidades del Jacobiano inverso en los puntos
  medios, como aquí.
- *Lo que falta.* En la variante AC, la Ref. [10] no usa todo el rango [−1, 1] de cada símbolo: primero resuelve
  un problema no lineal y no convexo por punto interior (método de escenarios de optimización) que entrega el
  rango [ε_min, ε_max] de cada elemento de ruido, y la varianza de cada término sale de ese rango. Además, si el
  número de términos efectivos (los de frontera mayor que 10⁻⁶) es menor que M = 6, no usa el intervalo 3σ sino
  una contracción con p′ = 0,9973^(1/m). Ninguno de los dos pasos está implementado.
- *Consecuencia.* La "CLTM Box" del manuscrito es el paso de TCL de la Ref. [10] aplicado a las sensibilidades de
  CD-HIM. **No es una reproducción de la Ref. [10]** y el manuscrito no debería llamarla así.
- *Escala.* La Ref. [10] prueba la variante AC solo en IEEE-30 (cerca de 30 s) y la DC en el sistema polaco de
  2.383 barras (cerca de 3 s). No reporta ningún caso AC del tamaño de PEGASE, y el "12,5 s" de la Tabla II del
  manuscrito no tiene código que lo respalde.

**Ref. [11].** Leído el texto completo. El conjunto híbrido caja ∩ elipsoide se define sobre las *entradas*
(las inyecciones), no sobre pares de voltajes, y los intervalos de salida se obtienen con el mismo método de
escenarios de optimización: dos problemas no lineales por variable de salida. El radio del elipsoide se calibra
para normales independientes (cuantil χ² con tantos grados de libertad como entradas) o con un programa con
restricciones de azar resuelto por enjambre de partículas. El único caso de prueba es IEEE-30 con ±20 %. No hay
código en el repositorio y no la implementé (sección 15).

**Ref. [12].** Sin código y sin PDF. Sigue bloqueada (sección 15).
"""
