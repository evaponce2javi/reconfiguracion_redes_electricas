# AUDIT.md — Auditoría del código CD-HIM (Parte A, solo lectura)

Código auditado: `finalFinalCode.py` (commit `4750ce2`, 1.021 líneas, único archivo de código del repositorio).
Todas las citas son `finalFinalCode.py:línea`. Las cifras marcadas "medido" se obtuvieron ejecutando el
código original **sin modificarlo** (pandapower 3.5.6, NumPy 2.4.6, SciPy 1.18.1, Python 3.12.10), con las
definiciones del archivo cargadas tal cual y, aparte, con el script completo corrido de principio a fin.

`caja-elipse.tex` no está en el repositorio, así que no pude ver el texto de los `[TODO]`; uso la tabla de
correspondencia del encargo.

---

## Resumen: lo que hay que saber antes de leer el detalle

| # | Hallazgo | Evidencia | Consecuencia |
|---|---|---|---|
| 1 | **El Jacobiano tiene un error de signo** en los términos fuera de la diagonal del bloque H = ∂P/∂θ. | `:115` `H = -ViVj * Gsin_Bcos`. Lo correcto es `+`. Comparado contra diferencias finitas: `max|H_código + H_exacto| = 6·10⁻⁹` fuera de la diagonal (es la negación exacta); los otros tres bloques coinciden a 10⁻⁸. | En IEEE-33 el 27,5 % de las sensibilidades ∂V/∂P sale con el signo cambiado y la matriz de sensibilidad de voltaje difiere 123 % (norma de Frobenius relativa) de la correcta. Todas las cifras del paper (ρ, G, CRI, reparto Elipse/Caja, Tabla I, Beta) se calcularon con esta matriz. |
| 2 | **El Jacobiano se evalúa en arranque plano**, no en el punto convergido. | `:104-107` fija V = 1, θ = 0 (V de consigna en PV). `to_ppc(..., init='flat')` en `:945` y `:979`. No se ejecuta ningún flujo de potencia base antes de linealizar. | Contradice la Sección II-A ("pre-converged AC-NR solution (V\*, θ\*)… Unlike flat-start approximations"). |
| 3 | **γ_th = 15 % no hace nada.** G = 1 − (π/4)·√(1−ρ²) ≥ 21,46 % para cualquier ρ. | `:197-200`. Medido: G mínimo = 46,6 % (IEEE-33) y 21,46 % (PEGASE). Con γ_th ∈ {5, 10, 15, 20} la selección es idéntica en ambos sistemas. | La decisión Elipse/Caja la toma **solo** la guarda `abs(rho) < 0.999` (`:204`). Las ramas en modo Caja son las de **mayor** correlación (\|ρ\| ≥ 0,999), no las "débilmente acopladas". |
| 4 | **La validación contra AC no existe en el código.** | Las muestras AC-NR (`delta_v_mcs`) solo se usan para dibujar la nube de la Fig. 4 (`:420-421`). Ninguna línea compara una cota con una muestra AC. | El "R_viol = 0,00 % en 20.000 muestras AC-NR" (IEEE-33) y el "0,06 % en el punto de diseño" (PEGASE) no tienen respaldo en este código. El 0,06 % es el valor de la Tabla I (modelo lineal). |
| 5 | **Los intervalos de salida del Algoritmo 1 no se calculan** en el flujo principal. | El radio `sqrt(chi2 * max(eigvals))` solo aparece en una función de gráfico (`:744`), con χ² = 11,82 escrito a mano y recortado a la caja afín (`:752-755`). | El CRI reportado es reducción de área de una elipse que el método no entrega como salida. |
| 6 | **La distancia de Mahalanobis lleva una regularización de 10⁻⁹** en la diagonal de Q. | `:472-473`. Medido: λ_min < 10⁻⁹ en el 42 % de las ramas de IEEE-33 y el 47 % de PEGASE; en las ramas del experimento Beta de PEGASE λ_min ≈ 10⁻¹², mil veces menor que la regularización. | Engorda la elipse y baja artificialmente las tasas de violación. Es lo que hace que la Tabla I dependa de δ y que la curva P-P de la Fig. 5 quede sobre la diagonal (ver A5). |
| 7 | **El álgebra lineal es densa**, no dispersa. | `scipy.linalg.lu_factor` sobre un Jacobiano denso (`:130-131`); `Ybus` denso (`:79`). | Contradice Sección II-A, III-A y la Contribución 2 (LU disperso, Gilbert-Peierls, O(N^1,2)). |
| 8 | **El experimento Beta se corre sobre ramas que el algoritmo manda a Caja.** | Selección con `abs(rho) < 0.9999` (`:560`) frente a la regla `< 0.999` (`:204`). Medido: las tres ramas elegidas tienen ρ = 0,99990 y decisión `BOX` en ambos sistemas. | El 0,19 % de PEGASE y el 0,40 % de IEEE-33 no miden ramas en modo Elipse. |

Los puntos 1, 2 y 7 son errores de implementación; 3 a 6 y 8 cambian la interpretación de los resultados del paper.

---

## A1. Qué caja define A_box y cuánto vale κ

**Dónde se calcula G** — `solve_hybrid_strategy`, `:186-200`:

```python
C = np.vstack([c1, c2])                                   # :186
Q = (C @ C.T) / 3.0                                       # :187
sigma_x = np.sqrt(Q[0, 0]); sigma_y = np.sqrt(Q[1, 1])    # :188-189
rho = Q[0, 1] / (sigma_x * sigma_y)                       # :193
eigvals = scipy.linalg.eigvalsh(Q); eigvals[eigvals < 0] = 0   # :194-195
area_ell = np.pi * np.sqrt(eigvals[0]) * np.sqrt(eigvals[1])   # :197
area_box = 4 * sigma_x * sigma_y                               # :198
ganancia = max(0, (1.0 - area_ell/area_box) * 100)             # :200
```

- **A_ell = π·√λ₁·√λ₂** (elipse de covarianza a escala 1σ, sin χ²).
- **A_box = 4·σ_x·σ_y** con σ = √Q_ii.

**Cuál de las tres cajas:** la **caja de covarianza** con semianchos k·√Q_ii y **k = 1**. No es la envolvente afín Σ|c_k| ni la caja 3σ. La elipse tampoco lleva el factor χ² = 11,83: ambas figuras están a la misma escala.

**κ:** como λ₁λ₂ = det Q = Q₁₁Q₂₂(1−ρ²),

  G = 1 − (π/4)·√(1−ρ²)  ⟹  **κ = π/4 ≈ 0,7854**, igual para todas las ramas.

Medido sobre las 31 ramas de IEEE-33: κ mínimo = mediana = máximo = 0,785398163397 (coincide con π/4 a 12 cifras). En PEGASE la mediana es π/4; las únicas desviaciones aparecen en ramas con |ρ| → 1, donde 1−ρ² se pierde por redondeo.

**Es el primer caso de referencia del encargo, con la contradicción resuelta así:**

- Con ρ = 0, G = 21,46 % ≥ 15 %. **Toda rama pasa el umbral de ganancia.** Medido: en IEEE-33 ninguna rama tiene G < 15 % (mínimo 46,6 %); en PEGASE tampoco (mínimo 21,46 %).
- El 51,6 % de IEEE-33 **no** lo produce γ_th. Lo produce la segunda condición de `:204`, `abs(rho) < 0.999`: 15 de las 31 ramas tienen |ρ| ≥ 0,999 y caen a Caja. En PEGASE, 234 de 3.218.

**ρ_th equivalente a γ_th = 15 %: no existe.** La condición G ≥ 15 % se cumple para todo ρ. La regla de selección efectiva del código es

  Elipse ⟺ |ρ| < 0,999.

El umbral de ganancia solo empieza a actuar para γ_th > 1 − π/4 = 21,46 %. Para γ_th = 25 % equivale a |ρ| ≥ 0,297.

Barrido medido de γ_th (ramas en Elipse / ramas analizadas):

| γ_th | IEEE-33 | PEGASE-2869 |
|---|---|---|
| 5 % | 16 / 31 | 2.984 / 3.218 |
| 10 % | 16 / 31 | 2.984 / 3.218 |
| 15 % | 16 / 31 | 2.984 / 3.218 |
| 20 % | 16 / 31 | 2.984 / 3.218 |
| 25 % | 16 / 31 | 2.649 / 3.218 |

En IEEE-33 los cinco valores del barrido dan exactamente la misma selección (|ρ| mínimo = 0,734). El estudio OAT de la Sección IV-A (γ_th < 10 % sube el tiempo 25 %, γ_th > 20 % devuelve demasiadas ramas a Caja, 15 % es el óptimo) **no se puede reproducir con este código**: ni la selección ni el tiempo cambian, porque además la descomposición espectral se calcula para todas las ramas antes de decidir (`:194`).

**Base del CRI:** la misma caja de covarianza con k = 1 (el CRI es el promedio de `ganancia`, ver A4).

**"CLTM Box" de la Fig. 4:** un rectángulo de semianchos **3·√Q_ii** (`:443-445`), o sea la caja de covarianza con k = 3, calculada con la misma Q = ⅓CCᵀ. No hay ninguna implementación del método CLTM de la Ref. [10] en el repositorio.

**"CLTM Box" de la Tabla II (12,5 s; ≈100×; 0 % Ref):** no hay código que la calcule ni la cronometre. Ver A4 y el punto bloqueado correspondiente.

Resumen de las tres cajas que el paper usa como si fueran una:

| Objeto | Semiancho | Dónde aparece en el código |
|---|---|---|
| Caja de covarianza k = 1 | √Q_ii | G y CRI (`:198`) |
| Caja de covarianza k = 3 ("CLTM Box") | 3·√Q_ii | Fig. 4 (`:443`) |
| Envolvente afín | Σ\|c_k\| | Solo en el gráfico del perfil de voltaje (`:728`) |

Y de las dos elipses:

| Objeto | Escala | Dónde aparece |
|---|---|---|
| Elipse de G | 1 (sin χ²) | `:197` |
| Elipse dibujada y de validación | χ²₂(0,9973) = 11,83 | `:409-412`, `:511`, `:563` |

Con la elipse a χ² = 11,83 y la caja a k = 3 el cociente sería κ = π·11,83/36 = 1,032 y γ_th = 15 % equivaldría a |ρ| ≥ 0,567. Eso **no** es lo que hace el código.

---

## A2. Modelo de incertidumbre

**Símbolos de P_d y Q_d: separados e independientes.** `:143`:

```python
noise_vec = np.concatenate([self.dP_inj[self.ns], self.dQ_inj[self.pq]])
```

Un símbolo por barra no-slack para P y otro, distinto, por barra PQ para Q. P y Q de una misma carga varían de forma independiente (el factor de potencia no se conserva). El Monte Carlo AC hace lo mismo: sortea `ruido_p` y `ruido_q` por separado (`:255-256`).

**N_u (medido):**

| | Longitud del vector | Símbolos P no nulos | Símbolos Q no nulos | N_u efectivo |
|---|---|---|---|---|
| IEEE-33 | 64 | 32 | 32 | **64** |
| PEGASE-2869 | 5.227 | 1.485 | 1.466 | **2.951** |

El vector tiene una posición por cada fila del Jacobiano; las barras sin carga quedan con semiancho 0. Se descartan cargas con |P| o |Q| ≤ 10⁻⁶ p.u. (`:53-54`).

**Cómo se aplica δ** — `:56-57`:

```python
self.dP_inj[mask_p] = np.abs(p_val[mask_p] * tol_p)
self.dQ_inj[mask_q] = np.abs(q_val[mask_q] * tol_q)
```

Sí: δ es el **semiancho relativo** de cada inyección (semiancho = δ·|P_d|, δ·|Q_d|), y es **el mismo para todas las fuentes** y para P y Q (`tol_p = tol_q` en `:946` y `:982`). El coeficiente afín es sensibilidad × semiancho (`:160`).

**Qué inyecciones son inciertas:** solo la demanda neta por barra, columnas PD y QD del `ppc` (`:46-47`). La generación despachada (P de los generadores, consigna de V) es determinista. Dos matices:

1. En PEGASE-2869 hay **180 barras con PD < 0** (inyecciones modeladas como carga negativa). `from_mpc` las convierte en `sgen`. El modelo lineal las trata como inciertas porque lee PD del `ppc`; el Monte Carlo AC **no las perturba**, porque solo modifica `net.load` (`:245-246`, `:265-266`). Los dos modelos no tienen el mismo conjunto de fuentes.
2. Las cargas en la barra slack y la Q de cargas en barras PV quedan fuera del modelo lineal. Es correcto: no afectan los voltajes de las barras PQ.

---

## A3. Álgebra lineal

**Rutina de LU:** `scipy.linalg.lu_factor` (`:131` para el Jacobiano AC; `:101` para la matriz DC). Es **LAPACK `getrf` denso** con pivoteo parcial por filas. No se usa `scipy.sparse` en ningún punto del archivo (el único import es `import scipy.linalg`, `:12`). No hay SuperLU ni Gilbert-Peierls.

La matriz que se factoriza también es densa: `Ybus` se arma con `np.zeros((nb, nb), dtype=complex)` (`:79`), los cuatro bloques del Jacobiano son matrices `nb × nb` completas (`:110-118`) y se ensamblan con `np.block` (`:130`). Para PEGASE-2869 el Jacobiano es de 5.227 × 5.227 (27,3 millones de entradas, 219 MB).

**¿Se forma J⁻¹ completa?** No. Se extraen filas: para cada barra pedida se resuelve Jᵀy = e_i con `scipy.linalg.lu_solve(..., trans=1)` (`:159`), dos sustituciones triangulares **densas** O(n²) por llamada.

**Cuántas resoluciones:** dos llamadas a `get_affine_coeffs` por rama PQ-PQ (`:181-182`), sin caché. Una barra que toca cuatro ramas se resuelve cuatro veces.

| | Ramas PQ-PQ | Llamadas a `lu_solve` | Barras PQ distintas |
|---|---|---|---|
| IEEE-33 | 31 | 62 | 32 |
| PEGASE-2869 | 3.218 | 6.436 | 2.359 |

En PEGASE se hacen 6.436 resoluciones para obtener como mucho 2.359 filas distintas. El experimento de la Tabla I vuelve a correr todo el lazo (`:503`).

**Ordenamiento de reducción de relleno:** ninguno. Con LU denso el concepto no aplica; el relleno es total.

**Factorización DC sin uso:** `solve_dc_sensitivities()` se llama en `:949` y `:985` y su resultado no se usa en ningún experimento (todo corre con `model='AC'`).

Para el `[TODO]` del manuscrito: lo correcto hoy es "LU denso (LAPACK `getrf` vía `scipy.linalg.lu_factor`)". La frase del paper sobre `scipy.linalg` era literal, y el revisor 3 (punto 8) tenía razón.

---

## A4. Cómo se obtienen las cifras del paper

### 58,18 % y 86,92 % (CRI)

`seccion_4_ejecucion_y_tabla`, `:838-841`:

```python
critical_lines = [r for r in results if r['decision'] == 'ELLIPSE']
avg_gain_critical = np.mean([r['ganancia'] for r in critical_lines])
```

Es el **promedio de G solo sobre las ramas en modo Elipse** (se imprime como "Reducción Promedio (Zona Crítica)", `:853`). G es reducción de **área**: elipse de covarianza a escala 1σ frente a caja de covarianza a escala 1σ, es decir 1 − (π/4)√(1−ρ²). No es reducción de ancho, no usa los intervalos de salida, no usa la caja afín y no es un promedio global.

| | Paper | Código actual (medido) | Promedio sobre todas las ramas |
|---|---|---|---|
| IEEE-33 | 86,92 % | **86,70 %** | 92,29 % |
| PEGASE-2869 | 58,18 % | **57,80 %** | 60,77 % |

El código del repositorio no reproduce el segundo decimal (diferencias de 0,22 y 0,38 puntos). Los conteos sí se reproducen exactos, y el título de la Fig. 4 también ("Branch 5--6, Reduction: 96.45 %, Rho: 0.9990"). No encontré la causa de la diferencia en el CRI; lo más probable es otra versión de pandapower o del propio script al generar la tabla del paper. Queda como punto abierto para los autores.

### 51,6 % y 92,7 %

`:851`: `len(critical_lines) / total_lines`. El **denominador** es `len(results)`: ramas con ambos extremos PQ (`:178`), con coeficientes no nulos (`:184`) y σ_xσ_y > 10⁻¹² (`:191`). El numerador son las que cumplen `ganancia >= 15 and abs(rho) < 0.999`, que por A1 equivale a `abs(rho) < 0.999`.

| | Paper | Medido |
|---|---|---|
| IEEE-33 | 16 / 31 = 51,6 % | 16 / 31 = 51,6 % ✔ |
| PEGASE-2869 | 92,7 % de 3.218 | 2.984 / 3.218 = 92,73 % ✔ |

Lectura correcta: el 48,4 % de IEEE-33 y el 7,3 % de PEGASE que van a Caja son las ramas con |ρ| ≥ 0,999.

### R_viol

Hay una sola definición en el código, y es **por rama, contra la elipse, dentro del modelo lineal**:

```python
points = C @ eps                                   # :530  (o c1 @ eps, c2 @ eps en :582, :625)
d2 = solve_mahalanobis_robust(Q, points)           # :533
violaciones = np.sum(d2 > chi2_val)                # :535
tasa = (violaciones / 5000) * 100                  # :536
```

- Se evalúa **una rama** (Tabla I, Fig. 5) o **tres** (Fig. 6).
- El denominador es el número de muestras: 5.000 (Tabla I), 10.000 (Beta, `:564`), 20.000 (P-P, `:624`).
- **No hay evaluación por barra contra el intervalo de salida.**
- **No hay evaluación contra AC.** Ver hallazgo 4.

### Tiempos

- **1.246,66 s y 1.263,48 s:** el tiempo del Monte Carlo AC lo imprime `generate_nonlinear_mcs_samples` (`:285-286`) cada vez que se ejecuta el script. Son dos valores de la misma medición en **dos ejecuciones distintas**; no hay otra fuente en el código. Ninguno se guarda en disco.
- **52,05 s de CD-HIM:** `:831-834` cronometra **solo** `solve_hybrid_strategy`: el lazo por ramas con las 6.436 sustituciones densas, Q, autovalores y decisión. **No incluye** la construcción de Ybus, la construcción del Jacobiano ni la factorización LU (`:984-986`), ni un flujo de potencia base (no existe). Es una sola medición; el "± 0,15 s, n = 4" no sale de este código. Medido aquí en una corrida aislada: construcción de Ybus 0,07 s, LU DC 0,29 s, Jacobiano + LU denso 2,1 s, lazo por ramas 66 s (con otra tarea en paralelo).
- **12,5 s de CLTM:** **no hay código**. No existe implementación ni cronometraje de CLTM en el repositorio.
- **0,0037 s de IEEE-33:** mismo cronómetro (`:831-834`); medido 0,004–0,005 s.

Los tiempos de Monte Carlo medidos en esta máquina con el script original están en la sección "Corrida completa del script original" al final.

### Extrapolación de IEEE-33 (377 s)

`estimar_tiempo_monte_carlo`, `:904-930`: corre 50 flujos `pandapower.runpp(..., numba=False)` (`:908-919`), promedia el tiempo de pared de cada uno y lo multiplica por 20.000 (`:923-924`). Observaciones:

1. Antes de cada corrida hace `net.load['p_mw'] *= 1.001` (`:909`), acumulativo y sobre la red original: al terminar la carga quedó 5,1 % más alta. No se restaura.
2. 377 s / 20.000 = 18,9 ms por flujo en una red de 33 barras. Ese tiempo es casi todo sobrecosto de pandapower (conversión de la red en cada llamada), no resolución numérica. El speed-up de 1,08·10⁵ compara un lazo de 4 ms con ese sobrecosto.
3. El mismo script **sí mide** el Monte Carlo completo de IEEE-33 (`:956`, 20.000 muestras). La extrapolación no era necesaria.

---

## A5. Experimentos existentes

### Tabla I

**Confirmado: las muestras se propagan por el mismo modelo lineal.** `:529-530`: `eps = rng.uniform(-1, 1, ...)`, `points = C @ eps`, con la misma C que define Q (`:524-525`). No interviene ningún flujo AC. La tabla mide la cobertura de la elipse χ² para una suma de uniformes, no el error de linealización.

**Confirmado: sin la regularización, la distancia de Mahalanobis no depende de δ.** C es proporcional a δ (`:517-518`), así que d² = (Cε)ᵀ(⅓CCᵀ)⁻¹(Cε) es invariante. Además la semilla se reinicia dentro del lazo (`rng = np.random.default_rng(99)`, `:528`), por lo que **las ε son idénticas para los cuatro valores de δ**. La tabla debería tener cuatro filas iguales.

**Por qué en el paper no son iguales:** `solve_mahalanobis_robust` suma 10⁻⁹ a la diagonal de Q (`:472-473`). Q escala con δ² y la regularización no, así que su peso relativo cambia con δ. A δ pequeño la elipse queda proporcionalmente más gruesa y la violación baja. La tendencia creciente de la Tabla I es ese artefacto.

Medido con el código actual (misma rama, misma semilla):

| δ | IEEE-33 paper | IEEE-33 medido | PEGASE paper | PEGASE medido |
|---|---|---|---|---|
| 5 % | 0,00 % | 0,00 % | 0,06 % | 0,04 % |
| 10 % | 0,06 % | 0,00 % | 0,08 % | 0,04 % |
| 20 % | 0,08 % | 0,02 % | 0,12 % | 0,10 % |
| 30 % | 0,10 % | 0,04 % | 0,12 % | 0,10 % |

Los valores no coinciden con los del paper (misma forma, distintas cifras). Con N = 5.000 cada muestra vale 0,02 %, así que las diferencias son de una a tres muestras.

### Fig. 5 (P-P): el mismo artefacto

Medido en PEGASE, rama de la figura, N = 20.000, semilla 101:

| Confianza teórica | Cobertura con regularización (código) | Cobertura exacta (sin regularización) |
|---|---|---|
| 0,50 | 0,7312 | 0,4699 |
| 0,70 | 0,8666 | 0,6874 |
| 0,80 | 0,9234 | 0,7988 |
| 0,90 | 0,9670 | 0,9131 |
| 0,95 | 0,9863 | 0,9608 |
| 0,99 | 0,9990 | 0,9953 |
| 0,9973 | 0,9998 | 0,9993 |

La columna central reproduce la Fig. 5 (arranca en 0,73). La columna derecha es la cobertura real: por debajo de la diagonal hasta 0,80 y por encima después. La afirmación "the curve lying above the diagonal confirms conservatively calibrated bounds" se apoya en la regularización, no en el método.

### Experimento Beta

Transformación del código, `:579-580`:

```python
raw_beta = rng.beta(2, 5, size=(len(c1), n_mc))
eps = (raw_beta - 0.285) * 4.0
```

- Constante 0,285 (el paper dice 0,286; la media exacta es 2/7 = 0,285714). Queda un sesgo de media de +0,0029.
- Soporte de ε: **[−1,140; 2,860]**. No es [−1, 1].
- Varianza de ε: 16 · 0,02551 = **0,408**. No es 1/3.
- N = 10.000 por rama (`:564`), semilla 42 (`:565`).
- Q sigue siendo ⅓CCᵀ, con la regularización de 10⁻⁹.

La crítica 5 del revisor 3 es correcta punto por punto.

Selección de ramas (`:560-561`): las tres de mayor `ganancia` entre las de |ρ| < **0,9999**. Como G crece con |ρ|, son las tres ramas con |ρ| más cercano a 0,9999, y todas están por encima del 0,999 de la regla de decisión:

| Sistema | Rama (índices) | ρ | Decisión del algoritmo | λ_min | Violación Beta medida |
|---|---|---|---|---|---|
| IEEE-33 | 28–29 | 0,999896 | BOX | 1,2·10⁻⁹ | 0,40 % |
| IEEE-33 | 5–25 | 0,999881 | BOX | 3,8·10⁻¹⁰ | 0,28 % |
| IEEE-33 | 14–15 | 0,999848 | BOX | 5,9·10⁻¹⁰ | 0,30 % |
| PEGASE | 1523–698 | 0,999899 | BOX | 1,0·10⁻¹² | 0,19 % |
| PEGASE | 83–2494 | 0,999899 | BOX | 1,3·10⁻¹² | 0,19 % |
| PEGASE | 1866–1926 | 0,999896 | BOX | 1,5·10⁻¹² | 0,23 % |

El 0,40 % de IEEE-33 y el 0,19 % de PEGASE se reproducen (el máximo de PEGASE con este código es 0,23 %, no 0,19 %). En PEGASE la regularización es mil veces mayor que λ_min: el eje menor de la elipse lo fija el 10⁻⁹, no la covarianza.

### Selección de la rama de "máximo estrés"

El paper describe un cribado de Monte Carlo en dos etapas con el cociente d²_max/χ². **Ese procedimiento no está en el código.** Lo que hay (`:503-508`, `:616-618`):

1. Volver a correr `solve_hybrid_strategy` (con umbral 10 % en `:503`; da lo mismo que 15 % por A1).
2. Descartar ramas con |ρ| ≥ 0,999.
3. Ordenar por `ganancia` descendente y tomar la primera.

Es la rama con el mayor |ρ| por debajo de 0,999. No hay muestreo ni cociente de estrés. Ramas resultantes: 5–6 en IEEE-33 (ρ = 0,99898, G = 96,45 %) y 2314–494 en PEGASE (ρ = 0,99899, G = 96,47 %).

### Muestras AC-NR de PEGASE-2869

**No están guardadas.** `generate_nonlinear_mcs_samples` devuelve las matrices en memoria (`:289`) y el script no escribe ningún `.npy`, `.npz` ni `.csv`. Hay que regenerarlas. La semilla es 42 (`:219`, `:248`), así que el sorteo es reproducible. Las muestras que no convergen se cuentan y se descartan (`:281-283`); el conteo se imprime pero no se guarda. El comentario de `:992` dice "Muestra reducida 100 flujos" y el código pide 20.000 (`:993`).

---

## Otros defectos del código (revisión completa)

### Lógica

| Línea | Defecto |
|---|---|
| `:115` | Signo de H fuera de la diagonal (hallazgo 1). |
| `:104-107` | Linealización en arranque plano (hallazgo 2). |
| `:120-123` | Las diagonales omiten los términos de inyección (H_ii = −Q_i − B_ii V², N_ii = P_i + G_ii V², etc.). Solo es exacto con V = 1 y sin shunts ni taps. |
| `:109-123` | El Jacobiano está en la forma V·∂/∂V, así que la incógnita es ΔV/V y no ΔV. En arranque plano coincide; en un punto convergido habría que multiplicar por V. |
| `:87` | El shunt de línea del lado "from" no se divide por \|a\|² (en MATPOWER, Y_ff = (y + jb/2)/\|a\|²). Afecta a transformadores con carga capacitiva. |
| `:65-73` | No se filtra BR_STATUS; `tap[tap == 0] = 1.0` modifica el `ppc` del llamador. |
| `:204` vs `:398`, `:505`, `:560`, `:616` | El umbral de estabilidad está escrito a mano cinco veces, con dos valores distintos (0,999 y 0,9999). |
| `:472-473` | Regularización absoluta de 10⁻⁹ (hallazgo 6). |
| `:580` | Transformación Beta (soporte, varianza, media). |
| `:245-266` | El Monte Carlo no perturba los `sgen` que el modelo lineal sí considera inciertos. |
| `:236`, `:277` | Las desviaciones AC se miden respecto del caso base AC convergido, mientras el modelo lineal está linealizado en arranque plano. Los centros no corresponden al mismo punto. |
| `:720` | El "voltaje base" del perfil es la columna VM del `ppc` con `init='flat'`: vale 1,0 en todas las barras, no es el perfil resuelto. |
| `:743` | χ² = 11,82 escrito a mano; en el resto se usa `chi2.ppf(0.9973, 2)` = 11,829. |
| `:752-755` | El radio de salida se recorta a la caja afín, cosa que el Algoritmo 1 del paper no hace. |
| `:909` | `net.load['p_mw'] *= 1.001` acumulativo dentro del cronometraje. |
| `:920` | `except:` desnudo; un fallo de convergencia desaparece del promedio sin aviso. |

### Estructura

- Sin `if __name__ == "__main__":`. Importar el archivo ejecuta todo, incluidos 40.000 flujos de potencia.
- `:4-9`: instala paquetes con `os.system('pip install ...')` y sigue sin volver a importar.
- Imports duplicados (`:16` y `:22`), bloque "MAIN" duplicado (`:932-937`).
- Parámetros repartidos por el archivo: δ en `:946`, `:956`, `:982`, `:993`; γ_th en `:503`, `:675`, `:833`, `:882`; semillas 42, 99, 101, 55; tamaños 5.000, 10.000, 20.000.
- Ningún resultado se guarda; todo sale por `print` o en el título de un gráfico (las tasas Beta solo existen como título de figura, `:602`).
- `plt.savefig` antes de poner las etiquetas de los ejes (`:386`, `:642`); la figura Beta lleva el título "P-P Plot" (`:605`).
- Descarga datos de internet dentro del flujo principal (`:974`), con `except Exception` que oculta cualquier error (`:998`).

### Rendimiento

| Qué | Costo actual | Cómo debería ser |
|---|---|---|
| Ybus | Denso, nb² complejos (`:79`) | CSR disperso |
| Jacobiano | Ocho temporales densos nb × nb (`:110-118`) + `np.block` | Fórmula dispersa dS/dV |
| Factorización | LU denso O(n³) (`:131`) | `scipy.sparse.linalg.splu` |
| Sensibilidades | 2 sustituciones densas O(n²) por rama, sin caché (`:181-182`) | Una resolución dispersa con múltiples lados derechos |
| Pertenencia | `idx in self.pq`, `np.where(self.ns == idx)` (`:136`, `:149`, `:153`, `:155`, `:178`): búsquedas O(n) por llamada | Tabla de índices precalculada |
| Autovalores | `scipy.linalg.eigvalsh` por rama (`:194`) | Fórmula cerrada 2×2 vectorizada |
| Monte Carlo AC | `pandapower.runpp` por muestra, con `numba=False` (`:271`): 20–60 ms por flujo, casi todo sobrecosto | Resolución por lotes con Ybus y LU reutilizados |

---

## Correspondencia entre la revisión y el código

| Crítica | ¿La confirma el código? |
|---|---|
| Editor; R1.1, R1.4; R2.5; R3.1 — la "robustez de linealización" no mide la linealización | **Sí.** Tabla I usa `C @ eps` (`:530`). No hay ninguna comparación AC en el código. |
| R2.1; R3.2 — G es función solo de ρ | **Sí**, y más fuerte: con κ = π/4 el umbral de 15 % es inoperante. |
| R3.2 — tres "cajas" distintas | **Sí** (tabla de A1). |
| R2.3; R3.3 — λ_max no es la proyección por coordenada; salida ≠ elipse | **Sí.** Además el código ni siquiera calcula los intervalos de salida. |
| R2.2; R3.4 — Q = ⅓CCᵀ supone fuentes independientes | **Sí** (`:187`). No hay Σ_ε. |
| R3.5 — transformación Beta | **Sí**, exacta (`:580`). |
| R2.4; R3.6 — garantías determinista y probabilística mezcladas; "no branch is under-covered" | **Sí.** No hay análisis por barra ni a nivel sistema. |
| R3.7 — 1.246,66 vs 1.263,48; 24,0 vs 24,8; "sub-second" con 52 s; "width" vs "area" | **Sí.** Dos corridas distintas; el CRI es de área y condicional a Elipse. |
| R3.8 — complejidad O(N^1,2), `scipy.linalg` | **Sí.** El LU es denso. No hay estudio de escalamiento. |
| R2.7 — γ_th calibrado en IEEE-33 | **Sí**, y la calibración no puede haber distinguido 5, 10, 15 y 20 %. |
| R1.3 — comparación con [11] y [12] | No hay código de CLTM, ni de [11], ni de [12]. |
| R1.5 — δ no está definido | En el código δ es el semiancho relativo por inyección (A2). |
| R1.6; R2.6 — datos reales, más distribuciones | No hay nada en el código. |

La revisión no detectó (no podía, sin el código) los hallazgos 1, 2, 3, 6 y 8 del resumen.

---

## Corrida completa del script original

`finalFinalCode.py` ejecutado sin modificar, de principio a fin, con `MPLBACKEND=Agg` (para que `plt.show()` no bloquee) y desde un directorio temporal. La máquina estaba ocupada con otros procesos durante la corrida, así que **los tiempos de pared de esta tabla son cotas superiores, no mediciones limpias**. Las mediciones aisladas están en RESULTS.md (C7 y C8).

| Salida del script | IEEE-33 | PEGASE-2869 |
|---|---|---|
| Total Líneas Analizadas (PQ-PQ) | 31 | 3.218 |
| Líneas Críticas (Usan Elipse) | 16 (51,6 %) | 2.984 (92,7 %) |
| Reducción Promedio (Zona Crítica) | 86,70 % | 57,80 % |
| Reducción Máxima Detectada | 99,31 % | 100,00 % |
| Tiempo de Ejecución (Solver) | 0,0047 s | 58,10 s |
| Monte Carlo AC, 20.000 flujos | 961,39 s; 20.000 éxitos, 0 fallos | no terminó: detenido tras más de dos horas (ver RESULTS.md, C8) |
| Tabla I (5 / 10 / 20 / 30 %) | 0,00 / 0,00 / 0,02 / 0,04 % | 0,04 / 0,04 / 0,10 / 0,10 % |

Todas las cantidades que no son tiempos coinciden con las que da el paquete nuevo en modo *legacy* (RESULTS.md, sección 2), que es la comprobación de que el refactor no cambió el comportamiento.
