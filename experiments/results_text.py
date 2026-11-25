"""Commentary of RESULTS.md. Tables are injected by build_report.py at the {{placeholders}}."""

TEXT = r"""# RESULTS.md — Resultados de los experimentos CD-HIM

Generado por `experiments/build_report.py` a partir de los CSV de `results/`. Ninguna cifra de las tablas
está escrita a mano.

- **Entorno:** {{versions}}.
- **Semilla base:** {{seed}} (`config.yaml`). Cada experimento deriva la suya y la guarda en la columna `seed`
  del CSV y en el `.meta.json` contiguo.
- **Regenerar todo:** `python experiments/run_all.py` y luego `python experiments/build_report.py`.
- **Dos linealizaciones.** *legacy* es lo que hace `finalFinalCode.py` (Jacobiano en arranque plano, con el
  error de signo en ∂P/∂θ): son las cotas del manuscrito enviado. *corrected* es lo que el manuscrito dice
  hacer (Jacobiano exacto en el punto AC convergido) y es el valor por defecto del paquete nuevo. Salvo que
  se indique, los experimentos C usan *corrected*.
- **Qué es "la cota CD-HIM" por barra.** El Algoritmo 1 asigna un radio por rama; una barra que toca varias
  ramas en modo Elipse recibe varios. Aquí se usa el menor (`bus_combine: min`), que es el caso más exigente
  para la violación. Las barras sin ninguna rama en modo Elipse llevan la caja afín. Todo intervalo se
  centra en la solución AC convergida.
- **Violación.** Una muestra viola una cota si |ΔV| supera el radio en más de 10⁻⁸ p.u. (precisión del solver
  AC). Las muestras AC que no convergen se cuentan aparte y se excluyen del denominador.
- **Barras IEEE-33.** Los números de barra son los índices de pandapower (parten en 0), igual que en la
  Fig. 4 del manuscrito.

---

## 1. Correspondencia con los `[TODO]` del manuscrito

| `[TODO]` del paper | Sale de | Valor para el paper |
|---|---|---|
| Caja de A_box, κ y ρ_th | A1 (AUDIT.md) | Caja de covarianza con semianchos √Q_ii, a la misma escala que la elipse (sin χ²). **κ = π/4 = 0,7854**. **No existe ρ_th**: G ≥ 21,46 % para todo ρ, así que γ_th = 15 % no filtra nada y la selección la hace solo \|ρ\| < 1 − δ_stab. Si se redefine G con la elipse χ² frente a la caja 3σ, κ = 1,032 y ρ_th = 0,567 (ver C4). |
| R_viol de PEGASE contra AC | B1 | Con el código corregido, δ = 5 %, N = 20.000: **por barra, media 0,0024 %, P99 0,047 %, máximo 0,41 %** (1 barra de 2.359 sobre 0,27 %; es una barra en modo Caja). **Por rama (elipse), media 0,103 %, P99 0,265 %, máximo 2,69 %** (22 de 2.543 sobre 0,27 %). Con el código original las cotas no contienen las muestras AC: media 3,69 %, máximo 89,4 %. |
| Símbolos de P_d y Q_d | A2 (AUDIT.md) | **Separados e independientes**: un símbolo por P_d y otro por Q_d de cada barra con carga. N_u = 64 (IEEE-33) y 2.951 (PEGASE-2869: 1.485 de P y 1.466 de Q). δ es el semiancho relativo, igual para todas las fuentes. |
| Rutina de LU | A3 (AUDIT.md) | El código original usa **LU denso** (`scipy.linalg.lu_factor`, LAPACK `getrf`). El paquete nuevo usa **`scipy.sparse.linalg.splu` (SuperLU)** con ordenamiento COLAMD; relleno medido en C7. |
| Reducción de ancho y sobreconservadurismo de λ_max | B2 | Intervalos de salida frente a la caja afín, barras en modo Elipse: **−42,8 % (IEEE-33) y −34,1 % (PEGASE)** con λ_max; −56,2 % y −64,5 % con la proyección exacta. Factor √(λ_max/Q_ii): **media 1,57 / mediana 1,41 / máx. 6,28 (IEEE-33); media 2,26 / mediana 1,40 / máx. 136 (PEGASE)**. Frente a la caja 3σ el intervalo de salida es *más ancho*: +50,6 % y +110 % con λ_max, +14,6 % con la proyección exacta. |

---

## 2. Antes y después: ¿se reproducen las cifras del paper?

La columna *legacy* corre el paquete nuevo emulando el código original (incluida la regularización de 10⁻⁹,
las semillas 99/42/101 y los tamaños 5.000/10.000/20.000). Coincide con lo que imprime `finalFinalCode.py`
sin modificar en todas las cantidades de la tabla (verificado contra la corrida del script original; ver AUDIT.md).

{{reproduction}}

*Script: `experiments/a_reproduce_paper.py` · CSV: `results/a_reproduction.csv`*

Lectura:

- El refactor no cambia el comportamiento: en modo *legacy* salen los mismos conteos (16/31, 2.984/3.218),
  la misma rama de la Fig. 4 y las mismas tasas de la Tabla I, Beta y P-P que con el script original.
- El código del repositorio da CRI = 86,70 % y 57,80 %, no 86,92 % y 58,18 %. La diferencia ya existe en el
  script original; no la introduce el refactor. Causa no identificada (ver AUDIT.md, A4).
- La columna central aísla los dos errores: corregir solo el signo cambia la selección; evaluar además en
  el punto convergido cambia sobre todo PEGASE.
- Sin la regularización de 10⁻⁹ la Tabla I es constante en δ y la prueba Beta del paper da 0,8–0,9 %, no 0,19 %.

---

## 3. B1 — Violación de las cotas contra AC no lineal (punto de diseño)

δ = 15 % (IEEE-33) y δ = 5 % (PEGASE-2869), N = 20.000 flujos AC por sistema, las mismas muestras para las
tres linealizaciones. Ninguna muestra dejó de converger.

**Por barra, contra el intervalo de salida (λ_max):**

{{b1_main}}

**Por rama, contra la elipse (distancia de Mahalanobis), solo ramas en modo Elipse:**

{{b1_branch}}

**Otras cotas, linealización corregida, contra AC:**

{{b1_bounds}}

*Script: `experiments/b1_ac_violation.py` · CSV: `results/b1_bus_violation.csv`, `results/b1_branch_violation.csv` · N = 20.000 · semilla en la columna `seed`*

Lectura:

- **IEEE-33: el 0,00 % del paper se confirma.** Con las cotas originales, 1 muestra de 20.000 sale del intervalo en
  una barra (0,005 %). Con las corregidas, el máximo por barra es 0,020 %.
- **PEGASE-2869: las cotas del paper no contienen la física AC.** Con el código original, 867 de 2.359 barras
  superan 0,27 %, el máximo es 89,4 % y todas las muestras tienen al menos una barra fuera. Incluso la caja afín
  "determinista" del código original se viola (media 1,5 %, máximo 80 %). La causa es el Jacobiano: corregir
  solo el signo no basta (fila *flat_exact*); el arranque plano por sí solo rompe la linealización en una red con
  voltajes entre 0,925 y 1,068 p.u.
- **Con el Jacobiano exacto en el punto convergido las cotas sí se sostienen por barra:** media 0,0024 %,
  P99 0,047 %. La única barra sobre 0,27 % (0,41 %) está en modo Caja: es una barra casi insensible, donde el
  error de segundo orden supera un radio afín muy pequeño.
- **La elipse es más frágil que el intervalo.** En PEGASE, 22 elipses superan 0,27 % contra AC y ninguna en el
  modelo lineal: son ramas con el eje menor muy fino, donde el error de linealización es comparable a √(χ²λ_min).
- **Nivel sistema:** en PEGASE el 3,1 % de las muestras tiene al menos una barra fuera y el 49,9 % al menos una
  elipse violada (48,3 % ya en el modelo lineal). Ver C5.
- La caja 3σ ("CLTM Box") viola más que CD-HIM, como corresponde a un intervalo más angosto: 80 barras sobre
  0,27 % en PEGASE.

---

## 4. B2 — Conservadurismo de la salida real

Determinista (sin muestreo). "Reducción de ancho" = promedio sobre barras de 1 − r_salida/r_referencia; un
valor negativo significa que la salida es más ancha que la referencia.

{{b2}}

*Script: `experiments/b2_conservatism.py` · CSV: `results/b2_conservatism.csv`*

Lectura:

- **El CRI del paper no mide lo que el método entrega.** El 58,18 % (ahora 71,1 % con el código corregido) es
  reducción de área de la elipse frente a la caja de covarianza. El intervalo de salida, que es lo que devuelve el
  Algoritmo 1, se reduce 34,1 % frente a la caja afín en PEGASE y 42,8 % en IEEE-33.
- **λ_max cuesta caro.** La mediana de √(λ_max/Q_ii) es 1,41 ≈ √2 en ambos sistemas: con dos barras muy
  correlacionadas y de varianza parecida, λ_max ≈ σ_i² + σ_j². Es decir, la regla infla el radio un 41 % justo en
  las ramas donde la elipse debía ayudar. Con varianzas desiguales el factor crece sin cota (136 en PEGASE).
- **En 140 barras de PEGASE (427 pares rama-extremo) el radio λ_max es más ancho que la caja afín**, que ya es
  una envolvente determinista del modelo lineal. Con la proyección exacta pasa en 50 barras (el radio 3,44σ
  supera Σ|c_k| cuando domina una sola fuente). Recomendación: activar `clip_to_affine_box`, que no cuesta
  nada y nunca pierde cobertura en el modelo lineal.
- **Frente a la caja 3σ el método siempre es más ancho.** La proyección exacta de la elipse conjunta al 99,73 %
  es √11,83 σ = 3,44σ, un 14,6 % más que 3σ; con λ_max, 50,6 % y 110 % más. La ventaja sobre "CLTM" que
  muestra la Fig. 4 existe solo en área bidimensional, no en los intervalos por barra.

---

## 5. B3 — Tabla I regenerada

Muestras uniformes independientes propagadas por el modelo lineal (ΔX = Cε), N = 200.000, mismas muestras para
todos los δ. La rama es la del paper (mayor G con |ρ| < 0,999).

{{b3}}

**Todas las ramas en modo Elipse (modelo lineal, entradas uniformes):**

{{b3_all}}

*Script: `experiments/b3_table1.py` · CSV: `results/b3_table1.csv`, `results/b3_all_branches.csv` · N = 200.000 · semilla en la columna `seed`*

Lectura:

- **La violación exacta no depende de δ**: las cuatro filas de cada bloque son idénticas y d² cambia menos de
  10⁻¹² entre niveles. La Tabla I no informa nada sobre la validez de la linealización.
- **La tendencia creciente del paper es la regularización**: con el 10⁻⁹ del código original la tasa sube con δ.
- Valor real de la celda única: 0,107 % (IEEE-33) y 0,070 % (PEGASE) con el código corregido.
- Con entradas uniformes independientes ninguna rama supera 0,27 % en el modelo lineal (máximo 0,148 % y
  0,261 %): la suma de uniformes tiene colas más livianas que la normal, así que la elipse χ² es conservadora.
  Ese margen es del modelo lineal con esas entradas; no cubre el error de linealización (C1) ni otras
  distribuciones (C3).

---

<<C_SECTIONS>>
"""
