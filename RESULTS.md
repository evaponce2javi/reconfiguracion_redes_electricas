# RESULTS.md — Resultados de los experimentos CD-HIM

Generado por `experiments/build_report.py` a partir de los CSV de `results/`. Ninguna cifra de las tablas
está escrita a mano.

- **Entorno:** Python 3.12.10, pandapower 3.5.6, NumPy 2.4.6, SciPy 1.18.1 (Windows-11-10.0.26200-SP0).
- **Semilla base:** 20261114 (`config.yaml`). Cada experimento deriva la suya y la guarda en la columna `seed`
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

**IEEE-33**

| Cantidad | Paper | Paquete nuevo, modo *legacy* | Signo corregido, arranque plano | Corregido |
|---|---|---|---|---|
| Ramas PQ-PQ analizadas | 31 | 31 | 31 | 31 |
| Ramas en modo Elipse | 16 | 16 | 21 | 21 |
| Ramas en modo Elipse (%) | 51,60 | 51,61 | 67,74 | 67,74 |
| CRI: G medio en ramas Elipse (%) | 86,92 | 86,70 | 85,88 | 86,37 |
| Ramas con \|ρ\| ≥ 0,999 (a Caja) | – | 15 | 10 | 10 |
| Rama de la Fig. 4 | – | 5-6 | 25-26 | 25-26 |
| ρ de esa rama | 0,9990 | 0,9990 | 0,9988 | 0,9990 |
| G de esa rama (%) | 96,45 | 96,45 | 96,15 | 96,40 |
| Tabla I, δ = 5 % (%) | 0,00 | 0,00 | 0,10 | 0,10 |
| Tabla I, δ = 10 % (%) | 0,06 | 0,00 | 0,10 | 0,10 |
| Tabla I, δ = 20 % (%) | 0,08 | 0,02 | 0,10 | 0,10 |
| Tabla I, δ = 30 % (%) | 0,10 | 0,04 | 0,10 | 0,10 |
| Beta: máximo de las 3 ramas (%) | 0,40 | 0,40 | 0,84 | 0,83 |
| Beta: ramas probadas que están en modo Caja | – | 3 | 3 | 3 |
| P-P: cobertura a confianza 0,50 | – | 0,4845 | 0,4878 | 0,4860 |

**PEGASE-2869**

| Cantidad | Paper | Paquete nuevo, modo *legacy* | Signo corregido, arranque plano | Corregido |
|---|---|---|---|---|
| Ramas PQ-PQ analizadas | 3218 | 3218 | 3216 | 3216 |
| Ramas en modo Elipse | – | 2984 | 2614 | 2543 |
| Ramas en modo Elipse (%) | 92,70 | 92,73 | 81,28 | 79,07 |
| CRI: G medio en ramas Elipse (%) | 58,18 | 57,80 | 67,07 | 71,08 |
| Ramas con \|ρ\| ≥ 0,999 (a Caja) | – | 234 | 602 | 673 |
| Rama de la Fig. 4 | – | 7484-1561 | 3928-9172 | 6150-6312 |
| ρ de esa rama | – | 0,9990 | 0,9990 | 0,9990 |
| G de esa rama (%) | – | 96,47 | 96,48 | 96,48 |
| Tabla I, δ = 5 % (%) | 0,06 | 0,04 | 0,02 | 0,10 |
| Tabla I, δ = 10 % (%) | 0,08 | 0,04 | 0,02 | 0,10 |
| Tabla I, δ = 20 % (%) | 0,12 | 0,10 | 0,02 | 0,10 |
| Tabla I, δ = 30 % (%) | 0,12 | 0,10 | 0,02 | 0,10 |
| Beta: máximo de las 3 ramas (%) | 0,19 | 0,23 | 0,98 | 0,91 |
| Beta: ramas probadas que están en modo Caja | – | 3 | 3 | 3 |
| P-P: cobertura a confianza 0,50 | – | 0,7312 | 0,4556 | 0,4469 |

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

| Sistema | Linealización | δ | Barras PQ | Media (%) | P99 (%) | Máx. (%) | Barras > 0,27 % | Muestras con ≥ 1 violación (%) | AC sin converger |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| IEEE-33 | legacy | 0,15 | 32 | 0,0002 | 0,003 | 0,005 | 0 | 0,01 | 0 |
| IEEE-33 | flat_exact | 0,15 | 32 | 0,0053 | 0,089 | 0,125 | 0 | 0,13 | 0 |
| IEEE-33 | corrected | 0,15 | 32 | 0,0006 | 0,014 | 0,020 | 0 | 0,02 | 0 |
| PEGASE-2869 | legacy | 0,05 | 2359 | 3,6932 | 61,641 | 89,420 | 867 | 100,00 | 0 |
| PEGASE-2869 | flat_exact | 0,05 | 2359 | 2,4169 | 49,413 | 88,525 | 634 | 100,00 | 0 |
| PEGASE-2869 | corrected | 0,05 | 2359 | 0,0024 | 0,047 | 0,410 | 1 | 3,15 | 0 |

**Por rama, contra la elipse (distancia de Mahalanobis), solo ramas en modo Elipse:**

| Sistema | Linealización | Muestras | Ramas Elipse | Media (%) | P99 (%) | Máx. (%) | Ramas > 0,27 % | Muestras con ≥ 1 elipse violada (%) |
|---|---|---|---:|---:|---:|---:|---:|---:|
| IEEE-33 | legacy | AC no lineal | 16 | 0,126 | 0,477 | 0,490 | 3 | 1,55 |
| IEEE-33 | legacy | modelo lineal | 16 | 0,024 | 0,054 | 0,055 | 0 | 0,27 |
| IEEE-33 | flat_exact | AC no lineal | 21 | 0,486 | 0,894 | 0,905 | 18 | 3,25 |
| IEEE-33 | flat_exact | modelo lineal | 21 | 0,097 | 0,172 | 0,175 | 0 | 0,98 |
| IEEE-33 | corrected | AC no lineal | 21 | 0,096 | 0,164 | 0,165 | 0 | 0,99 |
| IEEE-33 | corrected | modelo lineal | 21 | 0,097 | 0,165 | 0,165 | 0 | 0,99 |
| PEGASE-2869 | legacy | AC no lineal | 2984 | 13,104 | 87,910 | 97,515 | 2578 | 100,00 |
| PEGASE-2869 | legacy | modelo lineal | 2984 | 0,054 | 0,190 | 0,235 | 0 | 42,90 |
| PEGASE-2869 | flat_exact | AC no lineal | 2614 | 11,551 | 86,994 | 95,475 | 2253 | 100,00 |
| PEGASE-2869 | flat_exact | modelo lineal | 2614 | 0,075 | 0,215 | 0,250 | 0 | 47,71 |
| PEGASE-2869 | corrected | AC no lineal | 2543 | 0,103 | 0,265 | 2,685 | 22 | 49,89 |
| PEGASE-2869 | corrected | modelo lineal | 2543 | 0,098 | 0,245 | 0,270 | 0 | 48,27 |

**Otras cotas, linealización corregida, contra AC:**

| Sistema | Cota | Barras | n | Media (%) | P99 (%) | Máx. (%) | Barras > 0,27 % |
|---|---|---|---:|---:|---:|---:|---:|
| IEEE-33 | cdhim_lmax | todas las PQ | 32 | 0,0006 | 0,014 | 0,020 | 0 |
| IEEE-33 | cdhim_lmax | solo modo Elipse | 26 | 0,0008 | 0,015 | 0,020 | 0 |
| IEEE-33 | cdhim_coord | todas las PQ | 32 | 0,0116 | 0,045 | 0,050 | 0 |
| IEEE-33 | cdhim_coord | solo modo Elipse | 26 | 0,0142 | 0,046 | 0,050 | 0 |
| IEEE-33 | affine_box | todas las PQ | 32 | 0,0000 | 0,000 | 0,000 | 0 |
| IEEE-33 | affine_box | solo modo Elipse | 26 | 0,0000 | 0,000 | 0,000 | 0 |
| IEEE-33 | cltm_3sigma | todas las PQ | 32 | 0,1512 | 0,260 | 0,275 | 1 |
| IEEE-33 | cltm_3sigma | solo modo Elipse | 26 | 0,1471 | 0,263 | 0,275 | 1 |
| PEGASE-2869 | cdhim_lmax | todas las PQ | 2359 | 0,0024 | 0,047 | 0,410 | 1 |
| PEGASE-2869 | cdhim_lmax | solo modo Elipse | 1766 | 0,0030 | 0,050 | 0,105 | 0 |
| PEGASE-2869 | cdhim_coord | todas las PQ | 2359 | 0,0153 | 0,070 | 0,630 | 2 |
| PEGASE-2869 | cdhim_coord | solo modo Elipse | 1766 | 0,0202 | 0,070 | 0,630 | 1 |
| PEGASE-2869 | affine_box | todas las PQ | 2359 | 0,0002 | 0,000 | 0,410 | 1 |
| PEGASE-2869 | affine_box | solo modo Elipse | 1766 | 0,0000 | 0,000 | 0,000 | 0 |
| PEGASE-2869 | cltm_3sigma | todas las PQ | 2359 | 0,1333 | 0,295 | 1,190 | 80 |
| PEGASE-2869 | cltm_3sigma | solo modo Elipse | 1766 | 0,1374 | 0,292 | 1,150 | 61 |

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

| Cantidad | IEEE-33, legacy | IEEE-33, corrected | PEGASE-2869, legacy | PEGASE-2869, corrected |
|---|---|---|---|---|
| Barras en modo Elipse | 20 | 26 | 2081 | 1766 |
| Barras PQ | 32 | 32 | 2359 | 2359 |
| √(λ_max/Q_ii) por (rama, extremo): media | 1,578 | 1,574 | 1,935 | 2,262 |
| √(λ_max/Q_ii) por (rama, extremo): mediana | 1,414 | 1,414 | 1,385 | 1,402 |
| √(λ_max/Q_ii) por (rama, extremo): máximo | 6,31 | 6,28 | 132,56 | 136,02 |
| Ídem por barra (radio mínimo): media | 1,327 | 1,313 | 1,437 | 1,828 |
| Ídem por barra: mediana | 1,338 | 1,340 | 1,218 | 1,269 |
| Ídem por barra: máximo | 1,69 | 1,48 | 132,56 | 136,02 |
| Reducción de ancho, λ_max vs caja afín, barras Elipse (%) | 9,00 | 42,78 | 11,05 | 34,12 |
| Reducción de ancho, λ_max vs caja afín, todas las PQ (%) | 5,63 | 34,76 | 9,75 | 25,62 |
| Reducción de ancho, proyección exacta vs caja afín, barras Elipse (%) | 31,49 | 56,20 | 39,48 | 64,51 |
| Reducción de ancho, proyección exacta vs caja afín, todas las PQ (%) | 19,68 | 45,67 | 34,85 | 48,44 |
| Reducción de ancho, λ_max vs caja 3σ, barras Elipse (%) | -52,12 | -50,55 | -64,75 | -109,59 |
| Reducción de ancho, proyección exacta vs caja 3σ, barras Elipse (%) | -14,64 | -14,64 | -14,64 | -14,64 |
| Reducción de área del rectángulo de salida vs caja afín, λ_max (%) | 22,07 | 68,56 | 33,16 | 66,41 |
| CRI del paper (área de elipse vs caja de covarianza) (%) | 86,70 | 86,37 | 57,80 | 71,08 |
| Barras donde λ_max es más ancho que la caja afín | 4 | 0 | 338 | 140 |
| (rama, extremo) donde λ_max es más ancho que la caja afín | 13 | 2 | 1506 | 427 |
| Total de (rama, extremo) | 32 | 42 | 5968 | 5086 |
| Barras donde la proyección exacta es más ancha que la caja afín | 0 | 0 | 128 | 50 |
| Caja afín / caja 3σ (mediana) | 1,71 | 2,90 | 2,07 | 3,95 |

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

| Sistema | Linealización | Rama | ρ | δ | Violación exacta (%) | Con regularización 1e-9 (%) | Error estándar (%) | Cambio rel. máx. de d² entre δ |
|---|---|---|---:|---:|---:|---:|---:|---:|
| IEEE-33 | legacy | 5-6 | 0,99898 | 0,05 | 0,0380 | 0,0000 | 0,0044 | 1,0e-13 |
| IEEE-33 | legacy | 5-6 | 0,99898 | 0,10 | 0,0380 | 0,0000 | 0,0044 | 1,0e-13 |
| IEEE-33 | legacy | 5-6 | 0,99898 | 0,20 | 0,0380 | 0,0090 | 0,0044 | 1,0e-13 |
| IEEE-33 | legacy | 5-6 | 0,99898 | 0,30 | 0,0380 | 0,0205 | 0,0044 | 1,0e-13 |
| IEEE-33 | corrected | 25-26 | 0,99895 | 0,05 | 0,1070 | 0,0105 | 0,0073 | 1,8e-13 |
| IEEE-33 | corrected | 25-26 | 0,99895 | 0,10 | 0,1070 | 0,0175 | 0,0073 | 1,8e-13 |
| IEEE-33 | corrected | 25-26 | 0,99895 | 0,20 | 0,1070 | 0,0310 | 0,0073 | 1,8e-13 |
| IEEE-33 | corrected | 25-26 | 0,99895 | 0,30 | 0,1070 | 0,0480 | 0,0073 | 1,8e-13 |
| PEGASE-2869 | legacy | 7484-1561 | 0,99899 | 0,05 | 0,0925 | 0,0305 | 0,0068 | 9,3e-14 |
| PEGASE-2869 | legacy | 7484-1561 | 0,99899 | 0,10 | 0,0925 | 0,0350 | 0,0068 | 9,3e-14 |
| PEGASE-2869 | legacy | 7484-1561 | 0,99899 | 0,20 | 0,0925 | 0,0595 | 0,0068 | 9,3e-14 |
| PEGASE-2869 | legacy | 7484-1561 | 0,99899 | 0,30 | 0,0925 | 0,0710 | 0,0068 | 9,3e-14 |
| PEGASE-2869 | corrected | 6150-6312 | 0,99900 | 0,05 | 0,0695 | 0,0195 | 0,0059 | 1,6e-14 |
| PEGASE-2869 | corrected | 6150-6312 | 0,99900 | 0,10 | 0,0695 | 0,0260 | 0,0059 | 1,6e-14 |
| PEGASE-2869 | corrected | 6150-6312 | 0,99900 | 0,20 | 0,0695 | 0,0405 | 0,0059 | 1,6e-14 |
| PEGASE-2869 | corrected | 6150-6312 | 0,99900 | 0,30 | 0,0695 | 0,0555 | 0,0059 | 1,6e-14 |

**Todas las ramas en modo Elipse (modelo lineal, entradas uniformes):**

| Sistema | Linealización | Ramas Elipse | Media (%) | Mediana (%) | P99 (%) | Máx. (%) | Ramas > 0,27 % |
|---|---|---:|---:|---:|---:|---:|---:|
| IEEE-33 | legacy | 16 | 0,0225 | 0,0198 | 0,0618 | 0,0630 | 0 |
| IEEE-33 | corrected | 21 | 0,0948 | 0,0975 | 0,1470 | 0,1480 | 0 |
| PEGASE-2869 | legacy | 2984 | 0,0554 | 0,0445 | 0,1805 | 0,2225 | 0 |
| PEGASE-2869 | corrected | 2543 | 0,1018 | 0,1005 | 0,2410 | 0,2610 | 0 |

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

## 6. C1 — Validación contra AC no lineal: carga × incertidumbre

Carga: P_d, Q_d y la generación programada se multiplican por el factor de carga; la barra slack cubre las
pérdidas. Límite de cargabilidad por continuación: **3,62 (IEEE-33) y 1,59 (PEGASE-2869)**; el quinto nivel
es el 95 % de ese límite. N = 20.000 flujos AC por celda, las mismas muestras de ruido en todas las celdas
de un sistema. En cada celda: (a) error de linealización |V_AC − V_lineal| por barra, en p.u. y relativo al
radio del intervalo de salida; (b) miscobertura de la elipse dentro del modelo lineal; (c) violación total
de las muestras AC.

**IEEE-33**

| Carga | δ | No conv. NR (%) | (a) error máx. (p.u.) | (a) error P99 (p.u.) | (a) error/radio P99 | (a) error/radio máx. | (b) elipse lineal, media (%) | (b) máx. (%) | (c) barra AC, media (%) | (c) barra AC, máx. (%) | (c) barras > 0,27 % | (c) sistema (%) | (c) elipse AC, media (%) | (c) elipse AC, máx. (%) | (c) elipses > 0,27 % |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0,800 | 0,05 | 0,00 | 3,18e-06 | 1,22e-06 | 0,0007 | 0,002 | 0,078 | 0,130 | 0,0005 | 0,010 | 0 | 0,01 | 0,079 | 0,13 | 0 |
| 0,800 | 0,10 | 0,00 | 1,28e-05 | 4,84e-06 | 0,0014 | 0,003 | 0,078 | 0,130 | 0,0005 | 0,010 | 0 | 0,01 | 0,079 | 0,13 | 0 |
| 0,800 | 0,15 | 0,00 | 2,89e-05 | 1,10e-05 | 0,0021 | 0,005 | 0,078 | 0,130 | 0,0005 | 0,010 | 0 | 0,01 | 0,079 | 0,13 | 0 |
| 0,800 | 0,20 | 0,00 | 5,15e-05 | 1,95e-05 | 0,0028 | 0,007 | 0,078 | 0,130 | 0,0005 | 0,010 | 0 | 0,01 | 0,079 | 0,13 | 0 |
| 0,800 | 0,30 | 0,00 | 1,17e-04 | 4,37e-05 | 0,0042 | 0,010 | 0,078 | 0,130 | 0,0005 | 0,010 | 0 | 0,01 | 0,080 | 0,13 | 0 |
| 1,000 | 0,05 | 0,00 | 5,55e-06 | 2,09e-06 | 0,0009 | 0,002 | 0,078 | 0,130 | 0,0005 | 0,010 | 0 | 0,01 | 0,078 | 0,13 | 0 |
| 1,000 | 0,10 | 0,00 | 2,23e-05 | 8,41e-06 | 0,0018 | 0,005 | 0,078 | 0,130 | 0,0005 | 0,010 | 0 | 0,01 | 0,078 | 0,13 | 0 |
| 1,000 | 0,15 | 0,00 | 5,04e-05 | 1,88e-05 | 0,0028 | 0,007 | 0,078 | 0,130 | 0,0005 | 0,010 | 0 | 0,01 | 0,079 | 0,13 | 0 |
| 1,000 | 0,20 | 0,00 | 9,01e-05 | 3,35e-05 | 0,0037 | 0,009 | 0,078 | 0,130 | 0,0005 | 0,010 | 0 | 0,01 | 0,079 | 0,13 | 0 |
| 1,000 | 0,30 | 0,00 | 2,05e-04 | 7,59e-05 | 0,0055 | 0,014 | 0,078 | 0,130 | 0,0005 | 0,010 | 0 | 0,01 | 0,079 | 0,12 | 0 |
| 1,200 | 0,05 | 0,00 | 8,99e-06 | 3,35e-06 | 0,0012 | 0,003 | 0,078 | 0,130 | 0,0005 | 0,010 | 0 | 0,01 | 0,078 | 0,13 | 0 |
| 1,200 | 0,10 | 0,00 | 3,62e-05 | 1,35e-05 | 0,0023 | 0,006 | 0,078 | 0,130 | 0,0005 | 0,010 | 0 | 0,01 | 0,079 | 0,13 | 0 |
| 1,200 | 0,15 | 0,00 | 8,20e-05 | 3,02e-05 | 0,0035 | 0,009 | 0,078 | 0,130 | 0,0005 | 0,010 | 0 | 0,01 | 0,079 | 0,12 | 0 |
| 1,200 | 0,20 | 0,00 | 1,47e-04 | 5,37e-05 | 0,0047 | 0,012 | 0,078 | 0,130 | 0,0005 | 0,010 | 0 | 0,01 | 0,080 | 0,13 | 0 |
| 1,200 | 0,30 | 0,00 | 3,35e-04 | 1,22e-04 | 0,0071 | 0,018 | 0,078 | 0,130 | 0,0005 | 0,010 | 0 | 0,01 | 0,080 | 0,13 | 0 |
| 1,400 | 0,05 | 0,00 | 1,39e-05 | 5,13e-06 | 0,0014 | 0,004 | 0,075 | 0,130 | 0,0005 | 0,010 | 0 | 0,01 | 0,075 | 0,13 | 0 |
| 1,400 | 0,10 | 0,00 | 5,61e-05 | 2,07e-05 | 0,0029 | 0,007 | 0,075 | 0,130 | 0,0005 | 0,010 | 0 | 0,01 | 0,076 | 0,12 | 0 |
| 1,400 | 0,15 | 0,00 | 1,27e-04 | 4,62e-05 | 0,0044 | 0,011 | 0,075 | 0,130 | 0,0005 | 0,010 | 0 | 0,01 | 0,077 | 0,13 | 0 |
| 1,400 | 0,20 | 0,00 | 2,28e-04 | 8,22e-05 | 0,0058 | 0,015 | 0,075 | 0,130 | 0,0005 | 0,010 | 0 | 0,01 | 0,078 | 0,13 | 0 |
| 1,400 | 0,30 | 0,00 | 5,23e-04 | 1,86e-04 | 0,0087 | 0,023 | 0,075 | 0,130 | 0,0005 | 0,010 | 0 | 0,01 | 0,081 | 0,14 | 0 |
| 3,441 | 0,05 | 0,00 | 4,61e-03 | 1,17e-03 | 0,0380 | 0,161 | 0,074 | 0,120 | 0,0025 | 0,055 | 0 | 0,06 | 0,143 | 0,18 | 0 |
| 3,441 | 0,10 | 0,00 | 3,58e-02 | 5,31e-03 | 0,0881 | 0,618 | 0,074 | 0,120 | 0,0175 | 0,195 | 0 | 0,20 | 0,373 | 0,52 | 13 |
| 3,441 | 0,15 | 0,28 | 5,71e-02 | 1,46e-02 | 0,1660 | 0,649 | 0,074 | 0,120 | 0,0470 | 0,471 | 1 | 0,48 | 0,776 | 1,09 | 15 |
| 3,441 | 0,20 | 2,37 | 6,40e-02 | 2,34e-02 | 0,1928 | 0,536 | 0,074 | 0,120 | 0,0048 | 0,138 | 0 | 0,14 | 0,453 | 0,82 | 13 |
| 3,441 | 0,30 | 10,84 | 6,83e-02 | 3,16e-02 | 0,1679 | 0,354 | 0,074 | 0,120 | 0,0000 | 0,000 | 0 | 0,00 | 0,086 | 0,30 | 3 |

**PEGASE-2869**

| Carga | δ | No conv. NR (%) | (a) error máx. (p.u.) | (a) error P99 (p.u.) | (a) error/radio P99 | (a) error/radio máx. | (b) elipse lineal, media (%) | (b) máx. (%) | (c) barra AC, media (%) | (c) barra AC, máx. (%) | (c) barras > 0,27 % | (c) sistema (%) | (c) elipse AC, media (%) | (c) elipse AC, máx. (%) | (c) elipses > 0,27 % |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0,800 | 0,05 | 0,00 | 1,46e-04 | 1,50e-05 | 0,0123 | 0,479 | 0,095 | 0,330 | 0,0022 | 0,270 | 0 | 2,97 | 0,099 | 1,71 | 27 |
| 0,800 | 0,10 | 0,00 | 5,75e-04 | 5,96e-05 | 0,0245 | 0,955 | 0,095 | 0,330 | 0,0029 | 0,670 | 2 | 3,86 | 0,110 | 5,90 | 57 |
| 0,800 | 0,15 | 0,00 | 1,27e-03 | 1,35e-04 | 0,0367 | 1,430 | 0,095 | 0,330 | 0,0040 | 0,995 | 3 | 5,20 | 0,127 | 11,01 | 124 |
| 0,800 | 0,20 | 0,00 | 2,23e-03 | 2,40e-04 | 0,0490 | 1,903 | 0,095 | 0,330 | 0,0056 | 1,720 | 6 | 7,03 | 0,148 | 15,90 | 208 |
| 0,800 | 0,30 | 0,00 | 5,12e-03 | 5,37e-04 | 0,0741 | 2,845 | 0,095 | 0,330 | 0,0101 | 3,550 | 12 | 10,87 | 0,202 | 23,70 | 358 |
| 1,000 | 0,05 | 0,00 | 2,74e-04 | 2,79e-05 | 0,0140 | 0,530 | 0,099 | 0,330 | 0,0023 | 0,365 | 1 | 3,03 | 0,104 | 2,62 | 33 |
| 1,000 | 0,10 | 0,00 | 1,07e-03 | 1,11e-04 | 0,0279 | 1,058 | 0,099 | 0,330 | 0,0034 | 0,915 | 2 | 4,28 | 0,118 | 8,29 | 77 |
| 1,000 | 0,15 | 0,00 | 2,38e-03 | 2,51e-04 | 0,0417 | 1,583 | 0,099 | 0,330 | 0,0050 | 1,360 | 3 | 6,04 | 0,139 | 14,61 | 172 |
| 1,000 | 0,20 | 0,00 | 4,38e-03 | 4,47e-04 | 0,0562 | 2,107 | 0,099 | 0,330 | 0,0072 | 2,015 | 6 | 8,15 | 0,165 | 19,86 | 281 |
| 1,000 | 0,30 | 0,00 | 1,07e-02 | 1,00e-03 | 0,0841 | 3,160 | 0,099 | 0,330 | 0,0131 | 4,090 | 18 | 12,88 | 0,234 | 28,52 | 472 |
| 1,200 | 0,05 | 0,00 | 5,06e-04 | 5,25e-05 | 0,0160 | 0,521 | 0,103 | 0,330 | 0,0026 | 0,595 | 1 | 3,29 | 0,109 | 3,23 | 40 |
| 1,200 | 0,10 | 0,00 | 2,02e-03 | 2,11e-04 | 0,0320 | 1,039 | 0,103 | 0,330 | 0,0040 | 1,200 | 3 | 4,91 | 0,126 | 10,07 | 118 |
| 1,200 | 0,15 | 0,00 | 4,80e-03 | 4,73e-04 | 0,0479 | 1,555 | 0,103 | 0,330 | 0,0061 | 1,760 | 4 | 7,02 | 0,152 | 16,84 | 246 |
| 1,200 | 0,20 | 0,00 | 9,08e-03 | 8,41e-04 | 0,0646 | 2,068 | 0,103 | 0,330 | 0,0091 | 2,290 | 9 | 9,43 | 0,185 | 22,52 | 382 |
| 1,200 | 0,30 | 0,00 | 2,38e-02 | 1,93e-03 | 0,0966 | 3,125 | 0,103 | 0,330 | 0,0177 | 4,265 | 26 | 15,13 | 0,276 | 31,55 | 655 |
| 1,400 | 0,05 | 0,00 | 1,15e-03 | 1,22e-04 | 0,0193 | 0,437 | 0,107 | 0,325 | 0,0030 | 0,800 | 1 | 3,68 | 0,115 | 3,44 | 56 |
| 1,400 | 0,10 | 0,00 | 5,09e-03 | 4,90e-04 | 0,0389 | 0,869 | 0,107 | 0,325 | 0,0050 | 1,585 | 3 | 5,54 | 0,139 | 10,43 | 189 |
| 1,400 | 0,15 | 0,00 | 1,30e-02 | 1,11e-03 | 0,0589 | 1,298 | 0,107 | 0,325 | 0,0084 | 2,290 | 5 | 8,06 | 0,176 | 17,23 | 410 |
| 1,400 | 0,20 | 0,00 | 2,70e-02 | 1,97e-03 | 0,0794 | 1,724 | 0,107 | 0,325 | 0,0134 | 3,020 | 14 | 11,18 | 0,228 | 23,31 | 615 |
| 1,400 | 0,30 | 0,01 | 1,40e-01 | 4,57e-03 | 0,1230 | 2,656 | 0,107 | 0,325 | 0,0321 | 4,436 | 73 | 18,86 | 0,388 | 32,32 | 1033 |
| 1,511 | 0,05 | 0,00 | 3,94e-03 | 2,92e-04 | 0,0251 | 0,370 | 0,109 | 0,315 | 0,0035 | 0,940 | 1 | 3,88 | 0,122 | 3,19 | 76 |
| 1,511 | 0,10 | 0,00 | 2,39e-02 | 1,19e-03 | 0,0513 | 0,735 | 0,109 | 0,315 | 0,0075 | 1,825 | 3 | 6,16 | 0,166 | 9,71 | 359 |
| 1,511 | 0,15 | 0,12 | 8,88e-02 | 2,75e-03 | 0,0794 | 1,096 | 0,109 | 0,315 | 0,0170 | 2,723 | 30 | 9,93 | 0,246 | 16,25 | 637 |
| 1,511 | 0,20 | 2,58 | 9,97e-02 | 4,90e-03 | 0,1072 | 1,477 | 0,109 | 0,315 | 0,0211 | 3,572 | 41 | 13,45 | 0,304 | 22,08 | 882 |
| 1,511 | 0,30 | 17,91 | 1,04e-01 | 1,00e-02 | 0,1496 | 2,317 | 0,109 | 0,315 | 0,0376 | 5,153 | 88 | 19,73 | 0,431 | 31,12 | 1219 |

**Primera celda en que cada criterio supera 0,27 %:**

| Sistema | Criterio | Celdas > 0,27 % | Celdas | Primera celda que lo supera (orden: carga, luego δ) |
|---|---|---:|---:|---|
| IEEE-33 | bus mean (lambda_max interval) | 0 | 25 | nunca |
| IEEE-33 | bus max (lambda_max interval) | 1 | 25 | carga 3,441, δ 0,15 → 0,471 % |
| IEEE-33 | bus max (coordinate interval) | 3 | 25 | carga 3,441, δ 0,10 → 0,295 % |
| IEEE-33 | ellipse mean (Mahalanobis) | 3 | 25 | carga 3,441, δ 0,10 → 0,373 % |
| IEEE-33 | ellipse max (Mahalanobis) | 4 | 25 | carga 3,441, δ 0,10 → 0,520 % |
| IEEE-33 | system level (lambda_max interval) | 1 | 25 | carga 3,441, δ 0,15 → 0,481 % |
| PEGASE-2869 | bus mean (lambda_max interval) | 0 | 25 | nunca |
| PEGASE-2869 | bus max (lambda_max interval) | 24 | 25 | carga 0,800, δ 0,10 → 0,670 % |
| PEGASE-2869 | bus max (coordinate interval) | 25 | 25 | carga 0,800, δ 0,05 → 0,435 % |
| PEGASE-2869 | ellipse mean (Mahalanobis) | 4 | 25 | carga 1,200, δ 0,30 → 0,276 % |
| PEGASE-2869 | ellipse max (Mahalanobis) | 25 | 25 | carga 0,800, δ 0,05 → 1,705 % |
| PEGASE-2869 | system level (lambda_max interval) | 25 | 25 | carga 0,800, δ 0,05 → 2,970 % |

**Otras cotas contra AC (por barra), IEEE-33:**

| Carga | δ | Proy. exacta: media (%) | máx. (%) | > 0,27 % | Caja 3σ: media (%) | máx. (%) | > 0,27 % | Caja afín: media (%) | máx. (%) | > 0,27 % |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0,800 | 0,05 | 0,0098 | 0,025 | 0 | 0,1414 | 0,220 | 0 | 0,0000 | 0,000 | 0 |
| 0,800 | 0,10 | 0,0102 | 0,025 | 0 | 0,1412 | 0,220 | 0 | 0,0000 | 0,000 | 0 |
| 0,800 | 0,15 | 0,0102 | 0,025 | 0 | 0,1406 | 0,220 | 0 | 0,0000 | 0,000 | 0 |
| 0,800 | 0,20 | 0,0103 | 0,025 | 0 | 0,1411 | 0,220 | 0 | 0,0000 | 0,000 | 0 |
| 0,800 | 0,30 | 0,0106 | 0,025 | 0 | 0,1403 | 0,215 | 0 | 0,0000 | 0,000 | 0 |
| 1,000 | 0,05 | 0,0102 | 0,025 | 0 | 0,1419 | 0,225 | 0 | 0,0000 | 0,000 | 0 |
| 1,000 | 0,10 | 0,0102 | 0,025 | 0 | 0,1405 | 0,220 | 0 | 0,0000 | 0,000 | 0 |
| 1,000 | 0,15 | 0,0103 | 0,025 | 0 | 0,1409 | 0,220 | 0 | 0,0000 | 0,000 | 0 |
| 1,000 | 0,20 | 0,0105 | 0,025 | 0 | 0,1419 | 0,220 | 0 | 0,0000 | 0,000 | 0 |
| 1,000 | 0,30 | 0,0108 | 0,025 | 0 | 0,1414 | 0,210 | 0 | 0,0000 | 0,000 | 0 |
| 1,200 | 0,05 | 0,0103 | 0,025 | 0 | 0,1411 | 0,215 | 0 | 0,0000 | 0,000 | 0 |
| 1,200 | 0,10 | 0,0106 | 0,025 | 0 | 0,1408 | 0,210 | 0 | 0,0000 | 0,000 | 0 |
| 1,200 | 0,15 | 0,0103 | 0,025 | 0 | 0,1403 | 0,210 | 0 | 0,0000 | 0,000 | 0 |
| 1,200 | 0,20 | 0,0108 | 0,025 | 0 | 0,1405 | 0,210 | 0 | 0,0000 | 0,000 | 0 |
| 1,200 | 0,30 | 0,0109 | 0,025 | 0 | 0,1419 | 0,210 | 0 | 0,0000 | 0,000 | 0 |
| 1,400 | 0,05 | 0,0098 | 0,025 | 0 | 0,1400 | 0,210 | 0 | 0,0000 | 0,000 | 0 |
| 1,400 | 0,10 | 0,0100 | 0,025 | 0 | 0,1409 | 0,210 | 0 | 0,0000 | 0,000 | 0 |
| 1,400 | 0,15 | 0,0103 | 0,025 | 0 | 0,1397 | 0,200 | 0 | 0,0000 | 0,000 | 0 |
| 1,400 | 0,20 | 0,0105 | 0,030 | 0 | 0,1402 | 0,200 | 0 | 0,0000 | 0,000 | 0 |
| 1,400 | 0,30 | 0,0109 | 0,030 | 0 | 0,1452 | 0,210 | 0 | 0,0000 | 0,000 | 0 |
| 3,441 | 0,05 | 0,0339 | 0,075 | 0 | 0,2248 | 0,300 | 6 | 0,0000 | 0,000 | 0 |
| 3,441 | 0,10 | 0,1203 | 0,295 | 3 | 0,5294 | 0,655 | 29 | 0,0000 | 0,000 | 0 |
| 3,441 | 0,15 | 0,2880 | 0,672 | 15 | 1,0391 | 1,314 | 31 | 0,0000 | 0,000 | 0 |
| 3,441 | 0,20 | 0,0957 | 0,307 | 2 | 0,6768 | 0,942 | 28 | 0,0000 | 0,000 | 0 |
| 3,441 | 0,30 | 0,0002 | 0,006 | 0 | 0,0170 | 0,095 | 0 | 0,0000 | 0,000 | 0 |

**Otras cotas contra AC (por barra), PEGASE-2869:**

| Carga | δ | Proy. exacta: media (%) | máx. (%) | > 0,27 % | Caja 3σ: media (%) | máx. (%) | > 0,27 % | Caja afín: media (%) | máx. (%) | > 0,27 % |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0,800 | 0,05 | 0,0146 | 0,435 | 1 | 0,1253 | 0,980 | 83 | 0,0001 | 0,270 | 0 |
| 0,800 | 0,10 | 0,0165 | 1,395 | 2 | 0,1314 | 2,340 | 128 | 0,0003 | 0,670 | 1 |
| 0,800 | 0,15 | 0,0194 | 2,520 | 5 | 0,1404 | 4,025 | 183 | 0,0005 | 0,995 | 1 |
| 0,800 | 0,20 | 0,0231 | 3,760 | 10 | 0,1522 | 5,650 | 241 | 0,0006 | 1,310 | 1 |
| 0,800 | 0,30 | 0,0335 | 6,465 | 41 | 0,1832 | 10,145 | 356 | 0,0009 | 1,870 | 1 |
| 1,000 | 0,05 | 0,0155 | 0,545 | 2 | 0,1319 | 1,145 | 115 | 0,0002 | 0,365 | 1 |
| 1,000 | 0,10 | 0,0183 | 1,765 | 2 | 0,1398 | 2,790 | 175 | 0,0004 | 0,915 | 1 |
| 1,000 | 0,15 | 0,0223 | 3,100 | 4 | 0,1518 | 4,540 | 271 | 0,0006 | 1,360 | 1 |
| 1,000 | 0,20 | 0,0277 | 4,590 | 12 | 0,1672 | 6,255 | 328 | 0,0009 | 1,750 | 1 |
| 1,000 | 0,30 | 0,0426 | 7,905 | 64 | 0,2099 | 10,025 | 480 | 0,0013 | 2,555 | 1 |
| 1,200 | 0,05 | 0,0167 | 0,595 | 2 | 0,1385 | 1,070 | 152 | 0,0003 | 0,595 | 1 |
| 1,200 | 0,10 | 0,0204 | 1,565 | 3 | 0,1495 | 2,540 | 246 | 0,0005 | 1,200 | 1 |
| 1,200 | 0,15 | 0,0261 | 2,775 | 8 | 0,1668 | 4,040 | 354 | 0,0008 | 1,760 | 1 |
| 1,200 | 0,20 | 0,0339 | 4,080 | 24 | 0,1893 | 5,605 | 443 | 0,0011 | 2,290 | 1 |
| 1,200 | 0,30 | 0,0564 | 6,895 | 92 | 0,2535 | 8,920 | 742 | 0,0017 | 3,375 | 1 |
| 1,400 | 0,05 | 0,0181 | 0,800 | 2 | 0,1466 | 1,305 | 178 | 0,0003 | 0,800 | 1 |
| 1,400 | 0,10 | 0,0240 | 1,585 | 4 | 0,1654 | 2,935 | 355 | 0,0007 | 1,585 | 1 |
| 1,400 | 0,15 | 0,0338 | 2,290 | 9 | 0,1960 | 4,770 | 556 | 0,0011 | 2,290 | 1 |
| 1,400 | 0,20 | 0,0482 | 3,020 | 64 | 0,2384 | 6,270 | 746 | 0,0015 | 3,020 | 1 |
| 1,400 | 0,30 | 0,0989 | 4,936 | 252 | 0,3690 | 9,736 | 1157 | 0,0023 | 4,436 | 2 |
| 1,511 | 0,05 | 0,0200 | 0,940 | 1 | 0,1584 | 3,495 | 302 | 0,0004 | 0,940 | 1 |
| 1,511 | 0,10 | 0,0328 | 1,825 | 17 | 0,2004 | 7,305 | 551 | 0,0008 | 1,825 | 1 |
| 1,511 | 0,15 | 0,0594 | 2,723 | 129 | 0,2730 | 11,384 | 796 | 0,0013 | 2,723 | 1 |
| 1,511 | 0,20 | 0,0747 | 3,572 | 182 | 0,3275 | 15,711 | 1030 | 0,0018 | 3,572 | 1 |
| 1,511 | 0,30 | 0,1171 | 5,153 | 289 | 0,4276 | 23,838 | 1246 | 0,0029 | 5,153 | 3 |

*Script: `experiments/c1_ac_sweep.py` · CSV: `results/c1_ac_sweep.csv`, `results/c1_first_exceedance.csv` · N = 20.000 por celda · semilla en la columna `seed`*

Lectura:

- **(b) no cambia con la carga ni con δ:** la miscobertura de la elipse en el modelo lineal es 0,07–0,11 % en todas
  las celdas. Todo lo que varía en la tabla es error de linealización, que ahora queda medido por separado.
- **(a) El error de linealización crece como δ² en valor absoluto** (cada vez que δ se duplica, el error se
  cuadruplica) y como δ relativo al radio. A carga nominal y δ de diseño: 5,0·10⁻⁵ p.u. en IEEE-33 (0,7 % del
  radio en la peor barra) y 2,7·10⁻⁴ p.u. en PEGASE (53 % del radio en la peor barra; P99 de 1,4 %).
- **IEEE-33 se sostiene en todo el rango pedido salvo cerca del límite.** Para carga entre 0,8 y 1,4 y δ hasta
  30 %, la violación por barra no pasa de 0,01 % y la de la elipse de 0,14 %. Al 95 % del límite de cargabilidad
  (factor 3,44) la elipse supera 0,27 % en promedio desde δ = 10 % y la peor barra desde δ = 15 %.
- **PEGASE-2869: depende del criterio.**
  - *Promedio sobre barras:* nunca supera 0,27 % (máximo 0,038 %, en la celda más exigente).
  - *Peor barra:* supera 0,27 % en 24 de 25 celdas, incluido el punto de diseño (0,37 % a carga 1,0 y δ = 5 %).
    Son pocas barras: 1 a δ = 5 %, 2 a 10 %, 6 a 20 %, 18 a 30 % (carga nominal).
  - *Elipse:* el promedio sobre ramas pasa 0,27 % solo con δ = 30 % y carga ≥ 1,2, o cerca del límite con
    δ ≥ 20 %. Pero la peor elipse ya se viola en el 2,6 % de las muestras en el punto de diseño y en el 20 % con
    δ = 20 %. A carga nominal, 33 elipses superan 0,27 % con δ = 5 % y 281 con δ = 20 %.
  - *Sistema:* entre 3 % y 20 % de las muestras tiene al menos una barra fuera.
- **Respuesta a "desde qué combinación se supera 0,27 %":** en IEEE-33, solo al 95 % del límite de cargabilidad
  y con δ ≥ 10 %. En PEGASE-2869 no hay una frontera en el plano carga × δ: el promedio nunca la cruza y la peor
  barra la cruza casi siempre. El "δ ≤ 20 %" del manuscrito no se puede sostener como rango de validez general;
  hay que decir con qué criterio y para qué sistema.
- **No convergencia de Newton-Raphson:** prácticamente solo cerca del límite de cargabilidad: 0,3 %, 2,4 % y
  10,8 % en IEEE-33 (δ = 15, 20 y 30 %) y 0,1 %, 2,6 % y 17,9 % en PEGASE. Fuera de ahí, una única celda no es
  cero: PEGASE a carga 1,4 y δ = 30 % (0,01 %).
- **Ojo con las celdas de alta no convergencia.** En IEEE-33 al 95 % del límite, la violación entre las muestras
  convergidas *baja* de 0,47 % (δ = 15 %) a 0,00 % (δ = 30 %). No es que el método mejore: las muestras de mayor
  carga, que son las que violarían, son las que dejan de tener solución AC. Si se cuentan como violaciones, la
  tasa es al menos 2,4 % y 10,8 %. Lo correcto en esas celdas es reportar ambas cifras.
- La caja afín, que es determinista para el modelo lineal, también se viola contra AC en PEGASE (tablas "otras
  cotas"). Ninguna de las dos modalidades del método da una garantía determinista sobre el sistema no lineal.

---

## 7. C2 — Incertidumbre correlacionada entre fuentes

El código acepta ahora una covarianza de entrada: Q = C·Σ_ε·Cᵀ (`cdhim.SourceCovariance`, argumento
`sigma_eps` de `Study.run`). El valor por defecto es Σ_ε = I/3 y da exactamente los resultados anteriores
(filas r = 0, idénticas entre sí; además lo comprueba una prueba unitaria). Las muestras se generan con una
cópula gaussiana de un factor y marginales U[−1, 1]; la correlación de Pearson resultante es
ρ = (6/π)·asin(r/2). Punto de diseño, N = 20.000, modelo lineal y AC.

| Sistema | r cópula | ρ entre ε | Q | Ramas Elipse (%) | CRI área (%) | Red. ancho λ_max vs afín (%) | Elipse, lineal: media (%) | máx. (%) | > 0,27 % | Elipse, AC: media (%) | Barra λ_max, AC: media (%) | máx. (%) | > 0,27 % | Caja 3σ, AC: máx. (%) | Caja afín, AC: máx. (%) |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| IEEE-33 | 0,0 | 0,000 | misspecified_I/3 | 67,7 | 86,37 | 42,78 | 0,101 | 0,18 | 0 | 0,102 | 0,0013 | 0,04 | 0 | 0,30 | 0,000 |
| IEEE-33 | 0,0 | 0,000 | true_Sigma_eps | 67,7 | 86,37 | 42,78 | 0,101 | 0,18 | 0 | 0,102 | 0,0013 | 0,04 | 0 | 0,30 | 0,000 |
| IEEE-33 | 0,3 | 0,288 | misspecified_I/3 | 67,7 | 86,37 | 42,78 | 26,345 | 31,25 | 21 | 26,342 | 7,4744 | 23,25 | 26 | 36,37 | 0,000 |
| IEEE-33 | 0,3 | 0,288 | true_Sigma_eps | 35,5 | 92,58 | -51,59 | 0,048 | 0,09 | 0 | 0,047 | 0,0000 | 0,00 | 0 | 0,01 | 0,000 |
| IEEE-33 | 0,6 | 0,582 | misspecified_I/3 | 67,7 | 86,37 | 42,78 | 45,305 | 50,55 | 21 | 45,299 | 20,1409 | 43,28 | 26 | 55,49 | 0,010 |
| IEEE-33 | 0,6 | 0,582 | true_Sigma_eps | 16,1 | 93,81 | -108,11 | 0,108 | 0,14 | 0 | 0,106 | 0,0072 | 0,01 | 0 | 0,00 | 0,010 |
| IEEE-33 | 0,9 | 0,891 | misspecified_I/3 | 67,7 | 86,37 | 42,78 | 57,453 | 62,02 | 21 | 57,451 | 31,6256 | 56,55 | 32 | 66,64 | 0,320 |
| IEEE-33 | 0,9 | 0,891 | true_Sigma_eps | 3,2 | 94,55 | -214,06 | 0,295 | 0,29 | 1 | 0,305 | 0,2472 | 0,32 | 21 | 0,00 | 0,320 |
| PEGASE-2869 | 0,0 | 0,000 | misspecified_I/3 | 79,1 | 71,08 | 34,12 | 0,106 | 0,34 | 28 | 0,111 | 0,0027 | 0,46 | 1 | 1,19 | 0,455 |
| PEGASE-2869 | 0,0 | 0,000 | true_Sigma_eps | 79,1 | 71,08 | 34,12 | 0,106 | 0,34 | 28 | 0,111 | 0,0027 | 0,46 | 1 | 1,19 | 0,455 |
| PEGASE-2869 | 0,3 | 0,288 | misspecified_I/3 | 79,1 | 71,08 | 34,12 | 34,520 | 81,93 | 2434 | 34,338 | 12,4912 | 79,38 | 1407 | 82,95 | 7,035 |
| PEGASE-2869 | 0,3 | 0,288 | true_Sigma_eps | 59,3 | 82,57 | -52,29 | 0,052 | 0,23 | 0 | 1,139 | 0,0804 | 9,16 | 129 | 36,16 | 7,035 |
| PEGASE-2869 | 0,6 | 0,582 | misspecified_I/3 | 79,1 | 71,08 | 34,12 | 48,761 | 87,50 | 2446 | 48,533 | 21,3131 | 86,18 | 1527 | 88,46 | 21,950 |
| PEGASE-2869 | 0,6 | 0,582 | true_Sigma_eps | 44,7 | 85,12 | -106,98 | 0,100 | 0,28 | 1 | 3,519 | 0,3900 | 24,16 | 274 | 50,87 | 21,950 |
| PEGASE-2869 | 0,9 | 0,891 | misspecified_I/3 | 79,1 | 71,08 | 34,12 | 58,673 | 90,95 | 2448 | 58,374 | 28,5428 | 89,96 | 1600 | 91,71 | 36,680 |
| PEGASE-2869 | 0,9 | 0,891 | true_Sigma_eps | 23,0 | 86,63 | -180,91 | 0,235 | 0,52 | 99 | 10,385 | 1,3498 | 40,66 | 460 | 63,95 | 36,680 |

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

| Distribución | Soporte mín. | Soporte máx. | Media | Varianza | Dentro de [−1, 1] | Masa fuera de [−1, 1] (%) | Definición |
|---|---:|---:|---:|---:|---:|---:|---|
| uniform | -1,000 | 1,000 | -0,0000 | 0,3333 | sí | 0,00 | U[-1,1], the manuscript's assumption |
| beta25_support | -0,400 | 1,000 | 0,0000 | 0,0500 | sí | 0,00 | Beta(2,5): (u-mu)/max(mu,1-mu), support respected |
| beta25_var13 | -1,033 | 2,582 | 0,0000 | 0,3333 | no | 6,24 | Beta(2,5) rescaled to variance 1/3; leaves [-1,1] |
| beta25_paper | -1,140 | 2,860 | 0,0029 | 0,4082 | no | 9,66 | (u-0.285)*4 exactly as in finalFinalCode.py:580 |
| beta52_support | -1,000 | 0,400 | 0,0000 | 0,0500 | sí | 0,00 | Beta(5,2): opposite skew, support respected |
| beta0505_ushape | -1,000 | 1,000 | 0,0000 | 0,5000 | sí | 0,00 | Beta(0.5,0.5): U-shaped, variance 1/2 > 1/3 |
| triangular | -1,000 | 1,000 | -0,0000 | 0,1667 | sí | 0,00 | symmetric triangular on [-1,1] |
| truncnorm_3s | -1,000 | 1,000 | 0,0000 | 0,1081 | sí | 0,00 | N(0,1/3^2) truncated at +-3 sigma |
| student_t3_trunc | -1,000 | 1,000 | 0,0000 | 0,0271 | sí | 0,00 | Student t (3 dof) truncated to its central 99.73 %, rescaled into [-1,1] |
| student_t3 | -inf | inf | 0,0000 | 0,3333 | no | 5,77 | Student t (3 dof), variance 1/3, unbounded: outside the model |
| lognormal_trunc | -0,388 | 1,000 | 0,0000 | 0,0393 | sí | 0,00 | lognormal(s=0.5) truncated at its 99.73 % quantile, rescaled into [-1,1] |
| lognormal | -1,083 | inf | 0,0000 | 0,3333 | no | 5,96 | lognormal(s=0.5), variance 1/3, unbounded above: outside the model |

Para cada distribución, dos covarianzas: *true_variance* usa Q = Var(ε)·CCᵀ y *uniform_1/3* usa Q = ⅓CCᵀ
(el supuesto del paper, sea cual sea la distribución real). Violación sobre **todas** las ramas en modo
Elipse: modelo lineal con N = 200.000 y AC con N = 5.000.

**IEEE-33 (δ = 15 %, 21 ramas en modo Elipse)**

| Distribución | Varianza en Q | Var(ε) real | Elipse, lineal: media (%) | P99 (%) | máx. (%) | ramas > 0,27 % | Barra λ_max, lineal: máx. (%) | Barra proy., lineal: máx. (%) | Elipse, AC: media (%) | Barra λ_max, AC: máx. (%) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| uniform | true_variance | 0,3333 | 0,100 | 0,155 | 0,155 | 0 | 0,021 | 0,039 | 0,104 | 0,020 |
| uniform | uniform_1/3 | 0,3333 | 0,100 | 0,155 | 0,155 | 0 | 0,021 | 0,039 | 0,104 | 0,020 |
| beta25_support | true_variance | 0,0500 | 0,306 | 0,322 | 0,323 | 21 | 0,073 | 0,088 | 0,351 | 0,100 |
| beta25_support | uniform_1/3 | 0,0500 | 0,000 | 0,000 | 0,000 | 0 | 0,000 | 0,000 | 0,000 | 0,000 |
| beta25_var13 | true_variance | 0,3333 | 0,301 | 0,327 | 0,329 | 21 | 0,063 | 0,086 | 0,321 | 0,080 |
| beta25_var13 | uniform_1/3 | 0,3333 | 0,301 | 0,327 | 0,329 | 21 | 0,063 | 0,086 | 0,321 | 0,080 |
| beta25_paper | true_variance | 0,4082 | 0,312 | 0,344 | 0,344 | 21 | 0,064 | 0,083 | 0,340 | 0,100 |
| beta25_paper | uniform_1/3 | 0,4082 | 0,836 | 0,900 | 0,901 | 21 | 0,170 | 0,228 | 0,886 | 0,240 |
| beta52_support | true_variance | 0,0500 | 0,301 | 0,318 | 0,319 | 21 | 0,060 | 0,083 | 0,261 | 0,080 |
| beta52_support | uniform_1/3 | 0,0500 | 0,000 | 0,000 | 0,000 | 0 | 0,000 | 0,000 | 0,000 | 0,000 |
| beta0505_ushape | true_variance | 0,5000 | 0,074 | 0,131 | 0,132 | 0 | 0,013 | 0,036 | 0,080 | 0,000 |
| beta0505_ushape | uniform_1/3 | 0,5000 | 1,103 | 1,545 | 1,549 | 20 | 0,306 | 0,420 | 1,266 | 0,460 |
| triangular | true_variance | 0,1667 | 0,158 | 0,205 | 0,206 | 0 | 0,030 | 0,051 | 0,191 | 0,080 |
| triangular | uniform_1/3 | 0,1667 | 0,000 | 0,001 | 0,002 | 0 | 0,000 | 0,001 | 0,000 | 0,000 |
| truncnorm_3s | true_variance | 0,1081 | 0,237 | 0,269 | 0,269 | 0 | 0,046 | 0,068 | 0,272 | 0,060 |
| truncnorm_3s | uniform_1/3 | 0,1081 | 0,000 | 0,000 | 0,000 | 0 | 0,000 | 0,000 | 0,000 | 0,000 |
| student_t3_trunc | true_variance | 0,0271 | 1,055 | 1,847 | 1,902 | 21 | 0,258 | 0,563 | 1,010 | 0,300 |
| student_t3_trunc | uniform_1/3 | 0,0271 | 0,000 | 0,000 | 0,000 | 0 | 0,000 | 0,000 | 0,000 | 0,000 |
| student_t3 | true_variance | 0,3333 | 1,481 | 1,831 | 1,855 | 21 | 0,572 | 0,753 | 1,416 | 0,500 |
| student_t3 | uniform_1/3 | 0,3333 | 1,481 | 1,831 | 1,855 | 21 | 0,572 | 0,753 | 1,416 | 0,500 |
| lognormal_trunc | true_variance | 0,0393 | 0,851 | 1,552 | 1,629 | 21 | 0,208 | 0,420 | 0,739 | 0,140 |
| lognormal_trunc | uniform_1/3 | 0,0393 | 0,000 | 0,000 | 0,000 | 0 | 0,000 | 0,000 | 0,000 | 0,000 |
| lognormal | true_variance | 0,3333 | 1,120 | 1,766 | 1,806 | 21 | 0,364 | 0,607 | 1,129 | 0,360 |
| lognormal | uniform_1/3 | 0,3333 | 1,120 | 1,766 | 1,806 | 21 | 0,364 | 0,607 | 1,129 | 0,360 |

**PEGASE-2869 (δ = 5 %, 2.543 ramas en modo Elipse)**

| Distribución | Varianza en Q | Var(ε) real | Elipse, lineal: media (%) | P99 (%) | máx. (%) | ramas > 0,27 % | Barra λ_max, lineal: máx. (%) | Barra proy., lineal: máx. (%) | Elipse, AC: media (%) | Barra λ_max, AC: máx. (%) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| uniform | true_variance | 0,3333 | 0,102 | 0,242 | 0,264 | 0 | 0,059 | 0,068 | 0,103 | 0,320 |
| uniform | uniform_1/3 | 0,3333 | 0,102 | 0,242 | 0,264 | 0 | 0,059 | 0,068 | 0,103 | 0,320 |
| beta25_support | true_variance | 0,0500 | 0,291 | 0,343 | 0,370 | 2082 | 0,097 | 0,104 | 0,304 | 0,200 |
| beta25_support | uniform_1/3 | 0,0500 | 0,000 | 0,000 | 0,000 | 0 | 0,000 | 0,000 | 0,000 | 0,000 |
| beta25_var13 | true_variance | 0,3333 | 0,288 | 0,339 | 0,364 | 2024 | 6,261 | 6,261 | 0,319 | 6,500 |
| beta25_var13 | uniform_1/3 | 0,3333 | 0,288 | 0,339 | 0,364 | 2024 | 6,261 | 6,261 | 0,319 | 6,500 |
| beta25_paper | true_variance | 0,4082 | 0,300 | 0,350 | 0,382 | 2338 | 9,730 | 9,730 | 0,329 | 10,040 |
| beta25_paper | uniform_1/3 | 0,4082 | 0,828 | 0,912 | 0,952 | 2543 | 9,730 | 9,730 | 0,882 | 10,040 |
| beta52_support | true_variance | 0,0500 | 0,289 | 0,343 | 0,368 | 2080 | 0,117 | 0,117 | 0,297 | 0,180 |
| beta52_support | uniform_1/3 | 0,0500 | 0,000 | 0,000 | 0,000 | 0 | 0,000 | 0,000 | 0,000 | 0,000 |
| beta0505_ushape | true_variance | 0,5000 | 0,084 | 0,235 | 0,267 | 0 | 0,063 | 0,063 | 0,091 | 3,840 |
| beta0505_ushape | uniform_1/3 | 0,5000 | 1,067 | 1,864 | 1,942 | 2385 | 0,496 | 0,511 | 1,085 | 3,840 |
| triangular | true_variance | 0,1667 | 0,156 | 0,260 | 0,283 | 13 | 0,065 | 0,070 | 0,157 | 0,120 |
| triangular | uniform_1/3 | 0,1667 | 0,000 | 0,001 | 0,002 | 0 | 0,001 | 0,001 | 0,000 | 0,000 |
| truncnorm_3s | true_variance | 0,1081 | 0,225 | 0,275 | 0,300 | 57 | 0,068 | 0,070 | 0,226 | 0,140 |
| truncnorm_3s | uniform_1/3 | 0,1081 | 0,000 | 0,000 | 0,000 | 0 | 0,000 | 0,000 | 0,000 | 0,020 |
| student_t3_trunc | true_variance | 0,0271 | 1,092 | 2,086 | 2,522 | 2543 | 1,124 | 1,125 | 1,113 | 1,460 |
| student_t3_trunc | uniform_1/3 | 0,0271 | 0,000 | 0,000 | 0,000 | 0 | 0,000 | 0,000 | 0,000 | 0,000 |
| student_t3 | true_variance | 0,3333 | 1,399 | 1,887 | 2,087 | 2543 | 5,769 | 5,769 | 1,432 | 5,800 |
| student_t3 | uniform_1/3 | 0,3333 | 1,399 | 1,887 | 2,087 | 2543 | 5,769 | 5,769 | 1,432 | 5,800 |
| lognormal_trunc | true_variance | 0,0393 | 0,859 | 1,731 | 2,199 | 2543 | 0,910 | 0,950 | 0,868 | 1,040 |
| lognormal_trunc | uniform_1/3 | 0,0393 | 0,000 | 0,000 | 0,000 | 0 | 0,000 | 0,000 | 0,000 | 0,000 |
| lognormal | true_variance | 0,3333 | 1,114 | 1,907 | 2,257 | 2543 | 5,957 | 5,957 | 1,162 | 6,120 |
| lognormal | uniform_1/3 | 0,3333 | 1,114 | 1,907 | 2,257 | 2543 | 5,957 | 5,957 | 1,162 | 6,120 |

**Distribución sesgada con fuentes correlacionadas** (Beta(2,5) con soporte respetado, cópula r = 0,6, N = 20.000):

| Sistema | Q | ρ entre ε | CRI área (%) | Elipse, lineal: media (%) | máx. (%) | > 0,27 % | Elipse, AC: media (%) | Barra λ_max, lineal: máx. (%) | Barra λ_max, AC: máx. (%) | barras > 0,27 % |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| IEEE-33 | true (variance and correlation) | 0,592 | 93,93 | 0,451 | 0,52 | 5 | 0,473 | 0,030 | 0,030 | 0 |
| IEEE-33 | true variance, independent | 0,592 | 86,37 | 40,496 | 45,66 | 21 | 40,482 | 38,625 | 38,595 | 26 |
| IEEE-33 | uniform 1/3, independent | 0,592 | 86,37 | 2,605 | 3,58 | 21 | 2,644 | 2,145 | 2,205 | 14 |
| PEGASE-2869 | true (variance and correlation) | 0,592 | 85,20 | 0,449 | 0,60 | 1423 | 1,722 | 0,145 | 7,190 | 137 |
| PEGASE-2869 | true variance, independent | 0,592 | 71,08 | 46,316 | 86,76 | 2527 | 46,178 | 85,205 | 85,175 | 1594 |
| PEGASE-2869 | uniform 1/3, independent | 0,592 | 71,08 | 12,054 | 65,86 | 2107 | 11,964 | 62,895 | 62,400 | 965 |

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

| Sistema | Linealización | Caja en G | κ mín. | κ mediana | κ máx. | G mínimo (%) |
|---|---|---|---:|---:|---:|---:|
| IEEE-33 | legacy | covariance_1sigma | 0,7854 | 0,7854 | 0,7854 | 46,64 |
| IEEE-33 | legacy | covariance_3sigma | 1,0323 | 1,0323 | 1,0323 | 29,87 |
| IEEE-33 | legacy | affine | 0,2706 | 0,3468 | 0,5712 | 81,62 |
| IEEE-33 | corrected | covariance_1sigma | 0,7854 | 0,7854 | 0,7854 | 41,63 |
| IEEE-33 | corrected | covariance_3sigma | 1,0323 | 1,0323 | 1,0323 | 23,29 |
| IEEE-33 | corrected | affine | 0,1021 | 0,1197 | 0,2916 | 89,02 |
| PEGASE-2869 | legacy | covariance_1sigma | 0,7854 | 0,7854 | 0,7854 | 21,46 |
| PEGASE-2869 | legacy | covariance_3sigma | 1,0323 | 1,0323 | 1,0323 | 0,00 |
| PEGASE-2869 | legacy | affine | 0,0407 | 0,2219 | 1,8843 | 0,00 |
| PEGASE-2869 | corrected | covariance_1sigma | 0,7854 | 0,7854 | 0,7854 | 21,58 |
| PEGASE-2869 | corrected | covariance_3sigma | 1,0323 | 1,0323 | 1,0323 | 0,00 |
| PEGASE-2869 | corrected | affine | 0,0047 | 0,0628 | 2,7326 | 8,34 |

**Las dos reglas, con γ_th = 15 %:**

| Sistema | Linealización | Caja en G | ρ_th equivalente | Elipse (regla G) | Elipse (regla ρ) | Ramas que cambian | CRI regla G (%) | CRI regla ρ (%) | Tiempo G (ms) | Tiempo ρ (ms) |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| IEEE-33 | legacy | covariance_1sigma | 0,0000 | 16 | 16 | 0 | 86,70 | 86,70 | 0,55 | 0,35 |
| IEEE-33 | legacy | covariance_3sigma | 0,5674 | 16 | 16 | 0 | 82,52 | 82,52 | 0,33 | 0,39 |
| IEEE-33 | legacy | affine | 0,0000 | 16 | 16 | 0 | 94,39 | 94,39 | 0,58 | 0,34 |
| IEEE-33 | corrected | covariance_1sigma | 0,0000 | 21 | 21 | 0 | 86,37 | 86,37 | 0,35 | 0,35 |
| IEEE-33 | corrected | covariance_3sigma | 0,5674 | 21 | 21 | 0 | 82,09 | 82,09 | 0,34 | 0,33 |
| IEEE-33 | corrected | affine | 0,0000 | 21 | 21 | 0 | 97,37 | 97,37 | 0,35 | 0,32 |
| PEGASE-2869 | legacy | covariance_1sigma | 0,0000 | 2984 | 2984 | 0 | 57,80 | 57,80 | 95,88 | 103,52 |
| PEGASE-2869 | legacy | covariance_3sigma | 0,5674 | 2227 | 2227 | 0 | 58,55 | 58,55 | 101,77 | 99,59 |
| PEGASE-2869 | legacy | affine | 0,0000 | 2962 | 2984 | 22 | 84,29 | 83,69 | 93,72 | 102,76 |
| PEGASE-2869 | corrected | covariance_1sigma | 0,0000 | 2543 | 2543 | 0 | 71,08 | 71,08 | 103,01 | 98,75 |
| PEGASE-2869 | corrected | covariance_3sigma | 0,5674 | 2339 | 2339 | 0 | 66,90 | 66,90 | 106,36 | 105,33 |
| PEGASE-2869 | corrected | affine | 0,0000 | 2542 | 2543 | 1 | 96,01 | 95,98 | 96,09 | 102,10 |

**Efecto de las ramas que cambian (solo aparece con la caja afín):**

| Sistema | Caja en G | Ramas que cambian | Red. ancho regla G (%) | Red. ancho regla ρ (%) | Barra AC máx., G (%) | Barra AC máx., ρ (%) | Elipse AC media, G (%) | Elipse AC media, ρ (%) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| PEGASE-2869 | affine | 1 | 34,92 | 34,12 | 0,440 | 0,440 | 0,108 | 0,108 |

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

| Cantidad | IEEE-33, 99,73 % por rama | IEEE-33, Bonferroni | PEGASE-2869, 99,73 % por rama | PEGASE-2869, Bonferroni |
|---|---|---|---|---|
| m (ramas en modo Elipse) | 21 | 21 | 2543 | 2543 |
| Cuantil χ² | 11,83 | 17,92 | 11,83 | 27,51 |
| Factor sobre los radios | 1,000 | 1,231 | 1,000 | 1,525 |
| ≥ 1 barra fuera, intervalo λ_max, modelo lineal (%) | 0,02 | 0,00 | 2,52 | 0,00 |
| ≥ 1 barra fuera, intervalo λ_max, AC (%) | 0,02 | 0,00 | 3,15 | 0,42 |
| Barras fuera por muestra (esperado), λ_max, AC | 0,0002 | 0,0000 | 0,0571 | 0,0043 |
| ≥ 1 barra fuera, proyección exacta, lineal (%) | 0,10 | 0,00 | 8,20 | 0,00 |
| ≥ 1 barra fuera, proyección exacta, AC (%) | 0,10 | 0,00 | 9,31 | 0,44 |
| Barras fuera por muestra, proyección exacta, AC | 0,0037 | 0,0000 | 0,3616 | 0,0044 |
| ≥ 1 barra fuera, caja 3σ, AC (%) | 0,72 | 0,72 | 41,84 | 41,84 |
| Barras fuera por muestra, caja 3σ, AC | 0,048 | 0,048 | 3,144 | 3,144 |
| ≥ 1 barra fuera, caja afín, AC (%) | 0,00 | 0,00 | 0,42 | 0,42 |
| ≥ 1 elipse violada, modelo lineal (%) | 0,99 | 0,03 | 48,27 | 0,00 |
| ≥ 1 elipse violada, AC (%) | 0,99 | 0,03 | 49,89 | 0,32 |
| Elipses violadas por muestra, lineal | 0,020 | 0,000 | 2,494 | 0,000 |
| Elipses violadas por muestra, AC | 0,020 | 0,000 | 2,619 | 0,003 |
| CRI de área vs caja de covarianza (%) | 86,37 | 86,37 | 71,08 | 71,08 |
| CRI de área vs caja afín (%) | 97,37 | 96,01 | 95,98 | 90,92 |
| Reducción de ancho λ_max vs caja afín (%) | 42,78 | 29,58 | 34,12 | -0,47 |
| Reducción de ancho proyección exacta vs caja afín (%) | 56,20 | 46,10 | 64,51 | 45,88 |
| Barras con λ_max más ancho que la caja afín | 0 | 1 | 140 | 308 |
| Barras con proyección más ancha que la caja afín | 0 | 0 | 50 | 140 |

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

Caja en G: `covariance_1sigma` — ramas en modo Elipse (%) / CRI de área (%)

| Sistema | γ_th = 5 % | γ_th = 10 % | γ_th = 15 % | γ_th = 20 % | γ_th = 25 % |
|---|---|---|---|---|---|
| case33bw | 67,7 / 86,4 | 67,7 / 86,4 | 67,7 / 86,4 | 67,7 / 86,4 | 67,7 / 86,4 |
| case118 | 100,0 / 55,6 | 100,0 / 55,6 | 100,0 / 55,6 | 100,0 / 55,6 | 100,0 / 55,6 |
| case300 | 95,1 / 64,4 | 95,1 / 64,4 | 95,1 / 64,4 | 95,1 / 64,4 | 91,2 / 66,1 |
| case1354pegase | 76,8 / 70,0 | 76,8 / 70,0 | 76,8 / 70,0 | 76,8 / 70,0 | 75,7 / 70,6 |
| case2869pegase | 79,7 / 70,3 | 79,7 / 70,3 | 79,7 / 70,3 | 79,7 / 70,3 | 77,8 / 71,5 |

Caja en G: `covariance_3sigma` — ramas en modo Elipse (%) / CRI de área (%)

| Sistema | γ_th = 5 % | γ_th = 10 % | γ_th = 15 % | γ_th = 20 % | γ_th = 25 % |
|---|---|---|---|---|---|
| case33bw | 67,7 / 82,1 | 67,7 / 82,1 | 67,7 / 82,1 | 67,7 / 82,1 | 64,5 / 85,0 |
| case118 | 88,6 / 46,6 | 77,1 / 52,3 | 77,1 / 52,3 | 77,1 / 52,3 | 71,4 / 54,8 |
| case300 | 88,7 / 56,9 | 84,5 / 59,4 | 81,3 / 61,2 | 79,5 / 62,2 | 76,7 / 63,7 |
| case1354pegase | 75,0 / 61,9 | 73,0 / 63,4 | 70,4 / 65,4 | 68,1 / 66,9 | 66,3 / 68,2 |
| case2869pegase | 76,4 / 63,6 | 74,9 / 64,8 | 72,7 / 66,4 | 70,6 / 67,7 | 68,8 / 69,0 |

Caja en G: `affine` — ramas en modo Elipse (%) / CRI de área (%)

| Sistema | γ_th = 5 % | γ_th = 10 % | γ_th = 15 % | γ_th = 20 % | γ_th = 25 % |
|---|---|---|---|---|---|
| case33bw | 67,7 / 97,4 | 67,7 / 97,4 | 67,7 / 97,4 | 67,7 / 97,4 | 67,7 / 97,4 |
| case118 | 100,0 / 80,6 | 100,0 / 80,6 | 100,0 / 80,6 | 100,0 / 80,6 | 100,0 / 80,6 |
| case300 | 95,1 / 91,1 | 94,7 / 91,4 | 94,7 / 91,4 | 94,7 / 91,4 | 94,3 / 91,7 |
| case1354pegase | 76,8 / 94,9 | 76,7 / 95,0 | 76,7 / 95,0 | 76,7 / 95,0 | 76,7 / 95,0 |
| case2869pegase | 79,7 / 95,8 | 79,7 / 95,8 | 79,7 / 95,8 | 79,7 / 95,9 | 79,6 / 95,9 |

**Óptimo según el criterio del paper (CRI sobre tiempo normalizado) y estabilidad:**

| Sistema | Caja en G | Selecciones distintas en el barrido | Elipse mín. (%) | Elipse máx. (%) | γ óptimo (CRI/tiempo) | γ óptimo (solo CRI) | Variación de tiempo en el barrido (%) | Violación AC máx. por barra, rango (%) |
|---|---|---:|---:|---:|---:|---:|---:|---|
| case33bw | covariance_1sigma | 1 | 67,7 | 67,7 | 0,05 | 0,05 | 14,2 | 0,000 – 0,000 |
| case33bw | covariance_3sigma | 2 | 64,5 | 67,7 | 0,20 | 0,25 | 7,4 | 0,000 – 0,000 |
| case33bw | affine | 1 | 67,7 | 67,7 | 0,25 | 0,05 | 75,0 | 0,000 – 0,000 |
| case118 | covariance_1sigma | 1 | 100,0 | 100,0 | 0,20 | 0,05 | 9,6 | 0,080 – 0,080 |
| case118 | covariance_3sigma | 3 | 71,4 | 88,6 | 0,25 | 0,25 | 11,7 | 0,040 – 0,080 |
| case118 | affine | 1 | 100,0 | 100,0 | 0,10 | 0,05 | 9,8 | 0,080 – 0,080 |
| case300 | covariance_1sigma | 2 | 91,2 | 95,1 | 0,20 | 0,25 | 34,1 | 1,060 – 1,260 |
| case300 | covariance_3sigma | 5 | 76,7 | 88,7 | 0,15 | 0,25 | 24,9 | 0,980 – 1,060 |
| case300 | affine | 3 | 94,3 | 95,1 | 0,10 | 0,25 | 59,9 | 1,260 – 1,260 |
| case1354pegase | covariance_1sigma | 2 | 75,7 | 76,8 | 0,20 | 0,25 | 13,9 | 0,600 – 0,600 |
| case1354pegase | covariance_3sigma | 5 | 66,3 | 75,0 | 0,25 | 0,25 | 10,8 | 0,600 – 0,600 |
| case1354pegase | affine | 2 | 76,7 | 76,8 | 0,25 | 0,10 | 10,5 | 0,600 – 0,600 |
| case2869pegase | covariance_1sigma | 2 | 77,8 | 79,7 | 0,25 | 0,25 | 13,6 | 0,360 – 0,360 |
| case2869pegase | covariance_3sigma | 5 | 68,8 | 76,4 | 0,15 | 0,25 | 5,9 | 0,360 – 0,360 |
| case2869pegase | affine | 4 | 79,6 | 79,7 | 0,25 | 0,25 | 3,6 | 0,360 – 0,360 |

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

| Sistema | Barras | N_u | Ramas PQ-PQ | Jacobiano (ms) | LU (ms) | Resoluciones (ms) | Ramas (ms) | Total (ms) | Resoluciones / total (%) | Resoluciones con hilos (ms) | Flujo base pandapower (ms) | Bloque denso de sensibilidades (MB) |
|---|---:|---:|---:|---|---|---|---|---|---:|---:|---:|---:|
| case33bw | 33 | 64 | 31 | 2,49 ± 0,23 | 0,13 ± 0,03 | 0,16 ± 0,05 | 0,41 ± 0,12 | 3,19 ± 0,37 | 5,0 | 1,82 | 42 | 0,0 |
| case118 | 118 | 152 | 35 | 2,60 ± 0,09 | 0,35 ± 0,07 | 0,84 ± 0,08 | 0,47 ± 0,04 | 4,27 ± 0,18 | 19,7 | 2,75 | 52 | 0,1 |
| case300 | 300 | 346 | 283 | 3,20 ± 0,18 | 0,95 ± 0,20 | 7,21 ± 0,84 | 1,27 ± 0,28 | 12,63 ± 1,12 | 57,1 | 6,80 | 64 | 0,6 |
| case1354pegase | 1354 | 1332 | 1309 | 4,19 ± 0,40 | 4,41 ± 0,38 | 132,76 ± 3,33 | 21,69 ± 0,53 | 163,05 ± 3,43 | 81,4 | 86,52 | 88 | 11,1 |
| case2869pegase | 2869 | 2951 | 3216 | 5,74 ± 0,80 | 10,25 ± 0,69 | 765,35 ± 21,78 | 105,74 ± 1,07 | 887,08 ± 21,32 | 86,3 | 405,49 | 133 | 53,1 |
| case9241pegase | 9241 | 9687 | 11633 | 14,94 ± 0,84 | 43,40 ± 2,94 | 10433,41 ± 97,99 | 1189,38 ± 13,03 | 11681,13 ± 101,21 | 89,3 | 4630,65 | 372 | 576,2 |

**Relleno de la factorización:**

| Sistema | Ordenamiento | nnz(J) | nnz(L) | nnz(U) | (L+U)/J | LU (ms) |
|---|---|---:|---:|---:|---:|---:|
| case33bw | COLAMD | 376 | 271 | 323 | 1,58 | 0,12 |
| case33bw | MMD_AT_PLUS_A | 376 | 220 | 288 | 1,35 | 0,09 |
| case33bw | MMD_ATA | 376 | 220 | 298 | 1,38 | 0,15 |
| case118 | COLAMD | 1051 | 1013 | 1075 | 1,99 | 0,31 |
| case118 | MMD_AT_PLUS_A | 1051 | 758 | 805 | 1,49 | 0,37 |
| case118 | MMD_ATA | 1051 | 1037 | 1088 | 2,02 | 0,39 |
| case300 | COLAMD | 3736 | 4124 | 4490 | 2,31 | 1,00 |
| case300 | MMD_AT_PLUS_A | 3736 | 3144 | 3614 | 1,81 | 0,91 |
| case300 | MMD_ATA | 3736 | 3766 | 4287 | 2,16 | 1,32 |
| case1354pegase | COLAMD | 15803 | 16349 | 20373 | 2,32 | 3,60 |
| case1354pegase | MMD_AT_PLUS_A | 15803 | 12299 | 16839 | 1,84 | 4,03 |
| case1354pegase | MMD_ATA | 15803 | 15764 | 19957 | 2,26 | 5,89 |
| case2869pegase | COLAMD | 36591 | 43775 | 49265 | 2,54 | 9,42 |
| case2869pegase | MMD_AT_PLUS_A | 36591 | 30411 | 36346 | 1,82 | 9,73 |
| case2869pegase | MMD_ATA | 36591 | 39896 | 45613 | 2,34 | 14,60 |
| case9241pegase | COLAMD | 129412 | 184604 | 197696 | 2,95 | 39,60 |
| case9241pegase | MMD_AT_PLUS_A | 129412 | 118999 | 130669 | 1,93 | 51,83 |
| case9241pegase | MMD_ATA | 129412 | 176876 | 189782 | 2,83 | 59,59 |

**Exponentes en escala log-log (y ∝ N^k, N = número de barras):**

| Cantidad | Sistemas | Puntos | Exponente | IC 95 % | R² |
|---|---|---:|---:|---|---:|
| jacobian | all systems | 6 | 0,29 | [0,11; 0,48] | 0,828 |
| lu | all systems | 6 | 1,04 | [0,94; 1,13] | 0,996 |
| solves | all systems | 6 | 2,00 | [1,78; 2,22] | 0,994 |
| branch | all systems | 6 | 1,49 | [0,95; 2,04] | 0,935 |
| total | all systems | 6 | 1,51 | [0,98; 2,04] | 0,940 |
| nnz(L+U), COLAMD | all systems | 6 | 1,15 | [1,07; 1,23] | 0,997 |
| nnz(J) | all systems | 6 | 1,05 | [0,98; 1,12] | 0,997 |
| dense |V| block entries (n_pq x N_u) | all systems | 6 | 1,91 | [1,72; 2,11] | 0,995 |
| jacobian | n_bus >= 118 | 5 | 0,37 | [0,11; 0,62] | 0,874 |
| lu | n_bus >= 118 | 5 | 1,09 | [1,01; 1,18] | 0,998 |
| solves | n_bus >= 118 | 5 | 2,14 | [2,01; 2,26] | 0,999 |
| branch | n_bus >= 118 | 5 | 1,83 | [1,51; 2,15] | 0,991 |
| total | n_bus >= 118 | 5 | 1,82 | [1,45; 2,20] | 0,988 |
| nnz(L+U), COLAMD | n_bus >= 118 | 5 | 1,16 | [1,02; 1,30] | 0,996 |
| nnz(J) | n_bus >= 118 | 5 | 1,08 | [0,98; 1,18] | 0,997 |
| dense |V| block entries (n_pq x N_u) | n_bus >= 118 | 5 | 2,04 | [1,95; 2,12] | 1,000 |

*Script: `experiments/c7_scalability.py` · CSV: `results/c7_timing.csv`, `results/c7_fill.csv`, `results/c7_fit.csv`*

Lectura:

- **La factorización es barata y casi lineal.** El tiempo de LU crece como N^1,04 (IC 95 % [0,94; 1,13]) y el
  relleno, nnz(L+U), como N^1,15 (IC [1,07; 1,23]). El cociente (L+U)/J va de 1,6 a 3,0 con COLAMD y de 1,35 a
  1,93 con MMD_AT_PLUS_A, que deja cerca de un tercio menos de no-ceros en los sistemas grandes con un tiempo de
  factorización similar.
- **El costo del método no está en la LU sino en extraer las sensibilidades.** Las N_u resoluciones triangulares
  son el 57 % del total con 300 barras, el 86 % con 2.869 y el 89 % con 9.241. Crecen como
  N^2,00 (IC [1,78; 2,22]) porque el resultado es denso: el bloque de |V| tiene n_PQ × N_u entradas
  (576 MB con 9.241 barras). El revisor 3 tenía razón en el punto 8.
- **Exponente del método completo: 1,51 con los seis sistemas (IC 95 % [0,98; 2,04]) y 1,82 sin el de
  33 barras (IC [1,45; 2,20]).** Con cinco o seis puntos los intervalos son anchos. Lo que sí se puede afirmar:
  α ≈ 1,2 describe el relleno de la factorización, no el método, y el comportamiento asintótico es cuadrático
  mientras se necesite el bloque denso.
- En valores absolutos: 0,89 s para 2.869 barras y 11,7 s para 9.241 en un hilo. Con las resoluciones
  repartidas en hilos (SuperLU libera el GIL) esa etapa baja a 0,41 s y 4,6 s.
- El flujo de potencia base (pandapower, ya en caliente) cuesta entre 0,04 y 0,4 s y hay que sumarlo: el método
  parte de un punto AC convergido.
- No hice el ajuste por separado para el ordenamiento alternativo ni probé más allá de 9.241 barras.

---

## 13. C8 — Comparación con otros métodos

Todos los métodos de intervalo usan el mismo modelo lineal corregido. Punto de diseño, N = 20.000 flujos AC
(las mismas muestras que B1). "Rango AC 99,73 %" es el semiancho empírico entre los cuantiles 0,135 % y
99,865 % de las muestras AC.

| Sistema | Método | Tipo de garantía | Ancho medio (p.u.) | Ancho / rango AC 99,73 % (mediana) | Violación AC media (%) | máx. (%) | Barras > 0,27 % | Sistema (%) | Tiempo (s) |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| IEEE-33 | affine_box | deterministic (linear model) | 1,724e-02 | 3,09 | 0,0000 | 0,000 | 0 | 0,00 | 0,0040 |
| IEEE-33 | cltm_3sigma | 99.73 % per bus (marginal) | 6,328e-03 | 1,05 | 0,1512 | 0,275 | 1 | 0,72 | 0,0041 |
| IEEE-33 | cdhim_lmax | hybrid: 99.73 % joint per branch / deterministic | 1,194e-02 | 1,64 | 0,0006 | 0,020 | 0 | 0,02 | 0,0045 |
| IEEE-33 | cdhim_coord | hybrid: 99.73 % joint per branch / deterministic | 1,016e-02 | 1,23 | 0,0116 | 0,050 | 0 | 0,10 | 0,0045 |
| IEEE-33 | ac_mc_central_99.73 | empirical | 5,991e-03 | 1,00 | 0,2742 | 0,290 | 22 | 1,38 | 0,3462 |
| IEEE-33 | ac_mc_full_range | empirical | 8,036e-03 | 1,32 | 0,0000 | 0,000 | 0 | 0,00 | 0,3462 |
| PEGASE-2869 | affine_box | deterministic (linear model) | 8,061e-03 | 4,12 | 0,0002 | 0,410 | 1 | 0,42 | 0,8556 |
| PEGASE-2869 | cltm_3sigma | 99.73 % per bus (marginal) | 2,010e-03 | 1,06 | 0,1333 | 1,190 | 80 | 41,84 | 0,8356 |
| PEGASE-2869 | cdhim_lmax | hybrid: 99.73 % joint per branch / deterministic | 4,303e-03 | 1,61 | 0,0024 | 0,410 | 1 | 3,15 | 0,9231 |
| PEGASE-2869 | cdhim_coord | hybrid: 99.73 % joint per branch / deterministic | 3,847e-03 | 1,23 | 0,0153 | 0,630 | 2 | 9,31 | 0,9231 |
| PEGASE-2869 | ac_mc_central_99.73 | empirical | 1,870e-03 | 1,00 | 0,2733 | 1,250 | 1468 | 78,89 | 60,0622 |
| PEGASE-2869 | ac_mc_full_range | empirical | 2,447e-03 | 1,30 | 0,0000 | 0,000 | 0 | 0,00 | 60,0622 |

**Tiempos (un hilo):**

| Cantidad | IEEE-33 | PEGASE-2869 |
|---|---|---|
| Jacobiano + LU + resoluciones (s) | 0,0040 | 0,8304 |
| Etapa por ramas (s) | 0,00049 | 0,09268 |
| CD-HIM total, sin flujo base (s) | 0,0045 | 0,9231 |
| Flujo base, pandapower (s) | 0,040 | 0,145 |
| Monte Carlo AC por lotes, N muestras (s) | 0,3 | 60,1 |
| ídem, ms por muestra | 0,017 | 3,003 |
| pandapower.runpp, ms por flujo (protocolo original) | 52,6 | 189,0 |
| Lazo pandapower extrapolado a N muestras (s) | 1051 | 3780 |
| Aceleración frente al Monte Carlo por lotes | 76,3 | 65,1 |
| Aceleración frente al lazo pandapower (extrapolado) | 231796 | 4095 |

*Script: `experiments/c8_methods.py` · CSV: `results/c8_methods.csv`, `results/c8_timing.csv` · N = 20.000*

Lectura:

- **Ancho frente al rango empírico AC al 99,73 % (mediana sobre barras):** caja afín 3,1× y 4,1×; caja 3σ 1,05× y
  1,06×; CD-HIM con λ_max 1,64× y 1,61×; CD-HIM con proyección exacta 1,23× en ambos sistemas.
- **La caja 3σ reproduce casi exactamente el rango empírico**, con violación media 0,15 % y 0,13 %. Su problema
  es la cola: 80 barras sobre 0,27 % en PEGASE y 42 % de las muestras con alguna barra fuera.
- **CD-HIM es más ancho que la caja 3σ y viola menos.** No es que la elipse aproveche la correlación para estrechar
  los intervalos por barra: el cuantil conjunto (3,44σ en vez de 3σ) y λ_max los ensanchan. La correlación solo
  se aprovecha si lo que se entrega es el conjunto bidimensional (por ejemplo, para acotar una caída de voltaje o
  un flujo por la rama), no dos intervalos.
- **Todos los métodos lineales cuestan lo mismo**, porque comparten las sensibilidades: 0,92 s para CD-HIM en
  PEGASE, de los cuales la etapa propia (ramas) es el 10 %. La Tabla II del paper (CLTM ≈ 100× más
  rápido que el Monte Carlo y CD-HIM 4,2× más lento que CLTM) no se puede reproducir.
- **El speed-up depende del Monte Carlo con que se compare.** Frente al solver AC por lotes de este paquete:
  76× (IEEE-33) y 65× (PEGASE). Frente al lazo `pandapower.runpp(numba=False)` del script original,
  que aquí tarda 53 y 189 ms por flujo: 231.796× y 4.095×, extrapolando desde 200 flujos.
  El paper midió 18,9 y 62,3 ms por flujo; no pude reproducir esos tiempos con pandapower 3.5.6. Recomiendo
  reportar tiempos absolutos y decir qué implementación de Monte Carlo se usó.
- Las muestras AC por lotes son soluciones exactas del mismo problema (desajuste < 10⁻⁹ p.u.); una prueba unitaria
  las compara con pandapower.

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

---

## 14. Cifras del paper que cambian

"Código original" es lo que da `finalFinalCode.py` tal como está en el repositorio (reproducido). "Corregido" es
el paquete nuevo con el Jacobiano exacto en el punto AC convergido.

| Cifra | Paper | Código original | Corregido | Causa de la diferencia |
|---|---|---|---|---|
| Ramas en modo Elipse, IEEE-33 | 51,6 % (16/31) | 51,6 % (16/31) | **67,7 % (21/31)** | Jacobiano (signo de ∂P/∂θ y punto de linealización) |
| Ramas en modo Elipse, PEGASE | 92,7 % | 92,73 % (2.984/3.218) | **79,07 % (2.543/3.216)** | Ídem. El total baja a 3.216 porque el filtro absoluto σ_xσ_y > 10⁻¹² se reemplazó por un piso de ruido numérico en el radio |
| CRI, IEEE-33 | 86,92 % | 86,70 % | **86,37 %** | Jacobiano. El 86,92 % no se reproduce con el código del repositorio |
| CRI, PEGASE | 58,18 % | 57,80 % | **71,08 %** | Ídem |
| Qué mide el CRI | "reducción de volumen frente a las cotas afines independientes" | Área de elipse 1σ frente a caja de covarianza 1σ | La salida real se reduce **42,8 % (IEEE-33) y 34,1 % (PEGASE) en ancho** frente a la caja afín | Definición (B2). Frente a la caja 3σ la salida es más ancha |
| R_viol contra AC, IEEE-33 | 0,00 % | No se calcula. Con sus cotas: media 0,0002 %, máx. 0,005 % | **media 0,0006 %, máx. 0,020 %** | El paper acierta en este sistema |
| R_viol contra AC, PEGASE | 0,06 % "en el punto de diseño" | No se calcula; el 0,06 % es del modelo lineal. Con sus cotas: **media 3,69 %, máx. 89,4 %** | **media 0,0024 %, máx. 0,41 %** por barra; elipse: media 0,103 %, máx. 2,69 % | La validación AC no existía; el Jacobiano original no linealiza la red |
| Tabla I, IEEE-33 (5/10/20/30 %) | 0,00 / 0,06 / 0,08 / 0,10 % | 0,00 / 0,00 / 0,02 / 0,04 % | **0,107 % para todo δ** | Regularización de 10⁻⁹ en Mahalanobis. La tabla mide la elipse en el modelo lineal, no la linealización |
| Tabla I, PEGASE | 0,06 / 0,08 / 0,12 / 0,12 % | 0,04 / 0,04 / 0,10 / 0,10 % | **0,070 % para todo δ** | Ídem |
| Violación Beta, PEGASE | 0,19 % (máximo de 3 ramas) | 0,23 % (3 ramas en modo Caja, con regularización) | **0,83 % de media, 0,95 % máx.**, todas las ramas Elipse, con la transformación del paper; 0,29 % con Beta(2,5) dentro de [−1, 1] y su varianza real | Regularización, selección de ramas y transformación (C3) |
| Violación Beta, IEEE-33 | hasta 0,40 % | 0,40 % | **0,84 % de media, 0,90 % máx.**; 0,31 % con Beta(2,5) bien planteada | Ídem |
| Curva P-P (Fig. 5) | Sobre la diagonal ("conservadora") | 0,73 a confianza 0,50 (PEGASE) | **0,45 a confianza 0,50**: bajo la diagonal hasta ≈ 0,80, sobre ella después | Regularización de 10⁻⁹ |
| Tiempo CD-HIM, PEGASE | 52,05 s | 58–66 s medidos aquí; solo el lazo por ramas, LU denso | **0,92 s** en un hilo (Jacobiano + LU + resoluciones + ramas); 0,53 s con las resoluciones en hilos | Álgebra dispersa y vectorización |
| Tiempo CD-HIM, IEEE-33 | 0,0037 s | 0,0047 s | **0,0045 s** | Ídem (en 33 barras domina el armado del Jacobiano) |
| Monte Carlo AC, 20.000 flujos, PEGASE | 1.246,66 s (y 1.263,48 s) | No terminó en la corrida de referencia (máquina compartida) | **60 s** por lotes; lazo pandapower: 189 ms por flujo medidos, 3.780 s extrapolados | Solver por lotes; el tiempo de pandapower depende de la versión |
| Speed-up, PEGASE | 24,0× (24,8× en IV-D) | — | **65×** frente al Monte Carlo por lotes; 4.095× frente al lazo pandapower extrapolado | La cifra depende de cómo se implemente el Monte Carlo |
| Speed-up, IEEE-33 | 1,08·10⁵ (extrapolado) | — | **76×** frente al Monte Carlo por lotes | Ídem: el 10⁵ compara contra el sobrecosto de pandapower en una red de 33 barras |
| "Sub-second" | Afirmado con 52,05 s | — | 0,92 s en un hilo, 0,53 s con hilos, más 0,15 s del flujo base | — |
| CLTM: 12,5 s, ≈ 100×, sobrecosto 4,2× | Tabla II | Sin código | La caja 3σ cuesta lo mismo que CD-HIM menos la etapa por ramas (10 % del total en PEGASE) | No hay implementación de CLTM en el repositorio |
| Complejidad | O(N^α), α ≈ 1,2 | LU denso, O(N³) | LU: N^1,04; relleno: N^1,15; **método completo: N^1,51** (IC 95 % [0,98; 2,04]), dominado por las resoluciones, N^2,00 | El bloque de sensibilidades es denso (C7) |
| γ_th = 15 % como óptimo | Estudio OAT en IEEE-33 | El umbral no actúa: 5, 10, 15, 20 y 25 % dan la misma selección en IEEE-33 | Con la caja del paper sigue sin actuar; con la caja 3σ actúa y depende del sistema | G ≥ 21,46 % para todo ρ cuando κ = π/4 (A1, C4, C6) |
| Rango de validez | δ ≤ 20 % | Sin respaldo | Depende del sistema y del criterio (C1) | — |
| Robustez por escala ("rewards system scale") | Afirmada | Sin respaldo | Mismas tasas en 64 y en 2.951 fuentes (C3) | El TCL no actúa a escala de rama |

---

## 15. Puntos bloqueados y qué se necesita de los autores

| # | Punto | Estado | Qué se necesita |
|---|---|---|---|
| 1 | Datos reales de carga y renovables (C2, R1.6) | Función lista (`estimate_source_covariance`), experimento sin correr | Qué serie usar: fuente, resolución temporal, periodo y cómo se asignan las columnas a las barras de cada sistema |
| 2 | Ref. [11], conjuntos híbridos caja-elipsoide (C8, R1.3) | PDF recibido y leído (ver C8). No hay código en el repositorio y no la implementé | Decisión de implementarla. Requiere un flujo de potencia por intervalos basado en optimización (dos problemas no lineales y no convexos por variable de salida) y, para fuentes correlacionadas, un programa con restricciones de azar. Estimación: 2 a 3 semanas con validación; en PEGASE-2869 son miles de problemas de optimización, y la Ref. [11] solo lo muestra en IEEE-30 |
| 3 | Ref. [12], reflexión de escenarios (C8, R1.3) | No hay código; no pude acceder ni al resumen | El PDF. Sin leerlo no puedo estimar bien; orden de magnitud, 1 a 2 semanas |
| 4 | CLTM, Ref. [10], variante AC completa (C8) | PDF recibido y leído: multiplicador (3σ) y modelo de varianza (uniforme) confirmados. Faltan el paso de escenarios de optimización y la regla de M = 6 términos efectivos | Decisión de implementar esos dos pasos. Estimación: 1 a 2 semanas. Sin ellos, el manuscrito no puede llamar "CLTM" a la caja 3σ |
| 5 | CRI 86,92 % y 58,18 % | El repositorio da 86,70 % y 57,80 % | Con qué versión del script y de pandapower se generaron las tablas del paper |
| 6 | "12,5 s" de CLTM y "n = 4 stable runs" | Sin código que los produzca | Origen de esas mediciones, o retirarlas. La Ref. [10] no trae ningún tiempo sobre PEGASE: su caso grande es el sistema polaco de 2.383 barras en modelo DC (cerca de 3 s, frente a cerca de 1.250 s de Monte Carlo) |
| 7 | Selección de la rama de "máximo estrés" en dos etapas | El código toma la rama de mayor G con \|ρ\| < 0,999; el cribado descrito en el paper no existe | Implementarlo o corregir el texto. Los experimentos nuevos reportan todas las ramas en modo Elipse, así que la selección deja de importar |
| 8 | `caja-elipse.tex` | No está en el repositorio | No pude contrastar el texto de los `[TODO]`; usé la tabla del encargo |
| 9 | Monte Carlo AC de PEGASE con el script original | Detenido tras más de dos horas con la máquina compartida | Nada: la medición aislada del mismo protocolo está en C8 |
| 10 | Revisión bibliográfica (R1.7: procesos gaussianos, caos polinomial, bajo rango) y definición de δ en el texto (R1.5) | Fuera del alcance del código | Trabajo de redacción |

---

## 16. Estado de cada crítica con la evidencia nueva

| Crítica | Experimento | Qué muestra | Qué queda por hacer |
|---|---|---|---|
| Editor; R1.1; R2.5; R3.1: la linealización no está validada | B1, C1 | Error de linealización, cobertura de la elipse y violación total medidos por separado en 25 celdas por sistema. Las cotas del manuscrito no pasan en PEGASE; las corregidas pasan en promedio y fallan en la peor barra | Reescribir la Tabla I y las afirmaciones de validez con el criterio explícito |
| R1.4: condiciones de alta carga | C1 | Hasta el 95 % del límite de cargabilidad, con tasa de no convergencia | — |
| R2.1; R3.2: G es función de ρ; tres cajas distintas | A1, C4 | Confirmado. Con la caja del código el umbral no actúa. Las tres cajas tienen ahora nombre y función propios | Decidir entre regla de ρ o G frente a la caja afín |
| R2.3; R3.3: λ_max frente a la proyección exacta; objetos mezclados | B2 | Factor mediano √2, hasta 136. Opción `radius: coordinate` implementada | Elegir el radio del método; aclarar los tres objetos en el texto |
| R3.3; R3.7: el CRI no se calcula sobre la salida | B2 | CRI de ancho sobre los intervalos de salida implementado junto al de área | Usar un solo CRI en el paper |
| R2.2; R3.4: fuentes correlacionadas | C2 | Q = C·Σ_ε·Cᵀ implementado. Con Σ_ε = I/3 y fuentes correlacionadas la violación es 26–59 % | Datos reales (bloqueado 1) |
| R3.5; R2.6: prueba Beta y otras distribuciones | C3 | Transformación corregida, 12 distribuciones, todas las ramas. La robustez afirmada no se sostiene | Reescribir la sección IV-E |
| R2.4; R3.6: garantías mezcladas; nivel sistema | B1, C5 | 3 % (barras) y 50 % (elipses) de las muestras con alguna violación en PEGASE. Bonferroni cuantificado | Enunciar la garantía con precisión: por rama, en el modelo lineal |
| R3.7: inconsistencias numéricas | A4, sección 14 | Todas localizadas | Corregir el texto con la sección 14 |
| R3.8: complejidad y escalabilidad | A3, C7 | LU disperso real; seis sistemas; exponentes con intervalo de confianza | Retirar α ≈ 1,2 para el método completo |
| R2.7: transferibilidad de γ_th | C6 | Cinco sistemas. Con la caja del paper el umbral no actúa; con otra, depende del sistema | Redefinir o eliminar γ_th |
| R1.2: parámetros heurísticos | A1, C4, C6 | γ_th no tiene efecto; δ_stab = 10⁻³ es lo que decide | Justificar δ_stab o reemplazarlo |
| R1.3: comparación con [11] y [12] | C8 | Solo caja afín y caja 3σ | Bloqueados 2, 3 y 4 |
| R1.5: definición de δ | A2 | Semiancho relativo de cada inyección, igual para todas | Escribirlo en el paper |
| R1.6: datos reales | C2 | Función preparada | Bloqueado 1 |
| R1.7: revisión de flujo probabilístico | — | — | Redacción |

---

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

| Cantidad | IEEE-33 | PEGASE-2869 |
|---|---|---|
| Símbolos de ruido N_u | 64 | 2951 |
| Barras PQ con radio afín > 10⁻⁶ p.u. | 32 | 2312 |
| Error de linealización máx. / radio afín: mediana sobre barras (%) | 0,25 | 0,55 |
| ídem, peor barra (%) | 0,36 | 15,53 |
| Error de linealización máx. / radio 3σ: mediana (%) | 0,78 | 2,04 |
| ídem, peor barra (%) | 0,8 | 101,0 |
| Factor en que baja el error al sumar εᵀH_iε: mediana | 75 | 115 |
| ídem, mínimo | 69 | 35 |
| Barras con H_i ensamblada | 32 | 24 |
| Parte negativa del espectro de H_i: mediana (%) | 100,0 | 99,8 |
| Mayor valor positivo alcanzado / mayor valor negativo alcanzado: máximo (%) | 0,00 | 11,56 |
| Autovalores que juntan el 90 % del espectro: mediana | 5,0 | 11,0 |

*Script: `experiments/x1_aa_second_order.py` · CSV: `results/x1_summary.csv`, `results/x1_bus_remainder.csv` ·
N = 20.000 flujos AC por sistema · δ de diseño (15 % y 5 %)*

- Sumar εᵀH_iε a la forma afín baja el error contra el flujo AC en un factor mediano de 75 (IEEE-33) y
  115 (PEGASE). El error de linealización medido en C1 es, casi entero, este término.
- **El resto tiene un solo signo.** H_i es semidefinida negativa en las 32 barras de IEEE-33 y lo es casi por
  completo en PEGASE (parte negativa del espectro: 99,8 % en mediana). La forma afín sobreestima el
  voltaje: es una envolvente superior, y solo la cota inferior necesita corrección.
- Pocas direcciones concentran el espectro (5 y 11 autovalores para el 90 %, en mediana),
  con N_u = 64 y 2.951.

### 17.2 La caja afín no encierra la solución AC

Para cada barra se resuelve el flujo AC exacto en los dos vértices de la caja de ruido donde su forma afín es
extrema: ε = −sign(c_i) (todas las fuentes en el extremo que baja el voltaje de esa barra) y ε = +sign(c_i).
No hay muestreo ni semilla.

| Cantidad | IEEE-33 | PEGASE-2869 |
|---|---|---|
| Vértice inferior: barras con la solución AC fuera de la caja afín | 32 | 2270 |
| Exceso / radio afín: mediana (%) | 1,36 | 9,20 |
| ídem, percentil 90 (%) | 1,44 | 58,38 |
| ídem, máximo (%) | 1,45 | 536,55 |
| Barras con exceso > 10 % del radio | 0 | 1130 |
| Barras con exceso > 100 % del radio | 0 | 70 |
| Parte del exceso que εᵀH_iε no explica: mediana (%) | 2,8 | 3,0 |
| ídem, percentil 90 (%) | 2,9 | 16,9 |
| Vértice superior: barras con la solución AC fuera de la caja afín | 0 | 3 |
| Vértice superior: holgura / radio afín, mediana (%) | 1,29 | 8,75 |
| Vértices inferiores sin solución AC | 0 | 0 |
| Vértices superiores sin solución AC | 0 | 0 |

*CSV: `results/x1_bus_remainder.csv` (columnas `vertex_*`)*

- **En el vértice inferior, la solución AC queda por debajo de la caja afín en 32 de 32
  barras de IEEE-33 y en 2.270 de 2.312 de PEGASE.** El exceso es chico en IEEE-33
  (1,4 % del radio en mediana, 1,4 % en la peor barra) y grande en PEGASE:
  9,2 % en mediana, 58,4 % en el percentil 90 y 537 % en la peor barra.
  1.130 barras exceden el radio en más de 10 % y 70 en más de 100 %.
- El Monte Carlo de B1 y C1 no lo ve: su mayor error es 15,5 % del radio afín. Con 2.951 símbolos
  independientes los vértices son inalcanzables por muestreo. La caja afín "pasa" el Monte Carlo porque es entre
  tres y cuatro veces más ancha que el rango al 99,73 %, no porque sea una cota.
- La forma cuadrática predice ese exceso: en mediana deja sin explicar el 2,8 % en IEEE-33 y el
  3,0 % en PEGASE.
- En el vértice superior la solución AC queda dentro en todas las barras de IEEE-33 y en todas menos
  3 de PEGASE, como corresponde a un resto negativo.
- **Cuidado con la lectura física del vértice en PEGASE.** Es el escenario en que todas las cargas se mueven a la
  vez hacia el mismo lado y la barra *slack* absorbe todo el desbalance. El radio afín del ángulo es
  0,964 rad en mediana, contra 0,0024 p.u. de la magnitud. El caso es válido dentro de la caja que
  el propio modelo declara, pero es extremo; con *slack* distribuida sería bastante menor. No lo probé.

### 17.3 Cuánto cuesta acotar el resto

Sobre todas las barras, sin ensamblar H_i (costo N·N_u por barra):

| Cantidad | IEEE-33 | PEGASE-2869 |
|---|---|---|
| Producto afín estándar: ancho de la cota / radio afín, mediana (%) | 5,0 | 37520,4 |
| Binomios escalados: ídem (%) | 4,3 | 34051,1 |
| Binomios / producto estándar, mediana | 0,842 | 0,904 |
| De la cota estándar: productos dV·dI (%) | 93,2 | 19,2 |
| ídem: cuadrados en barras PV (%) | 0,0 | 41,9 |
| ídem: curvatura de \|V\| (%) | 6,8 | 38,8 |
| Radio afín del ángulo, mediana (rad) | 0,0085 | 0,9639 |
| Radio afín de la magnitud, mediana (p.u.) | 0,0100 | 0,0024 |
| Profundidad alcanzada / radio afín, barras de mayor error: mediana (%) | – | 123 |
| ídem, máximo (%) | – | 430 |
| Profundidad alcanzada / radio afín, barras al azar: mediana (%) | – | 10,3 |
| ídem, máximo (%) | – | 78,6 |
| Profundidad alcanzada / radio afín, todas las barras: mediana (%) | 1,33 | – |
| ídem, máximo (%) | 1,41 | – |

Sobre un subconjunto con H_i ensamblada (las 32 barras de IEEE-33; en PEGASE, 8 con el mayor error de
linealización relativo y 16 al azar):

| Cota de εᵀH_iε | Costo por barra | IEEE-33: ancho / alcanzado, mediana | máx. | lado positivo / profundidad (%) | PEGASE-2869: ancho / alcanzado, mediana | máx. | lado positivo / profundidad (%) |
|---|---|---|---|---|---|---|---|
| Producto afín estándar, término a término | N·N_u | 3,79 | 3,95 | 202 | 3620,28 | 82286,01 | 180405 |
| Los mismos términos como cuadrados de binomios escalados | N·N_u | 3,22 | 3,32 | 158 | 3585,23 | 81749,33 | 179331 |
| H_i ensamblada, monomio a monomio (k = 1, AF2) | N_u² | 1,88 | 1,91 | 88 | 2,00 | 2,33 | 100 |
| H_i ensamblada, bloques binomiales con presupuesto diagonal (k = 2) | N_u² log N_u | 1,78 | 1,82 | 78 | 1,99 | 2,20 | 99 |
| Teorema k-nomial: 1 bloque espectral + residuo monomial | N_u² + 1 autovector | 1,77 | 2,02 | 53 | 1,55 | 2,25 | 36 |
| Teorema k-nomial: 8 bloques espectrales + residuo monomial | N_u² + 8 autovectores | 1,61 | 1,99 | 2 | 1,30 | 2,50 | 6 |
| Teorema k-nomial: todos los bloques espectrales, sin residuo | N_u³ | 1,62 | 2,03 | 0 | 1,35 | 3,10 | 0 |

*CSV: `results/x1_quadratic_bounds.csv` · "ancho / alcanzado": ancho de la cota sobre el ancho del rango alcanzado
por búsqueda local · "lado positivo / profundidad": cota superior sobre el mayor valor negativo alcanzado; como
el resto es casi siempre negativo, lo que aparezca ahí es sobreestimación pura*

- **La evaluación afín estándar, producto a producto, no sirve a escala.** En IEEE-33 queda a 3,8 veces
  el rango alcanzado. En PEGASE, a 3.620 veces en mediana: la cota del resto es 375 veces el
  radio afín. Escribir cada producto como diferencia de cuadrados de binomios escalados (dos bloques k-nomiales
  por producto, nunca peor que el producto estándar) casi no cambia eso.
- **La causa son las cancelaciones entre términos, no cada producto.** Con ángulos que se mueven
  0,964 rad, la rotación de los fasores aparece en coordenadas rectangulares como términos de segundo
  orden grandes (cuadrados en barras PV y curvatura de |V|: 42 % y 39 % de la cota
  estándar, en mediana) que se anulan entre sí al sumar. Acotar término a término pierde esa cancelación. En IEEE-33, sin barras PV y con ángulos
  de 0,009 rad, el efecto no existe.
- **Con H_i ensamblada el problema es tratable.** Monomio a monomio (AF2): 1,88× y 2,00× el
  rango alcanzado. Bloques binomiales con presupuesto diagonal: 1,78× y 1,99×; ganan poco porque
  cada bloque solo puede absorber lo que le alcance la diagonal y aquí domina lo que está fuera de ella.
- **El teorema k-nomial con bloques grandes es lo que mejora.** Un solo bloque espectral más el residuo:
  1,77× y 1,55×. Ocho bloques: 1,61× y 1,30×. Y como el cuadrado de un k-nomio
  tiene rango de un solo lado, desaparece el lado positivo espurio: de 100 % de la profundidad
  con AF2 a 0,2 % con todos los bloques, en PEGASE.
- **Lo que no está resuelto es el costo.** Ensamblar H_i es N_u² por barra (70 MB en PEGASE) y los bloques
  espectrales necesitan autovectores. Para todas las barras de un sistema grande hace falta elegir los bloques
  sin formar H_i, o cambiar la descomposición para que las rotaciones no aparezcan (diferencias por rama, o
  coordenadas polares con resto acotado). Ese es el problema abierto.
- **Tamaño del resto en el peor caso.** En IEEE-33 la profundidad alcanzada es 1,3 % del radio afín. En
  PEGASE, entre las 8 barras de mayor error, 1,2 veces el radio en mediana y
  4,3 en la peor; entre las 16 al azar, 0,10 veces en mediana y
  0,79 en la mayor.

### 17.4 Capa probabilística sobre la forma afín

| Sistema | Intervalo | Violación AC media (%) | máx. (%) | Barras > 0,27 % | Muestras con ≥ 1 barra fuera (%) |
|---|---|---:|---:|---:|---:|
| IEEE-33 | Caja 3σ, centrada en V* | 0,1306 | 0,215 | 0 | 0,68 |
| IEEE-33 | Caja 3σ, centro corrido en tr(H_i)/3 | 0,1302 | 0,215 | 0 | 0,69 |
| IEEE-33 | Proyección exacta √(χ²Q_ii), centrada en V* | 0,0184 | 0,040 | 0 | 0,10 |
| IEEE-33 | Proyección exacta, centro corrido | 0,0183 | 0,040 | 0 | 0,10 |
| IEEE-33 | Radio de Hoeffding, centrado en V* | 0,0000 | 0,000 | 0 | 0,00 |
| IEEE-33 | Caja afín, centrada en V* | 0,0000 | 0,000 | 0 | 0,00 |
| IEEE-33 | Caja afín, centro corrido | 0,0000 | 0,000 | 0 | 0,00 |
| PEGASE-2869 | Caja 3σ, centrada en V* | 0,1419 | 1,150 | 158 | 42,27 |
| PEGASE-2869 | Caja 3σ, centro corrido en tr(H_i)/3 | 0,1407 | 0,910 | 149 | 41,77 |
| PEGASE-2869 | Proyección exacta √(χ²Q_ii), centrada en V* | 0,0241 | 0,570 | 2 | 9,93 |
| PEGASE-2869 | Proyección exacta, centro corrido | 0,0237 | 0,480 | 2 | 9,73 |
| PEGASE-2869 | Radio de Hoeffding, centrado en V* | 0,0000 | 0,005 | 0 | 0,01 |
| PEGASE-2869 | Caja afín, centrada en V* | 0,0002 | 0,380 | 1 | 0,38 |
| PEGASE-2869 | Caja afín, centro corrido | 0,0001 | 0,210 | 0 | 0,21 |

| Cantidad | IEEE-33 | PEGASE-2869 |
|---|---|---|
| Radio de Hoeffding / radio afín, mediana | 0,715 | 0,540 |
| Radio de Hoeffding / radio 3σ | 2,099 | 2,099 |
| Corrimiento del centro \|tr(H_i)/3\| / radio 3σ: mediana (%) | 0,14 | 0,37 |
| ídem, máximo (%) | 0,16 | 5,47 |
| Corrimiento predicho frente a la media del Monte Carlo: diferencia máxima / radio 3σ (%) | 0,0001 | 0,0323 |
| Barras con menos de 6 términos efectivos (regla M = 6 de la Ref. [10]) | 0 | 143 |

*CSV: `results/x1_summary.csv` · mismas 20.000 muestras AC que 17.1*

- **Correr el centro en la media del término cuadrático, tr(H_i)/3, ayuda poco.** La peor barra de la caja 3σ pasa
  de 1,15 % a 0,91 % en PEGASE y el nivel sistema casi no se mueve (42,3 % y
  41,8 %). El segundo orden no arregla la caja 3σ: su problema a nivel sistema es que el 99,73 % es
  por barra.
- **Radio de Hoeffding**, ‖c_i‖₂·√(2 ln(2/α)): vale para cualquier distribución de símbolos independientes, de
  media cero y soporte en [−1, 1], sin teorema central del límite. Es 2,10 veces el radio 3σ y
  0,71 (IEEE-33) y 0,54 (PEGASE) veces el radio afín. Contra AC: 0,01 % de las
  muestras de PEGASE con alguna barra fuera, frente a 42,3 % de la caja 3σ. No cubre símbolos
  correlacionados ni soportes no acotados.
- Regla M = 6 de la Ref. [10]: 143 barras de PEGASE tienen menos de 6 términos efectivos (ninguna en
  IEEE-33). En esas barras la Ref. [10] no usaría el intervalo 3σ.

### 17.5 Dependencia por rama: el equivalente afín de la ganancia G

| Cantidad | IEEE-33 | PEGASE-2869 |
|---|---|---|
| Ramas PQ-PQ | 31 | 3216 |
| Área del zonotopo frente a la caja afín: reducción mediana (%) | 95,0 | 89,8 |
| ídem, mínima (%) | 58,4 | 3,2 |
| Área del zonotopo / área de la elipse al 99,73 %: mediana | 3,79 | 7,60 |
| Ramas donde el zonotopo es más chico que la elipse | 0 | 119 |
| Radio de V_i − V_j: reducción de la forma afín frente a la resta de intervalos, mediana (%) | 95,9 | 87,8 |
| ídem, mínima (%) | 29,5 | 0,0 |
| ídem, máxima (%) | 99,8 | 100,0 |

*CSV: `results/x1_branch_dependency.csv`*

- El conjunto exacto de la pareja (V̂_i, V̂_j) es un zonotopo, no una caja ni una elipse. Ocupa 95 %
  (IEEE-33) y 90 % (PEGASE) menos área que la caja afín, en mediana. Es la versión determinista de lo
  que G quería medir, depende de todos los generadores y no solo de ρ, y no necesita Σ_ε ni independencia: la
  objeción de las fuentes correlacionadas (C2) no le aplica.
- A cambio es más grande que la elipse al 99,73 %: 3,8 y 7,6 veces en mediana (más chico
  en 119 de 3.216 ramas de PEGASE).
- **Donde la dependencia rinde es en cantidades de rama.** El radio de V_i − V_j calculado con la forma afín es
  96 % y 88 % menor, en mediana, que restando los dos intervalos. Los intervalos por barra
  no pueden mejorar con ningún conjunto por pares (B2); las caídas de voltaje y los flujos sí.

### 17.6 Lectura

1. Hay una orientación desde la aritmética afín que se sostiene con datos: el resto de la forma afín del flujo
   de potencia es una forma cuadrática en los símbolos de ruido, de un solo signo y dominada por pocas
   direcciones, y acotarla bien es justo el problema del teorema k-nomial.
2. La motivación es medible: la caja afín, que solo es determinista para el modelo lineal, no encierra la
   solución AC (17.2), y la aritmética afín estándar aplicada al resto sobreestima en tres órdenes de magnitud
   en un sistema de transmisión (17.3).
3. La ganancia medida del k-nomial sobre AF2, barra a barra, es un factor mediano de 1,17 (IEEE-33) y
   1,51 (PEGASE) en el ancho con ocho bloques, y de 1,07 y 1,29 con uno, más la
   eliminación del lado espurio. No es un orden de magnitud. El salto grande está entre acotar término a término
   y acotar la forma ensamblada, y ese salto hoy cuesta N_u² por barra.
4. La elipse y G no entran en esta orientación. Lo que se conserva de CD-HIM es el motor de sensibilidades
   dispersas y la observación de que las barras vecinas comparten casi todos sus símbolos.

---

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

| Sistema | Balance | K | Vértices con flujo AC | Radio afín, mediana (p.u.) | Mayor \|ΔV\| en un vértice (p.u.) | Radio afín del ángulo, mediana (rad) | Radio afín / radio 3σ, mediana | Barras donde 3σ es más ancho que la caja afín | Caja afín / rango MC 99,73 % | Caja 3σ / rango MC 99,73 % | Caja 3σ: muestras con ≥ 1 barra fuera (%) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| IEEE-33 | slack única | 4 | 16 | 0,0100 | 0,0144 | 0,002 | 1,04 | 11 | 1,21 | 1,17 | 0,0 |
| IEEE-33 | slack única | 8 | 256 | 0,0100 | 0,0144 | 0,002 | 1,39 | 4 | 1,53 | 1,10 | 0,2 |
| IEEE-33 | slack única | 16 | 65536 | 0,0100 | 0,0144 | 0,002 | 1,78 | 3 | 1,94 | 1,10 | 0,3 |
| PEGASE-2869 | participación | 8 | 256 | 0,0018 | 0,0533 | 0,346 | 0,77 | 1768 | 1,17 | 1,54 | 33,8 |
| PEGASE-2869 | participación | 16 | 65536 | 0,0019 | 0,0533 | 0,353 | 0,91 | 1479 | 1,29 | 1,43 | 26,8 |
| PEGASE-2869 | participación | 32 | 710 | 0,0020 | 0,0562 | 0,409 | 1,04 | 1069 | 1,43 | 1,37 | 21,4 |
| PEGASE-2869 | participación | 64 | 1172 | 0,0020 | 0,0583 | 0,429 | 1,16 | 902 | 1,56 | 1,32 | 19,5 |
| PEGASE-2869 | slack única | 8 | 256 | 0,0020 | 0,0423 | 0,810 | 0,73 | 1856 | 1,09 | 1,56 | 47,2 |
| PEGASE-2869 | slack única | 16 | 65536 | 0,0020 | 0,0428 | 0,810 | 0,85 | 1600 | 1,19 | 1,43 | 35,0 |
| PEGASE-2869 | slack única | 32 | 480 | 0,0021 | 0,0436 | 0,836 | 0,95 | 1288 | 1,32 | 1,35 | 26,8 |
| PEGASE-2869 | slack única | 64 | 742 | 0,0021 | 0,0443 | 0,836 | 1,05 | 1030 | 1,44 | 1,31 | 22,6 |

*Script: `experiments/x2_few_sources.py` · CSV: `results/x2_summary.csv`, `results/x2_bus.csv`*

- **Con pocas fuentes la caja 3σ deja de tener sentido.** En PEGASE con K = 8 es más ancha que la caja afín en
  1.768 de 2.322 barras, es decir, se sale del soporte, y aun así el 34 % de
  las muestras tiene alguna barra fuera, porque está centrada en V* y el voltaje real cae por debajo.
- La caja afín queda a entre 1,1 y 1,6 veces el rango de Monte Carlo al 99,73 %. Con
  2.951 símbolos independientes eran 3 a 4 veces (C8).
- **El coste baja a K resoluciones.** El primer orden tarda 0,03 s como máximo en PEGASE; el modelo de
  segundo orden para todas las barras, 3,1 s, y las cotas, hasta 25,6 s (tabla de 18.6).
- El modelo de participación deja el radio del ángulo en 0,43 a 0,51 veces el de la
  slack única, pero no cambia el cuadro: las barras que más se apartan son barras débiles, donde el voltaje cae
  hasta 0,058 p.u. en un vértice.

### 18.2 La caja afín de primer orden

| Sistema | Balance | K | Barras | Barras con algún vértice fuera | Exceso / radio: mediana (%) | P90 (%) | máx. (%) | Exceso máx. (p.u.) | Barras con exceso > 0,001 p.u. | Barras con alguna muestra MC fuera | Muestras MC con ≥ 1 barra fuera (%) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| IEEE-33 | slack única | 4 | 32 | 32 | 1,4 | 1,4 | 1,4 | 0,0002 | 0 | 0 | 0,00 |
| IEEE-33 | slack única | 8 | 32 | 32 | 1,4 | 1,4 | 1,4 | 0,0002 | 0 | 0 | 0,00 |
| IEEE-33 | slack única | 16 | 32 | 32 | 1,4 | 1,4 | 1,4 | 0,0002 | 0 | 0 | 0,00 |
| PEGASE-2869 | participación | 8 | 2322 | 2284 | 5,5 | 32,5 | 666,8 | 0,0216 | 177 | 737 | 24,34 |
| PEGASE-2869 | participación | 16 | 2322 | 2284 | 6,5 | 35,4 | 579,3 | 0,0216 | 188 | 164 | 8,28 |
| PEGASE-2869 | participación | 32 | 2322 | 2284 | 8,1 | 39,2 | 630,6 | 0,0235 | 263 | 82 | 3,92 |
| PEGASE-2869 | participación | 64 | 2322 | 2284 | 8,7 | 40,7 | 682,8 | 0,0248 | 303 | 54 | 2,63 |
| PEGASE-2869 | slack única | 8 | 2312 | 2273 | 6,2 | 53,9 | 723,3 | 0,0219 | 380 | 1569 | 51,16 |
| PEGASE-2869 | slack única | 16 | 2312 | 2273 | 7,0 | 53,5 | 723,3 | 0,0219 | 416 | 788 | 22,05 |
| PEGASE-2869 | slack única | 32 | 2312 | 2273 | 7,5 | 54,6 | 751,3 | 0,0225 | 435 | 164 | 8,02 |
| PEGASE-2869 | slack única | 64 | 2312 | 2273 | 8,4 | 55,1 | 751,6 | 0,0225 | 455 | 72 | 3,56 |

- **Con pocas fuentes la caja afín falla también por muestreo, no solo en los vértices.** En PEGASE con K = 8,
  el 24 % de las muestras (participación) y el 51 % (slack única) tiene alguna barra
  fuera. Con un símbolo por carga el Monte Carlo no lo veía (sección 17.2).
- En los vértices el exceso llega a 0,025 p.u. de voltaje, y entre 177 y
  455 barras de PEGASE exceden su caja en más de 0,001 p.u.
- En IEEE-33 el exceso es 1,4 % del radio en mediana y 1,4 % como máximo, y ninguna
  muestra sale.

### 18.3 Tamaño de cada orden

| Sistema | Balance | K | 2.º orden / radio: mediana (%) | P90 (%) | máx. (%) | 2.º orden máx. (p.u.) | 3.er orden / radio: mediana (%) | P90 (%) | máx. (%) | 3.er orden máx. (p.u.) | 3.er / 2.º orden: mediana (%) | P90 (%) | Barras con H_i semidefinida negativa (%) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| IEEE-33 | slack única | 4 | 1,3 | 1,4 | 1,4 | 0,0002 | 0,04 | 0,04 | 0,04 | 0,00001 | 2,9 | 2,9 | 100 |
| IEEE-33 | slack única | 8 | 1,3 | 1,4 | 1,4 | 0,0002 | 0,04 | 0,04 | 0,04 | 0,00001 | 2,9 | 2,9 | 100 |
| IEEE-33 | slack única | 16 | 1,3 | 1,4 | 1,4 | 0,0002 | 0,04 | 0,04 | 0,04 | 0,00001 | 2,9 | 2,9 | 100 |
| PEGASE-2869 | participación | 8 | 8,2 | 39,1 | 637,6 | 0,0161 | 0,35 | 4,01 | 120,59 | 0,00552 | 3,6 | 16,3 | 42 |
| PEGASE-2869 | participación | 16 | 9,3 | 39,3 | 541,4 | 0,0161 | 0,35 | 3,88 | 120,59 | 0,00552 | 3,4 | 16,6 | 33 |
| PEGASE-2869 | participación | 32 | 10,7 | 44,6 | 586,8 | 0,0172 | 0,40 | 4,45 | 134,17 | 0,00627 | 3,5 | 16,7 | 18 |
| PEGASE-2869 | participación | 64 | 11,2 | 46,9 | 633,5 | 0,0181 | 0,45 | 4,88 | 141,55 | 0,00683 | 3,7 | 16,4 | 11 |
| PEGASE-2869 | slack única | 8 | 9,0 | 51,3 | 641,5 | 0,0182 | 0,38 | 8,28 | 95,32 | 0,00375 | 3,7 | 18,8 | 41 |
| PEGASE-2869 | slack única | 16 | 9,3 | 50,1 | 628,0 | 0,0182 | 0,36 | 8,18 | 95,32 | 0,00375 | 3,6 | 18,6 | 18 |
| PEGASE-2869 | slack única | 32 | 9,9 | 51,3 | 657,8 | 0,0186 | 0,38 | 8,34 | 100,47 | 0,00389 | 3,5 | 18,6 | 15 |
| PEGASE-2869 | slack única | 64 | 10,5 | 52,2 | 658,2 | 0,0186 | 0,38 | 8,34 | 100,47 | 0,00389 | 3,4 | 18,6 | 7 |

- En la barra mediana de PEGASE el segundo orden es 8 a 11 % del radio y el tercero
  0,35 a 0,45 %. El tercer orden es el 3,4 a 3,7 % del segundo en
  mediana y el 16 a 19 % en el percentil 90.
- **En las peores barras la serie converge despacio.** El segundo orden supera varias veces el radio y el
  tercero llega a 0,0037–0,0068 p.u. Son barras de radio grande, no artefactos de radios
  diminutos: en la peor, el radio afín es 0,031 p.u., el segundo orden agrega hasta 0,016 y el tercero
  hasta 0,0055 (participación, K = 8; `results/x2_bus.csv`, barra 8853).
- **El resto ya no es de un solo signo en PEGASE.** H_i es semidefinida negativa en todas las barras de IEEE-33,
  pero solo en el 7 a 42 % de las de PEGASE, y la proporción baja al subir K.

### 18.4 Criterio 1: ¿la cota de segundo orden contiene la solución AC?

| Sistema | Balance | K | Cota del resto | Barras con algún vértice fuera | Déficit máx. (p.u.) | Déficit máx. / radio (%) | Barras con déficit > 10⁻⁴ p.u. | > 10⁻³ p.u. | Barras con alguna muestra MC fuera | Muestras MC con ≥ 1 barra fuera (%) | Criterio 1 |
|---|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---|
| IEEE-33 | slack única | 4 | AF2 (monomio a monomio) | 32 | 0,00001 | 0,04 | 0 | 0 | 0 | 0,00 | False |
| IEEE-33 | slack única | 4 | RAA-K (cliques de hasta K) | 32 | 0,00001 | 0,04 | 0 | 0 | 0 | 0,00 | False |
| IEEE-33 | slack única | 4 | Bloques espectrales | 0 | 0,00000 | 0,00 | 0 | 0 | 0 | 0,00 | True |
| IEEE-33 | slack única | 8 | AF2 (monomio a monomio) | 32 | 0,00001 | 0,04 | 0 | 0 | 0 | 0,00 | False |
| IEEE-33 | slack única | 8 | RAA-K (cliques de hasta K) | 32 | 0,00001 | 0,04 | 0 | 0 | 0 | 0,00 | False |
| IEEE-33 | slack única | 8 | Bloques espectrales | 0 | 0,00000 | 0,00 | 0 | 0 | 0 | 0,00 | True |
| IEEE-33 | slack única | 16 | AF2 (monomio a monomio) | 32 | 0,00001 | 0,04 | 0 | 0 | 0 | 0,00 | False |
| IEEE-33 | slack única | 16 | RAA-K (cliques de hasta K) | 32 | 0,00001 | 0,04 | 0 | 0 | 0 | 0,00 | False |
| IEEE-33 | slack única | 16 | Bloques espectrales | 0 | 0,00000 | 0,00 | 0 | 0 | 0 | 0,00 | True |
| PEGASE-2869 | participación | 8 | AF2 (monomio a monomio) | 480 | 0,00552 | 29,17 | 94 | 26 | 0 | 0,00 | False |
| PEGASE-2869 | participación | 8 | RAA-K (cliques de hasta K) | 480 | 0,00552 | 29,17 | 94 | 26 | 0 | 0,00 | False |
| PEGASE-2869 | participación | 8 | Bloques espectrales | 180 | 0,00532 | 26,92 | 76 | 24 | 0 | 0,00 | False |
| PEGASE-2869 | participación | 16 | AF2 (monomio a monomio) | 388 | 0,00552 | 37,96 | 107 | 26 | 0 | 0,00 | False |
| PEGASE-2869 | participación | 16 | RAA-K (cliques de hasta K) | 388 | 0,00552 | 37,96 | 107 | 26 | 0 | 0,00 | False |
| PEGASE-2869 | participación | 16 | Bloques espectrales | 150 | 0,00522 | 34,01 | 62 | 24 | 0 | 0,00 | False |
| PEGASE-2869 | participación | 32 | AF2 (monomio a monomio) | 329 | 0,00627 | 43,81 | 129 | 29 | 2 | 0,01 | False |
| PEGASE-2869 | participación | 32 | RAA-K (cliques de hasta K) | 329 | 0,00627 | 43,81 | 129 | 29 | 2 | 0,01 | False |
| PEGASE-2869 | participación | 32 | Bloques espectrales | 139 | 0,00583 | 39,73 | 60 | 26 | 2 | 0,01 | False |
| PEGASE-2869 | participación | 64 | AF2 (monomio a monomio) | 321 | 0,00683 | 49,30 | 124 | 29 | 1 | 0,01 | False |
| PEGASE-2869 | participación | 64 | RAA-K (cliques de hasta K) | 321 | 0,00683 | 49,30 | 124 | 29 | 1 | 0,01 | False |
| PEGASE-2869 | participación | 64 | Bloques espectrales | 119 | 0,00636 | 44,92 | 57 | 27 | 1 | 0,01 | False |
| PEGASE-2869 | slack única | 8 | AF2 (monomio a monomio) | 1050 | 0,00375 | 95,32 | 315 | 41 | 4 | 0,01 | False |
| PEGASE-2869 | slack única | 8 | RAA-K (cliques de hasta K) | 1050 | 0,00375 | 95,32 | 315 | 41 | 4 | 0,01 | False |
| PEGASE-2869 | slack única | 8 | Bloques espectrales | 418 | 0,00297 | 87,02 | 182 | 35 | 4 | 0,01 | False |
| PEGASE-2869 | slack única | 16 | AF2 (monomio a monomio) | 992 | 0,00375 | 95,32 | 280 | 41 | 4 | 0,01 | False |
| PEGASE-2869 | slack única | 16 | RAA-K (cliques de hasta K) | 992 | 0,00375 | 95,32 | 280 | 41 | 4 | 0,01 | False |
| PEGASE-2869 | slack única | 16 | Bloques espectrales | 286 | 0,00295 | 85,72 | 150 | 33 | 4 | 0,01 | False |
| PEGASE-2869 | slack única | 32 | AF2 (monomio a monomio) | 632 | 0,00389 | 93,55 | 272 | 42 | 1 | 0,01 | False |
| PEGASE-2869 | slack única | 32 | RAA-K (cliques de hasta K) | 632 | 0,00389 | 93,55 | 272 | 42 | 1 | 0,01 | False |
| PEGASE-2869 | slack única | 32 | Bloques espectrales | 237 | 0,00289 | 77,52 | 127 | 33 | 1 | 0,01 | False |
| PEGASE-2869 | slack única | 64 | AF2 (monomio a monomio) | 570 | 0,00389 | 93,37 | 255 | 42 | 4 | 0,01 | False |
| PEGASE-2869 | slack única | 64 | RAA-K (cliques de hasta K) | 570 | 0,00389 | 93,37 | 255 | 42 | 4 | 0,01 | False |
| PEGASE-2869 | slack única | 64 | Bloques espectrales | 222 | 0,00296 | 75,40 | 115 | 33 | 4 | 0,01 | False |

**No se cumple.**

- **Una cota de segundo orden ajustada falla exactamente por el tercer orden.** En IEEE-33 el lado inferior de
  AF2 y de RAA-K coincide con el mínimo exacto de la forma cuadrática, y la solución AC queda
  5,9·10⁻⁶ p.u. por debajo en las 32 barras. Los bloques espectrales "contienen" solo porque
  sobreestiman la profundidad (18.6); eso es holgura, no garantía.
- **En PEGASE fallan todas.** Con participación y K = 8, los bloques espectrales dejan 180 barras con
  algún vértice fuera, 76 por más de 10⁻⁴ p.u. y 24 por más de 10⁻³ p.u.; el déficit
  máximo es 0,0053 p.u. Con slack única: 418, 182, 35 y
  0,0030 p.u.
- **Por muestreo, en cambio, casi no falla.** Como máximo el 0,01 % de las muestras tiene alguna
  barra fuera de la cota de segundo orden (4 barras a lo más), contra 24 % y
  51 % de la caja afín. La corrección de segundo orden es lo que hace usable el intervalo; lo que
  no da es contención en el peor caso.

Barrido de δ con el menor K de cada sistema:

| Sistema | Balance | K | δ (%) | Exceso de la caja afín / radio: mediana (%) | máx. (%) | 2.º orden / radio, máx. (%) | 3.er orden / radio, máx. (%) | 3.er orden máx. (p.u.) | AF2: barras con algún vértice fuera | AF2: déficit máx. (p.u.) | Espectral: barras fuera | Espectral: déficit máx. (p.u.) | Espectral: barras con déficit > 10⁻⁴ p.u. |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| IEEE-33 | slack única | 4 | 5 | 0,45 | 0,5 | 0,5 | 0,00 | 2,13e-07 | 27 | 2,13e-07 | 0 | 0,00e+00 | 0 |
| IEEE-33 | slack única | 4 | 15 | 1,36 | 1,4 | 1,4 | 0,04 | 5,89e-06 | 32 | 5,89e-06 | 0 | 0,00e+00 | 0 |
| IEEE-33 | slack única | 4 | 30 | 2,81 | 3,0 | 2,8 | 0,17 | 4,89e-05 | 32 | 4,89e-05 | 0 | 0,00e+00 | 0 |
| PEGASE-2869 | participación | 8 | 1 | 1,03 | 105,0 | 127,5 | 3,25 | 3,00e-05 | 280 | 3,00e-05 | 33 | 2,22e-05 | 0 |
| PEGASE-2869 | participación | 8 | 2 | 2,06 | 231,5 | 255,0 | 14,15 | 2,60e-04 | 344 | 2,60e-04 | 106 | 2,29e-04 | 19 |
| PEGASE-2869 | participación | 8 | 5 | 5,53 | 666,8 | 637,6 | 120,59 | 5,52e-03 | 480 | 5,52e-03 | 180 | 5,32e-03 | 76 |
| PEGASE-2869 | slack única | 8 | 1 | 1,18 | 128,5 | 128,3 | 2,94 | 2,30e-05 | 637 | 2,30e-05 | 29 | 9,81e-06 | 0 |
| PEGASE-2869 | slack única | 8 | 2 | 2,37 | 263,7 | 256,6 | 12,48 | 1,96e-04 | 818 | 1,96e-04 | 182 | 1,23e-04 | 4 |
| PEGASE-2869 | slack única | 8 | 5 | 6,17 | 723,3 | 641,5 | 95,32 | 3,75e-03 | 1050 | 3,75e-03 | 418 | 2,97e-03 | 182 |

- El déficit crece como δ³: en PEGASE con participación, 2,2·10⁻⁵ p.u. con δ = 1 %, 2,3·10⁻⁴ con 2 % y
  5,3·10⁻³ con 5 %. Con δ = 1 % ninguna barra pasa de 10⁻⁴ p.u.; con 2 % ya son 19.

### 18.5 Criterio 2: ancho frente al rango de referencia

| Sistema | Balance | K | Rango AC en vértices / caja afín | AF2 / rango AC: mediana | máx. | RAA-K / rango AC: mediana | máx. | Bloques espectrales / rango AC: mediana | máx. | Producto estándar / rango AC: mediana | Bloques espectrales / rango MC 99,73 % | Criterio 2 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| IEEE-33 | slack única | 4 | 1,000 | 1,010 | 1,01 | 1,006 | 1,01 | 1,007 | 1,01 | 1,01 | 1,22 | True |
| IEEE-33 | slack única | 8 | 1,000 | 1,011 | 1,01 | 1,007 | 1,01 | 1,008 | 1,01 | 1,02 | 1,54 | True |
| IEEE-33 | slack única | 16 | 1,000 | 1,012 | 1,01 | 1,010 | 1,01 | 1,008 | 1,01 | 1,02 | 1,96 | True |
| PEGASE-2869 | participación | 8 | 1,001 | 1,063 | 2,07 | 1,041 | 1,57 | 1,043 | 1,51 | 26,28 | 1,26 | True |
| PEGASE-2869 | participación | 16 | 1,001 | 1,081 | 2,13 | 1,062 | 1,76 | 1,052 | 1,39 | 27,66 | 1,39 | True |
| PEGASE-2869 | participación | 32 | 1,003 | 1,099 | 2,17 | 1,093 | 2,03 | 1,067 | 1,38 | 31,70 | 1,56 | True |
| PEGASE-2869 | participación | 64 | 1,003 | 1,106 | 2,19 | 1,103 | 2,16 | 1,070 | 1,40 | 34,27 | 1,72 | True |
| PEGASE-2869 | slack única | 8 | 1,002 | 1,065 | 1,93 | 1,045 | 1,51 | 1,045 | 1,47 | 176,13 | 1,18 | True |
| PEGASE-2869 | slack única | 16 | 1,002 | 1,074 | 2,02 | 1,055 | 1,55 | 1,051 | 1,47 | 174,40 | 1,27 | True |
| PEGASE-2869 | slack única | 32 | 1,003 | 1,085 | 2,08 | 1,076 | 1,95 | 1,060 | 1,53 | 173,72 | 1,43 | True |
| PEGASE-2869 | slack única | 64 | 1,003 | 1,094 | 2,11 | 1,088 | 2,03 | 1,060 | 1,61 | 172,87 | 1,57 | True |

**Se cumple.** El intervalo con bloques espectrales queda entre 1,01 y 1,07 veces el rango AC
de los vértices en mediana, y 1,61 veces en la peor barra. Con AF2 la peor barra llega a
2,19 veces. El producto afín estándar, término a término, da 26 a 34
veces en PEGASE con participación y 173 a 176 con slack única.

### 18.6 Criterio 3: ¿cuánto mejoran los bloques a AF2?

> **Aviso (15-11-2026).** Las cifras de RAA-K de esta sección se calcularon sin normalizar las formas
> cuadráticas. El operador usa una tolerancia absoluta y con formas de este tamaño pierde bloques, así que
> estas cifras **subestiman a RAA-K** (sección 25.8, con la medición). Las cotas siguen siendo válidas; solo son
> más anchas de lo que el operador puede dar. Falta rehacer la sección con la forma normalizada.

Ancho de la cota sobre el rango exacto de εᵀH_iε, mediana; entre paréntesis, el lado positivo espurio como
porcentaje de la profundidad. Solo barras con H_i semidefinida negativa y K ≤ 16, donde el rango exacto es
[mínimo en los vértices, 0]:

| Sistema | Balance | K | Barras con rango exacto | Producto estándar | AF2 | RAA-K, k = 2 | RAA-K, k = 3 | RAA-K, k = 6 | RAA-K, k = K | Espectral, 1 bloque | Espectral, 8 bloques | Espectral, todos |
|---|---|---:|---:|---|---|---|---|---|---|---|---|---|
| IEEE-33 | slack única | 4 | 32 | 2,29 (115 %) | 1,52 (52 %) | 1,15 (15 %) | 1,00 (0 %) | 1,00 (0 %) | 1,00 (0 %) | 1,41 (17 %) | 1,29 (0 %) | 1,29 (0 %) |
| IEEE-33 | slack única | 8 | 32 | 2,34 (118 %) | 1,70 (70 %) | 1,50 (50 %) | 1,31 (31 %) | 1,07 (7 %) | 1,07 (7 %) | 1,43 (20 %) | 1,36 (0 %) | 1,36 (0 %) |
| IEEE-33 | slack única | 16 | 32 | 2,37 (120 %) | 1,78 (78 %) | 1,68 (68 %) | 1,59 (59 %) | 1,48 (48 %) | 1,51 (51 %) | 1,47 (25 %) | 1,38 (0 %) | 1,38 (0 %) |
| PEGASE-2869 | participación | 8 | 984 | 375,47 (18683 %) | 1,73 (71 %) | 1,54 (52 %) | 1,39 (34 %) | 1,15 (14 %) | 1,10 (9 %) | 1,23 (16 %) | 1,17 (0 %) | 1,17 (0 %) |
| PEGASE-2869 | participación | 16 | 757 | 402,48 (19126 %) | 1,84 (82 %) | 1,73 (67 %) | 1,65 (59 %) | 1,46 (44 %) | 1,39 (37 %) | 1,34 (24 %) | 1,21 (0 %) | 1,21 (0 %) |
| PEGASE-2869 | slack única | 8 | 945 | 2012,84 (98518 %) | 1,69 (67 %) | 1,43 (43 %) | 1,19 (19 %) | 1,07 (7 %) | 1,05 (5 %) | 1,10 (5 %) | 1,08 (0 %) | 1,08 (0 %) |
| PEGASE-2869 | slack única | 16 | 419 | 503,03 (23889 %) | 1,89 (87 %) | 1,76 (74 %) | 1,65 (62 %) | 1,42 (40 %) | 1,22 (17 %) | 1,21 (10 %) | 1,14 (0 %) | 1,15 (0 %) |

Sobre todas las barras, ancho de AF2 dividido por el de cada alternativa (mediana barra a barra), y tiempos:

| Sistema | Balance | K | AF2 / RAA-K k = 2 | k = 3 | k = 6 | k = K | AF2 / espectral 1 bloque | 8 bloques | todos | Resto AF2 / radio, mediana (%) | 1.er orden (s) | 2.º orden (s) | Cotas (s) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| IEEE-33 | slack única | 4 | 1,32 | 1,50 | 1,50 | 1,50 | 1,09 | 1,22 | 1,22 | 2,0 | 0,004 | 0,00 | 0,07 |
| IEEE-33 | slack única | 8 | 1,14 | 1,32 | 1,59 | 1,59 | 1,21 | 1,27 | 1,27 | 2,3 | 0,003 | 0,00 | 0,02 |
| IEEE-33 | slack única | 16 | 1,08 | 1,14 | 1,20 | 1,20 | 1,22 | 1,32 | 1,32 | 2,4 | 0,003 | 0,00 | 0,09 |
| PEGASE-2869 | participación | 8 | 1,11 | 1,20 | 1,36 | 1,37 | 1,20 | 1,33 | 1,33 | 13,8 | 0,013 | 0,94 | 0,98 |
| PEGASE-2869 | participación | 16 | 1,07 | 1,11 | 1,15 | 1,14 | 1,18 | 1,37 | 1,37 | 17,3 | 0,014 | 1,01 | 4,33 |
| PEGASE-2869 | participación | 32 | 1,04 | 1,06 | 1,07 | 1,05 | 1,13 | 1,37 | 1,36 | 20,8 | 0,017 | 1,25 | 10,45 |
| PEGASE-2869 | participación | 64 | 1,02 | 1,03 | 1,03 | 1,03 | 1,12 | 1,33 | 1,33 | 22,8 | 0,022 | 2,18 | 16,86 |
| PEGASE-2869 | slack única | 8 | 1,16 | 1,24 | 1,35 | 1,35 | 1,14 | 1,22 | 1,22 | 14,4 | 0,020 | 1,56 | 1,99 |
| PEGASE-2869 | slack única | 16 | 1,10 | 1,18 | 1,30 | 1,33 | 1,14 | 1,26 | 1,26 | 16,2 | 0,023 | 1,58 | 6,02 |
| PEGASE-2869 | slack única | 32 | 1,06 | 1,10 | 1,14 | 1,09 | 1,10 | 1,27 | 1,27 | 18,2 | 0,034 | 1,98 | 14,71 |
| PEGASE-2869 | slack única | 64 | 1,03 | 1,05 | 1,05 | 1,04 | 1,09 | 1,28 | 1,27 | 20,1 | 0,035 | 3,10 | 25,58 |

**Se cumple, y el operador del grupo es el mejor cuando las fuentes son pocas.**

- **Con hasta 8 zonas gana RAA-K.** En IEEE-33 con K = 4 da el rango exacto (1,00, contra
  1,29 de los bloques espectrales y 1,52 de AF2). Con K = 8 queda a
  1,05–1,10 veces el rango exacto en los tres casos; los bloques espectrales, a
  1,08–1,36, y AF2, a 1,69–1,73.
- **Con más zonas se invierte.** Con K = 16, RAA-K queda a 1,22–1,51 y los bloques
  espectrales a 1,15–1,38. Sobre todas las barras de PEGASE, la ventaja de RAA-K
  frente a AF2 baja de 1,35–1,37 veces con K = 8 a 1,03–1,04
  con K = 64, mientras la de los bloques espectrales se mantiene en 1,22–1,37. Con muchas variables
  el presupuesto diagonal se agota antes de absorber lo que está fuera de la diagonal, y la búsqueda solo une
  variables cuya diagonal tiene el mismo signo.
- **La ganancia está en el lado superior, no en el que manda.** En formas semidefinidas negativas, el lado
  inferior de AF2 ya es casi exacto (1,00 veces la profundidad) y RAA-K lo conserva
  (1,00). Lo que RAA-K quita es el lado positivo espurio, que con AF2 vale
  52 a 87 % de la profundidad y con RAA-K 0 a 51 %. Los
  bloques espectrales lo eliminan, pero aflojan el lado inferior (1,08 a 1,38 veces).
  Para la cota inferior de voltaje, que es la que importa, AF2 y RAA-K dan lo mismo; RAA-K evita inflar la
  superior.
- **El tamaño del clique importa.** Con pares (el método binomial) la ganancia es la menor: 1,02 a
  1,32 veces. Pero dejar crecer el clique sin límite no siempre ayuda: en 2 de los
  11 casos, k = K queda peor que k = 6, por el orden voraz y el tope de nodos.
- **Costo.** Las cotas de todas las barras tardan hasta 25,6 s (K = 64), casi todo en la búsqueda de
  cliques.
- **Lectura para el grupo.** El teorema vale para bloques cualesquiera. Una combinación (los autovectores
  dominantes primero, RAA-K sobre el residuo) podría juntar lo mejor de las dos; no está probada.

### 18.7 Veredicto de la Puerta 2

| Criterio | Resultado | Evidencia |
|---|---|---|
| 1. La cota contiene la solución AC en todos los vértices y muestras | **No** | Falla por el término de tercer orden: hasta 5,9·10⁻⁶ p.u. en IEEE-33 y 0,0064 p.u. en PEGASE (18.4) |
| 2. Ancho ≤ 1,5 veces el rango de referencia, en mediana | **Sí** | 1,01 a 1,07 (18.5) |
| 3. El operador por bloques mejora a AF2 | **Sí; cuánto depende del número de fuentes** | RAA-K: 1,20–1,59 veces en IEEE-33; en PEGASE, de 1,35–1,37 con K = 8 a 1,03–1,04 con K = 64. Bloques espectrales: 1,22–1,37. La ganancia es en el lado superior (18.6) |

Qué cambia en el plan:

1. **El punto 3.4 deja de ser opcional.** Sin control del tercer orden no se puede afirmar contención, ni siquiera
   en IEEE-33. O se hace el punto fijo r = −J⁻¹q(Cε + r), o el paper declara una cota de segundo orden con el
   déficit medido y su ley en δ.
2. **Hay dos afirmaciones que sí se sostienen hoy.** La corrección de segundo orden baja la fracción de muestras
   con alguna barra fuera de 24–51 % a 0,01 % como máximo, y cuesta
   segundos en 2.869 barras.
3. **El rango de δ importa.** Con zonas coherentes, δ = 5 % ya es un caso duro para un método de segundo orden en
   PEGASE; con δ = 1 % el déficit es de 2,2·10⁻⁵ p.u.
4. **El resto de un solo signo no es general.** Vale en IEEE-33; en PEGASE con zonas, no. La salida "solo se
   corrige la cota inferior" (punto 3.2) hay que plantearla con H_i indefinida.
5. **La elección de la slack pesa menos de lo que parecía** en estas medidas: el modelo de participación baja los
   ángulos y la sobreestimación del producto estándar, pero no el déficit en las barras débiles.

---

## 19. X3 — ¿El modelo de segundo orden encuentra el peor escenario de cada barra?

Para cada barra se proponen los vértices de la caja de ruido donde su voltaje sería extremo: con el modelo de
primer orden (los signos de c_i) y con el de segundo orden (el vértice que minimiza o maximiza c_iᵀx + xᵀH_ix).
Se resuelve el flujo AC exacto solo en esos vértices y se compara con el extremo verdadero sobre los 2^K
vértices, que X2 calculó por enumeración. No hay números aleatorios.

| Sistema | Balance | K | Barras | Vértices de la enumeración | Voltaje mínimo, 1.er orden: barras que aciertan (%) | brecha máx. (p.u.) | Voltaje mínimo, 2.º orden: aciertan (%) | brecha máx. (p.u.) | flujos AC | Voltaje máximo, 1.er orden: aciertan (%) | brecha máx. (p.u.) | Voltaje máximo, 2.º orden: aciertan (%) | brecha máx. (p.u.) | flujos AC |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| IEEE-33 | slack única | 4 | 32 | 16 | 100,0 | 0,00000 | 100,0 | 0,00000 | 1 | 100,0 | 0,00000 | 100,0 | 0,00000 | 1 |
| IEEE-33 | slack única | 8 | 32 | 256 | 100,0 | 0,00000 | 100,0 | 0,00000 | 1 | 100,0 | 0,00000 | 100,0 | 0,00000 | 1 |
| IEEE-33 | slack única | 16 | 32 | 65536 | 100,0 | 0,00000 | 100,0 | 0,00000 | 1 | 100,0 | 0,00000 | 100,0 | 0,00000 | 1 |
| PEGASE-2869 | slack única | 8 | 2312 | 256 | 85,2 | 0,00559 | 97,6 | 0,00123 | 63 | 71,2 | 0,00451 | 81,2 | 0,00027 | 116 |
| PEGASE-2869 | slack única | 16 | 2312 | 65536 | 77,5 | 0,00541 | 94,6 | 0,00500 | 122 | 60,1 | 0,00453 | 75,0 | 0,00029 | 372 |
| PEGASE-2869 | participación | 8 | 2322 | 256 | 83,4 | 0,00268 | 97,2 | 0,00151 | 65 | 66,4 | 0,00293 | 92,5 | 0,00045 | 115 |
| PEGASE-2869 | participación | 16 | 2322 | 65536 | 77,6 | 0,00183 | 96,9 | 0,00029 | 143 | 61,0 | 0,00293 | 86,9 | 0,00048 | 367 |

*Script: `experiments/x3_worst_vertex.py` · CSV: `results/x3_worst_vertex.csv` · "acierta": el voltaje alcanzado
coincide con el extremo de la enumeración dentro de 10⁻⁹ p.u. · "brecha": cuánto le falta al peor candidato*

- En IEEE-33 el peor vértice es el mismo para todas las barras (todas las zonas al mismo lado) y los dos modelos
  lo encuentran.
- **En PEGASE el segundo orden acierta el voltaje mínimo en 95 a 98 % de las barras; el
  primer orden, en 77 a 85 %.** Para el máximo: 75 a 92 % contra
  60 a 71 %.
- Bastan 63 a 143 flujos AC para todas las barras, porque muchas comparten vértice;
  la enumeración usa 256 o 65.536.
- **No es exacto.** Donde falla, al mínimo le faltan hasta 0,0003–0,0050 p.u. según el caso
  (con el primer orden, 0,0018–0,0056). Lo que entrega es un valor alcanzado, es decir, una
  cota interior del rango: sirve para medir cuán conservadora es una cota exterior, no para reemplazarla.

---

## 20. Auditoría del modelo y controles de la referencia (aportes 2 y 5)

**Auditoría previa.** Antes de probar nada se revisó el modelo de segundo orden contra su derivación y se le
pusieron pruebas que no dependen de ella (`tests/test_second_order.py`, en una red radial sin barras PV y en una
mallada con barras PV, con y sin participación de generadores):

- H_i coincide con la segunda derivada de la solución AC exacta, calculada por diferencias finitas: el error
  relativo medido es cercano a 10⁻⁶ y la prueba exige menos de 2·10⁻⁵.
- El error que queda tras sumar la forma cuadrática se divide por ocho cuando δ se divide por dos.
- La forma ensamblada, la expresión término a término y la traza coinciden; las sensibilidades rectangulares y
  polares también.

### 20.1 ¿Son reales las violaciones? Solver independiente

En el vértice donde la cota de segundo orden queda más corta, el mismo escenario se arma en pandapower (escalando
sus cargas y generadores estáticos zona por zona, y moviendo las consignas de los generadores en el modelo de
participación) y se resuelve con su Newton-Raphson.

| Cantidad | IEEE-33, slack única, K = 8 | PEGASE-2869, slack única, K = 8 | PEGASE-2869, participación, K = 8 |
|---|---|---|---|
| Barra | 17 | 7639 | 8853 |
| Desajuste máximo del solver por lotes (p.u.) | 4,8·10⁻¹⁰ | 8,2·10⁻¹⁰ | 5,6·10⁻¹⁰ |
| Mayor diferencia de \|V\| con pandapower, todas las barras (p.u.) | 1,7·10⁻¹⁰ | 1,1·10⁻¹¹ | 8,9·10⁻¹² |
| Radio afín de la barra (p.u.) | 0,0142 | 0,0177 | 0,0314 |
| ΔV, solver por lotes (p.u.) | -0,014365 | -0,039642 | -0,053003 |
| ΔV, pandapower (p.u.) | -0,014365 | -0,039642 | -0,053003 |
| Caja afín: violada por (p.u.) | 2,1·10⁻⁴ | 0,0219 | 0,0216 |
| 2.º orden con AF2: violada por (p.u.) | 5,9·10⁻⁶ | 0,0037 | 0,0055 |
| 2.º orden con RAA-K: violada por (p.u.) | 5,9·10⁻⁶ | 0,0037 | 0,0055 |
| 2.º orden con bloques espectrales: violada por (p.u.; negativo = contiene) | -6,3·10⁻⁵ | 0,0025 | 0,0053 |

*Script: `experiments/x4_reference_checks.py` · CSV: `results/x4_independent.csv`*

- Los dos solvers dan el mismo voltaje en todas las barras: la mayor diferencia es 1,7·10⁻¹⁰ p.u. Las
  violaciones de las secciones 17 a 19 no son un artefacto del solver por lotes ni del armado del escenario.
- Límite del modelo, que también vale para pandapower tal como se usó: los generadores no tienen límites de
  potencia reactiva. Liao et al. sí conmutan el tipo de barra al violarlos.

### 20.2 ¿Está el extremo en un vértice?

X2 tomó como referencia el rango de la solución AC en los 2^K vértices. Aquí se minimiza y maximiza el modelo de
segundo orden sobre toda la caja (puntos interiores incluidos), se resuelve el flujo AC exacto en ese punto, y en
las diez barras donde más gana se afina con una búsqueda acotada sobre el flujo AC exacto.

| Sistema | Balance | K | Voltaje | Barras | Óptimo del modelo en el interior | Barras que superan a todos los vértices | por hasta (p.u.) | > 10⁻⁵ p.u. | > 10⁻⁴ p.u. | Tras afinar con el flujo AC exacto (p.u.) |
|---|---|---:|---|---:|---:|---:|---|---:|---:|---|
| IEEE-33 | slack única | 4 | mínimo | 32 | 0 | 0 | 3,4·10⁻¹² | 0 | 0 | 3,4·10⁻¹² |
| IEEE-33 | slack única | 4 | máximo | 32 | 0 | 0 | 4,1·10⁻¹² | 0 | 0 | 4,1·10⁻¹² |
| IEEE-33 | slack única | 8 | mínimo | 32 | 0 | 0 | 3,4·10⁻¹² | 0 | 0 | 3,4·10⁻¹² |
| IEEE-33 | slack única | 8 | máximo | 32 | 0 | 0 | 4,1·10⁻¹² | 0 | 0 | 4,1·10⁻¹² |
| IEEE-33 | slack única | 16 | mínimo | 32 | 0 | 0 | 3,4·10⁻¹² | 0 | 0 | 3,4·10⁻¹² |
| IEEE-33 | slack única | 16 | máximo | 32 | 0 | 0 | 4,1·10⁻¹² | 0 | 0 | 4,1·10⁻¹² |
| PEGASE-2869 | slack única | 8 | mínimo | 2312 | 44 | 1 | 1,0·10⁻⁸ | 0 | 0 | 2,3·10⁻⁸ |
| PEGASE-2869 | slack única | 8 | máximo | 2312 | 931 | 469 | 1,2·10⁻⁴ | 111 | 3 | 1,2·10⁻⁴ |
| PEGASE-2869 | slack única | 16 | mínimo | 2312 | 45 | 0 | 4,9·10⁻¹² | 0 | 0 | 4,9·10⁻¹² |
| PEGASE-2869 | slack única | 16 | máximo | 2312 | 990 | 414 | 1,3·10⁻⁴ | 35 | 6 | 1,3·10⁻⁴ |
| PEGASE-2869 | participación | 8 | mínimo | 2322 | 25 | 0 | 5,0·10⁻¹² | 0 | 0 | 5,0·10⁻¹² |
| PEGASE-2869 | participación | 8 | máximo | 2322 | 780 | 543 | 1,7·10⁻⁴ | 59 | 1 | 1,7·10⁻⁴ |
| PEGASE-2869 | participación | 16 | mínimo | 2322 | 30 | 0 | 4,9·10⁻¹² | 0 | 0 | 4,9·10⁻¹² |
| PEGASE-2869 | participación | 16 | máximo | 2322 | 820 | 494 | 1,3·10⁻⁵ | 3 | 0 | 1,3·10⁻⁵ |

*CSV: `results/x4_interior.csv`*

- **El voltaje mínimo está en un vértice.** Ningún punto interior baja de la enumeración en más de
  1,0·10⁻⁸ p.u. La referencia de la cota inferior, que es la que manda, es correcta.
- **El máximo no siempre.** En PEGASE, entre 18 y 23 % de las barras tienen un punto
  interior por encima de todos los vértices, por hasta 1,7·10⁻⁴ p.u.; a lo más 6 barras
  pasan de 10⁻⁴ p.u. Los anchos "sobre el rango de referencia" de las secciones 18 y 19 están apenas
  sobreestimados por el lado superior.
- La búsqueda sobre el flujo exacto agrega como máximo 2,4·10⁻⁶ p.u. a lo que ya encuentra el modelo
  de segundo orden: para ubicar el extremo, el modelo basta.

---

## 21. X5 — El método publicado más cercano: Liao et al. (2017) (aporte 4)

Se reimplementó la Ec. (26) del artículo: desarrollo de segundo orden de cada componente del estado (e, f),
acotado monomio a monomio. El artículo no dice cómo pasa de los intervalos de e y f al de |V|; se usó aritmética
de intervalos sobre √([e]² + [f]²). Se compara con la forma directa de este paquete, que proyecta la misma
información de segundo orden sobre |V| antes de acotar.

| Caso | K | Radio afín del ángulo (rad) | Liao sobre e y f: ancho / rango (componentes que no contienen) | \|V\| desde [e], [f] con intervalos: ancho / rango | Forma directa, AF2 | Forma directa, RAA-K | Resoluciones lineales, vía de Liao | vía adjunta |
|---|---:|---:|---|---:|---:|---:|---:|---:|
| IEEE 9, ±10 % | 8 | 0,073 | 1,13 (7 de 16) | 1,83 | 1,118 | 1,072 | 44 | 14 |
| IEEE 57, ±20 % | 79 | 0,227 | 1,80 (0 de 112) | 5,92 | 1,218 | 1,194 | 3239 | 129 |
| PEGASE-2869, 8 zonas, slack única | 8 | 0,810 | – | 324,97 | 1,065 | 1,045 | 44 | 2320 |
| PEGASE-2869, 8 zonas, participación | 8 | 0,346 | – | 87,78 | 1,063 | 1,041 | 44 | 2330 |

*Script: `experiments/x5_liao_baseline.py` · CSV: `results/x5_liao_summary.csv` · "ancho / rango": mediana sobre
barras del ancho del intervalo sobre el rango de la solución AC (todos los vértices en IEEE 9 y PEGASE; vértices
extremos y Monte Carlo en IEEE 57)*

- **Las dos derivaciones coinciden.** Los coeficientes de segundo orden de Liao (una resolución por cada par de
  variables), proyectados sobre |V|, son iguales a los de la vía adjunta (una resolución traspuesta por barra).
  El script lo comprueba en cada caso.
- **La Ec. (26) tampoco contiene la solución exacta.** En IEEE 9 con ±10 %, 7 de las
  16 componentes de e y f tienen algún vértice fuera de su intervalo, por hasta 6,7·10⁻⁵
  p.u. Es el tercer orden, en el sistema de prueba del propio artículo.
- **Pasar de e y f a |V| con intervalos no escala.** El intervalo queda a 1,8 veces el rango en 9
  barras, 5,9 en 57 y 88 a 325 en PEGASE, porque se pierde que e y f
  rotan juntas. La forma directa se mantiene entre 1,04 y 1,19.
- **La vía adjunta no es más barata con pocas fuentes.** Con 8 zonas en PEGASE, la vía de Liao necesita
  44 resoluciones lineales y la adjunta 2.320. La adjunta solo gana cuando hay más
  fuentes que barras a vigilar. Lo que marca la diferencia es proyectar sobre |V| antes de acotar, y eso se puede
  hacer con cualquiera de las dos.

**¿Se reproduce la Tabla 1 del artículo?** No.

| Barra | Nominal | Liao, Tabla 1 del artículo | Liao reimplementado | Forma directa, AF2 | Forma directa, RAA-K | Exacto en los 256 vértices | Monte Carlo del artículo | Monte Carlo aquí (extremos de 20.000) |
|---:|---:|---|---|---|---|---|---|---|
| 4 | 0,9870 | [0,9717; 0,9923] | [0,9780; 0,9966] | [0,9790; 0,9946] | [0,9790; 0,9937] | [0,9799; 0,9933] | [0,9793; 0,9910] | [0,9811; 0,9919] |
| 5 | 0,9755 | [0,9672; 0,9839] | [0,9610; 0,9923] | [0,9640; 0,9865] | [0,9640; 0,9855] | [0,9649; 0,9849] | [0,9711; 0,9805] | [0,9666; 0,9825] |
| 6 | 1,0034 | [1,0001; 1,0066] | [0,9938; 1,0180] | [0,9996; 1,0070] | [0,9996; 1,0066] | [1,0000; 1,0065] | [1,0010; 1,0062] | [1,0006; 1,0058] |
| 7 | 0,9856 | [0,9793; 0,9926] | [0,9743; 1,0010] | [0,9786; 0,9924] | [0,9786; 0,9920] | [0,9791; 0,9918] | [0,9823; 0,9896] | [0,9802; 0,9908] |
| 8 | 0,9962 | [0,9917; 1,0088] | [0,9799; 1,0191] | [0,9910; 1,0011] | [0,9910; 1,0007] | [0,9915; 1,0005] | [0,9931; 1,0058] | [0,9923; 0,9997] |
| 9 | 0,9576 | [0,9470; 0,9682] | [0,9398; 0,9779] | [0,9433; 0,9715] | [0,9433; 0,9704] | [0,9439; 0,9696] | [0,9485; 0,9649] | [0,9471; 0,9672] |

| Variables inciertas | K | Rango exacto aquí frente al Monte Carlo del artículo: mayor diferencia (p.u.) | Liao reimplementado frente a su Tabla 1 (p.u.) | Barras (de 6) donde el intervalo del artículo contiene el rango exacto |
|---|---:|---:|---:|---:|
| P de generadores PV; P y Q de cargas | 8 | 0,0062 | 0,0118 | 0 |
| P y Q de cargas | 6 | 0,0065 | 0,0094 | 6 |
| P de generadores PV; P de cargas | 5 | 0,0079 | 0,0095 | 6 |
| P de cargas | 3 | 0,0091 | 0,0126 | 6 |
| P de generadores PV | 2 | 0,0086 | 0,0128 | 6 |

*CSV: `results/x5_liao_table1.csv`, `results/x5_liao_variants.csv`*

- Los voltajes nominales coinciden con cuatro decimales: es el mismo sistema.
- El artículo habla de 12 variables de perturbación sin listarlas. Con cinco lecturas distintas del conjunto
  incierto, el rango exacto queda a no menos de 0,006 p.u. de su columna de Monte Carlo, que no depende
  del método. El escenario del artículo no se puede identificar a partir de su texto.
- Sus intervalos son bastante más angostos que los de mi lectura de "|V| con intervalos". O usan otro paso, que
  no describen, u otro conjunto incierto. **No se puede afirmar que el método de Liao se desborde a escala**:
  lo que se desborda es la recombinación con intervalos.

---

## 22. X6 — Punto fijo sobre el resto (aporte 1)

**Construcción** (`cdhim/enclosure.py`). Con D = Cε + r, el resto r es punto fijo de
T_ε(r) = −J⁻¹[q(Cε) + ΔJ(Cε)·r + q(r)]. Para |r| ≤ ρ y todo ε de la caja, |T_ε(r)| ≤ b₂ + Gρ + N(ρ), con b₂ la
cota de las formas cuadráticas del estado, G = Σ_k |J⁻¹ΔJ(C_k)| y N(ρ) la parte cuadrática, las dos con la
inversa aplicada antes de tomar valores absolutos. Si existe ρ* ≥ b₂ + Gρ* + N(ρ*), el teorema de Brouwer da una
solución con |D − Cε| ≤ ρ* para todo ε, y lo que queda por encima del segundo orden está acotado por
σ = Gρ* + N(ρ*). El intervalo de |V_i| es la forma de segundo orden más una holgura que junta σ y los términos
de orden superior del módulo. Vale en aritmética exacta; no hay redondeo dirigido.

**No es una construcción nueva.** Es, en lo esencial, la condición de auto-mapeo del preprint de Z. Zhu,
*Convex Inner Approximations for Safe AC Dispatch With Device Dynamics* (arXiv:2610.26634, octubre de 2026), que
a su vez sigue a Dvijotham, Nguyen y Turitsyn (2018) y a Lee et al. (2019, 2021). Esos trabajos la usan para
certificar factibilidad, no entregan intervalos de |V| por barra y no comparan el encierro con el rango real.
Lo que aquí se agrega es restar el segundo orden exacto antes de aplicar la holgura y medir el resultado contra
todos los vértices.

Ocho zonas en todos los casos:

| Sistema | Balance | δ (%) | Radio de G, rectangulares | Radio de G, polares (1.er orden) | Variación real del Jacobiano, rectangulares | ídem, polares | Punto fijo | Holgura añadida / radio, mediana (%) | Encierro / rango exacto, mediana | Barras con algún vértice fuera | ídem, solo 2.º orden |
|---|---|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|
| IEEE-33 | slack única | 5 | 0,023 | 0,024 | 0,009 | 0,011 | sí | 0,013 | 1,002 | 0 | 27 |
| IEEE-33 | slack única | 15 | 0,069 | 0,073 | 0,028 | 0,031 | sí | 0,250 | 1,009 | 0 | 32 |
| IEEE-33 | slack única | 30 | 0,139 | 0,145 | 0,057 | 0,062 | no | – | – | – | – |
| IEEE-118 | slack única | 2 | 0,086 | 0,041 | 0,044 | 0,011 | sí | 29,254 | 1,306 | 0 | 5 |
| IEEE-118 | slack única | 5 | 0,215 | 0,102 | 0,110 | 0,027 | no | – | – | – | – |
| IEEE-118 | slack única | 10 | 0,430 | 0,205 | 0,220 | 0,054 | no | – | – | – | – |
| IEEE-118 | participación | 2 | 0,074 | 0,046 | 0,027 | 0,011 | sí | 6,857 | 1,079 | 0 | 1 |
| IEEE-118 | participación | 5 | 0,186 | 0,115 | 0,068 | 0,026 | no | – | – | – | – |
| IEEE-118 | participación | 10 | 0,372 | 0,231 | 0,136 | 0,053 | no | – | – | – | – |
| IEEE-300 | slack única | 2 | 1,417 | 0,983 | 0,499 | 0,312 | no | – | – | – | – |
| IEEE-300 | slack única | 5 | 3,543 | 2,457 | 1,247 | 0,733 | no | – | – | – | – |
| IEEE-300 | slack única | 10 | 7,085 | 4,914 | 2,493 | 1,481 | no | – | – | – | – |
| IEEE-300 | participación | 2 | 0,863 | 0,731 | 0,150 | 0,096 | no | – | – | – | – |
| IEEE-300 | participación | 5 | 2,156 | 1,827 | 0,374 | 0,247 | no | – | – | – | – |
| IEEE-300 | participación | 10 | 4,313 | 3,654 | 0,748 | 0,517 | no | – | – | – | – |
| PEGASE-2869 | slack única | 1 | 0,625 | 0,330 | 0,243 | 0,054 | no | – | – | – | – |
| PEGASE-2869 | slack única | 2 | 1,251 | 0,661 | 0,486 | 0,109 | no | – | – | – | – |
| PEGASE-2869 | slack única | 5 | 3,127 | 1,652 | 1,214 | 0,280 | no | – | – | – | – |
| PEGASE-2869 | participación | 1 | 0,613 | 0,385 | 0,154 | 0,054 | no | – | – | – | – |
| PEGASE-2869 | participación | 2 | 1,225 | 0,770 | 0,307 | 0,110 | no | – | – | – | – |
| PEGASE-2869 | participación | 5 | 3,064 | 1,925 | 0,768 | 0,283 | no | – | – | – | – |

*Script: `experiments/x6_fixed_point.py` · CSV: `results/x6_fixed_point.csv` · "variación real del Jacobiano":
mayor radio espectral de J⁻¹(J(x* + Cε) − J(x*)) en 64 vértices, sin ninguna cota*

- **Funciona en el alimentador radial.** En IEEE-33 con δ de 5 y 15 % hay punto fijo; la holgura añadida es
  0,01 a 0,25 % del radio y el encierro queda a 1,002–1,009 veces el
  rango exacto.
- **Donde hay punto fijo, el encierro contiene.** En los 4 casos con punto fijo (de 21),
  el número de violaciones en los 256 vértices, sumando el resto del estado, la parte de tercer orden y |V|, es
  0. La cota de segundo orden sola fallaba en esos mismos casos.
- **En 118 barras solo cierra con δ = 2 %**, y ya suelta: la holgura es 7 a
  29 % del radio y el encierro queda a 1,08–1,31 veces el rango. Es el mismo
  orden de caja que usa Zhu (1,5 a 2 % en 39 y 118 barras).
- **En 300 y 2.869 barras no cierra.** El radio de G en PEGASE con δ = 5 % es 3,1; en coordenadas
  polares baja a 1,7 a 1,9, que sigue sobre 1.
- **El obstáculo es la cota componente a componente, no la física.** La variación real del Jacobiano en PEGASE
  con δ = 5 %, en polares, tiene radio espectral 0,28: hay margen para un argumento de contracción,
  pero los valores absolutos elemento a elemento lo multiplican por 6 o más.
- IEEE-33 con δ = 30 % tampoco cierra, aunque el radio de G es 0,14: ahí manda la parte cuadrática N(ρ), que
  trata como independientes restos de barras vecinas que en realidad se mueven juntos.

---

## 23. X7 — RAA-K frente a la relajación semidefinida (aporte 3)

> **Aviso (15-11-2026).** Las cifras de RAA-K de esta sección se calcularon sin normalizar las formas
> cuadráticas. El operador usa una tolerancia absoluta y con formas de este tamaño pierde bloques, así que
> estas cifras **subestiman a RAA-K** (sección 25.8, con la medición). Las cotas siguen siendo válidas; solo son
> más anchas de lo que el operador puede dar. Falta rehacer la sección con la forma normalizada.

Según la literatura de formas cuadráticas sobre una caja, hay una escalera de cotas: término a término (Messine
y Touhami), relajación semidefinida de Shor, y el rango exacto, que es NP-duro en general. Se agregó la
relajación semidefinida y una cota híbrida (los autovectores dominantes como bloques y RAA-K sobre lo que
queda). Muestra de hasta 300 barras por caso.

| Sistema | Balance | K | Cota | Formas | Ancho / rango en vértices: mediana | máx. | Lado inferior / profundidad (formas semidefinidas negativas) | Lado superior / profundidad | ms por forma |
|---|---|---:|---|---:|---:|---:|---:|---:|---:|
| IEEE-33 | slack única | 8 | AF2 | 32 | 1,696 | 1,74 | 1,000 | 0,696 | 0,07 |
| IEEE-33 | slack única | 8 | RAA-K | 32 | 1,069 | 1,68 | 1,000 | 0,069 | 2,18 |
| IEEE-33 | slack única | 8 | Bloques espectrales | 32 | 1,362 | 1,63 | 1,362 | 0,000 | 0,12 |
| IEEE-33 | slack única | 8 | Híbrida, 1 autovector + RAA-K | 32 | 1,328 | 1,74 | 1,243 | 0,093 | 0,24 |
| IEEE-33 | slack única | 8 | Híbrida, 2 autovectores | 32 | 1,336 | 1,70 | 1,295 | 0,038 | 0,20 |
| IEEE-33 | slack única | 8 | Híbrida, 4 autovectores | 32 | 1,363 | 1,62 | 1,345 | 0,019 | 0,18 |
| IEEE-33 | slack única | 8 | Relajación semidefinida | 32 | 1,000 | 1,00 | 1,000 | -0,000 | 13,88 |
| IEEE-33 | slack única | 16 | AF2 | 32 | 1,782 | 1,83 | 1,000 | 0,782 | 0,07 |
| IEEE-33 | slack única | 16 | RAA-K | 32 | 1,510 | 1,77 | 1,000 | 0,510 | 3,23 |
| IEEE-33 | slack única | 16 | Bloques espectrales | 32 | 1,376 | 1,77 | 1,376 | 0,000 | 0,12 |
| IEEE-33 | slack única | 16 | Híbrida, 1 autovector + RAA-K | 32 | 1,428 | 2,10 | 1,202 | 0,215 | 1,54 |
| IEEE-33 | slack única | 16 | Híbrida, 2 autovectores | 32 | 1,367 | 1,86 | 1,294 | 0,060 | 0,21 |
| IEEE-33 | slack única | 16 | Híbrida, 4 autovectores | 32 | 1,382 | 1,77 | 1,352 | 0,030 | 0,20 |
| IEEE-33 | slack única | 16 | Relajación semidefinida | 32 | 1,000 | 1,00 | 1,000 | -0,000 | 118,21 |
| PEGASE-2869 | slack única | 8 | AF2 | 300 | 1,572 | 2,13 | 1,000 | 0,596 | 0,06 |
| PEGASE-2869 | slack única | 8 | RAA-K | 300 | 1,072 | 1,91 | 1,000 | 0,035 | 0,17 |
| PEGASE-2869 | slack única | 8 | Bloques espectrales | 300 | 1,097 | 1,90 | 1,090 | 0,000 | 0,09 |
| PEGASE-2869 | slack única | 8 | Híbrida, 1 autovector + RAA-K | 300 | 1,105 | 2,00 | 1,059 | 0,017 | 0,18 |
| PEGASE-2869 | slack única | 8 | Híbrida, 2 autovectores | 300 | 1,101 | 2,06 | 1,086 | 0,009 | 0,16 |
| PEGASE-2869 | slack única | 8 | Híbrida, 4 autovectores | 300 | 1,097 | 1,90 | 1,090 | 0,001 | 0,15 |
| PEGASE-2869 | slack única | 8 | Relajación semidefinida | 300 | 1,000 | 1,07 | 1,000 | -0,000 | 11,67 |
| PEGASE-2869 | slack única | 16 | AF2 | 300 | 1,721 | 2,14 | 1,000 | 0,881 | 0,07 |
| PEGASE-2869 | slack única | 16 | RAA-K | 300 | 1,207 | 2,05 | 1,000 | 0,197 | 1,81 |
| PEGASE-2869 | slack única | 16 | Bloques espectrales | 300 | 1,168 | 2,08 | 1,159 | 0,000 | 0,09 |
| PEGASE-2869 | slack única | 16 | Híbrida, 1 autovector + RAA-K | 300 | 1,187 | 2,25 | 1,104 | 0,054 | 0,89 |
| PEGASE-2869 | slack única | 16 | Híbrida, 2 autovectores | 300 | 1,178 | 2,12 | 1,134 | 0,029 | 0,53 |
| PEGASE-2869 | slack única | 16 | Híbrida, 4 autovectores | 300 | 1,165 | 2,08 | 1,151 | 0,006 | 0,17 |
| PEGASE-2869 | slack única | 16 | Relajación semidefinida | 300 | 1,000 | 1,17 | 1,000 | -0,000 | 126,13 |
| PEGASE-2869 | participación | 8 | AF2 | 300 | 1,727 | 2,01 | 1,000 | 0,714 | 0,05 |
| PEGASE-2869 | participación | 8 | RAA-K | 300 | 1,170 | 1,98 | 1,000 | 0,093 | 0,18 |
| PEGASE-2869 | participación | 8 | Bloques espectrales | 300 | 1,158 | 1,94 | 1,162 | 0,000 | 0,09 |
| PEGASE-2869 | participación | 8 | Híbrida, 1 autovector + RAA-K | 300 | 1,178 | 2,17 | 1,080 | 0,051 | 0,19 |
| PEGASE-2869 | participación | 8 | Híbrida, 2 autovectores | 300 | 1,171 | 2,06 | 1,124 | 0,026 | 0,16 |
| PEGASE-2869 | participación | 8 | Híbrida, 4 autovectores | 300 | 1,159 | 1,95 | 1,162 | 0,002 | 0,15 |
| PEGASE-2869 | participación | 8 | Relajación semidefinida | 300 | 1,000 | 1,16 | 1,000 | -0,000 | 11,37 |
| PEGASE-2869 | participación | 16 | AF2 | 300 | 1,850 | 2,29 | 1,003 | 0,824 | 0,06 |
| PEGASE-2869 | participación | 16 | RAA-K | 300 | 1,450 | 2,02 | 1,003 | 0,374 | 2,97 |
| PEGASE-2869 | participación | 16 | Bloques espectrales | 300 | 1,223 | 2,18 | 1,213 | 0,000 | 0,10 |
| PEGASE-2869 | participación | 16 | Híbrida, 1 autovector + RAA-K | 300 | 1,334 | 2,17 | 1,098 | 0,162 | 1,66 |
| PEGASE-2869 | participación | 16 | Híbrida, 2 autovectores | 300 | 1,271 | 2,01 | 1,123 | 0,084 | 0,59 |
| PEGASE-2869 | participación | 16 | Híbrida, 4 autovectores | 300 | 1,251 | 2,06 | 1,210 | 0,025 | 0,20 |
| PEGASE-2869 | participación | 16 | Relajación semidefinida | 300 | 1,001 | 1,10 | 1,000 | -0,000 | 137,90 |

*Script: `experiments/x7_block_choices.py` · CSV: `results/x7_block_choices.csv` · necesita `cvxpy`*

- **La relajación semidefinida da prácticamente el rango exacto**: 1,000 a 1,001 veces en mediana
  (1,17 en la peor forma). RAA-K queda entre 1,07 y 1,51, los bloques espectrales
  entre 1,10 y 1,38 y AF2 entre 1,57 y 1,85.
- **RAA-K es mucho más rápida**: 0,2 a 3,2 ms por forma contra
  11 a 138 ms. Con 16 símbolos y 2.300 barras, la relajación semidefinida son unos
  cinco minutos; con miles de símbolos no es viable.
- **RAA-K se puede hacer rigurosa con redondeo dirigido** (el grupo ya lo tiene en C++). Una cota semidefinida
  sale de un solver numérico y necesita un certificado dual verificado para ser rigurosa.
- **La cota híbrida no ayuda:** supera a la mejor de sus dos partes solo en 1 de los
  6 casos, y por menos de 1 %.

---

## 24. Veredicto de los cinco aportes posibles

> **Aviso (15-11-2026).** El veredicto del aporte 3 (RAA-K) se apoya en cifras que subestiman al operador
> (sección 25.8) y el del aporte 1 (punto fijo) quedó superado para redes sin barras PV por la cota de la
> sección 25.6. La tabla se deja como estaba hasta rehacer las secciones 18 a 23.

| # | Aporte | Lo más cercano en la literatura | Lo que mostró la prueba | Veredicto |
|---|---|---|---|---|
| 1 | Encierro verificado por punto fijo | Zhu (2026), Lee et al. (2019, 2021), Dvijotham et al. (2018): la misma construcción, para certificar factibilidad; hasta 118 barras con cajas de 1,5 a 2 % | Cierra en IEEE-33 hasta δ = 15 % y queda a 1,002–1,009 veces el rango exacto; en 118 barras solo con δ = 2 %; en 300 y 2.869 no cierra | **El mecanismo no es nuevo y no escala.** Aporte posible, acotado a alimentadores: intervalos de \|V\| por barra, ajustados y contrastados contra todos los vértices |
| 2 | Ni el primer ni el segundo orden encierran la solución | Vaccaro (2010) lo dice del primer orden sin medirlo. Liao (2017) afirma que no hay órdenes superiores; no encontré quien lo corrija ni quien valide en vértices | Confirmado con un solver independiente; la Ec. (26) de Liao falla en su propio IEEE 9 | **Se sostiene.** Es el resultado más sólido, pero es negativo |
| 3 | RAA-K para acotar el término cuadrático | Messine y Touhami (término a término); relajación semidefinida; ningún uso en flujo de potencia | Mejor que AF2 y, con pocas fuentes, que los bloques espectrales. La relajación semidefinida es más ajustada en todos los casos | **Parcial.** No es la cota más ajustada; su valor es el costo y que admite redondeo dirigido |
| 4 | Formulación adjunta y forma ensamblada | Liao (2017) obtiene los mismos coeficientes con una resolución por par; adjunto de segundo orden conocido (Pacaud et al., 2022) | Coinciden numéricamente. La adjunta no ahorra con pocas fuentes. Proyectar sobre \|V\| antes de acotar sí cambia el resultado | **Parcial.** El aporte es la proyección antes de acotar, no el adjunto. No se pudo comparar con el procedimiento real de Liao |
| 5 | Peor escenario guiado por el segundo orden | Escenarios optimizados (Zhang et al., 2018); flujo de frontera (Dimitrovski y Tomsovic, 2004); en redes radiales, un solo flujo da todos los mínimos (Sou et al., 2025) | Acierta el mínimo en 95 a 98 % de las barras de PEGASE; el mínimo está en un vértice; el máximo no siempre | **Parcial.** Útil como cota interior en redes malladas con barras PV; en radiales ya está resuelto |

Lo que esto deja:

1. **No hay un método nuevo que encierre la solución AC a escala de transmisión.** El punto fijo, que era la
   apuesta fuerte, ya está publicado en lo esencial y deja de cerrar entre 118 y 300 barras.
2. **Sí hay un resultado empírico sólido y no publicado:** el primer y el segundo orden no encierran, medido con
   todos los vértices y un solver independiente, con leyes δ² y δ³, hasta 2.869 barras.
3. **Sí hay una combinación que funciona en alimentadores de distribución:** segundo orden exacto, RAA-K y punto
   fijo dan intervalos de |V| que contienen la solución y quedan a menos de 1 % del rango exacto. Falta el
   redondeo dirigido para llamarlo verificado.
4. **A escala de transmisión, lo defendible hoy es una cota de segundo orden con el déficit medido**, no una
   garantía. El problema abierto está identificado: un argumento de contracción que no tome valores absolutos
   elemento a elemento.

---

## 25. Opción A, paso 1: encierros de voltaje en redes de distribución frente a lo publicado

La sección 24 dejó una combinación que funcionaba en un alimentador de 33 barras. Esta sección contesta si eso
alcanza para un aporte, midiéndola contra los dos métodos publicados que dan cotas garantizadas en redes de
distribución: el de Sou, Khong y Sandberg (2026) para redes radiales y el de Wang, Bernstein, Le Boudec y Paolone
(2018) para cualquier topología. Los dos están reimplementados en este repositorio a partir de sus artículos.

Dos advertencias que valen para toda la sección:

- **Las cotas de aquí son válidas en aritmética exacta. No hay redondeo dirigido**, así que todavía no se puede
  decir "verificado" en sentido estricto.
- **No hay un máximo ni un mínimo exactos contra los cuales comparar.** La referencia es el mejor valor
  *alcanzado* en cada barra: una búsqueda local (maximización y minimización lineal repetida sobre los
  semidiscos), los nueve flujos de frontera del método de Sou et al. y 4000 inyecciones aleatorias por
  caso. Un valor alcanzado acota el extremo real por dentro; por eso "exceso" y "holgura" son cotas superiores de
  lo floja que es cada cota, y una cota del lado equivocado de un valor alcanzado sería un error.

### 25.1 Montaje

- **Incertidumbre.** La de Sou et al.: la parte incierta de cada fuente con inversor, en un semidisco
  (Δp ≥ 0, |Δp + jΔq| ≤ C). La demanda es determinista.
- **Fuentes.** Su regla (Sec. VI-I): una fuente en la barra 1 + 10k, todas de igual potencia, con potencia total
  igual a la carga aparente total. Sus dos escenarios: flujo directo (ρD = 1, ρG = 0,5, fp = 0,8) e inverso
  (ρD = 0,5, ρG = 1, fp = 0,8). ρU es la fracción incierta de la potencia instalada; ellos usan 0,8 y 0,5.
- **Redes.** Cuatro alimentadores radiales de MATPOWER 7.1 (33, 69, 85 y 141 barras), los mismos de sus Tablas
  VI y VII, y cuatro redes donde su método no aplica: la de 33 barras con sus cinco enlaces cerrados, la de 33
  barras con capacidad de líneas y dos bancos de condensadores (modificación propia), la red CIGRE de media
  tensión y la red Oberrhein de pandapower, ambas con todos los interruptores cerrados, cables con capacidad y
  transformadores 110/20 kV.
- **Métodos.** *Encierro*: el de este repositorio (modelo afín de segundo orden, cota RAA-K del término
  cuadrático y cota del resto por punto fijo; sección 25.6). *WCLF*: Sou et al. *Discos*: Wang et al., Teorema 1,
  y su refinamiento por barra.

### 25.2 La reimplementación de WCLF frente a su artículo

| Tabla del artículo | Red | ρU (%) | ε medio publicado (%) | ε medio obtenido (%) | ε máximo publicado (%) | ε máximo obtenido (%) | Converge | Iteraciones |
|---|---|---:|---|---:|---|---:|---|---:|
| I | 33 | 0 | 7,5e-09 | 5,75e-09 | 1,3e-08 | 1,62e-08 | sí | 17 |
| I | 33 | 50 | 0,086 | 0,086 | 0,17 | 0,172 | sí | 17 |
| I | 33 | 100 | 0,18 | 0,176 | 0,34 | 0,344 | sí | 17 |
| VI | 33 | 80 | 0,065 | 0,0737 | 0,14 | 0,129 | sí | 17 |
| VII | 33 | 50 | 0,074 | 0,0812 | 0,14 | 0,128 | sí | 17 |
| VI | 69 | 80 | 0,012 | 0,00985 | 0,026 | 0,0193 | sí | 18 |
| VII | 69 | 50 | 0,020 | 0,0155 | 0,049 | 0,0364 | sí | 18 |
| VI | 85 | 80 | 0,14 | 0,146 | 0,19 | 0,188 | sí | 22 |
| VII | 85 | 50 | 0,17 | 0,166 | 0,20 | 0,197 | sí | 22 |
| VI | 141 | 80 | 0,051 | 0,0514 | 0,084 | 0,0843 | sí | 9 |
| VII | 141 | 50 | 0,058 | 0,0582 | 0,074 | 0,074 | sí | 8 |

*Script: `experiments/x8_feeders.py` · CSV: `results/x8_wclf_validation.csv`*

- **La Tabla I se reproduce** a los dígitos publicados (sí). Es el único caso con el montaje
  completamente especificado: las cinco fuentes y sus potencias están dibujadas en su Fig. 10.
- **Las Tablas VI y VII se reproducen en 3 de 8 filas** (redes: 85 barras, Tabla VII; 141 barras, Tabla VI; 141 barras, Tabla VII). En las
  demás la diferencia no la explica el método, que es el mismo de la Tabla I, sino el montaje: la regla "barra
  1 + 10k" y "carga aparente total" admite varias lecturas. Probé otras (ubicación de las fuentes, potencia
  total, numeración por recorrido del árbol) y ninguna reproduce las cuatro redes a la vez; esas pruebas fueron
  exploratorias y no están en los CSV. **Causa no identificada.**
- **Valor inicial.** El artículo propone partir de una cota de corriente "grande, p. ej. 5.000 A". Con ese valor
  las cotas de las líneas aguas arriba crecen como el cuadrado de las anteriores hasta desbordar; como la
  actualización de una línea solo usa las de sus descendientes, el resultado no depende del valor inicial una
  vez que se limita ese crecimiento. La reimplementación lo limita y una prueba unitaria comprueba que el
  resultado es el mismo partiendo de 1 a 1.000 veces la corriente muestreada, de 5.000 A y de 10⁹ A.

### 25.3 Redes radiales: frente a WCLF

| Red | Flujo | ρU (%) | Rango alcanzado, máx. (p.u.) | Encierro / rango alcanzado, mediana | Exceso de la cota superior: encierro (%) | ídem, WCLF (%) | Error certificado por WCLF (%) | Holgura de la cota inferior: encierro (%) | Discos / rango alcanzado, mediana | Encierro (ms) | WCLF (ms) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 33 | directo | 5 | 0,0063 | 1,001 | 0,0004 | 0,0013 | 0,0050 | 0,0000 | 1,43 | 66 | 9 |
| 33 | directo | 10 | 0,0125 | 1,002 | 0,0016 | 0,0023 | 0,0097 | 0,0001 | 1,44 | 25 | 8 |
| 33 | directo | 20 | 0,0253 | 1,005 | 0,0069 | 0,0038 | 0,0186 | 0,0009 | 1,44 | 24 | 8 |
| 33 | directo | 50 | 0,0645 | 1,019 | 0,0584 | 0,0087 | 0,0452 | 0,0158 | 1,46 | 24 | 9 |
| 33 | directo | 80 | 0,1059 | 1,044 | 0,2001 | 0,0158 | 0,0737 | 0,0776 | 1,50 | 24 | 8 |
| 33 | inverso | 5 | 0,0055 | 1,001 | 0,0002 | 0,0122 | 0,0153 | 0,0000 | 1,54 | 26 | 10 |
| 33 | inverso | 10 | 0,0110 | 1,002 | 0,0010 | 0,0220 | 0,0282 | 0,0001 | 1,53 | 26 | 9 |
| 33 | inverso | 20 | 0,0222 | 1,004 | 0,0046 | 0,0353 | 0,0478 | 0,0007 | 1,52 | 27 | 10 |
| 33 | inverso | 50 | 0,0577 | 1,016 | 0,0408 | 0,0494 | 0,0812 | 0,0118 | 1,47 | 27 | 10 |
| 33 | inverso | 80 | 0,0962 | 1,036 | 0,1440 | 0,0516 | 0,1032 | 0,0562 | 1,46 | 26 | 10 |
| 69 | directo | 5 | 0,0037 | 1,001 | 0,0001 | 0,0002 | 0,0005 | 0,0000 | 1,57 | 215 | 12 |
| 69 | directo | 10 | 0,0074 | 1,001 | 0,0004 | 0,0004 | 0,0010 | 0,0000 | 1,57 | 215 | 10 |
| 69 | directo | 20 | 0,0148 | 1,002 | 0,0015 | 0,0007 | 0,0019 | 0,0002 | 1,58 | 220 | 16 |
| 69 | directo | 50 | 0,0375 | 1,008 | 0,0118 | 0,0022 | 0,0053 | 0,0027 | 1,60 | 269 | 10 |
| 69 | directo | 80 | 0,0606 | 1,015 | 0,0370 | 0,0049 | 0,0098 | 0,0123 | 1,63 | 220 | 11 |
| 69 | inverso | 5 | 0,0034 | 1,000 | 0,0001 | 0,0030 | 0,0032 | 0,0000 | 1,53 | 128 | 11 |
| 69 | inverso | 10 | 0,0068 | 1,001 | 0,0003 | 0,0054 | 0,0060 | 0,0000 | 1,53 | 130 | 16 |
| 69 | inverso | 20 | 0,0137 | 1,002 | 0,0011 | 0,0089 | 0,0100 | 0,0001 | 1,53 | 159 | 11 |
| 69 | inverso | 50 | 0,0351 | 1,007 | 0,0092 | 0,0128 | 0,0155 | 0,0022 | 1,53 | 147 | 10 |
| 69 | inverso | 80 | 0,0575 | 1,013 | 0,0295 | 0,0139 | 0,0183 | 0,0099 | 1,54 | 163 | 11 |
| 85 | directo | 5 | 0,0100 | 1,003 | 0,0021 | 0,0024 | 0,0099 | 0,0001 | 1,52 | 1220 | 13 |
| 85 | directo | 10 | 0,0201 | 1,006 | 0,0090 | 0,0045 | 0,0194 | 0,0006 | 1,53 | 1199 | 17 |
| 85 | directo | 20 | 0,0405 | 1,015 | 0,0398 | 0,0082 | 0,0378 | 0,0051 | 1,55 | 1245 | 13 |
| 85 | directo | 50 | 0,1050 | 1,056 | 0,3453 | 0,0184 | 0,0915 | 0,1041 | 1,62 | 1276 | 22 |
| 85 | directo | 80 | 0,1758 | 1,139 | 1,2895 | 0,0303 | 0,1459 | 0,6203 | 1,79 | 1274 | 15 |
| 85 | inverso | 5 | 0,0081 | 1,002 | 0,0012 | 0,0269 | 0,0321 | 0,0000 | 1,70 | 1195 | 15 |
| 85 | inverso | 10 | 0,0164 | 1,005 | 0,0049 | 0,0479 | 0,0585 | 0,0004 | 1,69 | 1145 | 17 |
| 85 | inverso | 20 | 0,0334 | 1,011 | 0,0225 | 0,0747 | 0,0962 | 0,0032 | 1,67 | 1206 | 15 |
| 85 | inverso | 50 | 0,0880 | 1,042 | 0,2053 | 0,1103 | 0,1663 | 0,0617 | 1,61 | 1153 | 16 |
| 85 | inverso | 80 | 0,1502 | 1,100 | 0,7675 | 0,1157 | 0,2092 | 0,3374 | 1,61 | 1219 | 16 |
| 141 | directo | 5 | 0,0053 | 1,002 | 0,0007 | 0,0002 | 0,0031 | 0,0000 | 1,36 | 6503 | 12 |
| 141 | directo | 10 | 0,0106 | 1,004 | 0,0030 | 0,0004 | 0,0061 | 0,0001 | 1,36 | 6406 | 12 |
| 141 | directo | 20 | 0,0212 | 1,008 | 0,0125 | 0,0008 | 0,0121 | 0,0006 | 1,37 | 6435 | 11 |
| 141 | directo | 50 | 0,0538 | 1,022 | 0,0870 | 0,0030 | 0,0310 | 0,0096 | 1,39 | 6371 | 13 |
| 141 | directo | 80 | 0,0876 | 1,043 | 0,2509 | 0,0070 | 0,0514 | 0,0448 | 1,41 | 6398 | 13 |
| 141 | inverso | 5 | 0,0048 | 1,002 | 0,0005 | 0,0092 | 0,0116 | 0,0000 | 1,43 | 5905 | 20 |
| 141 | inverso | 10 | 0,0097 | 1,003 | 0,0022 | 0,0165 | 0,0213 | 0,0001 | 1,42 | 6740 | 18 |
| 141 | inverso | 20 | 0,0195 | 1,007 | 0,0094 | 0,0258 | 0,0353 | 0,0005 | 1,41 | 6941 | 17 |
| 141 | inverso | 50 | 0,0501 | 1,021 | 0,0687 | 0,0339 | 0,0582 | 0,0078 | 1,39 | 7331 | 19 |
| 141 | inverso | 80 | 0,0826 | 1,040 | 0,2055 | 0,0351 | 0,0744 | 0,0357 | 1,38 | 7323 | 28 |

*Script: `experiments/x8_feeders.py` · CSV: `results/x8_feeders.csv` · semilla por fila en el CSV. "Exceso" =
(cota superior − máximo alcanzado) / máximo alcanzado; "holgura" = (mínimo alcanzado − cota inferior) / mínimo
alcanzado; ambos promediados sobre las barras. Los tiempos son de una sola corrida en Python.*

- **El punto fijo cierra en 40 de los 40 casos radiales**, incluido ρU = 80 %, donde el
  rango de voltaje alcanzado llega a 0,18 p.u. La cota de la sección 22 dejaba de cerrar mucho antes (sección 25.6).
- **Ninguna cota queda del lado equivocado de un valor alcanzado**: 0 violaciones entre los tres
  métodos en toda la sección.
- **Cota superior.** La del encierro es más ajustada que la de WCLF en 19 de 40 casos:
  14 de 20 en flujo inverso y 5 de 20 en flujo directo. Con
  incertidumbre grande (ρU ≥ 50 %) en flujo directo gana WCLF en 8 de 8
  casos; con incertidumbre chica (ρU ≤ 10 %) en flujo inverso gana el encierro en 8 de
  8. El peor caso para el encierro es la red de 85 barras con
  ρU = 80 %: exceso de 1,29 % contra 0,030 % de WCLF, una razón de
  43. En el otro extremo la razón a favor del encierro llega a 51.
- **El "error certificado" de WCLF lo domina su cota inferior.** Esa cota sale de nueve flujos con la misma
  dirección para todas las fuentes; la búsqueda local encuentra máximos mayores, y esa diferencia explica entre
  8 y 93 % del error que el método certifica. Su cota superior sola es bastante
  mejor que lo que su propia métrica sugiere.
- **Cota inferior.** WCLF da el mínimo con un solo flujo en (Δp, Δq) = (0, −C), suponiendo monotonía. En estos
  cuatro alimentadores la suposición se cumple: su valor coincide con el mínimo alcanzado en todas las barras
  (diferencia máxima 4,8·10⁻¹⁰ p.u.). El encierro da una cota inferior que no depende de esa suposición, con
  holgura de hasta 0,0010 % para ρU ≤ 10 % y de hasta 1,12 % en el peor caso. En otras
  redes la suposición no se cumple (sección 25.4).
- **Discos de Wang et al.** Aplican en todos los casos salvo 2. El disco uniforme es entre
  1,4 y 4,6 veces el rango alcanzado y el refinado por barra entre 1,24 y
  2,17 (todas las redes de la sección).
- **Costo.** El encierro tarda entre 23 ms y 16,5 s (todas las redes de la sección) y
  la reimplementación de WCLF entre 8 y 28 ms. Sus autores reportan 3 a 15 ms en estas redes y 0,53 s en 10.476
  barras; el encierro usa matrices densas y aquí no pasa de 179 barras.

**Lectura.** En redes radiales con este tipo de incertidumbre, WCLF sigue siendo la referencia: es más ajustado
cuando la incertidumbre es grande, da el mínimo exacto y escala a diez mil barras. Donde el encierro gana
(incertidumbre chica, flujo inverso) las diferencias son de milésimas de punto porcentual de voltaje, sin
consecuencia práctica. **Aquí no hay aporte.**

### 25.4 Redes malladas y con elementos en derivación: donde WCLF no aplica

| Red | Barras | Mallas | Símbolos | Flujo | ρU (%) | Rango alcanzado, máx. (p.u.) | Punto fijo | Encierro / rango alcanzado, mediana | máx. | Exceso superior medio (%) | Holgura inferior media (%) | Barras con un valor alcanzado fuera | ídem, solo 2.º orden | Discos aplicables | Discos / rango alcanzado, mediana | Encierro (ms) |
|---|---:|---:|---:|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---|---:|---:|
| 33 + ties | 33 | 5 | 8 | directo | 5 | 0,0036 | sí | 1,001 | 1,001 | 0,0001 | 0,0000 | 0 | 26 | sí | 1,32 | 50 |
| 33 + ties | 33 | 5 | 8 | directo | 10 | 0,0072 | sí | 1,001 | 1,001 | 0,0006 | 0,0000 | 0 | 32 | sí | 1,32 | 43 |
| 33 + ties | 33 | 5 | 8 | directo | 20 | 0,0144 | sí | 1,003 | 1,003 | 0,0024 | 0,0002 | 0 | 32 | sí | 1,32 | 42 |
| 33 + ties | 33 | 5 | 8 | directo | 50 | 0,0364 | sí | 1,008 | 1,009 | 0,0175 | 0,0028 | 0 | 32 | sí | 1,33 | 23 |
| 33 + ties | 33 | 5 | 8 | directo | 80 | 0,0590 | sí | 1,016 | 1,017 | 0,0526 | 0,0123 | 0 | 32 | sí | 1,34 | 39 |
| 33 + ties | 33 | 5 | 8 | inverso | 5 | 0,0034 | sí | 1,001 | 1,001 | 0,0001 | 0,0000 | 0 | 20 | sí | 1,37 | 44 |
| 33 + ties | 33 | 5 | 8 | inverso | 10 | 0,0067 | sí | 1,001 | 1,001 | 0,0004 | 0,0000 | 0 | 29 | sí | 1,37 | 50 |
| 33 + ties | 33 | 5 | 8 | inverso | 20 | 0,0136 | sí | 1,002 | 1,002 | 0,0019 | 0,0001 | 0 | 30 | sí | 1,36 | 48 |
| 33 + ties | 33 | 5 | 8 | inverso | 50 | 0,0347 | sí | 1,008 | 1,008 | 0,0147 | 0,0025 | 0 | 32 | sí | 1,34 | 31 |
| 33 + ties | 33 | 5 | 8 | inverso | 80 | 0,0567 | sí | 1,015 | 1,016 | 0,0456 | 0,0109 | 0 | 32 | sí | 1,33 | 23 |
| 33 + shunts | 33 | 0 | 8 | directo | 5 | 0,0062 | sí | 1,001 | 1,001 | 0,0003 | 0,0000 | 0 | 28 | sí | 1,45 | 37 |
| 33 + shunts | 33 | 0 | 8 | directo | 10 | 0,0124 | sí | 1,002 | 1,002 | 0,0015 | 0,0001 | 0 | 32 | sí | 1,45 | 27 |
| 33 + shunts | 33 | 0 | 8 | directo | 20 | 0,0250 | sí | 1,005 | 1,005 | 0,0065 | 0,0008 | 0 | 32 | sí | 1,45 | 30 |
| 33 + shunts | 33 | 0 | 8 | directo | 50 | 0,0638 | sí | 1,019 | 1,021 | 0,0552 | 0,0151 | 0 | 32 | sí | 1,48 | 52 |
| 33 + shunts | 33 | 0 | 8 | directo | 80 | 0,1046 | sí | 1,043 | 1,047 | 0,1887 | 0,0738 | 0 | 32 | sí | 1,51 | 45 |
| 33 + shunts | 33 | 0 | 8 | inverso | 5 | 0,0054 | sí | 1,001 | 1,001 | 0,0002 | 0,0000 | 0 | 21 | sí | 1,55 | 33 |
| 33 + shunts | 33 | 0 | 8 | inverso | 10 | 0,0109 | sí | 1,002 | 1,002 | 0,0010 | 0,0001 | 0 | 22 | sí | 1,54 | 55 |
| 33 + shunts | 33 | 0 | 8 | inverso | 20 | 0,0220 | sí | 1,004 | 1,004 | 0,0044 | 0,0007 | 0 | 26 | sí | 1,53 | 48 |
| 33 + shunts | 33 | 0 | 8 | inverso | 50 | 0,0571 | sí | 1,016 | 1,017 | 0,0387 | 0,0113 | 0 | 32 | sí | 1,49 | 52 |
| 33 + shunts | 33 | 0 | 8 | inverso | 80 | 0,0951 | sí | 1,036 | 1,039 | 0,1364 | 0,0537 | 0 | 32 | sí | 1,47 | 25 |
| CIGRE MV | 15 | 3 | 4 | directo | 5 | 0,0250 | sí | 1,005 | 1,009 | 0,0084 | 0,0010 | 0 | 12 | sí | 1,64 | 43 |
| CIGRE MV | 15 | 3 | 4 | directo | 10 | 0,0505 | sí | 1,012 | 1,020 | 0,0387 | 0,0080 | 0 | 12 | sí | 1,66 | 34 |
| CIGRE MV | 15 | 3 | 4 | directo | 20 | 0,1032 | sí | 1,036 | 1,051 | 0,2043 | 0,0746 | 0 | 13 | sí | 1,72 | 30 |
| CIGRE MV | 15 | 3 | 4 | directo | 50 | 0,2812 | sí | 1,305 | 1,335 | 3,5434 | 2,9439 | 0 | 14 | sí | 2,16 | 47 |
| CIGRE MV | 15 | 3 | 4 | directo | 80 | 0,5427 | no | – | – | – | – | – | 14 | no | – | 28 |
| CIGRE MV | 15 | 3 | 4 | inverso | 5 | 0,0210 | sí | 1,004 | 1,009 | 0,0056 | 0,0012 | 0 | 4 | sí | 2,06 | 35 |
| CIGRE MV | 15 | 3 | 4 | inverso | 10 | 0,0426 | sí | 1,011 | 1,020 | 0,0261 | 0,0078 | 0 | 9 | sí | 2,06 | 28 |
| CIGRE MV | 15 | 3 | 4 | inverso | 20 | 0,0875 | sí | 1,033 | 1,048 | 0,1402 | 0,0632 | 0 | 10 | sí | 2,06 | 31 |
| CIGRE MV | 15 | 3 | 4 | inverso | 50 | 0,2421 | sí | 1,253 | 1,264 | 2,3585 | 1,9648 | 0 | 13 | sí | 2,17 | 25 |
| CIGRE MV | 15 | 3 | 4 | inverso | 80 | 0,4606 | no | – | – | – | – | – | 14 | no | – | 31 |
| Oberrhein | 179 | 5 | 36 | directo | 5 | 0,0132 | sí | 1,004 | 1,004 | 0,0045 | 0,0001 | 0 | 173 | sí | 1,24 | 16483 |
| Oberrhein | 179 | 5 | 36 | directo | 10 | 0,0266 | sí | 1,009 | 1,009 | 0,0187 | 0,0007 | 0 | 177 | sí | 1,25 | 14868 |
| Oberrhein | 179 | 5 | 36 | directo | 20 | 0,0538 | sí | 1,019 | 1,020 | 0,0799 | 0,0064 | 0 | 177 | sí | 1,26 | 14768 |
| Oberrhein | 179 | 5 | 36 | directo | 50 | 0,1393 | sí | 1,066 | 1,069 | 0,6274 | 0,1361 | 0 | 177 | sí | 1,31 | 14558 |
| Oberrhein | 179 | 5 | 36 | directo | 80 | 0,2339 | sí | 1,156 | 1,164 | 2,1924 | 0,8731 | 0 | 177 | sí | 1,40 | 14777 |
| Oberrhein | 179 | 5 | 36 | inverso | 5 | 0,0117 | sí | 1,004 | 1,004 | 0,0032 | 0,0001 | 0 | 41 | sí | 1,31 | 9399 |
| Oberrhein | 179 | 5 | 36 | inverso | 10 | 0,0236 | sí | 1,007 | 1,009 | 0,0133 | 0,0009 | 0 | 75 | sí | 1,31 | 9876 |
| Oberrhein | 179 | 5 | 36 | inverso | 20 | 0,0479 | sí | 1,016 | 1,019 | 0,0571 | 0,0059 | 0 | 147 | sí | 1,30 | 10222 |
| Oberrhein | 179 | 5 | 36 | inverso | 50 | 0,1257 | sí | 1,056 | 1,060 | 0,4582 | 0,1028 | 0 | 175 | sí | 1,28 | 11341 |
| Oberrhein | 179 | 5 | 36 | inverso | 80 | 0,2140 | sí | 1,124 | 1,130 | 1,5920 | 0,5746 | 0 | 177 | sí | 1,27 | 12268 |

*Script: `experiments/x8_feeders.py` · CSV: `results/x8_feeders.csv`*

- **El encierro cierra en 38 de los 40 casos.** No cierra en: CIGRE MV con ρU = 80 % (flujo directo); CIGRE MV con ρU = 80 % (flujo inverso). Donde
  cierra queda entre 1,001 y 1,305 veces el rango alcanzado (mediana por caso;
  1,335 en la peor barra).
- **Se degrada cuando el rango de voltaje es grande.** Con ρU ≤ 20 % el encierro queda a lo más a
  1,036 veces el rango en las cuatro redes. El peor caso es la red CIGRE, donde la regla de
  montaje pone toda la potencia en dos fuentes: con ρU = 50 % el rango alcanzado llega a 0,28 p.u.,
  el encierro queda a 1,31 veces ese rango, y con ρU = 80 % ni el encierro ni los discos cierran.
- **La cota de segundo orden sola no alcanza:** sin la holgura del punto fijo, la cota inferior deja fuera algún
  valor alcanzado en 80 de los 80 casos de la sección (por hasta 0,0741 p.u.), y
  la superior en 0. Es el mismo hallazgo de las secciones 18 y 20, ahora en distribución.
- **La garantía publicada para estas redes es bastante más ancha.** Los discos de Wang et al. contienen todo lo
  alcanzado, pero en estas cuatro redes su ancho es de 1,24 a 2,17 veces el rango
  (refinados por barra) y de 1,4 a 2,4 veces sin refinar.
- **Tamaño.** La red Oberrhein tiene 179 barras y 36 símbolos; el encierro tarda hasta
  16,5 s. En un perfil aparte (sin CSV) casi todo ese tiempo se va en la cota RAA-K de las formas
  cuadráticas.

**La monotonía no se cumple en todas las redes.** El mínimo de voltaje de cada barra, ¿está en (0, −C)?

| Red | Casos | Barras PQ | Casos con alguna barra cuyo mínimo no está en (0, −C) | Barras en esa situación, máx. | Mayor diferencia (p.u.) | ídem, como % del rango alcanzado de la barra |
|---|---:|---:|---:|---:|---|---:|
| 33 | 10 | 32 | 0 | 0 | – | – |
| 69 | 10 | 68 | 0 | 0 | – | – |
| 85 | 10 | 84 | 0 | 0 | – | – |
| 141 | 10 | 140 | 0 | 0 | – | – |
| 33 + ties | 10 | 32 | 0 | 0 | – | – |
| 33 + shunts | 10 | 32 | 0 | 0 | – | – |
| CIGRE MV | 10 | 14 | 10 | 3 | 0,0017 | 4,2 |
| Oberrhein | 10 | 177 | 10 | 122 | 1,7·10⁻⁴ | 0,5 |

*Script: `experiments/x8_feeders.py monotone` · CSV: `results/x8_monotonicity.csv`. Se compara el voltaje en
(0, −C) con el menor valor que encuentra la búsqueda local; uno menor lo alcanza una inyección factible.*

- En 33, 69, 85, 141, 33 + ties, 33 + shunts el mínimo está en (0, −C) en todas las barras.
- En CIGRE MV y Oberrhein no: en 20 de los 80 casos hay barras con un mínimo menor, hasta
  122 barras a la vez y por hasta 0,0017 p.u. (4,2 % del rango de esa
  barra). Son las dos redes con transformadores y cables con capacidad. La diferencia es chica, pero un método
  que tome el mínimo de un solo flujo deja de ser exacto ahí. La tesis de C. Wang reporta lo mismo en el IEEE 37.

### 25.5 Cajas de demanda: las dos cotas del resto

El montaje de la sección 22 (demanda incierta ±δ en 8 zonas, caja simétrica), ahora con las dos cotas del resto
y contra los 256 vértices.

| Red | K | δ (%) | Rango exacto en vértices, máx. (p.u.) | Punto fijo, cota de la sección 22 | Encierro / rango, mediana | Punto fijo, cota por inyección de corriente | Encierro / rango, mediana | máx. | Barras con algún vértice fuera |
|---|---:|---:|---:|---|---:|---|---:|---:|---:|
| 33 | 8 | 5 | 0,0094 | sí | 1,0024 | sí | 1,0023 | 1,0024 | 0 |
| 33 | 8 | 15 | 0,0283 | sí | 1,0090 | sí | 1,0068 | 1,0071 | 0 |
| 33 | 8 | 30 | 0,0567 | no | – | sí | 1,0138 | 1,0143 | 0 |
| 33 | 8 | 50 | 0,0948 | no | – | sí | 1,0234 | 1,0244 | 0 |
| 33 + ties | 8 | 5 | 0,0049 | sí | 1,0010 | sí | 1,0010 | 1,0011 | 0 |
| 33 + ties | 8 | 15 | 0,0146 | sí | 1,0033 | sí | 1,0031 | 1,0033 | 0 |
| 33 + ties | 8 | 30 | 0,0292 | sí | 1,0073 | sí | 1,0062 | 1,0067 | 0 |
| 33 + ties | 8 | 50 | 0,0487 | sí | 1,0156 | sí | 1,0104 | 1,0113 | 0 |

*Script: `experiments/x8_feeders.py` · CSV: `results/x8_load_boxes.csv`*

- La cota de la sección 22 cierra en 6 de 8 casos; la nueva, en 8, hasta
  δ = 50 %, con el encierro entre 1,001 y 1,023 veces el rango exacto en
  vértices y 0 barras con algún vértice fuera.

### 25.6 Qué cambió en el método

**Por qué fallaba la cota de la sección 22.** Escribía el resto en coordenadas rectangulares y acotaba el término
cuadrático con valores absolutos de J⁻¹ΔJ. Cuando dos barras están unidas por una rama muy corta, esa cota
trata sus voltajes como independientes y la corriente entre ellas se dispara. Era un defecto de la cota, no del
problema: deja de cerrar aunque la parte lineal contraiga con holgura (radio de G desde 0,02 en
los casos en que falla).

| Red | Casos | Rama más corta, \|z\| (p.u.) | Cierra la cota de la sección 22 | Cierra la cota por inyección de corriente | Cota del tercer orden: sección 22 / nueva, mediana donde ambas cierran |
|---|---:|---|---:|---:|---:|
| 33 | 10 | 0,0065 | 6 | 10 | 2,2 |
| 69 | 10 | 8,1·10⁻⁵ | 2 | 10 | 8,7 |
| 85 | 10 | 8,1·10⁻⁴ | 2 | 10 | 6,1 |
| 141 | 10 | 6,4·10⁻⁷ | 0 | 10 | – |
| 33 + ties | 10 | 0,0065 | 8 | 10 | 2,0 |
| 33 + shunts | 10 | 0,0065 | 6 | 10 | 2,2 |
| CIGRE MV | 10 | 5,2·10⁻⁴ | 3 | 8 | 2,9 |
| Oberrhein | 10 | 3,5·10⁻⁵ | 0 | 10 | – |

*Script: `experiments/x8_feeders.py remainder` · CSV: `results/x8_remainder_bounds.csv`. Mismos casos que las
tablas de 25.3 y 25.4.*

- La cota de la sección 22 cierra en 27 de 80 casos y la nueva en 78. La de la
  sección 22 no cierra en ningún caso de: 141, Oberrhein.
- Donde cierran las dos, la cota del tercer orden de la sección 22 es entre 1,6 y
  9 veces la nueva.

**La cota nueva** (`cdhim/enclosure.py`, clase `InjectionEnclosure`). Las ecuaciones exactas son
J D = ΔS − D·conj(ΔI). En una red con barras PQ la corriente inyectada es conj(S/V), así que conj(ΔI) depende solo
del voltaje de la misma barra. Con d = D/V*, σ = S*/V*² y τ = ΔS/V*²:

    J D = ΔS + V*² [ −τ d + (σ + τ) d² / (1 + d) ]

Los términos de primer y segundo orden dan el modelo afín y la forma cuadrática de siempre; lo que queda es de
tercer orden y se acota barra por barra con módulos, pesado por las normas de los bloques 2×2 de J⁻¹. La matriz
de admitancias no aparece en la cota y las barras sin inyección no contribuyen. Si la iteración de esas cotas
tiene un punto fijo, el teorema de Brouwer da una solución dentro del producto de discos para toda inyección
del conjunto. El programa infla el punto fijo calculado en 10⁻⁹ y comprueba la desigualdad de auto-mapeo.

**Comprobado contra el flujo AC exacto** (`tests/test_injection_enclosure.py`, en red radial, mallada y con
elementos en derivación): la identidad de arriba, a 10⁻⁸; cada cota parcial y la del resto en los 256 vértices;
que el resto escala como el cubo del tamaño de la caja; y que el encierro contiene los vértices, 2.000 muestras
y los extremos de la búsqueda local.

**Semidiscos.** El modelo afín trabaja con cajas, así que cada semidisco se cubre con su caja. El término de
primer orden, que domina, se maximiza sobre los semidiscos mismos con una fórmula cerrada; los órdenes
superiores se acotan sobre las cajas.

**Parentesco con lo publicado.** La forma de inyección de corriente y el punto fijo son los de Bolognani y
Zampieri (2016) y Wang et al. (2018), que acotan con un radio el error respecto de un punto. Lo que se agrega
aquí es restar el primer y el segundo orden exactos, acotar el segundo con RAA-K y llevar la cota del resto
barra por barra. Jabr (2019) usa órdenes 2 y 3 en distribución, mallada incluida, sin cota del resto.

### 25.7 Lectura crítica: qué queda de la opción A

1. **En redes radiales no hay aporte** (sección 25.3).
2. **En redes malladas o con elementos en derivación y barras PQ sí hay algo medible:** un encierro por barra a
   1,001–1,305 veces el rango alcanzado, donde WCLF no aplica y la garantía publicada es de
   1,24 a 2,17 veces el rango. En lo que se pudo leer, nadie reporta esa comparación en una red de distribución mallada. La ventaja
   se achica cuando el rango de voltaje es grande, y el encierro deja de cerrar antes que el flujo de carga
   deje de tener solución (red CIGRE, sección 25.4). En dos de esas redes, además, el mínimo no está donde lo
   pondría la monotonía, así que una cota inferior que no la suponga tiene sentido.
3. **El mecanismo no es nuevo; la combinación y la medición sí lo parecen.** Un revisor puede leerlo como una
   mejora incremental de la familia de Wang et al. El argumento tendría que ser la tabla de la sección 25.4.
4. **Lo que falta para sostenerlo:**
   - redondeo dirigido, sin el cual no es un encierro verificado;
   - una referencia exacta, o al menos una cota interior certificada, en vez de valores alcanzados;
   - la comparación con el método de admisibilidad de Wang, Le Boudec y Paolone (2019), que no se pudo leer
     completo; aquí está solo su cota por barra en lectura monofásica;
   - redes desbalanceadas: la distribución real lo es y la familia de Wang et al. ya es multifásica;
   - barras de voltaje controlado, que esta cota no admite;
   - escala: todo es denso y se probó hasta 179 barras.
5. **El conjunto de incertidumbre favorece a los métodos ajenos.** Semidiscos y demanda determinista son el caso
   de Sou et al. Con demanda incierta correlacionada por zonas (sección 25.5), que es el caso natural de la
   aritmética afín, WCLF no está definido y el encierro queda aún más cerca del rango exacto.
### 25.8 Hallazgo lateral: el operador RAA-K depende de la escala de la forma

Al perfilar la red de 179 barras apareció una inconsistencia entre la teoría y el código. La cota de
εᵀ(sQ)ε debe ser s veces la de εᵀQε, pero la búsqueda de bloques usa **umbrales absolutos de 10⁻¹²**. El port
de este repositorio los hereda del código C++ del grupo, al que reproduce bit a bit. El que actúa primero es el
del producto de dos entradas de la diagonal (`cp <= 1e-12` en `_benefit`): descarta todo par cuyas diagonales
tengan media geométrica menor o igual que 10⁻⁶. Las formas cuadráticas del flujo de potencia en p.u. están en
esa zona (columna «Mayor coeficiente» de la tabla): el operador deja de formar bloques y la cota se acerca a la
de AF2, término a término. Basta normalizar la forma por su mayor coeficiente antes de reescribirla y reescalar
la cota.

La tabla repite el montaje de las secciones 18 a 23 (demanda incierta por zonas) con y sin normalizar:

| Sistema | K | δ (%) | δ de diseño | Formas | Mayor coeficiente, mediana | RAA-K sin normalizar / AF2 | RAA-K normalizada / AF2 | Sin normalizar / rango en vértices | Normalizada / rango en vértices | Formas que mejoran más de 0,1 % | Mayor mejora (%) |
|---|---:|---:|---|---:|---|---:|---:|---:|---:|---:|---:|
| IEEE-33 | 8 | 5 | no | 32 | 2,0·10⁻⁶ | 0,966 | 0,599 | 1,779 | 1,091 | 32 | 41 |
| IEEE-33 | 8 | 15 | sí | 32 | 1,8·10⁻⁵ | 0,631 | 0,599 | 1,154 | 1,091 | 32 | 40 |
| IEEE-33 | 8 | 30 | no | 32 | 7,1·10⁻⁵ | 0,604 | 0,599 | 1,099 | 1,091 | 31 | 36 |
| IEEE-33 | 16 | 5 | no | 32 | 1,3·10⁻⁶ | 1,000 | 0,683 | – | – | 32 | 37 |
| IEEE-33 | 16 | 15 | sí | 32 | 1,2·10⁻⁵ | 0,834 | 0,683 | – | – | 32 | 31 |
| IEEE-33 | 16 | 30 | no | 32 | 4,9·10⁻⁵ | 0,703 | 0,683 | – | – | 31 | 30 |
| PEGASE-2869 | 8 | 1 | no | 2296 | 2,2·10⁻⁶ | 0,957 | 0,691 | 1,855 | 1,338 | 2218 | 45 |
| PEGASE-2869 | 8 | 2 | no | 2308 | 8,6·10⁻⁶ | 0,821 | 0,694 | 1,706 | 1,341 | 2139 | 45 |
| PEGASE-2869 | 8 | 5 | sí | 2312 | 5,3·10⁻⁵ | 0,741 | 0,694 | 1,460 | 1,343 | 2033 | 44 |
| PEGASE-2869 | 16 | 1 | no | 2296 | 1,2·10⁻⁶ | 1,000 | 0,666 | – | – | 2218 | 45 |
| PEGASE-2869 | 16 | 2 | no | 2309 | 4,8·10⁻⁶ | 0,887 | 0,667 | – | – | 2217 | 44 |
| PEGASE-2869 | 16 | 5 | sí | 2312 | 3,0·10⁻⁵ | 0,753 | 0,667 | – | – | 2166 | 43 |

*Script: `experiments/x8_feeders.py` · CSV: `results/x8_raak_scale.csv`*

- **Normalizada, la cota es la misma fracción de AF2 para todo δ**, como dice la teoría: entre 0,60 y
  0,69. Sin normalizar va de 0,60 a 1,00, y con δ chico no mejora nada a AF2.
- **En los δ de diseño**, que son los que usan las secciones 18 a 23, la cota sin normalizar es
  0,63 a 0,83 veces AF2 y la normalizada 0,60 a
  0,69. Frente al rango exacto en vértices, con 8 zonas: 1,15 a
  1,46 sin normalizar y 1,09 a 1,34 normalizada. Hay formas que
  se estrechan hasta 45 %.
- **No es un problema de validez.** Las dos cotas son válidas; la normalizada contiene los valores en los vértices
  en todas las formas con 8 zonas (0 excepciones). La búsqueda es voraz y no es monótona: en
  31 de las 14.025 formas la cota normalizada sale más ancha que la otra.
- **Qué se llama AF2 en esta sección.** Es la cota término a término de la forma ya simetrizada,
  Σ H_kk·[0, 1] + Σ_{k<l} 2H_kl·[−1, 1] (`bound_monomial`). En el código del grupo esa cota se llama RAA con
  k = 0 (`knomial_quad_bound` con `max_k < 2`), y «fAF2» es otra: la del producto de dos formas afines sin
  simetrizar, de radio Sx·Sy − Sxy/2, igual o más ancha. Las razones de la tabla se comparan con las razones
  «k / k = 0» del grupo, no con «k / fAF2».

Consecuencias:

1. **Las secciones 18, 21, 22 y 23 subestiman a RAA-K.** En particular la comparación con los bloques
   espectrales con 16 zonas (sección 18.6) y con la relajación semidefinida (sección 23) hay que rehacerlas; el
   veredicto "parcial" del aporte 3 puede cambiar. Quedan marcadas con un aviso. No las rehice aquí: son varias
   horas de cómputo y cambian texto ya entregado.
2. **Esta sección 25 ya usa la forma normalizada.**
3. **Para el grupo:** el operador conviene llamarlo siempre con la forma normalizada, o cambiar su tolerancia a
   relativa. Vale igual para la versión C++ con redondeo dirigido. No toqué `ibex-lib`.

---

## 26. Las comprobaciones que faltaban para sostener la opción A

La sección 25.7 listó lo que faltaba. Esta sección lo hace, o dice por qué no se pudo.

| Pendiente de 25.7 | Estado | Resultado en una línea | Dónde |
|---|---|---|---|
| Comparar con Krawczyk | Hecho | Un Krawczyk bien armado cierra en los mismos casos y queda cerca; el encierro es más ajustado en todas las barras, por un orden de exactitud, y más caro | 26.1 |
| Leer Wang, Le Boudec y Paolone (2019) | **Bloqueado: falta el PDF** | Por fuentes indirectas, su cota es la ya implementada y su vía al ajuste es subdividir; eso se midió | 26.2 |
| Redondeo dirigido | Hecho | El encierro está verificado en punto flotante, con unicidad, sin perder ancho | 26.3 |
| Referencia exacta o certificada | Hecho | Cota interior certificada en todos los casos; rango exacto certificado en cajas de pocos símbolos | 26.4 |
| Escala | Hecho, con un límite medido | Verifica hasta 2.940 barras con 64 símbolos y hasta 533 con 108; lo que limita es el número de símbolos | 26.5 |
| Barras de voltaje controlado | Hecho en distribución; límite en transmisión | Verifica con fuentes de voltaje controlado en la red de 33 barras y en la de 30; no en las de 118 y 300 | 26.6 |
| Redes desbalanceadas | Hecho, con supuestos | Funciona sin cambios con cargas en estrella; delta y transformadores quedan fuera | 26.7 |

Un resultado no previsto: al comparar con Krawczyk apareció que casi todo el exceso del encierro venía de acotar
por separado el término de primer orden y el de segundo. Se agregó una cota conjunta (26.3) y **todas las cifras
de esta sección son con ella**; las de la sección 25 son anteriores y más anchas.

### 26.1 Krawczyk: el baseline que podía mostrar que el aporte ya existía

**Qué se armó** (`cdhim/krawczyk.py`). El operador de Krawczyk paramétrico sobre las ecuaciones del flujo de
potencia, K(X) = x* − C·F(x*, P) + (I − C·S(X, P))(X − x*), con las incógnitas en rectangulares y las inyecciones
afines en los mismos símbolos que usa el encierro. El término de primer orden se evalúa exacto por componente,
sin perder la dependencia entre símbolos. Variantes: formulación por desbalance de potencia o por inyección de
corriente; matriz de derivadas (da unicidad) o de pendientes (solo existencia); ejes globales o rotados barra
por barra para perder menos al pasar de la caja a |V|. No redondea hacia fuera.

No es más flojo que lo publicado. Rodrigues Junior, Melo y Nepomuceno (2022), que usan Krawczyk por inyección de
corriente con INTLAB, reportan intervalos de 1,8 a 2,5 veces el rango que da Monte Carlo en IEEE 13, 33 y 69 con
cargas de ±5 a ±15 %. La comparación no es directa, porque Monte Carlo no cubre el rango exacto.

| Formulación | Matriz | Ejes | Casos | Cierra (semidiscos) | Cierra (cajas) | Violaciones | Ancho / rango alcanzado, mediana | Tiempo, mediana (ms) |
|---|---|---|---:|---:|---:|---:|---|---|
| desbalance de potencia | derivadas | globales | 80 | 0 | 0 | 0 | – | 173 |
| desbalance de potencia | derivadas | rotados | 80 | 0 | 0 | 0 | – | 176 |
| desbalance de potencia | pendientes | globales | 80 | 0 | 0 | 0 | – | 166 |
| desbalance de potencia | pendientes | rotados | 80 | 0 | 0 | 0 | – | 161 |
| inyección de corriente | derivadas | globales | 80 | 75 | 75 | 0 | 1,004 a 1,706 | 26 |
| inyección de corriente | derivadas | rotados | 80 | 75 | 75 | 0 | 1,004 a 1,494 | 27 |
| inyección de corriente | pendientes | globales | 80 | 78 | 78 | 0 | 1,004 a 1,917 | 25 |
| inyección de corriente | pendientes | rotados | 80 | 78 | 78 | 0 | 1,003 a 1,501 | 26 |

*Script: `experiments/x9_krawczyk.py` · CSV: `results/x9_krawczyk.csv`. Los 80 casos de la sección 25.*

- **Por desbalance de potencia no cierra en ningún caso** (0 de 320 corridas). Trata como independientes
  los voltajes de barras vecinas; es el mismo defecto de la cota de la sección 22.
- **Por inyección de corriente cierra:** en 78 de 80 casos con pendientes y en 75
  con derivadas. Con pendientes falla en los mismos dos casos que el encierro.
- **Ninguna violación** en ninguna variante (0).
- Con demanda incierta en 8 zonas (`results/x9_krawczyk_load_boxes.csv`) cierra en 8 de
  8 casos y mide entre 1,002 y 1,047 veces el rango exacto en vértices.

**Frente al encierro verificado**, con la mejor variante (corriente, pendientes, ejes rotados) y contra la misma
referencia certificada de 26.4:

| Red | Casos | Krawczyk cierra | Encierro verifica | Krawczyk / rango exacto, a lo más (mediana) | Encierro / rango exacto, a lo más (mediana) | Ancho Krawczyk / ancho encierro | Exceso Krawczyk / exceso encierro | Barras en que Krawczyk es más ajustado, máx. (%) | Krawczyk, mediana (ms) | Encierro, mediana (ms) |
|---|---:|---:|---:|---|---|---|---|---|---|---|
| 33 | 10 | 10 | 10 | 1,006 a 1,137 | 1,000 a 1,017 | 1,006 a 1,119 | 8 a 141 | 0 | 6 | 19 |
| 69 | 10 | 10 | 10 | 1,003 a 1,061 | 1,000 a 1,004 | 1,003 a 1,056 | 14 a 224 | 0 | 16 | 196 |
| 85 | 10 | 10 | 10 | 1,010 a 1,263 | 1,000 a 1,065 | 1,010 a 1,185 | 4 a 75 | 0 | 19 | 1189 |
| 141 | 10 | 10 | 10 | 1,005 a 1,099 | 1,000 a 1,012 | 1,005 a 1,085 | 8 a 132 | 0 | 50 | 6097 |
| 33 + ties | 10 | 10 | 10 | 1,004 a 1,067 | 1,000 a 1,005 | 1,004 a 1,062 | 14 a 225 | 0 | 6 | 19 |
| 33 + shunts | 10 | 10 | 10 | 1,007 a 1,137 | 1,000 a 1,016 | 1,006 a 1,118 | 8 a 143 | 0 | 5 | 18 |
| CIGRE MV | 10 | 8 | 8 | 1,023 a 1,501 | 1,001 a 1,120 | 1,022 a 1,340 | 4 a 43 | 0 | 4 | 11 |
| Oberrhein | 10 | 10 | 10 | 1,009 a 1,212 | 1,000 a 1,047 | 1,008 a 1,158 | 4 a 83 | 0 | 91 | 13563 |
| European LV, trifásica | 10 | 10 | 10 | 1,005 a 1,089 | 1,000 a 1,018 | 1,004 a 1,070 | 5 a 10 | 0 | 1393 | 3297 |

*Script: `experiments/x12_verified.py verified` · CSV: `results/x12_verified.csv`.*

- **Krawczyk no es más ajustado que el encierro en ninguna cota de ninguna barra** (0 barras en los 78 casos en que ambos cierran).
- **Su intervalo es entre 1,003 y 1,340 veces el del encierro** (mediana por caso). En
  exceso sobre el rango exacto, que es lo que distingue a dos métodos ajustados, el de Krawczyk es entre
  4 y 225 veces el del encierro.
- **El encierro es más caro:** 81 ms de mediana frente a 12 ms, y la diferencia crece con el
  número de símbolos (26.5).
- **Qué es cada uno.** Krawczyk es una forma afín de primer orden más una cota del resto, que es de segundo
  orden. El encierro es lo mismo un orden más arriba: modelo de segundo orden más un resto de tercer orden. La
  ventaja medida es esa diferencia de orden y se paga con las formas cuadráticas.

### 26.2 Wang, Le Boudec y Paolone (2019): la subdivisión del conjunto

**El artículo no se pudo leer completo:** no tiene preprint y está tras un muro de pago (RELATED_WORK.md, cuarta
pasada). Por sus láminas, su resumen y la tesis de C. Wang, su cota por nodo es la que ya estaba implementada
(`cdhim/discs.py`), y su vía al ajuste es **subdividir el conjunto de incertidumbre** y empalmar los dominios.
Eso sí se puede medir. La caja se bisecta donde sobra hasta que ninguna hoja exceda los valores alcanzados en
más de 1 % del ancho, con tres encierros distintos en cada hoja: el verificado, Krawczyk y el disco por barra.

| Conjunto | Red | Casos | Sin subdividir / rango exacto: encierro | Krawczyk | disco | Encierros hasta 1 %: encierro | Krawczyk | no llega | disco | no llega |
|---|---|---:|---|---|---|---|---|---:|---|---:|
| Demanda en 8 zonas | 33 | 4 | 1,000 a 1,001 | 1,005 a 1,047 | 1,124 a 1,192 | 1 | 1 a 15 | 0 | 1917 a 2999 | 1 |
| Demanda en 8 zonas | 33 + ties | 4 | 1,000 | 1,002 a 1,022 | 1,083 a 1,113 | 1 | 1 a 9 | 0 | 59 a 65 | 0 |
| Cajas envolventes de los semidiscos | 33 | 10 | 1,000 a 1,010 | 1,006 a 1,135 | 1,727 a 1,964 | 1 | 1 a 2999 | 2 | 2999 | 10 |
| Cajas envolventes de los semidiscos | 33 + ties | 10 | 1,000 a 1,003 | 1,004 a 1,066 | 1,589 a 1,674 | 1 | 1 a 17 | 0 | 2999 | 10 |
| Cajas envolventes de los semidiscos | 33 + shunts | 10 | 1,000 a 1,010 | 1,006 a 1,134 | 1,744 a 1,974 | 1 | 1 a 2999 | 2 | 2999 | 10 |

*Script: `experiments/x12_verified.py subdivision` · CSV: `results/x12_subdivision.csv`. A lo más 3.000 encierros
por caso y método. Krawczyk y el disco no controlan el redondeo: sus columnas valen en aritmética exacta.*

- **Encierros necesarios para llegar a 1 %:**
  - el verificado: 1, es decir, ya llega sin subdividir;
  - Krawczyk: 1 a 2989 en los casos en que llega; en 4 de 38 no llega dentro del
    presupuesto y queda a 2,3 %;
  - el disco: 59 a 2761 en los casos en que llega; en 31 de 38
    no llega.
- La razón es el orden del exceso en el tamaño de la hoja: primero para el disco, segundo para Krawczyk y tercero
  para el encierro.
- Sin subdividir: el disco mide entre 1,08 y 1,97 veces el rango exacto,
  Krawczyk entre 1,002 y 1,135, y el encierro entre 1,0000 y
  1,0104.
- Qué no es esta comparación: no es el algoritmo del artículo, que decide admisibilidad con empalme de dominios y
  recentrado; es su cota por nodo con una subdivisión hecha aquí. Si el artículo trae una cota más fina, no la vi.

### 26.3 Redondeo dirigido: el encierro verificado

**Qué se demuestra.** Para los datos guardados de la red (matriz de admitancias, inyecciones del centro y matriz
de inyecciones), y para todo punto del conjunto de incertidumbre:

1. existe una solución del flujo de potencia a distancia a lo más z_j, en cada barra, del modelo de segundo orden;
2. es la única solución dentro de ese tubo;
3. el módulo de su voltaje está entre las cotas que se entregan.

**Cómo.** Es una verificación a posteriori (`cdhim/verified.py`). El centro, las sensibilidades de primer y
segundo orden y la inversa aproximada del jacobiano se calculan en punto flotante corriente y se toman como
datos; ninguno tiene que ser exacto. Sobre ellos se evalúa un operador de punto fijo en forma de inyección de
corriente y se comprueba que manda un producto de discos dentro de sí mismo (Brouwer) y que contrae (unicidad).

- El **defecto del modelo** es un polinomio de grado 4 en los símbolos. Sus coeficientes de grado 0, 1 y 2, que
  serían cero con datos exactos, se evalúan con cota de error y se conservan: ahí entra el residuo del flujo de
  carga y el de cada resolución lineal. Los de grado 3 y 4 se acotan con las normas de los cambios de corriente
  calculados, así que las admitancias grandes no aparecen.
- **Redondeo.** Tras cada operación elemental el resultado se mueve al número de punto flotante siguiente en la
  dirección segura. Las sumas y los productos de matrices se corrigen con la cota a priori del error de una suma
  de n términos (Higham, 2002), que vale para cualquier orden de evaluación. Las pruebas comparan estas cotas
  con aritmética racional exacta.
- **Supuestos.** Aritmética IEEE-754 de doble precisión con redondeo al más cercano, operaciones básicas y raíz
  cuadrada correctamente redondeadas, sin desbordamiento por abajo. Los datos de la red se toman como exactos.
  No se usa una biblioteca de intervalos certificada.

**La cota conjunta.** El módulo es |V*| + cᵀε + εᵀHε + resto. Acotar cᵀε y εᵀHε cada uno por su lado suma dos
máximos que no ocurren en el mismo punto, y esa diferencia es de segundo orden: era casi todo el exceso de la
cota superior. Ahora se usa que −(uᵀε)² ≤ a² − 2a·uᵀε para cualquier vector u y cualquier número a (una
tangente queda por encima de una función cóncava). Con u los autovectores de los autovalores negativos de H y a
su valor en el máximo de la parte lineal, la parte cóncava de la forma se vuelve lineal y se maximiza junto con
cᵀε, de forma exacta sobre los semidiscos. La cota inferior es la imagen especular con los autovalores positivos.
Ni los autovectores ni los valores a tienen que ser exactos. En cada barra y cada lado se conserva la mejor de
las dos cotas, así que nunca es peor que la separada.

| Red | Barras PQ | K | Casos | Verificado | Con unicidad | Ancho / rango exacto, a lo más: cota separada | con la cota conjunta | peor barra | Barras en contradicción con un punto certificado | Puntos sin certificar | Verificación, mediana (ms) |
|---|---:|---:|---:|---:|---:|---|---|---|---:|---:|---|
| 33 | 32 | 8 | 10 | 10 | 10 | 1,001 a 1,032 | 1,000 a 1,017 | 1,020 | 0 | 0 | 19 |
| 69 | 68 | 14 | 10 | 10 | 10 | 1,000 a 1,012 | 1,000 a 1,004 | 1,009 | 0 | 0 | 196 |
| 85 | 84 | 18 | 10 | 10 | 10 | 1,002 a 1,110 | 1,000 a 1,065 | 1,087 | 0 | 0 | 1189 |
| 141 | 140 | 28 | 10 | 10 | 10 | 1,002 a 1,039 | 1,000 a 1,012 | 1,014 | 0 | 0 | 6097 |
| 33 + ties | 32 | 8 | 10 | 10 | 10 | 1,001 a 1,014 | 1,000 a 1,005 | 1,006 | 0 | 0 | 19 |
| 33 + shunts | 32 | 8 | 10 | 10 | 10 | 1,001 a 1,031 | 1,000 a 1,016 | 1,020 | 0 | 0 | 18 |
| CIGRE MV | 14 | 4 | 10 | 8 | 8 | 1,004 a 1,160 | 1,001 a 1,120 | 1,127 | 0 | 0 | 11 |
| Oberrhein | 177 | 36 | 10 | 10 | 10 | 1,003 a 1,114 | 1,000 a 1,047 | 1,053 | 0 | 0 | 13563 |
| European LV, trifásica | 612 | 22 | 10 | 10 | 10 | 1,001 a 1,031 | 1,000 a 1,018 | 1,020 | 0 | 0 | 3297 |

*Script: `experiments/x12_verified.py verified` · CSV: `results/x12_verified.csv`. Los 80 casos de la sección 25
y los 10 de la red trifásica de 26.7.*

- **Cierra en 78 de 80 casos de la sección 25 y prueba unicidad en 78.** El encierro de
  esa sección cerraba en 78. No verifica: CIGRE MV, F, ρU = 80 % (el tubo alcanza voltaje cero); CIGRE MV, R, ρU = 80 % (el tubo alcanza voltaje cero).
- **El rigor no cuesta ancho.** Con la cota separada, el encierro verificado mide entre 0,889 y
  1,000 veces el de la sección 25: sale igual o más angosto porque la cota del resto es otra (el radio
  del tubo es entre 0,49 y 0,84 veces la cota del resto de esa sección).
- **La cota conjunta lo estrecha más.** En las redes de la sección 25 la cota superior sale de ella en
  todas las barras y la inferior en ninguna; en la red trifásica, las dos. El ancho sobre el rango exacto
  baja de 1,0004–1,160 a 1,0000–1,120 (mediana por caso).
- **El redondeo no pesa.** El radio del tubo, redondeo incluido, es a lo más 7,0 % del mayor ancho
  alcanzado del caso; el piso numérico se mide aparte en 26.5.
- **Costo.** La verificación tarda entre 5 ms y 14,9 s por caso, en una corrida a solas. En
  los casos lentos casi todo es la cota RAA-K de las formas cuadráticas, que resulta prescindible (26.5).

### 26.4 Referencia certificada

Hasta aquí el ancho se comparaba con valores alcanzados por búsqueda. Ahora hay dos referencias con garantía
(`cdhim/certified.py`).

**Cota interior certificada.** El mejor punto encontrado para cada barra se verifica: se comprueba en aritmética
racional que pertenece al conjunto y se encierra su solución del flujo de potencia con el mismo verificador,
sin incertidumbre. El extremo inferior de ese encierro es una cota inferior demostrada del máximo exacto, y lo
análogo para el mínimo. Entonces el cociente entre el ancho del encierro y el ancho entre esas dos cotas es
una **cota superior demostrada** de cuánto sobreestima el encierro.

- Puntos que no se pudieron certificar: 0. El valor certificado difiere del alcanzado en a lo más
  2,9·10⁻⁹ p.u.
- **Sobreestimación demostrada:** el encierro mide a lo más entre 1,0000 y 1,120 veces el rango
  exacto (mediana por caso; la peor barra, 1,127). Con ρU de hasta 20 %, a lo más 1,0109.
- **Dos demostraciones independientes no se contradicen:** en 0 barras el encierro deja fuera un
  punto certificado.

**Rango exacto.** En cajas con pocos símbolos se puede acotar el rango exacto por los dos lados: por dentro con
puntos certificados y por fuera con la unión de los encierros verificados de una partición de la caja, que se
refina por bisección donde sobra. El exceso de cada encierro es de tercer orden en el tamaño de la caja, así que
bastan pocas cajas.

| Conjunto | Red | Flujo | Tamaños: δ o ρU (%) | Encierros usados | Error del rango exacto, máx. (%) | Encierro / rango exacto, mediana, por tamaño | peor barra | Contradicciones |
|---|---|---|---|---|---|---|---|---:|
| Demanda en 8 zonas | 33 | – | 5 · 15 · 30 · 50 | 1 a 13 | 0,015 | 1,0000 · 1,0000 · 1,0002 · 1,0008 | 1,0010 | 0 |
| Demanda en 8 zonas | 33 + ties | – | 5 · 15 · 30 · 50 | 1 a 3 | 0,010 | 1,0000 · 1,0000 · 1,0000 · 1,0001 | 1,0002 | 0 |
| Cajas envolventes de los semidiscos | 33 | F | 5 · 10 · 20 · 50 · 80 | 1 a 69 | 0,015 | 1,0000 · 1,0001 · 1,0005 · 1,0034 · 1,0104 | 1,0132 | 0 |
| Cajas envolventes de los semidiscos | 33 | R | 5 · 10 · 20 · 50 · 80 | 3 a 67 | 0,012 | 1,0000 · 1,0001 · 1,0005 · 1,0031 · 1,0087 | 1,0112 | 0 |
| Cajas envolventes de los semidiscos | 33 + ties | F | 5 · 10 · 20 · 50 · 80 | 1 a 31 | 0,019 | 1,0000 · 1,0000 · 1,0002 · 1,0010 · 1,0029 | 1,0036 | 0 |
| Cajas envolventes de los semidiscos | 33 + ties | R | 5 · 10 · 20 · 50 · 80 | 1 a 31 | 0,018 | 1,0000 · 1,0000 · 1,0002 · 1,0010 · 1,0028 | 1,0034 | 0 |
| Cajas envolventes de los semidiscos | 33 + shunts | F | 5 · 10 · 20 · 50 · 80 | 1 a 69 | 0,014 | 1,0000 · 1,0001 · 1,0005 · 1,0033 · 1,0102 | 1,0128 | 0 |
| Cajas envolventes de los semidiscos | 33 + shunts | R | 5 · 10 · 20 · 50 · 80 | 3 a 67 | 0,012 | 1,0000 · 1,0001 · 1,0005 · 1,0030 · 1,0086 | 1,0109 | 0 |

*Script: `experiments/x12_verified.py exact` · CSV: `results/x12_exact_range.csv`. La bisección se detiene cuando
ninguna hoja excede los valores alcanzados en más de 10⁻⁴ del ancho, o tras 3.000 encierros.*

- En los 38 casos el rango exacto queda determinado con un error relativo de a lo más
  0,019 %, con entre 1 y 69 encierros y a lo más 3 s.
- **Demanda incierta en 8 zonas:** el encierro mide entre 1,0000 y 1,0008 veces el rango
  exacto.
- **Cajas envolventes de los semidiscos:** entre 1,0000 y 1,0104 veces.
- Límite: la bisección necesita cajas, así que no sirve tal cual para los semidiscos, y su costo crece con el
  número de símbolos. Para las redes grandes queda solo la cota interior certificada.

### 26.5 Escala

**Qué se midió** (`experiments/x13_scale.py`). El encierro verificado en redes PQ de distribución más grandes:
los casos de MATPOWER 7.1 de 118, 136 y 533 barras con sus enlaces cerrados, la de 141 barras, Oberrhein y la red
de baja tensión de Schutterwald de pandapower. Tres diseños de fuentes: una cada 10 barras, como en la sección
25 (K crece con la red), y 8 o 32 fuentes (K fijo), para separar el efecto del tamaño de la red del de K.
Flujos F y R con ρU de 10 y 50 %. La cota de los módulos se calcula de tres maneras: con la forma cuadrática
acotada término a término, con bloques de hasta 3 símbolos y con RAA-K completa.

| Red | Barras | Fuentes | K | Casos | Verificado / con unicidad | Modelo (s) | Tubo (s) | Módulos, término a término (s) | Módulos, RAA-K completa (s) | Ancho / rango alcanzado, peor caso: término a término | bloques de hasta 3 | RAA-K completa | Memoria pico (GB) |
|---|---:|---|---:|---:|---|---|---|---|---|---|---|---|---|
| 118 + ties | 118 | una cada 10 barras | 24 | 4 | 4 / 4 | 0,03 | 0,06 | 0,11 | 1,7 | 1,011 | 1,011 | 1,011 | 0,3 |
| 118 + ties | 118 | 8 | 16 | 4 | 4 / 4 | 0,03 | 0,04 | 0,08 | 0,5 | 1,011 | 1,011 | 1,011 | 0,3 |
| 118 + ties | 118 | 32 | 64 | 4 | 4 / 4 | 0,08 | 0,31 | 0,55 | 16,1 | 1,010 | 1,010 | 1,010 | 0,4 |
| 136 + ties | 136 | una cada 10 barras | 28 | 4 | 4 / 4 | 0,04 | 0,08 | 0,18 | 4,3 | 1,002 | 1,002 | 1,002 | 0,3 |
| 136 + ties | 136 | 8 | 16 | 4 | 4 / 4 | 0,03 | 0,09 | 0,11 | 1,5 | 1,002 | 1,002 | 1,002 | 0,3 |
| 136 + ties | 136 | 32 | 64 | 4 | 4 / 4 | 0,12 | 0,61 | 0,81 | 41,2 | 1,002 | 1,002 | 1,002 | 0,4 |
| 141 | 141 | una cada 10 barras | 28 | 4 | 4 / 4 | 0,05 | 0,12 | 0,23 | 7,5 | 1,010 | 1,010 | 1,010 | 0,3 |
| 141 | 141 | 8 | 16 | 4 | 4 / 4 | 0,03 | 0,07 | 0,11 | 1,5 | 1,011 | 1,011 | 1,011 | 0,3 |
| 141 | 141 | 32 | 64 | 4 | 4 / 4 | 0,09 | 0,37 | 0,64 | 46,3 | 1,010 | 1,010 | 1,010 | 0,4 |
| Oberrhein | 179 | una cada 10 barras | 36 | 4 | 4 / 4 | 0,06 | 0,20 | 0,30 | 13,1 | 1,021 | 1,021 | 1,021 | 0,3 |
| Oberrhein | 179 | 8 | 16 | 4 | 4 / 4 | 0,04 | 0,08 | 0,10 | 1,5 | 1,022 | 1,022 | 1,022 | 0,3 |
| Oberrhein | 179 | 32 | 64 | 4 | 4 / 4 | 0,12 | 0,48 | 0,81 | 49,2 | 1,022 | 1,022 | 1,022 | 0,4 |
| 533 | 533 | una cada 10 barras | 108 | 4 | 4 / 4 | 0,63 | 8,26 | 8,06 | 97,6 | 1,012 | 1,012 | 1,012 | 1,3 |
| 533 | 533 | 8 | 16 | 4 | 4 / 4 | 0,09 | 1,06 | 0,26 | 0,4 | 1,017 | 1,017 | 1,017 | 0,4 |
| 533 | 533 | 32 | 64 | 4 | 4 / 4 | 0,28 | 3,37 | 2,58 | 28,4 | 1,010 | 1,010 | 1,010 | 0,6 |
| 533 + ties | 533 | una cada 10 barras | 108 | 4 | 4 / 4 | 0,65 | 4,16 | 7,67 | 146,6 | 1,011 | 1,011 | 1,011 | 1,3 |
| 533 + ties | 533 | 8 | 16 | 4 | 4 / 4 | 0,10 | 0,55 | 0,23 | 1,1 | 1,015 | 1,015 | 1,015 | 0,4 |
| 533 + ties | 533 | 32 | 64 | 4 | 4 / 4 | 0,29 | 1,67 | 2,35 | 57,7 | 1,012 | 1,012 | 1,012 | 0,6 |
| Schutterwald | 2940 | una cada 10 barras | 586 | 4 | no se intentó | – | – | – | – | – | – | – | – |
| Schutterwald | 2940 | 8 | 16 | 4 | 4 / 4 | 1,89 | 18,08 | 0,98 | 2,0 | 1,063 | 1,063 | 1,063 | 3,8 |
| Schutterwald | 2940 | 32 | 64 | 4 | 4 / 4 | 2,97 | 26,16 | 11,79 | 198,9 | 1,053 | 1,053 | 1,053 | 5,0 |

*Script: `experiments/x13_scale.py` · CSV: `results/x13_scale.csv`. Tiempos de una corrida a solas. El rango
alcanzado es el de flujos de frontera y muestras, sin búsqueda local, para que el tiempo no se lo lleve la
referencia.*

- **Verifica, con unicidad, en 80 de los 80 casos corridos:** hasta 2.940 barras con
  64 símbolos y hasta 533 barras con 108. Ninguna violación
  (0).
- **Lo que limita es el número de símbolos, no el de barras.** El modelo de segundo orden guarda un tensor de
  barras × K × K. Con una fuente cada 10 barras en la red de 2.940 barras son 586 símbolos
  y un solo tensor ocupa 15 GB: no se intentó (4 casos).
- **Costo en la red más grande:** verificar el tubo toma de 18 a 28 s, por las
  matrices densas de 2m × 2m, y la memoria pico es de 5,0 GB.
- **RAA-K es lo más caro y, con la cota conjunta, aquí no aporta nada.** La versión completa llega a
  216 s donde la cota término a término toma 12,7 s y la de bloques de hasta 3,
  17,6 s. El ancho que entregan, sobre el de la cota término a término: 1,0000 la
  completa y 1,0000 la de bloques de hasta 3, en los 80 casos. Una explicación coherente con
  esto, que no se midió por separado: RAA-K ayudaba en la cota superior, que ahora sale de la cota conjunta, y
  para una forma cóncava su cota inferior coincide con la de término a término.
- **Sin RAA-K el encierro es mucho más barato.** En Oberrhein con 36 símbolos, modelo, tubo y módulos
  suman 0,6 s con la cota término a término y 13,3 s con RAA-K completa.
- **Piso numérico.** Sin incertidumbre, el radio del tubo es a lo más 6,1·10⁻⁹ p.u., que es a lo más
  0,2 % del mayor radio del caso. El redondeo no limita el tamaño de la red.
- **Redes que hubo que arreglar para cargarlas:** el caso de 533 barras no está en MATPOWER 7.1 y su archivo de la
  versión 8.0 trae expresiones que el lector no evalúa (se usa una copia con las expresiones evaluadas); el de 136
  barras solo se lee en modo UTF-8; y la red de Schutterwald exige ordenar la tabla de barras. Tal como se
  publican, las redes de pandapower con interruptores abiertos siguen sin cargar.

### 26.6 Barras de voltaje controlado

**Qué se hizo** (`cdhim/verified_vi.py`). Una segunda forma del operador. Las ecuaciones se escriben con la
corriente inyectada en cada barra: conj(V_j)·I_j = conj(S_j) en las barras PQ, y Re(conj(V_j)·I_j) = P_j junto con
|V_j|² = v_j² en las PV. La no linealidad es un producto del voltaje y la corriente de la misma barra, así que se
acotan a la vez el resto de voltaje y el de corriente, y ninguna admitancia multiplica una cota del resto. Misma
disciplina de redondeo que en 26.3; la cota de los módulos es la misma.

**Distribución con fuentes de voltaje controlado.** Las 4 barras con generación de la regla de la
sección 25 pasan a ser PV, con su potencia activa incierta (4 símbolos):

| Red | Flujo | Barras PV | K | Verificado / con unicidad | Radio de la parte lineal | Encierro / rango exacto en vértices, mediana, con ρU = 5 · 10 · 20 · 50 · 80 % | peor barra | Violaciones |
|---|---|---:|---:|---|---|---|---|---:|
| 33 | F | 4 | 4 | 5 / 5 de 5 | 0,02 a 0,27 | 1,002 · 1,007 · 1,029 · 1,237 · 1,993 | 1,01 · 1,03 · 1,13 · 2,25 · 7,82 | 0 |
| 33 | R | 4 | 4 | 5 / 5 de 5 | 0,02 a 0,26 | 1,002 · 1,007 · 1,027 · 1,216 · 1,895 | 1,01 · 1,03 · 1,13 · 2,22 · 7,43 | 0 |
| 33 + ties | F | 4 | 4 | 5 / 5 de 5 | 0,01 a 0,13 | 1,001 · 1,003 · 1,013 · 1,094 · 1,392 | 1,03 · 1,10 · 1,32 · 2,97 · 7,66 | 0 |
| 33 + ties | R | 4 | 4 | 5 / 5 de 5 | 0,01 a 0,13 | 1,001 · 1,003 · 1,012 · 1,092 · 1,378 | 1,03 · 1,10 · 1,31 · 2,96 · 7,53 | 0 |
| 33 + shunts | F | 4 | 4 | 5 / 5 de 5 | 0,02 a 0,28 | 1,002 · 1,007 · 1,026 · 1,210 · 1,824 | 1,01 · 1,03 · 1,15 · 2,42 · 9,64 | 0 |
| 33 + shunts | R | 4 | 4 | 5 / 5 de 5 | 0,02 a 0,27 | 1,001 · 1,006 · 1,024 · 1,192 · 1,741 | 1,01 · 1,03 · 1,15 · 2,39 · 9,16 | 0 |

*Script: `experiments/x14_pv.py` · CSV: `results/x14_pv_distribution.csv`.*

- **Verifica, con unicidad, en 30 de 30 casos**, sin violaciones (0).
- **Es ajustado con incertidumbre moderada y deja de serlo con mucha.** Con ρU de hasta 20 % mide entre
  1,001 y 1,029 veces el rango exacto en vértices; con 50 y 80 %, entre 1,09 y
  1,99, y hasta 9,6 veces en la peor barra. Es peor que con barras PQ a igual ρU. No se
  aisló cuánto se debe a que el rango de las barras vecinas a una PV es muy chico y cuánto a la cota.

**Transmisión**, con demanda incierta en 8 zonas:

| Red | Barras PV | Barras PQ | δ (%) | Verificado | Radio de la parte lineal | Radio del mayorante de la forma V-I | Encierro / rango exacto en vértices, mediana | peor barra |
|---|---:|---:|---:|---|---:|---:|---:|---:|
| case30 | 5 | 24 | 2 | sí (única) | 0,016 | 0,05 | 1,002 | 1,01 |
| case30 | 5 | 24 | 5 | sí (única) | 0,040 | 0,14 | 1,011 | 1,05 |
| case30 | 5 | 24 | 10 | sí (única) | 0,081 | 0,27 | 1,047 | 1,25 |
| case30 | 5 | 24 | 20 | sí (única) | 0,161 | 0,56 | 1,286 | 2,64 |
| case118 | 53 | 64 | 2 | no | 0,086 | 3,96 | – | – |
| case118 | 53 | 64 | 5 | no | 0,215 | 10,09 | – | – |
| case118 | 53 | 64 | 10 | no | 0,430 | 20,86 | – | – |
| case118 | 53 | 64 | 20 | no | 0,861 | 44,41 | – | – |
| case300 | 68 | 231 | 2 | no | 1,417 | 148,20 | – | – |
| case300 | 68 | 231 | 5 | no | 3,543 | 479,82 | – | – |
| case300 | 68 | 231 | 10 | no | 7,085 | 1324,04 | – | – |
| case300 | 68 | 231 | 20 | no | 14,170 | 4105,65 | – | – |

*Script: `experiments/x14_pv.py` · CSV: `results/x14_pv_transmission.csv`. «Radio de la parte lineal»: radio
espectral de la suma de |J⁻¹ΔJ| de la sección 22; por encima de 1, ningún punto fijo con mayorante por entradas
puede cerrar.*

- **Verifica en case30** (4 de 12 casos), entre 1,002 y 1,29
  veces el rango exacto en vértices. **No verifica en case118 y case300, ni con δ = 2 %.**
- **En la red de 118 barras el límite es de esta cota, no del problema.** La parte lineal contrae (radio de
  0,09 a 0,86), pero el mayorante de la forma voltaje-corriente, que suma normas por bloque
  sobre todas las barras, tiene radio de 4,0 para arriba. Una cota que conserve las cancelaciones de
  la parte lineal podría cerrar ahí; no se implementó.
- **En la de 300 barras ni la parte lineal contrae** (radio desde 1,42 con δ = 2 %). Con esta familia
  de cotas no hay encierro verificado a esa escala.

**Las dos formas en las redes PQ de la sección 25:**

| Red | Casos | Rama más corta, \|z\| (p.u.) | Verifica la forma de corriente | Verifica la forma V-I | Radio del mayorante V-I, máx. | ídem sin los términos de redondeo | Radio del tubo V-I / corriente, mediana |
|---|---:|---|---:|---:|---|---|---|
| 33 | 10 | 0,0040 | 10 | 10 | 0,21 | 0,21 | 1,00 |
| 69 | 10 | 4,1·10⁻⁵ | 10 | 10 | 0,12 | 0,12 | 1,00 |
| 85 | 10 | 4,1·10⁻⁴ | 10 | 10 | 0,37 | 0,37 | 0,99 a 1,00 |
| 141 | 10 | 6,4·10⁻⁷ | 10 | 6 | 0,94 | 0,16 | 1,00 a 1,09 |
| 33 + ties | 10 | 0,0040 | 10 | 10 | 0,11 | 0,11 | 1,00 |
| 33 + shunts | 10 | 0,0040 | 10 | 10 | 0,21 | 0,21 | 1,00 |
| CIGRE MV | 10 | 4,0·10⁻⁴ | 8 | 8 | 0,88 | 0,88 | 0,93 a 1,00 |
| Oberrhein | 10 | 2,3·10⁻⁵ | 10 | 8 | 0,39 | 0,39 | 1,00 |

*Script: `experiments/x14_pv.py` · CSV: `results/x14_forms.csv`.*

- La forma de corriente verifica en 78 de 80 casos y la de voltaje-corriente en
  72. Esta falla donde la otra no en: 141 (4), Oberrhein (2). Al revés, en 0 casos.
- En la red de 141 barras la causa está aislada: el bloque que acota Y·(I − C·J) crece con las admitancias de las
  ramas muy cortas, y sin esos términos de redondeo el radio del mayorante baja a la sexta parte. Haría falta
  evaluar ese residuo con más precisión.
- Donde verifican las dos, el radio del tubo es el mismo (razón de 0,93 a 1,09).
- Consecuencia: en redes PQ se usa la forma de corriente; la de voltaje-corriente queda para cuando hay barras PV.

### 26.7 Redes desbalanceadas

**Observación.** En el dominio de fases, con el neutro a potencial de tierra en cada barra y cargas de potencia
constante entre fase y neutro, cada nodo-fase cumple (YV)ₙ + cₙ = conj(Sₙ/Vₙ): la misma forma que una barra
monofásica. El encierro se aplica entonces sin cambios, con la matriz de admitancias de fases (bloques 3 × 3 con
acoplamiento mutuo) y la fuente en 1∠0°, 1∠−120° y 1∠120°. Se verificó en vez de darla por buena
(`cdhim/threephase.py`, `tests/test_threephase.py`).

**Red.** El lado de baja tensión del IEEE European LV Test Feeder, tal como viene en pandapower: 906 barras a
416 V y 55 cargas monofásicas en estrella. Se usa completa (2715 nodos-fase PQ) y reducida de forma
exacta a 205 barras (612 nodos-fase). Los alimentadores IEEE 13, 37 y 123 se descartaron:
traen reguladores, transformadores y cargas en delta o de corriente constante.

| Instantánea de carga | Modelo | Barras | Nodos-fase | Carga (kVA) | Diferencia máx. con pandapower (p.u.) | ídem, descontado el residuo de pandapower | Diferencia entre fases, máx. (p.u.) | Secuencia negativa / positiva, máx. (%) |
|---|---|---:|---:|---:|---|---|---:|---:|
| on_peak_566 | full | 906 | 2718 | 57,6 | 4,4·10⁻¹⁰ | 2,3·10⁻¹¹ | 0,070 | 0,69 |
| on_peak_566 | merged | 205 | 615 | 57,6 | 4,3·10⁻¹⁰ | 4,0·10⁻¹¹ | 0,070 | 0,69 |
| off_peak_1 | full | 906 | 2718 | 2,8 | 4,7·10⁻⁹ | 9,0·10⁻¹² | 0,001 | 0,01 |
| off_peak_1 | merged | 205 | 615 | 2,8 | 4,7·10⁻⁹ | 1,7·10⁻¹² | 0,001 | 0,01 |
| off_peak_1440 | full | 906 | 2718 | 9,8 | 3,3·10⁻¹¹ | 1,7·10⁻¹¹ | 0,004 | 0,05 |
| off_peak_1440 | merged | 205 | 615 | 9,8 | 2,0·10⁻¹¹ | 2,6·10⁻¹² | 0,004 | 0,05 |

*Script: `experiments/x10_three_phase.py validation` · CSV: `results/x10_model_validation.csv`.*

- **El modelo coincide con `runpp_3ph` de pandapower** a 4,7·10⁻⁹ p.u., y a 4,0·10⁻¹¹ p.u. si se
  descuenta el residuo de convergencia de pandapower.
- **Es desbalanceada de verdad:** hasta 0,070 p.u. de diferencia de módulo entre fases, secuencia negativa
  de hasta 0,69 % y secuencia cero de hasta 3,38 %.
- Dejar fuera el transformador Dyn cambia los módulos en a lo más 7,8·10⁻⁴ p.u.

**Encierro sin redondeo dirigido** (el de la sección 25), 11 inversores monofásicos repartidos
5/5/1 entre las fases, 22 símbolos:

| Flujo | ρU (%) | Rango alcanzado, máx. (p.u.) | Cierra | Encierro / rango, mediana | máx. | Violaciones | Discos: cierra | Discos por nodo / rango, mediana | Nodos-fase con el mínimo fuera de (0, −C) |
|---|---:|---:|---|---:|---:|---:|---|---:|---:|
| F | 5 | 0,0039 | sí | 1,0014 | 1,0015 | 0 | sí | 1,61 | 612 |
| F | 10 | 0,0079 | sí | 1,0028 | 1,0031 | 0 | sí | 1,62 | 612 |
| F | 20 | 0,0158 | sí | 1,0060 | 1,0066 | 0 | sí | 1,63 | 612 |
| F | 50 | 0,0395 | sí | 1,0185 | 1,0202 | 0 | sí | 1,65 | 612 |
| F | 80 | 0,0634 | sí | 1,0353 | 1,0398 | 0 | sí | 1,69 | 612 |
| R | 5 | 0,0038 | sí | 1,0013 | 1,0014 | 0 | sí | 1,60 | 612 |
| R | 10 | 0,0075 | sí | 1,0026 | 1,0029 | 0 | sí | 1,60 | 612 |
| R | 20 | 0,0151 | sí | 1,0057 | 1,0061 | 0 | sí | 1,59 | 612 |
| R | 50 | 0,0380 | sí | 1,0172 | 1,0187 | 0 | sí | 1,61 | 612 |
| R | 80 | 0,0612 | sí | 1,0328 | 1,0368 | 0 | sí | 1,62 | 612 |

*Script: `experiments/x10_three_phase.py cases` · CSV: `results/x10_three_phase.csv`.*

- **Cierra en 10 de 10 casos** y queda entre 1,001 y 1,035 veces el rango alcanzado;
  los discos por nodo, entre 1,59 y 1,69. Ninguna violación (0).
- **Con redondeo dirigido y cota conjunta** (`results/x12_verified.csv`): verifica en 10 de
  10 casos, con unicidad en 10, y mide a lo más entre 1,0005 y 1,018
  veces el rango exacto. Krawczyk, entre 1,005 y 1,089.
- **La monotonía no se cumple nunca:** en los 10 casos el mínimo de todos los nodos-fase está fuera de
  (0, −C), a una distancia de 50 % del rango. Una fuente en una fase mueve
  las otras dos en sentido contrario, así que un mínimo tomado de un solo flujo no sirve en trifásico.
- **Tamaño completo** (`results/x10_full_model.csv`): con 2715 nodos-fase cierra en 10
  de 10 casos en 8,4 a 11,1 s, con las mismas cotas que la red reducida
  (diferencia máxima 1,7·10⁻¹¹ p.u.) y 0 violaciones en las muestras.

**Hasta dónde llega** (capacidad total multiplicada, flujo inverso, ρU = 80 %):

| Capacidad total / carga | \|V\| alcanzado (p.u.) | Encierro: cierra | Encierro / rango, mediana | Discos: cierra | Discos por nodo / rango, mediana | Muestras sin solución AC |
|---:|---|---|---:|---|---:|---:|
| 1 | 1,012 a 1,090 | sí | 1,033 | sí | 1,62 | 0 |
| 2 | 1,000 a 1,133 | sí | 1,104 | sí | 1,75 | 0 |
| 3 | 0,977 a 1,175 | sí | 1,251 | sí | 1,92 | 0 |
| 4 | 0,953 a 1,216 | sí | 1,634 | sí | 2,18 | 0 |
| 6 | 0,906 a 1,301 | no | – | no | – | 0 |
| 8 | 0,862 a 1,385 | no | – | no | – | 0 |

*Script: `experiments/x10_three_phase.py limits` · CSV: `results/x10_limits.csv`.*

- El encierro y los discos dejan de cerrar en el mismo escalón (capacidad 6 veces la carga),
  aunque el flujo AC sigue teniendo solución en todas las muestras. Con 4 veces, los voltajes ya
  están fuera de cualquier rango de operación.

**Qué queda fuera:** cargas en delta, transformadores, la impedancia de cortocircuito de la fuente y el neutro
explícito. Las líneas no transpuestas y las fases faltantes están probadas solo en una red sintética, sin
contraste con un solver externo.

### 26.8 Lectura crítica: qué cambia

1. **El aporte resistió la comparación que podía tumbarlo, y quedó mejor delimitado.** Frente a Krawczyk por
   inyección de corriente, que es el mejor baseline intervalar que se pudo armar y más fuerte que lo publicado,
   el encierro es más ajustado en todas las barras: su exceso sobre el rango exacto es entre 4 y
   225 veces menor. Lo que ya no se puede decir es que la garantía publicada queda a 1,24–2,17 veces
   el rango: esa es la de los discos, y Krawczyk queda a 1,003–1,501.
2. **Qué es el método, dicho sin adorno:** Krawczyk un orden más arriba. Modelo de segundo orden, resto de
   tercer orden acotado por un punto fijo en forma de inyección de corriente, y una cota conjunta de los
   términos de primer y segundo orden. Un revisor puede leerlo como incremental; el argumento es que ese orden
   de más cambia el exceso en un factor de 4 a 225 y la subdivisión necesaria de
   decenas de cajas a una.
3. **Ahora es un encierro verificado de verdad:** existencia y unicidad con el redondeo controlado, y la
   sobreestimación demostrada contra puntos certificados, no medida contra una búsqueda. En la revisión de
   literatura no apareció nada así a tamaño realista; es un resultado negativo de búsqueda, no una prueba.
4. **Alcance comprobado:** redes PQ radiales, malladas y con elementos en derivación; trifásicas desbalanceadas
   con cargas en estrella; fuentes de voltaje controlado en distribución, con menos ajuste; hasta
   2.940 barras con 64 símbolos.
5. **Límites comprobados:** no verifica en transmisión de 118 barras o más; el número de símbolos está acotado
   por la memoria (cientos, no miles); no cierra en CIGRE MV con ρU = 80 %; es más caro que Krawczyk.
6. **Lo que sigue sin hacerse:**
   - leer el artículo de Wang, Le Boudec y Paolone (2019): falta el PDF;
   - cargas en delta, transformadores y neutro explícito;
   - una segunda implementación del redondeo con una biblioteca de intervalos certificada, que es lo que
     pediría un revisor de computación verificada;
   - revisar si la cota conjunta por tangentes ya existe con otro nombre: es una idea corriente en optimización
     global y no hay que presentarla como nueva sin esa revisión;
   - rehacer las secciones 18 a 23 con RAA-K normalizada (punto 9.14), que sigue pendiente.
7. **RAA-K deja de aportar en este encierro.** Con la cota conjunta, en los 80 casos de la escala el ancho
   es el mismo con RAA-K completa que acotando la forma término a término (razón 1,0000), a un costo de
   hasta 216 s. Quitarlo abarata el encierro (en Oberrhein, de 13,3 a 0,6 s) y
   acerca su costo al de Krawczyk. Para un artículo cuyo interés es la aritmética afín del grupo, eso hay que
   decirlo tal cual: hoy el resultado no se sostiene en el operador. Lo que sí queda de esa línea es la estructura: la cota
   conjunta es una descomposición en bloques (uᵀε)², como la del teorema k-nomial, acotados junto con el término
   lineal en vez de uno por uno. Usar ahí los bloques de RAA-K en lugar de autovectores no se probó (punto 9.23).

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

| Sistema | Cota | Por | Violación máx. contra AC (%) | Sobre 0,27 % | Máx. contra el modelo lineal (%) | Máx. contra el sustituto de 2.º orden (%) | Ancho (área) / sin corregir, mediana | máximo |
|---|---|---|---:|---:|---:|---:|---:|---:|
| IEEE-33 | Intervalo CD-HIM (λ_max), sin corregir | barra | 0,020 | 0 | 0,020 | 0,020 | 1,000 | 1,00 |
| IEEE-33 |   + media y varianza de q | barra | 0,020 | 0 | 0,020 | 0,020 | 1,000 | 1,01 |
| IEEE-33 |   + rango de q (Caja: q en el vértice) | barra | 0,020 | 0 | 0,020 | 0,020 | 1,003 | 1,01 |
| IEEE-33 | Intervalo CD-HIM (proyección exacta), sin corregir | barra | 0,050 | 0 | 0,050 | 0,050 | 1,000 | 1,00 |
| IEEE-33 |   + rango de q | barra | 0,050 | 0 | 0,050 | 0,050 | 1,003 | 1,01 |
| IEEE-33 | Caja 3σ, sin corregir | barra | 0,275 | 1 | 0,275 | 0,275 | 1,000 | 1,00 |
| IEEE-33 |   + rango de q | barra | 0,270 | 0 | 0,275 | 0,270 | 1,004 | 1,00 |
| IEEE-33 | Caja afín, sin corregir | barra | 0,000 | 0 | 0,000 | 0,000 | 1,000 | 1,00 |
| IEEE-33 |   + cota rigurosa de q (término a término) | barra | 0,000 | 0 | 0,000 | 0,000 | 1,021 | 1,02 |
| IEEE-33 |   + q en el vértice extremo | barra | 0,000 | 0 | 0,000 | 0,000 | 1,007 | 1,01 |
| IEEE-33 | Elipse, sin corregir | rama | 0,165 | 0 | 0,165 | 0,165 | 1,000 | 1,00 |
| IEEE-33 |   + media de q | rama | 0,165 | 0 | 0,165 | 0,165 | 1,000 | 1,00 |
| IEEE-33 |   + media y covarianza de q | rama | 0,165 | 0 | 0,165 | 0,165 | 1,000 | 1,00 |
| IEEE-33 |   + rango de q en los ejes principales | rama | 0,160 | 0 | 0,165 | 0,160 | 1,008 | 1,01 |
| PEGASE-2869 | Intervalo CD-HIM (λ_max), sin corregir | barra | 0,410 | 1 | 0,065 | 0,410 | 1,000 | 1,00 |
| PEGASE-2869 |   + media y varianza de q | barra | 0,090 | 0 | 0,065 | 0,090 | 1,000 | 3,79 |
| PEGASE-2869 |   + rango de q (Caja: q en el vértice) | barra | 0,050 | 0 | 0,055 | 0,050 | 1,008 | 3,79 |
| PEGASE-2869 | Intervalo CD-HIM (proyección exacta), sin corregir | barra | 0,630 | 2 | 0,085 | 0,630 | 1,000 | 1,00 |
| PEGASE-2869 |   + rango de q | barra | 0,080 | 0 | 0,080 | 0,080 | 1,010 | 3,79 |
| PEGASE-2869 | Caja 3σ, sin corregir | barra | 1,190 | 80 | 0,315 | 1,205 | 1,000 | 1,00 |
| PEGASE-2869 |   + rango de q | barra | 0,290 | 3 | 0,290 | 0,290 | 1,009 | 1,45 |
| PEGASE-2869 | Caja afín, sin corregir | barra | 0,410 | 1 | 0,000 | 0,410 | 1,000 | 1,00 |
| PEGASE-2869 |   + cota rigurosa de q (término a término) | barra | 0,000 | 0 | 0,000 | 0,000 | 172,934 | 2332451,91 |
| PEGASE-2869 |   + q en el vértice extremo | barra | 0,000 | 0 | 0,000 | 0,000 | 1,044 | 3,79 |
| PEGASE-2869 | Elipse, sin corregir | rama | 2,685 | 22 | 0,270 | 2,700 | 1,000 | 1,00 |
| PEGASE-2869 |   + media de q | rama | 1,635 | 16 | 0,275 | 1,625 | 1,000 | 1,00 |
| PEGASE-2869 |   + media y covarianza de q | rama | 1,080 | 13 | 0,275 | 1,085 | 1,000 | 1,18 |
| PEGASE-2869 |   + rango de q en los ejes principales | rama | 0,235 | 0 | 0,235 | 0,235 | 1,024 | 3,07 |

*Script: `experiments/b4_second_order.py` · CSV: `results/b4_design_bus.csv`, `results/b4_design_branch.csv` ·
4.000 muestras del término cuadrático (semilla de B1 + 1)*

- **El intervalo CD-HIM por barra queda resuelto:** en PEGASE la peor barra pasa de 0,41 % a
  0,05 %, con un ancho 0,8 % mayor en mediana.
- **La elipse necesita el rango, no los momentos.** Con media y covarianza de q la peor elipse baja de
  2,69 % (22 ramas sobre 0,27 %) a 1,08 % (13 ramas). A lo largo del eje menor de
  las elipses delgadas el término cuadrático es casi de un solo signo y de cola pesada (asimetría cercana a ±2,7 y
  curtosis cercana a 13 en las ramas peores; diagnóstico en esta sesión, sin CSV), lo que un ajuste de dos momentos no representa. Con el rango
  proyectado la peor elipse queda en 0,235 %, con un área 2,4 % mayor en mediana y hasta
  3,1 veces en la peor rama.
- **El tercer orden no pesa en el punto de diseño:** el sustituto de segundo orden reproduce las violaciones AC.
- **La cota rigurosa de q término a término no sirve:** es válida para el modelo de segundo orden pero multiplica
  el ancho de la caja afín por más de cien en mediana en PEGASE.

### 27.2 Barrido carga × δ (las 50 celdas de C1, mismas muestras)

| Sistema | Cota | Celdas | Celdas con alguna barra o rama sobre 0,27 % | Máx. barras o ramas sobre 0,27 % en una celda | Peor violación (%) | Ancho / sin corregir, mediana de las celdas | máximo |
|---|---|---:|---:|---:|---:|---:|---:|
| IEEE-33 | Intervalo CD-HIM, sin corregir | 25 | 1 | 1 | 0,471 | 1,000 | 1,000 |
| IEEE-33 |   + media y varianza de q | 25 | 1 | 1 | 0,386 | 1,000 | 1,030 |
| IEEE-33 |   + rango de q | 25 | 0 | 0 | 0,040 | 1,004 | 1,285 |
| IEEE-33 | Caja 3σ, sin corregir | 25 | 4 | 31 | 1,314 | 1,000 | 1,000 |
| IEEE-33 |   + media y varianza de q | 25 | 3 | 31 | 1,048 | 1,000 | 1,028 |
| IEEE-33 | Caja afín, sin corregir | 25 | 0 | 0 | 0,000 | 1,000 | 1,000 |
| IEEE-33 |   + q en el vértice extremo | 25 | 0 | 0 | 0,000 | 1,009 | 1,675 |
| IEEE-33 | Elipse, sin corregir | 25 | 4 | 15 | 1,088 | – | – |
| IEEE-33 |   + media y covarianza de q | 25 | 3 | 15 | 0,933 | – | – |
| IEEE-33 |   + rango de q en los ejes principales | 25 | 0 | 0 | 0,130 | – | – |
| PEGASE-2869 | Intervalo CD-HIM, sin corregir | 25 | 24 | 88 | 5,153 | 1,000 | 1,000 |
| PEGASE-2869 |   + media y varianza de q | 25 | 15 | 43 | 1,102 | 1,000 | 1,005 |
| PEGASE-2869 |   + rango de q | 25 | 0 | 0 | 0,110 | 1,028 | 1,097 |
| PEGASE-2869 | Caja 3σ, sin corregir | 25 | 25 | 1246 | 23,838 | 1,000 | 1,000 |
| PEGASE-2869 |   + media y varianza de q | 25 | 25 | 922 | 2,420 | 1,000 | 1,003 |
| PEGASE-2869 | Caja afín, sin corregir | 25 | 24 | 3 | 5,153 | 1,000 | 1,000 |
| PEGASE-2869 |   + q en el vértice extremo | 25 | 0 | 0 | 0,006 | 1,148 | 1,420 |
| PEGASE-2869 | Elipse, sin corregir | 25 | 25 | 1219 | 32,325 | – | – |
| PEGASE-2869 |   + media y covarianza de q | 25 | 25 | 837 | 2,850 | – | – |
| PEGASE-2869 |   + rango de q en los ejes principales | 25 | 0 | 0 | 0,260 | – | – |

*Script: `experiments/b4_second_order.py` · CSV: `results/b4_sweep.csv`. El ancho de la elipse no está en el
barrido.*

- Con el rango de q, ninguna de las 50 celdas tiene una barra o una elipse sobre el nominal. Con solo
  media y varianza siguen fallando 15 de las 25 celdas de PEGASE.
- La caja 3σ corregida sigue fallando: lo que falla ahí es la normalidad, no la linealización.
- "Sobre 0,27 %" compara una frecuencia estimada con 20.000 muestras con el nivel nominal: una cota exacta al
  nivel nominal supera el umbral por azar en una parte de los objetivos. La elipse corregida queda bajo el umbral
  porque es algo más ancha.

### 27.3 Costo

| Sistema | Método lineal (s) | Modelo de 2.º orden (s) | Media exacta de q (s) | q en los vértices de las barras Caja (s) | 1.000 muestras de q (s) | 4.000 muestras de q (s) |
|---|---:|---:|---:|---:|---:|---:|
| IEEE-33 | 0,002 | 0,006 | 0,000 | 0,000 | 0,009 | 0,005 |
| PEGASE-2869 | 1,161 | 3,910 | 0,275 | 1,236 | 1,069 | 3,568 |

*Script: `experiments/b4_second_order.py` (modo `timing`) · CSV: `results/b4_timing.csv` · un proceso*

En PEGASE la corrección completa suma 9,0 s a los 1,2 s del método lineal. La sección 28.5 muestra
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

| Sistema | Objetivo | Cota | Violación máx. contra AC (%) | Sobre 0,27 % | Muestras con alguna fuera (%) | Máx. contra el modelo lineal (%) | Ancho / caja afín, mediana | Ancho / 3σ, mediana | Ancho mediano (p.u.) |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| IEEE-33 | barra | Caja afín (peor caso, 1.er orden) | 0,000 | 0 | 0,00 | 0,000 | 1,00 | 2,92 | 0,02005 |
| IEEE-33 | barra |   + q en el vértice extremo | 0,000 | 0 | 0,00 | 0,000 | 1,01 | 2,94 | 0,02018 |
| IEEE-33 | barra |   + extremo del 2.º orden hallado por búsqueda | 0,000 | 0 | 0,00 | 0,000 | 1,01 | 2,94 | 0,02018 |
| IEEE-33 | barra | Bernstein + rango de q | 0,000 | 0 | 0,00 | 0,000 | 0,70 | 1,97 | 0,01375 |
| IEEE-33 | barra | 3σ + rango de q (supone normalidad) | 0,270 | 0 | 0,70 | 0,275 | 0,34 | 1,00 | 0,00659 |
| IEEE-33 | barra | Intervalo CD-HIM + rango de q (B4) | 0,020 | 0 | 0,02 | 0,020 | 0,56 | 1,56 | 0,01030 |
| IEEE-33 | caída | Forma afín de la caída, peor caso | 0,000 | 0 | 0,00 | 0,000 | 1,00 | 1,49 | 0,00108 |
| IEEE-33 | caída |   + q en el vértice extremo | 0,000 | 0 | 0,00 | 0,000 | 1,01 | 1,51 | 0,00108 |
| IEEE-33 | caída |   + extremo hallado por búsqueda | 0,000 | 0 | 0,00 | 0,000 | 1,01 | 1,51 | 0,00108 |
| IEEE-33 | caída | Forma afín de la caída, Bernstein + rango de q | 0,000 | 0 | 0,00 | 0,000 | 1,00 | 1,50 | 0,00103 |
| IEEE-33 | caída | Forma afín de la caída, 3σ + rango de q | 0,195 | 0 | 1,31 | 0,195 | 0,67 | 1,00 | 0,00078 |
| IEEE-33 | caída | Resta de las cajas afines corregidas de las dos barras | 0,000 | 0 | 0,00 | 0,000 | 24,73 | 46,80 | 0,03957 |
| IEEE-33 | caída | Resta de los intervalos Bernstein de las dos barras | 0,000 | 0 | 0,00 | 0,000 | 20,96 | 31,05 | 0,02696 |
| IEEE-33 | caída | Resta de los intervalos CD-HIM corregidos | 0,000 | 0 | 0,00 | 0,000 | 13,97 | 24,48 | 0,02047 |
| PEGASE-2869 | barra | Caja afín (peor caso, 1.er orden) | 0,410 | 1 | 0,42 | 0,000 | 1,00 | 3,82 | 0,00478 |
| PEGASE-2869 | barra |   + q en el vértice extremo | 0,000 | 0 | 0,00 | 0,000 | 1,04 | 4,10 | 0,00527 |
| PEGASE-2869 | barra |   + extremo del 2.º orden hallado por búsqueda | 0,000 | 0 | 0,00 | 0,000 | 1,05 | 4,12 | 0,00528 |
| PEGASE-2869 | barra | Bernstein + rango de q | 0,005 | 0 | 0,01 | 0,005 | 0,54 | 1,95 | 0,00232 |
| PEGASE-2869 | barra | 3σ + rango de q (supone normalidad) | 0,290 | 3 | 35,55 | 0,290 | 0,26 | 1,00 | 0,00120 |
| PEGASE-2869 | barra | Intervalo CD-HIM + rango de q (B4) | 0,050 | 0 | 1,93 | 0,055 | 0,48 | 1,51 | 0,00207 |
| PEGASE-2869 | caída | Forma afín de la caída, peor caso | 0,000 | 0 | 0,00 | 0,000 | 1,00 | 2,74 | 0,00117 |
| PEGASE-2869 | caída |   + q en el vértice extremo | 0,000 | 0 | 0,00 | 0,000 | 1,02 | 2,87 | 0,00124 |
| PEGASE-2869 | caída |   + extremo hallado por búsqueda | 0,000 | 0 | 0,00 | 0,000 | 1,02 | 2,89 | 0,00124 |
| PEGASE-2869 | caída | Forma afín de la caída, Bernstein + rango de q | 0,000 | 0 | 0,00 | 0,000 | 0,80 | 1,88 | 0,00071 |
| PEGASE-2869 | caída | Forma afín de la caída, 3σ + rango de q | 0,275 | 1 | 66,05 | 0,275 | 0,36 | 1,00 | 0,00041 |
| PEGASE-2869 | caída | Resta de las cajas afines corregidas de las dos barras | 0,000 | 0 | 0,00 | 0,000 | 8,70 | 27,38 | 0,01367 |
| PEGASE-2869 | caída | Resta de los intervalos Bernstein de las dos barras | 0,000 | 0 | 0,00 | 0,000 | 4,30 | 11,17 | 0,00611 |
| PEGASE-2869 | caída | Resta de los intervalos CD-HIM corregidos | 0,010 | 0 | 0,04 | 0,030 | 3,19 | 8,20 | 0,00465 |

*Script: `experiments/b5_affine_forms.py` · CSV: `results/b5_mc.csv`, `results/b5_width.csv`*

- **Barras.** Bernstein con el rango de q tiene su peor barra en 0,005 % en PEGASE, con un ancho de
  0,54 veces la caja afín y 1,95 veces 3σ en mediana. El intervalo CD-HIM corregido mide
  0,81 veces el de Bernstein en mediana, sin una justificación de cobertura propia.
- **Caídas por rama.** La forma afín de la caída con Bernstein no supera el nominal (peor rama 0,000 %).
  Restar los intervalos de las dos barras da caídas 4,6 veces más anchas en mediana con CD-HIM,
  6,2 con Bernstein y 15,6 con las cajas afines (PEGASE; en IEEE-33, 14,0 veces con
  CD-HIM).

### 28.2 El peor caso en los vértices: la caja de segundo orden no encierra la solución AC

| Sistema | Objetivo | Vértice | Objetivos | AC sin converger | Fuera de la caja afín | Fuera de la caja + q en el vértice | Fuera de la caja + búsqueda | Exceso máx. sobre esta última (p.u.) | Exceso máx. / ancho |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| IEEE-33 | barra | 1.er orden, máximo | 32 | 0 | 0 | 0 | 0 | 0,000000 | 0,000 |
| IEEE-33 | barra | 1.er orden, mínimo | 32 | 0 | 32 | 32 | 32 | 0,000006 | 0,000 |
| IEEE-33 | barra | búsqueda de 2.º orden, máximo | 32 | 0 | 0 | 0 | 0 | 0,000000 | 0,000 |
| IEEE-33 | barra | búsqueda de 2.º orden, mínimo | 32 | 0 | 32 | 32 | 32 | 0,000006 | 0,000 |
| IEEE-33 | caída | 1.er orden, máximo | 31 | 0 | 31 | 27 | 27 | 0,000001 | 0,000 |
| IEEE-33 | caída | 1.er orden, mínimo | 31 | 0 | 0 | 0 | 0 | 0,000000 | 0,000 |
| IEEE-33 | caída | búsqueda de 2.º orden, máximo | 31 | 0 | 31 | 27 | 27 | 0,000001 | 0,000 |
| IEEE-33 | caída | búsqueda de 2.º orden, mínimo | 31 | 0 | 0 | 0 | 0 | 0,000000 | 0,000 |
| PEGASE-2869 | barra | 1.er orden, máximo | 2359 | 0 | 3 | 3 | 3 | 0,000003 | 0,006 |
| PEGASE-2869 | barra | 1.er orden, mínimo | 2359 | 0 | 2271 | 1920 | 1456 | 0,005151 | 0,444 |
| PEGASE-2869 | barra | búsqueda de 2.º orden, máximo | 2359 | 0 | 3 | 3 | 3 | 0,000003 | 0,006 |
| PEGASE-2869 | barra | búsqueda de 2.º orden, mínimo | 2359 | 0 | 2271 | 2065 | 1921 | 0,005155 | 0,633 |
| PEGASE-2869 | caída | 1.er orden, máximo | 3216 | 0 | 1485 | 1210 | 880 | 0,002721 | 0,112 |
| PEGASE-2869 | caída | 1.er orden, mínimo | 3216 | 0 | 1573 | 1269 | 941 | 0,001993 | 0,084 |
| PEGASE-2869 | caída | búsqueda de 2.º orden, máximo | 3216 | 0 | 1509 | 1337 | 1228 | 0,002726 | 0,112 |
| PEGASE-2869 | caída | búsqueda de 2.º orden, mínimo | 3216 | 0 | 1601 | 1390 | 1280 | 0,001999 | 0,132 |

*Script: `experiments/b5_affine_forms.py` · CSV: `results/b5_vertex.csv` · flujo AC exacto en el vértice propuesto
para cada barra o rama*

- **PEGASE-2869:** con todas las fuentes en un extremo, hasta 1921 de 2359 barras y 1280 de
  3216 caídas quedan fuera incluso de la cota que toma el extremo de segundo orden hallado por búsqueda; el
  exceso llega a 0,0052 p.u., el 63 % del ancho de la cota. Pesan los órdenes superiores al segundo.
- **IEEE-33:** también quedan fuera, por 5,9·10⁻⁶ p.u. como máximo (el tercer orden de la sección 18).
- El Monte Carlo no lo detecta: las muestras no llegan a los vértices.
- **Consecuencia:** sin una cota rigurosa del resto, que la sección 22 mostró que no escala a transmisión, no se
  puede afirmar una cota determinista a esta escala. La caja afín queda como referencia del peor caso del modelo
  lineal, no como garantía.

### 28.3 Barrido carga × δ (las 50 celdas de C1, mismas muestras)

| Sistema | Cota | Celdas | Celdas con alguna barra o rama sobre 0,27 % | Máx. barras o ramas sobre 0,27 % en una celda | Peor violación (%) | Ancho / caja afín, mediana de las celdas | máximo |
|---|---|---:|---:|---:|---:|---:|---:|
| IEEE-33 | Barra: caja afín | 25 | 0 | 0 | 0,000 | 1,000 | 1,000 |
| IEEE-33 | Barra: caja afín + q en el vértice | 25 | 0 | 0 | 0,000 | 1,009 | 1,675 |
| IEEE-33 | Barra: Bernstein + rango de q | 25 | 0 | 0 | 0,000 | 0,701 | 0,927 |
| IEEE-33 | Barra: 3σ + rango de q | 25 | 0 | 0 | 0,225 | 0,342 | 0,506 |
| IEEE-33 | Barra: CD-HIM + rango de q | 25 | 0 | 0 | 0,040 | 0,563 | 0,732 |
| IEEE-33 | Caída: forma afín, peor caso + q en el vértice | 25 | 0 | 0 | 0,000 | 1,010 | 1,724 |
| IEEE-33 | Caída: forma afín, Bernstein + rango de q | 25 | 0 | 0 | 0,010 | 1,002 | 1,005 |
| IEEE-33 | Caída: forma afín, 3σ + rango de q | 25 | 1 | 2 | 0,276 | 0,659 | 0,681 |
| IEEE-33 | Caída: resta de intervalos Bernstein | 25 | 0 | 0 | 0,000 | 21,073 | 27,437 |
| IEEE-33 | Caída: resta de intervalos CD-HIM | 25 | 0 | 0 | 0,000 | 14,003 | 21,096 |
| PEGASE-2869 | Barra: caja afín | 25 | 24 | 3 | 5,153 | 1,000 | 1,000 |
| PEGASE-2869 | Barra: caja afín + q en el vértice | 25 | 0 | 0 | 0,006 | 1,148 | 1,420 |
| PEGASE-2869 | Barra: Bernstein + rango de q | 25 | 0 | 0 | 0,018 | 0,517 | 0,586 |
| PEGASE-2869 | Barra: 3σ + rango de q | 25 | 7 | 18 | 0,341 | 0,263 | 0,285 |
| PEGASE-2869 | Barra: CD-HIM + rango de q | 25 | 0 | 0 | 0,125 | 0,473 | 0,504 |
| PEGASE-2869 | Caída: forma afín, peor caso + q en el vértice | 25 | 0 | 0 | 0,000 | 1,062 | 1,258 |
| PEGASE-2869 | Caída: forma afín, Bernstein + rango de q | 25 | 0 | 0 | 0,060 | 0,767 | 0,855 |
| PEGASE-2869 | Caída: forma afín, 3σ + rango de q | 25 | 8 | 10 | 0,310 | 0,359 | 0,388 |
| PEGASE-2869 | Caída: resta de intervalos Bernstein | 25 | 0 | 0 | 0,000 | 4,486 | 4,935 |
| PEGASE-2869 | Caída: resta de intervalos CD-HIM | 25 | 0 | 0 | 0,030 | 3,382 | 3,813 |

*Script: `experiments/b5_affine_forms.py` (modo `sweep`) · CSV: `results/b5_sweep.csv` · sin búsqueda de vértices*

En las 50 celdas, Bernstein con el rango de q no supera el nominal en ninguna barra (peor
0,018 %) ni en ninguna caída (peor 0,060 %). La versión 3σ corregida falla en 7
celdas.

### 28.4 Qué radio: lo que aporta suponer cada cosa

| Sistema | Objetivo | Radio (qué supone) | Semiancho / σ, mediana | Ancho / caja afín, mediana | Topado por la caja afín (% de objetivos) | Violación máx. AC sin corregir (%) | sobre 0,27 % | con el rango de q (%) | sobre 0,27 % |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| IEEE-33 | barra | 3σ (normal) | 3,01 | 0,34 | 0 | 0,275 | 1 | 0,270 | 0 |
| IEEE-33 | barra | Bernstein (soporte y varianza) | 5,92 | 0,70 | 28 | 0,000 | 0 | 0,000 | 0 |
| IEEE-33 | barra | Hoeffding (soporte) | 6,31 | 0,72 | 12 | 0,000 | 0 | 0,000 | 0 |
| IEEE-33 | barra | Gauss, unimodal (varianza) | 8,82 | 1,00 | 100 | 0,000 | 0 | 0,000 | 0 |
| IEEE-33 | barra | Chebyshev (varianza) | 8,82 | 1,00 | 100 | 0,000 | 0 | 0,000 | 0 |
| IEEE-33 | barra | Caja afín (peor caso) | 8,82 | 1,00 | 100 | 0,000 | 0 | 0,000 | 0 |
| IEEE-33 | caída | 3σ (normal) | 3,01 | 0,67 | 13 | 0,195 | 0 | 0,195 | 0 |
| IEEE-33 | caída | Bernstein (soporte y varianza) | 4,51 | 1,00 | 77 | 0,000 | 0 | 0,000 | 0 |
| IEEE-33 | caída | Hoeffding (soporte) | 4,51 | 1,00 | 77 | 0,000 | 0 | 0,000 | 0 |
| IEEE-33 | caída | Gauss, unimodal (varianza) | 4,51 | 1,00 | 100 | 0,000 | 0 | 0,000 | 0 |
| IEEE-33 | caída | Chebyshev (varianza) | 4,51 | 1,00 | 100 | 0,000 | 0 | 0,000 | 0 |
| IEEE-33 | caída | Caja afín (peor caso) | 4,51 | 1,00 | 100 | 0,000 | 0 | 0,000 | 0 |
| PEGASE-2869 | barra | 3σ (normal) | 3,02 | 0,26 | 4 | 1,190 | 81 | 0,295 | 4 |
| PEGASE-2869 | barra | Bernstein (soporte y varianza) | 5,90 | 0,54 | 20 | 0,410 | 1 | 0,005 | 0 |
| PEGASE-2869 | barra | Hoeffding (soporte) | 6,32 | 0,55 | 16 | 0,410 | 1 | 0,000 | 0 |
| PEGASE-2869 | barra | Gauss, unimodal (varianza) | 11,65 | 1,00 | 58 | 0,410 | 1 | 0,000 | 0 |
| PEGASE-2869 | barra | Chebyshev (varianza) | 11,65 | 1,00 | 85 | 0,410 | 1 | 0,000 | 0 |
| PEGASE-2869 | barra | Caja afín (peor caso) | 11,65 | 1,00 | 100 | 0,410 | 1 | 0,000 | 0 |
| PEGASE-2869 | caída | 3σ (normal) | 3,01 | 0,36 | 14 | 0,380 | 69 | 0,275 | 1 |
| PEGASE-2869 | caída | Bernstein (soporte y varianza) | 5,70 | 0,80 | 40 | 0,005 | 0 | 0,000 | 0 |
| PEGASE-2869 | caída | Hoeffding (soporte) | 6,31 | 0,76 | 36 | 0,000 | 0 | 0,000 | 0 |
| PEGASE-2869 | caída | Gauss, unimodal (varianza) | 8,32 | 1,00 | 73 | 0,000 | 0 | 0,000 | 0 |
| PEGASE-2869 | caída | Chebyshev (varianza) | 8,32 | 1,00 | 88 | 0,000 | 0 | 0,000 | 0 |
| PEGASE-2869 | caída | Caja afín (peor caso) | 8,32 | 1,00 | 100 | 0,000 | 0 | 0,000 | 0 |

*Script: `experiments/b5_affine_forms.py` (modo `radii`) · CSV: `results/b5_radii.csv` · 1.000 muestras de q*

- **Los radios que usan solo la varianza no sirven:** Chebyshev pide 11,7σ en PEGASE y queda topado por la
  caja afín en el 85 % de las barras; la desigualdad de Gauss para distribuciones unimodales, en el
  58 %. Equivalen al peor caso.
- **Los que usan el soporte acotado sí:** Hoeffding (6,3σ) y Bernstein (5,9σ) quedan en la mitad de la
  caja afín. La varianza que agrega Bernstein gana solo 7 % sobre Hoeffding.
- Todos necesitan la corrección de segundo orden en las barras con una sola fuente dominante, donde el radio queda
  topado por la caja afín y el error de linealización lo cruza.

### 28.5 Cuántas muestras del término cuadrático hacen falta, y el costo del método

| Sistema | Objetivo | Muestras de q | Peor violación, 3 semillas (%) | Sobre 0,27 %, peor semilla | Ancho / sin corregir, mediana | Tiempo del muestreo (s) |
|---|---|---|---:|---:|---:|---:|
| IEEE-33 | barra | 0 (sin corregir) | 0,000 | 0 | 1,0000 | 0,00 |
| IEEE-33 | barra | 250 | 0,000 | 0 | 1,0014 | 0,00 |
| IEEE-33 | barra | 500 | 0,000 | 0 | 1,0014 | 0,00 |
| IEEE-33 | barra | 1000 | 0,000 | 0 | 1,0014 | 0,00 |
| IEEE-33 | barra | 2000 | 0,000 | 0 | 1,0016 | 0,00 |
| IEEE-33 | barra | 4000 | 0,000 | 0 | 1,0019 | 0,01 |
| IEEE-33 | caída | 0 (sin corregir) | 0,000 | 0 | 1,0000 | 0,00 |
| IEEE-33 | caída | 250 | 0,000 | 0 | 1,0013 | 0,00 |
| IEEE-33 | caída | 500 | 0,000 | 0 | 1,0015 | 0,00 |
| IEEE-33 | caída | 1000 | 0,000 | 0 | 1,0015 | 0,00 |
| IEEE-33 | caída | 2000 | 0,000 | 0 | 1,0018 | 0,00 |
| IEEE-33 | caída | 4000 | 0,000 | 0 | 1,0020 | 0,01 |
| PEGASE-2869 | barra | 0 (sin corregir) | 0,410 | 1 | 1,0000 | 0,00 |
| PEGASE-2869 | barra | 250 | 0,010 | 0 | 1,0032 | 0,26 |
| PEGASE-2869 | barra | 500 | 0,010 | 0 | 1,0036 | 0,40 |
| PEGASE-2869 | barra | 1000 | 0,010 | 0 | 1,0038 | 0,87 |
| PEGASE-2869 | barra | 2000 | 0,005 | 0 | 1,0040 | 1,96 |
| PEGASE-2869 | barra | 4000 | 0,005 | 0 | 1,0044 | 3,76 |
| PEGASE-2869 | caída | 0 (sin corregir) | 0,005 | 0 | 1,0000 | 0,00 |
| PEGASE-2869 | caída | 250 | 0,005 | 0 | 1,0020 | 0,26 |
| PEGASE-2869 | caída | 500 | 0,005 | 0 | 1,0022 | 0,40 |
| PEGASE-2869 | caída | 1000 | 0,005 | 0 | 1,0023 | 0,87 |
| PEGASE-2869 | caída | 2000 | 0,000 | 0 | 1,0025 | 1,96 |
| PEGASE-2869 | caída | 4000 | 0,000 | 0 | 1,0027 | 3,76 |

*Script: `experiments/b5_affine_forms.py` (modo `surrogate`) · CSV: `results/b5_surrogate.csv` · tres semillas por
tamaño · un proceso*

Con 250 muestras la cobertura no cambia: el margen del radio absorbe el error de estimar la cola del rango. El
método completo en PEGASE (estudio lineal, modelo de segundo orden y 250 muestras) toma 4,1 s, unas
15 veces menos que el Monte Carlo AC por lotes de C8 (60 s). El barrido 28.3 usó 4.000 muestras;
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

| Sistema | Correlación de las fuentes | Balance | Bernstein suponiendo independencia: violación máx. (%) | Bernstein con el modelo de factores (cópula: factor con el mismo ρ): violación máx. (%) | sobre 0,27 % | su ancho / caja afín, mediana | 3σ con la Σ_ε verdadera: violación máx. (%) | Caja afín (peor caso lineal): violación máx. (%) | AC sin converger |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| IEEE-33 | independientes | generadores | 0,000 | 0,000 | 0 | 0,70 | 0,200 | 0,000 | 0 |
| IEEE-33 | independientes | slack | 0,000 | 0,000 | 0 | 0,70 | 0,290 | 0,000 | 0 |
| IEEE-33 | factor global, ρ = 0,155 | generadores | 0,000 | 0,000 | 0 | 1,00 | 0,005 | 0,000 | 0 |
| IEEE-33 | factor global, ρ = 0,155 | slack | 0,000 | 0,000 | 0 | 1,00 | 0,005 | 0,000 | 0 |
| IEEE-33 | factor global, ρ = 0,5 | generadores | 0,645 | 0,000 | 0 | 1,00 | 0,000 | 0,000 | 0 |
| IEEE-33 | factor global, ρ = 0,5 | slack | 0,625 | 0,000 | 0 | 1,00 | 0,000 | 0,000 | 0 |
| IEEE-33 | factor global, ρ = 0,845 | generadores | 19,415 | 0,000 | 0 | 1,00 | 0,000 | 0,000 | 0 |
| IEEE-33 | factor global, ρ = 0,845 | slack | 20,125 | 0,000 | 0 | 1,00 | 0,000 | 0,000 | 0 |
| IEEE-33 | 8 zonas, ρ = 0,155 | generadores | 0,000 | 0,000 | 0 | 0,60 | 0,225 | 0,000 | 0 |
| IEEE-33 | 8 zonas, ρ = 0,155 | slack | 0,000 | 0,000 | 0 | 0,59 | 0,255 | 0,000 | 0 |
| IEEE-33 | 8 zonas, ρ = 0,5 | generadores | 0,000 | 0,000 | 0 | 0,83 | 0,200 | 0,000 | 0 |
| IEEE-33 | 8 zonas, ρ = 0,5 | slack | 0,000 | 0,000 | 0 | 0,83 | 0,160 | 0,000 | 0 |
| IEEE-33 | 8 zonas, ρ = 0,845 | generadores | 0,005 | 0,000 | 0 | 1,00 | 0,160 | 0,000 | 0 |
| IEEE-33 | 8 zonas, ρ = 0,845 | slack | 0,005 | 0,000 | 0 | 1,00 | 0,125 | 0,000 | 0 |
| IEEE-33 | cópula de C2, r = 0,3 | generadores | 7,795 | 0,000 | 0 | 1,00 | 0,000 | 0,000 | 0 |
| IEEE-33 | cópula de C2, r = 0,3 | slack | 7,325 | 0,000 | 0 | 1,01 | 0,000 | 0,000 | 0 |
| IEEE-33 | cópula de C2, r = 0,6 | generadores | 24,715 | 0,000 | 0 | 1,01 | 0,000 | 0,000 | 0 |
| IEEE-33 | cópula de C2, r = 0,6 | slack | 24,405 | 0,000 | 0 | 1,01 | 0,000 | 0,000 | 0 |
| IEEE-33 | cópula de C2, r = 0,9 | generadores | 39,345 | 0,015 | 0 | 1,01 | 0,000 | 0,015 | 0 |
| IEEE-33 | cópula de C2, r = 0,9 | slack | 39,740 | 0,015 | 0 | 1,01 | 0,000 | 0,015 | 0 |
| PEGASE-2869 | independientes | generadores | 0,005 | 0,005 | 0 | 0,52 | 0,310 | 0,000 | 0 |
| PEGASE-2869 | independientes | slack | 0,015 | 0,015 | 0 | 0,54 | 0,295 | 0,000 | 0 |
| PEGASE-2869 | factor global, ρ = 0,155 | generadores | 65,590 | 0,005 | 0 | 0,92 | 0,260 | 0,000 | 0 |
| PEGASE-2869 | factor global, ρ = 0,155 | slack | 56,555 | 0,000 | 0 | 0,95 | 0,170 | 0,000 | 0 |
| PEGASE-2869 | factor global, ρ = 0,5 | generadores | 79,440 | 0,000 | 0 | 0,87 | 0,135 | 0,000 | 0 |
| PEGASE-2869 | factor global, ρ = 0,5 | slack | 68,985 | 0,000 | 0 | 0,92 | 0,110 | 0,000 | 0 |
| PEGASE-2869 | factor global, ρ = 0,845 | generadores | 85,415 | 0,000 | 0 | 0,82 | 0,075 | 0,000 | 0 |
| PEGASE-2869 | factor global, ρ = 0,845 | slack | 73,900 | 0,000 | 0 | 0,89 | 0,135 | 0,000 | 0 |
| PEGASE-2869 | 8 zonas, ρ = 0,155 | generadores | 11,255 | 0,000 | 0 | 0,93 | 0,215 | 0,000 | 0 |
| PEGASE-2869 | 8 zonas, ρ = 0,155 | slack | 10,540 | 0,000 | 0 | 0,94 | 0,185 | 0,000 | 0 |
| PEGASE-2869 | 8 zonas, ρ = 0,5 | generadores | 39,540 | 0,000 | 0 | 0,92 | 0,180 | 0,000 | 0 |
| PEGASE-2869 | 8 zonas, ρ = 0,5 | slack | 40,240 | 0,000 | 0 | 0,93 | 0,155 | 0,000 | 0 |
| PEGASE-2869 | 8 zonas, ρ = 0,845 | generadores | 56,315 | 0,000 | 0 | 0,90 | 0,165 | 0,000 | 0 |
| PEGASE-2869 | 8 zonas, ρ = 0,845 | slack | 57,000 | 0,000 | 0 | 0,92 | 0,080 | 0,000 | 0 |
| PEGASE-2869 | cópula de C2, r = 0,3 | generadores | 74,895 | 3,310 | 4 | 0,90 | 0,165 | 0,000 | 0 |
| PEGASE-2869 | cópula de C2, r = 0,3 | slack | 54,975 | 3,440 | 4 | 0,95 | 0,180 | 0,040 | 0 |
| PEGASE-2869 | cópula de C2, r = 0,6 | generadores | 83,615 | 3,195 | 4 | 0,86 | 0,275 | 0,000 | 0 |
| PEGASE-2869 | cópula de C2, r = 0,6 | slack | 67,670 | 3,015 | 25 | 0,93 | 0,380 | 0,425 | 0 |
| PEGASE-2869 | cópula de C2, r = 0,9 | generadores | 88,505 | 0,635 | 1 | 0,81 | 0,540 | 0,000 | 0 |
| PEGASE-2869 | cópula de C2, r = 0,9 | slack | 76,025 | 2,290 | 293 | 0,90 | 1,265 | 1,870 | 0 |

**Caídas por rama**

| Sistema | Correlación de las fuentes | Balance | Bernstein suponiendo independencia: violación máx. (%) | Bernstein con el modelo de factores (cópula: factor con el mismo ρ): violación máx. (%) | sobre 0,27 % | su ancho / caja afín, mediana | 3σ con la Σ_ε verdadera: violación máx. (%) | Caja afín (peor caso lineal): violación máx. (%) | AC sin converger |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| IEEE-33 | independientes | generadores | 0,000 | 0,000 | 0 | 1,00 | 0,185 | 0,000 | 0 |
| IEEE-33 | independientes | slack | 0,000 | 0,000 | 0 | 1,00 | 0,210 | 0,000 | 0 |
| IEEE-33 | factor global, ρ = 0,155 | generadores | 0,000 | 0,000 | 0 | 1,00 | 0,045 | 0,000 | 0 |
| IEEE-33 | factor global, ρ = 0,155 | slack | 0,000 | 0,000 | 0 | 1,00 | 0,055 | 0,000 | 0 |
| IEEE-33 | factor global, ρ = 0,5 | generadores | 0,000 | 0,000 | 0 | 1,00 | 0,000 | 0,000 | 0 |
| IEEE-33 | factor global, ρ = 0,5 | slack | 0,005 | 0,000 | 0 | 1,00 | 0,000 | 0,000 | 0 |
| IEEE-33 | factor global, ρ = 0,845 | generadores | 3,965 | 0,000 | 0 | 1,00 | 0,000 | 0,000 | 0 |
| IEEE-33 | factor global, ρ = 0,845 | slack | 4,125 | 0,000 | 0 | 1,00 | 0,000 | 0,000 | 0 |
| IEEE-33 | 8 zonas, ρ = 0,155 | generadores | 0,000 | 0,000 | 0 | 1,00 | 0,215 | 0,000 | 0 |
| IEEE-33 | 8 zonas, ρ = 0,155 | slack | 0,000 | 0,000 | 0 | 1,00 | 0,180 | 0,000 | 0 |
| IEEE-33 | 8 zonas, ρ = 0,5 | generadores | 0,000 | 0,000 | 0 | 1,00 | 0,125 | 0,000 | 0 |
| IEEE-33 | 8 zonas, ρ = 0,5 | slack | 0,000 | 0,000 | 0 | 1,00 | 0,135 | 0,000 | 0 |
| IEEE-33 | 8 zonas, ρ = 0,845 | generadores | 0,000 | 0,000 | 0 | 1,00 | 0,110 | 0,000 | 0 |
| IEEE-33 | 8 zonas, ρ = 0,845 | slack | 0,000 | 0,000 | 0 | 1,00 | 0,085 | 0,000 | 0 |
| IEEE-33 | cópula de C2, r = 0,3 | generadores | 2,665 | 0,000 | 0 | 1,01 | 0,000 | 0,000 | 0 |
| IEEE-33 | cópula de C2, r = 0,3 | slack | 2,420 | 0,000 | 0 | 1,01 | 0,000 | 0,000 | 0 |
| IEEE-33 | cópula de C2, r = 0,6 | generadores | 13,985 | 0,000 | 0 | 1,01 | 0,000 | 0,000 | 0 |
| IEEE-33 | cópula de C2, r = 0,6 | slack | 13,760 | 0,000 | 0 | 1,01 | 0,000 | 0,000 | 0 |
| IEEE-33 | cópula de C2, r = 0,9 | generadores | 27,895 | 0,025 | 0 | 1,01 | 0,000 | 0,025 | 0 |
| IEEE-33 | cópula de C2, r = 0,9 | slack | 27,595 | 0,025 | 0 | 1,01 | 0,000 | 0,025 | 0 |
| PEGASE-2869 | independientes | generadores | 0,005 | 0,005 | 0 | 0,78 | 0,310 | 0,000 | 0 |
| PEGASE-2869 | independientes | slack | 0,000 | 0,000 | 0 | 0,80 | 0,300 | 0,000 | 0 |
| PEGASE-2869 | factor global, ρ = 0,155 | generadores | 36,135 | 0,005 | 0 | 0,85 | 0,305 | 0,000 | 0 |
| PEGASE-2869 | factor global, ρ = 0,155 | slack | 57,105 | 0,000 | 0 | 0,90 | 0,245 | 0,000 | 0 |
| PEGASE-2869 | factor global, ρ = 0,5 | generadores | 61,600 | 0,005 | 0 | 0,76 | 0,320 | 0,000 | 0 |
| PEGASE-2869 | factor global, ρ = 0,5 | slack | 69,395 | 0,000 | 0 | 0,83 | 0,215 | 0,000 | 0 |
| PEGASE-2869 | factor global, ρ = 0,845 | generadores | 72,790 | 0,000 | 0 | 0,66 | 0,240 | 0,000 | 0 |
| PEGASE-2869 | factor global, ρ = 0,845 | slack | 74,330 | 0,000 | 0 | 0,77 | 0,205 | 0,000 | 0 |
| PEGASE-2869 | 8 zonas, ρ = 0,155 | generadores | 9,255 | 0,000 | 0 | 0,83 | 0,290 | 0,000 | 0 |
| PEGASE-2869 | 8 zonas, ρ = 0,155 | slack | 10,680 | 0,000 | 0 | 0,85 | 0,255 | 0,000 | 0 |
| PEGASE-2869 | 8 zonas, ρ = 0,5 | generadores | 38,735 | 0,000 | 0 | 0,83 | 0,215 | 0,000 | 0 |
| PEGASE-2869 | 8 zonas, ρ = 0,5 | slack | 39,340 | 0,000 | 0 | 0,85 | 0,245 | 0,000 | 0 |
| PEGASE-2869 | 8 zonas, ρ = 0,845 | generadores | 55,785 | 0,000 | 0 | 0,80 | 0,185 | 0,000 | 0 |
| PEGASE-2869 | 8 zonas, ρ = 0,845 | slack | 56,480 | 0,000 | 0 | 0,82 | 0,190 | 0,000 | 0 |
| PEGASE-2869 | cópula de C2, r = 0,3 | generadores | 55,325 | 5,990 | 22 | 0,81 | 0,295 | 0,000 | 0 |
| PEGASE-2869 | cópula de C2, r = 0,3 | slack | 55,815 | 5,870 | 24 | 0,89 | 0,290 | 0,030 | 0 |
| PEGASE-2869 | cópula de C2, r = 0,6 | generadores | 69,535 | 5,930 | 28 | 0,74 | 0,465 | 0,000 | 0 |
| PEGASE-2869 | cópula de C2, r = 0,6 | slack | 68,940 | 5,905 | 26 | 0,84 | 0,385 | 0,195 | 0 |
| PEGASE-2869 | cópula de C2, r = 0,9 | generadores | 77,740 | 2,895 | 6 | 0,64 | 0,720 | 0,000 | 0 |
| PEGASE-2869 | cópula de C2, r = 0,9 | slack | 77,020 | 2,890 | 226 | 0,77 | 0,820 | 1,090 | 0 |

*Script: `experiments/b6_correlated.py` · CSV: `results/b6_correlated.csv` · semilla por caso en la columna `seed`*

- **Suponer independencia con fuentes correlacionadas no sirve:** en PEGASE la peor barra queda fuera entre el
  10,5 % y el 89 % de las muestras.
- **Con la correlación escrita como símbolos compartidos el método cubre** (peor barra 0,005 % en PEGASE), pero
  la concentración se pierde: el ancho pasa de 0,54 veces la caja afín con fuentes independientes a entre
  0,82 y 0,95 veces. Un factor común no se promedia. Con un factor global en IEEE-33 el radio queda
  topado por la caja afín.
- **Un modelo de factores mal especificado falla:** el de un factor con la misma correlación que la cópula deja fuera
  hasta el 3,4 % de las muestras en una barra de PEGASE. Igualar la correlación no basta: la cópula tiene otra
  dependencia en las colas. 3σ con la Σ_ε verdadera llega a 1,26 %.
- **El peor caso lineal es lo único robusto a la especificación**, con una condición: con la slack sola y la cópula
  fuerte, en PEGASE incluso la caja afín con el rango de q queda fuera hasta el 1,87 % (el error de
  linealización de C2, que el segundo orden no alcanza a cubrir); con los generadores repartiendo el desbalance,
  0,000 %.
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

| Sistema | δ | Objetivo | Nivel | Radio | Peor objetivo (%) | Muestras con alguno fuera (%) | con alguna barra o caída fuera (%) | Ancho / caja afín, mediana | Topado por la caja afín (%) |
|---|---:|---|---|---|---:|---:|---:|---:|---:|
| IEEE-33 | 0,15 | barra | simultáneo, Bonferroni | Bernstein | 0,000 | 0,000 | 0,000 | 0,96 | 47 |
| IEEE-33 | 0,15 | barra | simultáneo, Bonferroni | Gauss | 0,005 | 0,005 | 0,010 | 0,45 | 0 |
| IEEE-33 | 0,15 | barra | simultáneo, Bonferroni | Chebyshev | 0,000 | 0,000 | 0,000 | 1,00 | 100 |
| IEEE-33 | 0,15 | barra | por objetivo (α) | Bernstein | 0,000 | 0,000 | 0,000 | 0,70 | 28 |
| IEEE-33 | 0,15 | barra | por objetivo (α) | Gauss | 0,215 | 0,670 | 1,430 | 0,34 | 0 |
| IEEE-33 | 0,15 | barra | por objetivo (α) | Chebyshev | 0,000 | 0,000 | 0,000 | 1,00 | 100 |
| IEEE-33 | 0,15 | barra | por objetivo (α) | caja afín | 0,000 | 0,000 | 0,000 | 1,00 | 100 |
| IEEE-33 | 0,15 | caída | simultáneo, Bonferroni | Bernstein | 0,000 | 0,000 | 0,000 | 1,00 | 97 |
| IEEE-33 | 0,15 | caída | simultáneo, Bonferroni | Gauss | 0,005 | 0,010 | 0,010 | 0,87 | 35 |
| IEEE-33 | 0,15 | caída | simultáneo, Bonferroni | Chebyshev | 0,000 | 0,000 | 0,000 | 1,00 | 100 |
| IEEE-33 | 0,15 | caída | por objetivo (α) | Bernstein | 0,000 | 0,000 | 0,000 | 1,00 | 77 |
| IEEE-33 | 0,15 | caída | por objetivo (α) | Gauss | 0,200 | 1,175 | 1,430 | 0,67 | 13 |
| IEEE-33 | 0,15 | caída | por objetivo (α) | Chebyshev | 0,000 | 0,000 | 0,000 | 1,00 | 100 |
| IEEE-33 | 0,15 | caída | por objetivo (α) | caja afín | 0,000 | 0,000 | 0,000 | 1,00 | 100 |
| PEGASE-2869 | 0,05 | barra | simultáneo, Bonferroni | Bernstein | 0,000 | 0,000 | 0,000 | 0,99 | 49 |
| PEGASE-2869 | 0,05 | barra | simultáneo, Bonferroni | Gauss | 0,005 | 0,005 | 0,005 | 0,42 | 9 |
| PEGASE-2869 | 0,05 | barra | simultáneo, Bonferroni | Chebyshev | 0,000 | 0,000 | 0,000 | 1,00 | 100 |
| PEGASE-2869 | 0,05 | barra | por objetivo (α) | Bernstein | 0,015 | 0,020 | 0,025 | 0,54 | 20 |
| PEGASE-2869 | 0,05 | barra | por objetivo (α) | Gauss | 0,295 | 37,325 | 72,550 | 0,26 | 4 |
| PEGASE-2869 | 0,05 | barra | por objetivo (α) | Chebyshev | 0,000 | 0,000 | 0,000 | 1,00 | 85 |
| PEGASE-2869 | 0,05 | barra | por objetivo (α) | caja afín | 0,000 | 0,000 | 0,000 | 1,00 | 100 |
| PEGASE-2869 | 0,05 | caída | simultáneo, Bonferroni | Bernstein | 0,000 | 0,000 | 0,000 | 1,00 | 68 |
| PEGASE-2869 | 0,05 | caída | simultáneo, Bonferroni | Gauss | 0,000 | 0,000 | 0,005 | 0,60 | 26 |
| PEGASE-2869 | 0,05 | caída | simultáneo, Bonferroni | Chebyshev | 0,000 | 0,000 | 0,000 | 1,00 | 100 |
| PEGASE-2869 | 0,05 | caída | por objetivo (α) | Bernstein | 0,005 | 0,010 | 0,025 | 0,80 | 40 |
| PEGASE-2869 | 0,05 | caída | por objetivo (α) | Gauss | 0,290 | 66,975 | 72,550 | 0,36 | 14 |
| PEGASE-2869 | 0,05 | caída | por objetivo (α) | Chebyshev | 0,000 | 0,000 | 0,000 | 1,00 | 88 |
| PEGASE-2869 | 0,05 | caída | por objetivo (α) | caja afín | 0,000 | 0,000 | 0,000 | 1,00 | 100 |
| IEEE-118 | 0,05 | barra | simultáneo, Bonferroni | Bernstein | 0,000 | 0,000 | 0,000 | 1,00 | 87 |
| IEEE-118 | 0,05 | barra | simultáneo, Bonferroni | Gauss | 0,005 | 0,020 | 0,025 | 0,89 | 43 |
| IEEE-118 | 0,05 | barra | simultáneo, Bonferroni | Chebyshev | 0,000 | 0,000 | 0,000 | 1,00 | 100 |
| IEEE-118 | 0,05 | barra | por objetivo (α) | Bernstein | 0,000 | 0,000 | 0,000 | 1,00 | 76 |
| IEEE-118 | 0,05 | barra | por objetivo (α) | Gauss | 0,185 | 2,130 | 3,205 | 0,66 | 21 |
| IEEE-118 | 0,05 | barra | por objetivo (α) | Chebyshev | 0,000 | 0,000 | 0,000 | 1,00 | 100 |
| IEEE-118 | 0,05 | barra | por objetivo (α) | caja afín | 0,000 | 0,000 | 0,000 | 1,00 | 100 |
| IEEE-118 | 0,05 | caída | simultáneo, Bonferroni | Bernstein | 0,000 | 0,000 | 0,000 | 1,00 | 91 |
| IEEE-118 | 0,05 | caída | simultáneo, Bonferroni | Gauss | 0,005 | 0,010 | 0,025 | 0,80 | 37 |
| IEEE-118 | 0,05 | caída | simultáneo, Bonferroni | Chebyshev | 0,000 | 0,000 | 0,000 | 1,00 | 100 |
| IEEE-118 | 0,05 | caída | por objetivo (α) | Bernstein | 0,000 | 0,000 | 0,000 | 1,00 | 74 |
| IEEE-118 | 0,05 | caída | por objetivo (α) | Gauss | 0,200 | 1,605 | 3,205 | 0,60 | 11 |
| IEEE-118 | 0,05 | caída | por objetivo (α) | Chebyshev | 0,000 | 0,000 | 0,000 | 1,00 | 100 |
| IEEE-118 | 0,05 | caída | por objetivo (α) | caja afín | 0,000 | 0,000 | 0,000 | 1,00 | 100 |
| PEGASE-1354 | 0,05 | barra | simultáneo, Bonferroni | Bernstein | 0,000 | 0,000 | 0,000 | 1,00 | 61 |
| PEGASE-1354 | 0,05 | barra | simultáneo, Bonferroni | Gauss | 0,000 | 0,000 | 0,000 | 0,53 | 16 |
| PEGASE-1354 | 0,05 | barra | simultáneo, Bonferroni | Chebyshev | 0,000 | 0,000 | 0,000 | 1,00 | 100 |
| PEGASE-1354 | 0,05 | barra | por objetivo (α) | Bernstein | 0,000 | 0,000 | 0,000 | 0,74 | 33 |
| PEGASE-1354 | 0,05 | barra | por objetivo (α) | Gauss | 0,305 | 20,645 | 39,990 | 0,34 | 8 |
| PEGASE-1354 | 0,05 | barra | por objetivo (α) | Chebyshev | 0,000 | 0,000 | 0,000 | 1,00 | 94 |
| PEGASE-1354 | 0,05 | barra | por objetivo (α) | caja afín | 0,000 | 0,000 | 0,000 | 1,00 | 100 |
| PEGASE-1354 | 0,05 | caída | simultáneo, Bonferroni | Bernstein | 0,000 | 0,000 | 0,000 | 1,00 | 73 |
| PEGASE-1354 | 0,05 | caída | simultáneo, Bonferroni | Gauss | 0,000 | 0,000 | 0,000 | 0,69 | 32 |
| PEGASE-1354 | 0,05 | caída | simultáneo, Bonferroni | Chebyshev | 0,000 | 0,000 | 0,000 | 1,00 | 100 |
| PEGASE-1354 | 0,05 | caída | por objetivo (α) | Bernstein | 0,000 | 0,000 | 0,000 | 0,98 | 49 |
| PEGASE-1354 | 0,05 | caída | por objetivo (α) | Gauss | 0,320 | 33,495 | 39,990 | 0,43 | 20 |
| PEGASE-1354 | 0,05 | caída | por objetivo (α) | Chebyshev | 0,000 | 0,000 | 0,000 | 1,00 | 93 |
| PEGASE-1354 | 0,05 | caída | por objetivo (α) | caja afín | 0,000 | 0,000 | 0,000 | 1,00 | 100 |

**Grilla de IEEE-118 y PEGASE-1354**

| Sistema | Objetivo | Nivel | Radio | Celdas | Celdas con algún objetivo sobre 0,27 % | Peor objetivo (%) | Peor celda: muestras con alguna barra o caída fuera (%) | Ancho / caja afín, mediana de las celdas | Máx. muestras AC sin converger en una celda |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
| IEEE-118 | barra | por objetivo (α) | Bernstein | 12 | 0 | 0,005 | 0,005 | 1,01 | 4413 |
| IEEE-118 | barra | por objetivo (α) | Gauss | 12 | 0 | 0,220 | 3,435 | 0,66 | 4413 |
| IEEE-118 | barra | por objetivo (α) | Chebyshev | 12 | 0 | 0,000 | 0,000 | 1,01 | 4413 |
| IEEE-118 | barra | por objetivo (α) | caja afín | 12 | 0 | 0,000 | 0,000 | 1,01 | 4413 |
| IEEE-118 | barra | simultáneo, Bonferroni | Bernstein | 12 | 0 | 0,000 | 0,000 | 1,01 | 4413 |
| IEEE-118 | barra | simultáneo, Bonferroni | Gauss | 12 | 0 | 0,060 | 0,080 | 0,90 | 4413 |
| IEEE-118 | barra | simultáneo, Bonferroni | Chebyshev | 12 | 0 | 0,000 | 0,000 | 1,01 | 4413 |
| IEEE-118 | caída | por objetivo (α) | Bernstein | 12 | 0 | 0,000 | 0,005 | 1,01 | 4413 |
| IEEE-118 | caída | por objetivo (α) | Gauss | 12 | 0 | 0,210 | 3,435 | 0,59 | 4413 |
| IEEE-118 | caída | por objetivo (α) | Chebyshev | 12 | 0 | 0,000 | 0,000 | 1,01 | 4413 |
| IEEE-118 | caída | por objetivo (α) | caja afín | 12 | 0 | 0,000 | 0,000 | 1,01 | 4413 |
| IEEE-118 | caída | simultáneo, Bonferroni | Bernstein | 12 | 0 | 0,000 | 0,000 | 1,01 | 4413 |
| IEEE-118 | caída | simultáneo, Bonferroni | Gauss | 12 | 0 | 0,015 | 0,080 | 0,78 | 4413 |
| IEEE-118 | caída | simultáneo, Bonferroni | Chebyshev | 12 | 0 | 0,000 | 0,000 | 1,01 | 4413 |
| PEGASE-1354 | barra | por objetivo (α) | Bernstein | 12 | 0 | 0,045 | 0,085 | 0,72 | 6756 |
| PEGASE-1354 | barra | por objetivo (α) | Gauss | 12 | 6 | 0,340 | 40,605 | 0,33 | 6756 |
| PEGASE-1354 | barra | por objetivo (α) | Chebyshev | 12 | 0 | 0,045 | 0,080 | 1,01 | 6756 |
| PEGASE-1354 | barra | por objetivo (α) | caja afín | 12 | 0 | 0,045 | 0,080 | 1,01 | 6756 |
| PEGASE-1354 | barra | simultáneo, Bonferroni | Bernstein | 12 | 0 | 0,045 | 0,080 | 1,00 | 6756 |
| PEGASE-1354 | barra | simultáneo, Bonferroni | Gauss | 12 | 0 | 0,045 | 0,110 | 0,52 | 6756 |
| PEGASE-1354 | barra | simultáneo, Bonferroni | Chebyshev | 12 | 0 | 0,045 | 0,080 | 1,01 | 6756 |
| PEGASE-1354 | caída | por objetivo (α) | Bernstein | 12 | 0 | 0,065 | 0,085 | 0,97 | 6756 |
| PEGASE-1354 | caída | por objetivo (α) | Gauss | 12 | 6 | 0,405 | 40,605 | 0,44 | 6756 |
| PEGASE-1354 | caída | por objetivo (α) | Chebyshev | 12 | 0 | 0,065 | 0,080 | 1,00 | 6756 |
| PEGASE-1354 | caída | por objetivo (α) | caja afín | 12 | 0 | 0,065 | 0,080 | 1,00 | 6756 |
| PEGASE-1354 | caída | simultáneo, Bonferroni | Bernstein | 12 | 0 | 0,065 | 0,080 | 1,00 | 6756 |
| PEGASE-1354 | caída | simultáneo, Bonferroni | Gauss | 12 | 0 | 0,065 | 0,110 | 0,68 | 6756 |
| PEGASE-1354 | caída | simultáneo, Bonferroni | Chebyshev | 12 | 0 | 0,065 | 0,080 | 1,00 | 6756 |

*Script: `experiments/b7_generality.py` · CSV: `results/b7_generality.csv` · las muestras AC sin converger se excluyen,
como en C1*

- **La cobertura simultánea con un radio sin distribución cuesta lo mismo que el peor caso.** En PEGASE-2869 el radio
  de Bernstein pasa de 0,54 veces la caja afín por barra a 0,99 con Bonferroni. Lo mismo ocurre en
  los cuatro sistemas. El crecimiento en √ln(m/α) es lento, pero el radio por barra ya estaba cerca del tope.
- **Gauss con Bonferroni cubre en la práctica** con 0,42 veces la caja afín (muestras con alguna barra fuera:
  0,005 %), pero sin garantía: por barra, Gauss deja 37 % de las muestras con alguna
  barra fuera, y en PEGASE-1354 falla por barra en 6 de las 12 celdas.
- **Bernstein por barra ya da una cobertura del sistema cercana al nominal, sin estar diseñado para eso:** muestras con
  alguna barra o caída fuera, 0,025 % en PEGASE-2869 y como máximo 0,085 % en la grilla.
  Es una observación empírica, no una garantía.
- **La ventaja de la concentración depende de cuántas fuentes comparables mueven cada barra.** Bernstein por barra mide
  0,74 veces la caja afín en PEGASE-1354 y 1,00 en IEEE-118, donde queda topado en el 76 % de las
  barras: allí pocas cargas mueven cada barra y no hay qué promediar.
- **Generalidad de la corrección:** en los dos sistemas nuevos Bernstein con el rango de q no supera el nominal en
  ninguna celda. Cerca del límite de cargabilidad con δ alto no converge hasta el 34 % de las muestras AC; esas
  celdas hay que leerlas con la advertencia de C1.

---

## 29. Lista de afirmaciones para el paper de formas afines

Cada afirmación con su respaldo y su límite. Sin fila de respaldo, no entra (CHECKLIST_AA.md, 6.1).

| # | Afirmación | Respaldo | Límite que hay que declarar |
|---|---|---|---|
| 1 | Las formas afines de |V| se obtienen de una factorización LU dispersa en el punto AC y escalan a 9.241 barras | C7, `results/c7_timing.csv` | El bloque de sensibilidades es denso: el método completo escala como N^1,5, no N^1,2 |
| 2 | El error de linealización es lo único que separa las cotas del 0,27 % nominal en el punto de diseño | B1 y 27.1 (columnas lineal y AC) | En la peor barra de PEGASE |
| 3 | El rango del término cuadrático corrige ese error en todo el rango carga × δ; un ajuste de media y varianza no basta | 27.2 y 28.3 | Validado por Monte Carlo (20.000 muestras por celda), en dos sistemas; el rango es muestreado, no acotado |
| 4 | Las caídas por rama desde la forma afín de la diferencia son varias veces más angostas que restar intervalos de barra, con la misma cobertura | 28.1 y 28.3 | Solo caídas de magnitud; corrientes y flujos no están medidos |
| 5 | Un radio que usa el soporte acotado (Bernstein o Hoeffding) cubre sin suponer normalidad; los que usan solo la varianza equivalen al peor caso | 28.4 | Supone símbolos independientes con soporte [−1, 1] |
| 5b | Fuentes correlacionadas: escritas como símbolos compartidos, el mismo método cubre; suponer independencia o un modelo mal especificado no | 28.6 | La ganancia sobre el peor caso cae con un factor común fuerte; la estructura de correlación es un supuesto del modelo; con la slack sola y correlación fuerte falla la linealización en PEGASE |
| 6 | 3σ, la reformulación gaussiana habitual, no cubre por objetivo aun corregida | 27.2, 28.3, 28.7 | Con Bonferroni sí cubre en la práctica, sin garantía |
| 6b | La cobertura simultánea con un radio sin distribución cuesta lo mismo que el peor caso; la ventaja de la concentración es por objetivo y depende del número de fuentes por barra | 28.7 | En IEEE-118 la concentración no aporta |
| 7 | La cota de peor caso de segundo orden no encierra la solución AC en los vértices de una red de transmisión | 28.2 | Es un resultado negativo; en IEEE-33 el exceso es despreciable |
| 8 | El método completo es unas 15 veces más rápido que un Monte Carlo AC por lotes en 2.869 barras | 28.5 y C8 | Depende de la implementación del Monte Carlo; 250 muestras probadas solo en el punto de diseño |

Lo que este informe no respalda todavía: una garantía simultánea más angosta que el peor caso, correlación estimada de datos
reales de carga y generación, la corrección en otros sistemas con el barrido completo de C1 (B7 usa una grilla de 12 celdas), y una revisión
bibliográfica en bases indexadas (RELATED_WORK.md, quinta pasada).

