"""Commentary that closes RESULTS.md: readings of C1/C7/C8, changed figures, blocked items, status of each critique."""

C1_LECTURA = r"""- **(b) no cambia con la carga ni con δ:** la miscobertura de la elipse en el modelo lineal es 0,07–0,11 % en todas
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
"""

FINAL_TEXT = r"""---

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
| Tiempo CD-HIM, PEGASE | 52,05 s | 58–66 s medidos aquí; solo el lazo por ramas, LU denso | **<<T_PEG>> s** en un hilo (Jacobiano + LU + resoluciones + ramas); <<T_PEG_THR>> s con las resoluciones en hilos | Álgebra dispersa y vectorización |
| Tiempo CD-HIM, IEEE-33 | 0,0037 s | 0,0047 s | **<<T_33>> s** | Ídem (en 33 barras domina el armado del Jacobiano) |
| Monte Carlo AC, 20.000 flujos, PEGASE | 1.246,66 s (y 1.263,48 s) | No terminó en la corrida de referencia (máquina compartida) | **<<T_MC_PEG>> s** por lotes; lazo pandapower: <<PP_PEG>> ms por flujo medidos, <<PP_PEG_TOT>> s extrapolados | Solver por lotes; el tiempo de pandapower depende de la versión |
| Speed-up, PEGASE | 24,0× (24,8× en IV-D) | — | **<<SP_PEG>>×** frente al Monte Carlo por lotes; <<SP_PEG_PP>>× frente al lazo pandapower extrapolado | La cifra depende de cómo se implemente el Monte Carlo |
| Speed-up, IEEE-33 | 1,08·10⁵ (extrapolado) | — | **<<SP_33>>×** frente al Monte Carlo por lotes | Ídem: el 10⁵ compara contra el sobrecosto de pandapower en una red de 33 barras |
| "Sub-second" | Afirmado con 52,05 s | — | <<T_PEG>> s en un hilo, <<T_PEG_THR>> s con hilos, más <<T_BASE_PEG>> s del flujo base | — |
| CLTM: 12,5 s, ≈ 100×, sobrecosto 4,2× | Tabla II | Sin código | La caja 3σ cuesta lo mismo que CD-HIM menos la etapa por ramas (<<BR_SHARE>> % del total en PEGASE) | No hay implementación de CLTM en el repositorio |
| Complejidad | O(N^α), α ≈ 1,2 | LU denso, O(N³) | LU: N^<<K_LU>>; relleno: N^<<K_NNZ>>; **método completo: N^<<K_TOT>>** (IC 95 % <<K_TOT_CI>>), dominado por las resoluciones, N^<<K_SOLVE>> | El bloque de sensibilidades es denso (C7) |
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
"""

C7_LECTURA = r"""- **La factorización es barata y casi lineal.** El tiempo de LU crece como N^<<K_LU>> (IC 95 % <<K_LU_CI>>) y el
  relleno, nnz(L+U), como N^<<K_NNZ>> (IC <<K_NNZ_CI>>). El cociente (L+U)/J va de 1,6 a 3,0 con COLAMD y de 1,35 a
  1,93 con MMD_AT_PLUS_A, que deja cerca de un tercio menos de no-ceros en los sistemas grandes con un tiempo de
  factorización similar.
- **El costo del método no está en la LU sino en extraer las sensibilidades.** Las N_u resoluciones triangulares
  son el <<SH_300>> % del total con 300 barras, el <<SH_2869>> % con 2.869 y el <<SH_9241>> % con 9.241. Crecen como
  N^<<K_SOLVE>> (IC <<K_SOLVE_CI>>) porque el resultado es denso: el bloque de |V| tiene n_PQ × N_u entradas
  (576 MB con 9.241 barras). El revisor 3 tenía razón en el punto 8.
- **Exponente del método completo: <<K_TOT>> con los seis sistemas (IC 95 % <<K_TOT_CI>>) y <<K_TOT5>> sin el de
  33 barras (IC <<K_TOT5_CI>>).** Con cinco o seis puntos los intervalos son anchos. Lo que sí se puede afirmar:
  α ≈ 1,2 describe el relleno de la factorización, no el método, y el comportamiento asintótico es cuadrático
  mientras se necesite el bloque denso.
- En valores absolutos: <<T7_2869>> s para 2.869 barras y <<T7_9241>> s para 9.241 en un hilo. Con las resoluciones
  repartidas en hilos (SuperLU libera el GIL) esa etapa baja a <<T7_2869_THR>> s y <<T7_9241_THR>> s.
- El flujo de potencia base (pandapower, ya en caliente) cuesta entre 0,04 y 0,4 s y hay que sumarlo: el método
  parte de un punto AC convergido.
- No hice el ajuste por separado para el ordenamiento alternativo ni probé más allá de 9.241 barras.
"""

C8_LECTURA = r"""- **Ancho frente al rango empírico AC al 99,73 % (mediana sobre barras):** caja afín 3,1× y 4,1×; caja 3σ 1,05× y
  1,06×; CD-HIM con λ_max 1,64× y 1,61×; CD-HIM con proyección exacta 1,23× en ambos sistemas.
- **La caja 3σ reproduce casi exactamente el rango empírico**, con violación media 0,15 % y 0,13 %. Su problema
  es la cola: 80 barras sobre 0,27 % en PEGASE y 42 % de las muestras con alguna barra fuera.
- **CD-HIM es más ancho que la caja 3σ y viola menos.** No es que la elipse aproveche la correlación para estrechar
  los intervalos por barra: el cuantil conjunto (3,44σ en vez de 3σ) y λ_max los ensanchan. La correlación solo
  se aprovecha si lo que se entrega es el conjunto bidimensional (por ejemplo, para acotar una caída de voltaje o
  un flujo por la rama), no dos intervalos.
- **Todos los métodos lineales cuestan lo mismo**, porque comparten las sensibilidades: <<T_PEG>> s para CD-HIM en
  PEGASE, de los cuales la etapa propia (ramas) es el <<BR_SHARE>> %. La Tabla II del paper (CLTM ≈ 100× más
  rápido que el Monte Carlo y CD-HIM 4,2× más lento que CLTM) no se puede reproducir.
- **El speed-up depende del Monte Carlo con que se compare.** Frente al solver AC por lotes de este paquete:
  <<SP_33>>× (IEEE-33) y <<SP_PEG>>× (PEGASE). Frente al lazo `pandapower.runpp(numba=False)` del script original,
  que aquí tarda <<PP_33>> y <<PP_PEG>> ms por flujo: <<SP_33_PP>>× y <<SP_PEG_PP>>×, extrapolando desde 200 flujos.
  El paper midió 18,9 y 62,3 ms por flujo; no pude reproducir esos tiempos con pandapower 3.5.6. Recomiendo
  reportar tiempos absolutos y decir qué implementación de Monte Carlo se usó.
- Las muestras AC por lotes son soluciones exactas del mismo problema (desajuste < 10⁻⁹ p.u.); una prueba unitaria
  las compara con pandapower.
"""
