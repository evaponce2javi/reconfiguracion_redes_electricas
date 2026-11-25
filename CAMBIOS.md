# CAMBIOS.md — Qué se cambió en CD-HIM y por qué

Este documento cubre todo lo que hay entre el commit `4750ce2` ("10k samples"), cuando el repositorio tenía solo
`README.md` y `finalFinalCode.py`, y el commit que lo acompaña. Para cada cambio dice qué se hizo y por qué.
Las cifras completas están en [RESULTS.md](RESULTS.md); aquí van solo las necesarias para entender la razón.
Lo que falta está al final, en la [sección 11](#11-lo-que-falta).

El método al que se llegó está explicado aparte, en [METODO.md](METODO.md).

---

## 1. Resumen

El repositorio pasó de un script de 1.021 líneas a:

- un paquete, `cdhim/`, con 22 módulos;
- 25 scripts de experimento que regeneran 56 tablas CSV, cada una con su archivo de metadatos (semilla, versiones,
  script);
- 315 pruebas unitarias;
- cinco documentos: `AUDIT.md`, `RESULTS.md`, `RELATED_WORK.md`, `CHECKLIST_AA.md` y el `README.md` ampliado.

El trabajo tuvo cinco etapas. Cada una se abrió por una razón distinta:

| Etapa | Qué se hizo | Por qué | Dónde quedó |
|---|---|---|---|
| 1. Auditoría | Lectura completa de `finalFinalCode.py`, sin modificarlo | Saber qué hace de verdad el código antes de tocarlo; era la primera regla del encargo | [AUDIT.md](AUDIT.md) |
| 2. Paquete nuevo y experimentos de la revisión | Reescritura corregida del método y los experimentos B1–B3 y C1–C8 | Responder con evidencia cada crítica del editor y de los tres revisores | RESULTS.md, secciones 1 a 16 |
| 3. Exploración de aritmética afín | Modelo de segundo orden, RAA-K, cinco aportes posibles puestos a prueba (X1–X7) | El método original no sostenía sus afirmaciones; el grupo trabaja en aritmética afín y se buscó un aporte defendible por ahí | Secciones 17 a 24 |
| 4. Opción A: encierros en distribución | Encierro de voltajes frente a los dos métodos publicados que dan garantía (X8) | De los cinco aportes posibles, era el único con una ruta clara | Sección 25 |
| 5. Las comprobaciones que faltaban | Krawczyk, redondeo dirigido, referencia certificada, escala, barras PV, redes trifásicas (X9–X14) | Eran las cinco cosas que podían tumbar o sostener la opción A | Sección 26 |

---

## 2. Lo que no se tocó, y por qué

| Qué | Por qué |
|---|---|
| `finalFinalCode.py` | Es el código del manuscrito enviado. Se dejó intacto para poder reproducir sus cifras; `python -m cdhim ieee33 --legacy` hace lo mismo que él desde el paquete nuevo. |
| `caja-elipse.tex` | No está en el repositorio y la regla era no editarlo: las cifras nuevas van a `RESULTS.md` y los autores las pasan al paper. |
| `ibex-lib` (código RAA-K del grupo) | Solo se leyó. Lo que se encontró ahí (sección 6.4) quedó anotado para que el grupo decida. |
| Las cifras del paper | Cuando un resultado nuevo contradice una cifra del paper no se ajustó nada para que calzara: se reportan las dos y la causa (RESULTS.md, sección 14). |
| El repositorio remoto | El commit es local. No se hizo `push`. |

---

## 3. Etapa 1: auditoría del código original

**Por qué:** las reglas del encargo pedían auditar antes de modificar, y la revisión del paper cuestionaba cosas
que solo se podían contestar mirando el código.

Los ocho hallazgos principales, y qué se hizo con cada uno en la etapa 2:

| # | Hallazgo | Consecuencia | Qué se hizo |
|---|---|---|---|
| 1 | El jacobiano tiene un error de signo en ∂P/∂θ fuera de la diagonal | En IEEE-33 el 27,5 % de las sensibilidades ∂V/∂P sale con el signo cambiado; todas las cifras del paper se calculan sobre esa matriz | Jacobiano exacto, comprobado contra diferencias finitas |
| 2 | El jacobiano se evalúa en arranque plano, no en el punto convergido | Contradice lo que dice el manuscrito | Se linealiza en la solución AC convergida |
| 3 | γ_th = 15 % no hace nada: G ≥ 21,46 % para cualquier ρ | La decisión Elipse/Caja la toma solo la guarda \|ρ\| < 0,999 | Documentado; regla alternativa estudiada en C4 y C6 |
| 4 | La validación contra AC no existe en el código | El "R_viol" del paper no tiene respaldo | Validación Monte Carlo contra AC, por barra, por rama y por sistema |
| 5 | Los intervalos de salida del Algoritmo 1 no se calculan | El CRI reportado es de un objeto que el método no entrega | Intervalos de salida implementados, con dos radios posibles |
| 6 | La distancia de Mahalanobis lleva una regularización de 10⁻⁹ | Engorda la elipse y baja las tasas de violación | Eliminada; piso de ruido numérico en su lugar |
| 7 | El álgebra lineal es densa | Contradice la contribución de LU disperso | Jacobiano y LU dispersos (SuperLU) |
| 8 | El experimento Beta se corre sobre ramas que el algoritmo manda a Caja | El 0,19 % del paper no mide ramas en modo Elipse | Rehecho sobre todas las ramas en modo Elipse |

La auditoría lista además 16 defectos de lógica, 7 de estructura y 7 de rendimiento, con la línea de cada uno, y
contrasta cada crítica de la revisión con el código: todas se confirman. Los hallazgos 1, 2, 3, 6 y 8 no los
podía ver la revisión, porque no tenía el código.

---

## 4. Etapa 2: paquete nuevo y experimentos de la revisión

### 4.1 Cambios de código

| Cambio | Por qué | Dónde |
|---|---|---|
| Paquete importable en vez de un script que ejecuta todo al importarlo | El original no tenía `if __name__ == "__main__"`, instalaba paquetes con `os.system` y no guardaba ningún resultado | `cdhim/` |
| Dos linealizaciones seleccionables: *corrected* (por defecto) y *legacy* | Poder reproducir las cifras del manuscrito y, a la vez, medir lo que el manuscrito dice hacer | `cdhim/study.py`, `cdhim/powerflow.py` |
| Jacobiano disperso exacto, LU disperso, sensibilidades con una resolución de varios lados derechos | Hallazgos 1, 2 y 7; crítica R3.8 | `cdhim/powerflow.py`, `cdhim/study.py` |
| Flujo AC por lotes reutilizando la matriz de admitancias y la factorización | El Monte Carlo original llamaba a pandapower por muestra (20 a 60 ms por flujo, casi todo sobrecosto) | `cdhim/powerflow.py` |
| Las tres "cajas" con nombre propio (covarianza, 3σ, afín), κ y G como funciones puras | Crítica R3.2: el manuscrito mezclaba tres cajas distintas | `cdhim/geometry.py` |
| Intervalos de salida del Algoritmo 1, con radio por λ_max o por proyección exacta | Hallazgo 5; críticas R2.3 y R3.3 | `cdhim/core.py` |
| Covarianza general de las fuentes, Q = C·Σ_ε·Cᵀ, y cópula gaussiana | Críticas R2.2 y R3.4: el original supone fuentes independientes | `cdhim/uncertainty.py` |
| Validación contra AC por barra, por rama y a nivel sistema | Hallazgo 4; críticas del editor, R1.1, R2.4, R2.5, R3.1 y R3.6 | `cdhim/validation.py` |
| Carga de casos de pandapower, PGLib y MATPOWER | Sistemas de prueba más allá de los dos del manuscrito | `cdhim/network.py` |
| Todos los parámetros en un archivo | En el original δ, γ_th, las semillas y los tamaños de muestra estaban repartidos y repetidos | `config.yaml`, `cdhim/config.py` |

### 4.2 Experimentos

| Experimento | Crítica que responde | Qué mostró |
|---|---|---|
| A, reproducción | — | El código del repositorio da CRI de 86,70 % y 57,80 %; el paper dice 86,92 % y 58,18 %. No se reproduce exacto |
| B1, violación contra AC | Editor, R1.1, R2.5, R3.1 | En PEGASE las cotas del manuscrito no contienen las muestras AC (media 3,69 %, máximo 89,4 %); las corregidas sí en promedio (0,0024 %) y fallan en la peor barra (0,41 %) |
| B2, conservadurismo | R2.3, R3.3, R3.7 | La salida real se reduce 42,8 % y 34,1 % en ancho frente a la caja afín; frente a la caja 3σ es más ancha |
| B3, Tabla I | R3.1 | La tabla mide la elipse en el modelo lineal, no la linealización; sin la regularización da lo mismo para todo δ |
| C1, carga × incertidumbre | R1.4 | 25 celdas por sistema, hasta el 95 % del límite de cargabilidad |
| C2, fuentes correlacionadas | R2.2, R3.4 | Con fuentes correlacionadas y la covarianza del manuscrito la violación es de 26 a 59 % |
| C3, distribuciones no gaussianas | R3.5, R2.6 | Doce distribuciones, todas las ramas: la robustez afirmada no se sostiene |
| C4, regla de ρ frente a G | R2.1, R3.2 | Con la caja del manuscrito el umbral no actúa |
| C5, nivel sistema | R2.4, R3.6 | En PEGASE, 3 % de las muestras viola alguna barra y 50 % alguna elipse; Bonferroni cuantificado |
| C6, transferibilidad de γ_th | R2.7 | Cinco sistemas: o no actúa, o depende del sistema |
| C7, escalabilidad | R3.8 | Seis sistemas; el método completo escala como N^1,51, no N^1,2 |
| C8, otros métodos | R1.3 | Solo caja afín y caja 3σ; las referencias [11] y [12] no están implementadas |

Las cifras del paper que cambian, con la causa de cada diferencia, están en la sección 14 de RESULTS.md.

### 4.3 Reproducibilidad

| Cambio | Por qué |
|---|---|
| Cada experimento escribe un CSV en `results/` y un `.meta.json` con semilla, tamaño de muestra, script y versiones | Regla del encargo; el original solo imprimía en pantalla o en el título de una figura |
| Semilla base en `config.yaml`; cada experimento deriva la suya | Regla del encargo |
| `experiments/run_all.py` regenera todos los CSV y `experiments/build_report.py` arma `RESULTS.md` desde ellos | Ninguna cifra del informe está escrita a mano |
| `requirements.txt` con versiones fijas | Regla del encargo: registrar Python, pandapower, NumPy y SciPy |
| `tests/` con pruebas contra diferencias finitas, fuerza bruta y casos publicados | Que un error como el del signo del jacobiano no pueda volver a pasar inadvertido |

---

## 5. Etapa 3: exploración de aritmética afín

**Por qué:** después de la etapa 2 el método original quedaba sin sus afirmaciones principales. El grupo investiga
aritmética afín (RAA, NPAA, teorema k-nomial), así que se buscó un aporte por ese lado, con una lista de
comprobaciones propia ([CHECKLIST_AA.md](CHECKLIST_AA.md)) y una revisión bibliográfica
([RELATED_WORK.md](RELATED_WORK.md)).

| Cambio | Por qué | Dónde |
|---|---|---|
| Modelo de segundo orden de \|V\|: una forma cuadrática exacta por barra | El resto del modelo lineal es de segundo orden; sin él no hay cómo acotarlo | `cdhim/second_order.py`, X1 |
| Fuentes agregadas por zonas | Con miles de símbolos independientes ninguna cota cierra; pocas fuentes correlacionadas es el caso natural de la aritmética afín | `cdhim/zones.py`, X2 |
| Port a Python del operador RAA-K del grupo | Poder usarlo dentro del método; validado bit a bit contra los 80 casos de prueba de `ibex-lib` | `cdhim/knomial.py` |
| Comprobación del peor escenario por barra | Saber si el segundo orden encuentra los extremos | X3 |
| Controles independientes de la referencia | Descartar que el contraste AC fuera el que estaba mal | X4 |
| Reimplementación de Liao et al. (2017) | Es el método publicado más cercano | X5 |
| Punto fijo sobre el resto del modelo | Convertir una cota aproximada en un encierro | `cdhim/enclosure.py`, X6 |
| RAA-K frente a la relajación semidefinida | Medir al operador contra la cota conocida más ajustada | X7 |

El veredicto de los cinco aportes posibles (sección 24):

| Aporte | Veredicto |
|---|---|
| 1. Encierro verificado por punto fijo | El mecanismo no es nuevo y no escalaba; posible solo en alimentadores |
| 2. Ni el primer ni el segundo orden encierran la solución | Se sostiene, pero es un resultado negativo |
| 3. RAA-K para acotar el término cuadrático | Parcial: no es la cota más ajustada |
| 4. Formulación adjunta | Parcial |
| 5. Peor escenario guiado por el segundo orden | Parcial; en redes radiales ya está resuelto |

De ahí salió la opción A: llevar el aporte 1 a redes de distribución.

---

## 6. Etapa 4: opción A, encierros en redes de distribución

**Por qué:** para saber si el encierro era un aporte había que medirlo contra lo publicado que da garantía, en
las mismas redes y con el mismo conjunto de incertidumbre.

### 6.1 Cambios

| Cambio | Por qué | Dónde |
|---|---|---|
| Reimplementación de Sou, Khong y Sandberg (2026) | Es el método publicado para redes radiales; se validó contra su Tabla I | `cdhim/distflow.py` |
| Discos de Wang, Bernstein, Le Boudec y Paolone (2018) | Es la garantía publicada para cualquier topología; reproduce el ejemplo numérico del artículo | `cdhim/discs.py` |
| Incertidumbre en semidiscos (generación con inversor) | Es el conjunto de Sou et al.; el primer orden se maximiza exacto sobre el semidisco | `cdhim/halfdisk.py` |
| Cota del resto en forma de inyección de corriente | La cota anterior, en rectangulares, dejaba de cerrar con ramas muy cortas: cerraba en 27 de 80 casos y la nueva en 78 | `cdhim/enclosure.py` |
| Los solucionadores de punto fijo inflan el resultado y comprueban el auto-mapeo | Aceptaban un iterado que se acercaba desde abajo, lo que no prueba nada | `cdhim/enclosure.py` |
| Conversión de unidades de los alimentadores de MATPOWER | Sus archivos traen ohms y kW y los convierten con sentencias MATLAB al final, que un lector simple se salta | `cdhim/network.py` |
| Tope en el método de Sou et al. | Con el arranque de 5.000 A hay desbordamiento y un máximo que se tragaba un NaN reportaba convergencia falsa | `cdhim/distflow.py` |

### 6.2 Lo que mostró

- **En redes radiales no hay aporte:** el método de Sou et al. es más ajustado con incertidumbre grande, da el
  mínimo exacto y es mucho más barato.
- **En redes malladas o con elementos en derivación sí:** el encierro cerraba en 38 de 40 casos a 1,001–1,305
  veces el rango alcanzado, donde el método de Sou et al. no está definido y los discos quedan a 1,24–2,17.
- **La cota de segundo orden sola no contiene la solución** en ninguno de los 80 casos: el resto hay que acotarlo.

### 6.3 Hallazgo: la monotonía no se cumple en todas las redes

El método de Sou et al. toma el mínimo de voltaje de un solo flujo de carga, suponiendo monotonía. En las redes
CIGRE MV y Oberrhein eso falla en 20 de 80 casos. Se agregó una tabla que lo mide (`results/x8_monotonicity.csv`)
porque es un argumento a favor de una cota inferior que no dependa de esa suposición.

### 6.4 Hallazgo: el operador RAA-K depende de la escala de la forma

Al perfilar la red de 179 barras apareció que la búsqueda de bloques usa umbrales absolutos de 10⁻¹². El que
actúa es el del producto de dos entradas de la diagonal: con formas de coeficientes chicos, como las del flujo de
potencia en p.u., el operador deja de formar bloques y entrega casi la cota término a término. La cota sigue
siendo válida, pero más ancha de lo que el operador puede dar.

| Cambio | Por qué |
|---|---|
| `knomial_bound(..., normalise=True)`: divide la forma por su mayor coeficiente antes de buscar bloques | Restablece la invariancia de escala que la teoría exige |
| Aviso en las secciones 18 a 23 de RESULTS.md | Esas secciones se corrieron sin normalizar y subestiman a RAA-K; no se rehicieron (punto 9.14 del checklist) |
| Nota de vocabulario en la sección 25.8 | Lo que este informe llama "AF2" es lo que el grupo llama RAA con k = 0; su "fAF2" es otra cota |

El mismo umbral está en el código C++ del grupo (`ibex_KNomial.h`, línea 156). No se modificó.

---

## 7. Etapa 5: las comprobaciones que faltaban

**Por qué:** la sección 25.7 dejó una lista de lo que faltaba para sostener la opción A. Se hicieron todas salvo
una, que quedó bloqueada.

| Comprobación | Qué se cambió | Por qué | Resultado |
|---|---|---|---|
| Krawczyk | Operador de Krawczyk paramétrico, en dos formulaciones y cuatro variantes (`cdhim/krawczyk.py`, X9) | Era la comparación que podía mostrar que el aporte ya existía; se armó lo más fuerte posible | Por inyección de corriente cierra en los mismos 78 de 80 casos y no es más ajustado en ninguna barra; por desbalance de potencia no cierra en ninguno |
| Leer Wang, Le Boudec y Paolone (2019) | Cuarta pasada de literatura en RELATED_WORK.md; subdivisión del conjunto con tres encierros distintos (X12) | Su método era el vecino más cercano en redes malladas | **Bloqueado: falta el PDF.** Por fuentes indirectas, su cota es la que ya estaba y su vía al ajuste es subdividir; eso se midió |
| Redondeo dirigido | Verificación a posteriori con existencia y unicidad (`cdhim/verified.py`) | Sin controlar el redondeo no se puede escribir "verificado" | Verifica en 78 de 80 casos sin perder ancho |
| Referencia certificada | Puntos interiores verificados y rango exacto por bisección (`cdhim/certified.py`) | Hasta entonces el ancho se comparaba con valores hallados por búsqueda, que no prueban nada | Cota interior certificada en todos los casos; rango exacto en 38 cajas |
| Escala | Redes de 118 a 2.940 barras con tres diseños de fuentes (X13) | Todo estaba probado hasta 179 barras | Hasta 2.940 barras con 64 símbolos; el límite es el número de símbolos |
| Barras PV | Segunda forma del operador, con voltaje y corriente (`cdhim/verified_vi.py`, X14) | La forma de inyección de corriente no admite barras de voltaje controlado | Verifica en distribución y en la red de 30 barras; no en las de 118 y 300 |
| Redes desbalanceadas | Modelo en el dominio de fases (`cdhim/threephase.py`, X10) | La distribución real es desbalanceada | Funciona sin cambios con cargas en estrella; validado contra pandapower |

### 7.1 Cambios que salieron de esas comprobaciones

| Cambio | Por qué | Efecto |
|---|---|---|
| **Cota conjunta del primer y el segundo orden** | Al comparar con Krawczyk se vio que casi todo el exceso venía de sumar dos máximos que no ocurren en el mismo punto | El ancho sobre el rango exacto bajó de 1,0004–1,160 a 1,0000–1,120; la cota superior sale de ella en todas las barras |
| Normas del modelo más finas (partes real e imaginaria, ejes a lo largo y a través del voltaje) | Sumar módulos sobreestimaba; por ejemplo 1,5 C en vez de 1,118 C para un par P/Q | Tubo entre 0,49 y 0,84 veces la cota del resto anterior |
| La misma norma en la forma voltaje-corriente | Lo afiné en una forma y no en la otra; una prueba lo detectó | Las dos formas vuelven a dar el mismo radio |
| Punto de partición compartido al certificar muchos puntos | Certificar mil puntos por caso en la red trifásica tomaba horas; el punto de partición puede ser cualquiera, así que la parte matricial se calcula una vez | De 875 a 50 ms por punto, con la misma precisión (medido aparte, sin CSV) |
| `attained_points`: la búsqueda local devuelve también el punto que alcanza cada extremo | Para certificar un punto hay que saber cuál es | — |
| El experimento X12 reemplaza al encierro de la sección 25 como método | Es el mismo modelo, verificado y más ajustado | Las cifras de la sección 25 quedan como históricas |

### 7.2 Lo que cambió en la lectura del aporte

- **Frente a Krawczyk**, el intervalo de este es entre 1,003 y 1,34 veces el del encierro, y su exceso sobre el
  rango exacto entre 4 y 225 veces mayor. El encierro es más caro: 81 ms de mediana contra 12 ms.
- **Ya no vale la comparación con los discos como argumento principal.** "La garantía publicada queda a 1,24–2,17
  veces el rango" es cierto de los discos; un Krawczyk bien armado queda a 1,003–1,50.
- **RAA-K deja de aportar.** Con la cota conjunta, en los 80 casos de escala el ancho es idéntico con RAA-K
  completa que acotando término a término, a un costo de hasta 216 s. Sin él, el caso de Oberrhein baja de 13,3 s
  a 0,6 s. Esto pesa en la decisión de la sección 11.

No existe un experimento X11: la numeración saltó de X10 a X12.

---

## 8. Cambios en los documentos

| Documento | Qué cambió | Por qué |
|---|---|---|
| `README.md` | De dos líneas a la guía de uso, la tabla de módulos y los comandos para regenerar todo | El original no decía cómo correr nada |
| `AUDIT.md` (nuevo) | La auditoría de la etapa 1 | Entregable de la Parte A |
| `RESULTS.md` (nuevo) | 26 secciones generadas desde los CSV | Entregable principal: de aquí salen las cifras para el paper |
| `RELATED_WORK.md` (nuevo) | Cuatro pasadas de literatura, con el nivel de lectura de cada fuente | Saber qué existe antes de reclamar novedad; cada pasada respondió a una pregunta distinta |
| `CHECKLIST_AA.md` (nuevo) | Lista de comprobaciones de la versión de aritmética afín, con el estado de cada punto | Separar lo hecho, lo bloqueado y lo que decide el grupo |
| `CAMBIOS.md` y `METODO.md` (nuevos) | Este documento y la explicación del método | Pedidos al cerrar esta etapa |

---

## 9. Inventario

**Paquete `cdhim/` (22 módulos)**

| Etapa | Módulos |
|---|---|
| 2 | `__init__`, `__main__`, `config`, `core`, `geometry`, `powerflow`, `study`, `uncertainty`, `validation`, `network` |
| 3 | `zones`, `second_order`, `knomial`, `enclosure` |
| 4 | `halfdisk`, `distflow`, `discs` |
| 5 | `verified`, `verified_vi`, `certified`, `krawczyk`, `threephase` |

**Experimentos (`experiments/`)**

| Grupo | Scripts |
|---|---|
| Revisión del paper | `a_reproduce_paper`, `b1` a `b3`, `c1` a `c8` |
| Aritmética afín | `x1_aa_second_order` a `x7_block_choices` |
| Opción A | `x8_feeders` |
| Comprobaciones | `x9_krawczyk`, `x10_three_phase`, `x12_verified`, `x13_scale`, `x14_pv` |
| Informe | `run_all`, `build_report`, `report_x12`, `figures`, `common` y los `results_text_*` |

**Pruebas (`tests/`):** 16 archivos, 315 pruebas. **Resultados (`results/`):** 56 CSV con sus 56 `.meta.json` y 5
figuras. Los casos de MATPOWER y PGLib se descargan a `data/` al correr y no se versionan.

---

## 10. Cifras del paper que cambian

El detalle está en la sección 14 de RESULTS.md. Las principales:

| Cifra | Paper | Código corregido | Causa |
|---|---|---|---|
| CRI, PEGASE | 58,18 % | 71,08 % (57,80 % con el código original) | Jacobiano; el 58,18 % no se reproduce |
| CRI, IEEE-33 | 86,92 % | 86,37 % (86,70 % con el original) | Ídem |
| Ramas en modo Elipse | 51,6 % y 92,7 % | 67,7 % y 79,07 % | Jacobiano |
| Violación Beta, PEGASE | 0,19 % | 0,83 % de media | Regularización, selección de ramas y transformación |
| Tiempo, PEGASE | 52,05 s | 0,92 s | Álgebra dispersa |
| Speed-up, PEGASE | 24,0× | 65× frente a un Monte Carlo por lotes | Depende de cómo se implemente el Monte Carlo |

---

## 11. Lo que falta

### 11.1 Lo que depende de los autores

1. **Mandar el PDF del artículo de Wang, Le Boudec y Paolone:** "Controlling the Electrical State via Uncertain
   Power Injections in Three-Phase Distribution Networks", IEEE Transactions on Smart Grid 10(2), 2019. No tiene
   preprint y está tras un muro de pago. Es el vecino más cercano en redes malladas y hoy se conoce solo por su
   resumen, unas láminas y la tesis de C. Wang. Sin leerlo no se puede descartar que traiga una cota más fina que
   la que se usó como referencia.
2. **Decidir si el paper va por los encierros o por el operador.**
   - *Por los encierros verificados en redes de distribución:* el caso es fuerte. Hay garantía con el redondeo
     controlado, sobreestimación demostrada y ventaja medida sobre Krawczyk. El costo es que el operador
     k-nomial del grupo no sostiene el resultado.
   - *Por el operador RAA-K:* hoy es débil, porque con la cota conjunta el operador no mejora el ancho. Depende
     de probar si sus bloques pueden usarse dentro de la cota conjunta en lugar de los autovectores (punto 9.23
     del checklist), que no se ha hecho.
3. **Decisión 0.1 del checklist, tipo de garantía.** En la práctica ya está tomada al elegir la opción A:
   encierro verificado, restringido a distribución. Falta escribir en una frase qué afirma el paper y qué no.
4. **Qué conjunto de incertidumbre usa el paper** (punto 9.13): semidiscos con demanda determinista, que es el
   caso de Sou et al., o demanda incierta por zonas, que es el caso natural de la aritmética afín.
5. **Avisar al grupo de la tolerancia absoluta en `ibex_KNomial.h`** (punto 9.15). Afecta a cualquiera que use
   RAA-K con cajas chicas, no solo a este proyecto.
6. **Dar el visto bueno para rehacer las secciones 18 a 23** con RAA-K normalizada (punto 9.14): son varias horas
   de cómputo.
7. **Lo bloqueado de la revisión original** (sección 15 de RESULTS.md): qué serie de datos reales usar, el PDF de
   la Ref. [12], y la decisión de implementar o no las referencias [10] y [11].

### 11.2 Trabajo técnico pendiente

| Pendiente | Por qué importa |
|---|---|
| Segunda implementación del redondeo con una biblioteca de intervalos certificada (9.19) | Hoy es disciplina escrita a mano, contrastada con aritmética racional exacta en las pruebas; es lo que pediría un revisor de computación verificada |
| Revisar si la cota conjunta ya existe con otro nombre (9.20) | Es una idea corriente en optimización global; no hay que presentarla como nueva sin esa revisión |
| Bloques de RAA-K dentro de la cota conjunta (9.23) | Es la forma de que el teorema k-nomial siga dentro del método |
| Cargas en delta, transformadores y neutro explícito | El modelo trifásico cubre solo cargas en estrella |
| Una cota que cierre con barras PV a mayor escala (9.18) | En la red de 118 barras la parte lineal contrae y lo que falla es la cota; en la de 300, ni la parte lineal |
| Más símbolos | Con 586 símbolos un solo tensor ocupa 15 GB |
| El caso CIGRE MV con ρU = 80 % | Es el único de los 80 en que el encierro no verifica |
| Detalles de carga de redes (9.22) | `load_grid` supone la tabla de barras ordenada; el lector depende de la codificación del sistema |
| Tablas VI y VII de Sou et al. para 33 y 69 barras | Se reproducen 3 de 8 filas; la causa de las otras no se identificó |
