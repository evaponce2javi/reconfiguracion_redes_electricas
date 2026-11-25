"""Commentary of sections 20 to 24 of RESULTS.md: critical tests of the five candidate contributions."""

X4_TEXT = r"""---

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

{{x4_independent}}

*Script: `experiments/x4_reference_checks.py` · CSV: `results/x4_independent.csv`*

- Los dos solvers dan el mismo voltaje en todas las barras: la mayor diferencia es <<X4_PP_MAX>> p.u. Las
  violaciones de las secciones 17 a 19 no son un artefacto del solver por lotes ni del armado del escenario.
- Límite del modelo, que también vale para pandapower tal como se usó: los generadores no tienen límites de
  potencia reactiva. Liao et al. sí conmutan el tipo de barra al violarlos.

### 20.2 ¿Está el extremo en un vértice?

X2 tomó como referencia el rango de la solución AC en los 2^K vértices. Aquí se minimiza y maximiza el modelo de
segundo orden sobre toda la caja (puntos interiores incluidos), se resuelve el flujo AC exacto en ese punto, y en
las diez barras donde más gana se afina con una búsqueda acotada sobre el flujo AC exacto.

{{x4_interior}}

*CSV: `results/x4_interior.csv`*

- **El voltaje mínimo está en un vértice.** Ningún punto interior baja de la enumeración en más de
  <<X4_LO_GAIN_MAX>> p.u. La referencia de la cota inferior, que es la que manda, es correcta.
- **El máximo no siempre.** En PEGASE, entre <<X4_HI_PCT_LO>> y <<X4_HI_PCT_HI>> % de las barras tienen un punto
  interior por encima de todos los vértices, por hasta <<X4_HI_GAIN_MAX>> p.u.; a lo más <<X4_HI_N4_MAX>> barras
  pasan de 10⁻⁴ p.u. Los anchos "sobre el rango de referencia" de las secciones 18 y 19 están apenas
  sobreestimados por el lado superior.
- La búsqueda sobre el flujo exacto agrega como máximo <<X4_POLISH_EXTRA>> p.u. a lo que ya encuentra el modelo
  de segundo orden: para ubicar el extremo, el modelo basta.

---

## 21. X5 — El método publicado más cercano: Liao et al. (2017) (aporte 4)

Se reimplementó la Ec. (26) del artículo: desarrollo de segundo orden de cada componente del estado (e, f),
acotado monomio a monomio. El artículo no dice cómo pasa de los intervalos de e y f al de |V|; se usó aritmética
de intervalos sobre √([e]² + [f]²). Se compara con la forma directa de este paquete, que proyecta la misma
información de segundo orden sobre |V| antes de acotar.

{{x5_summary}}

*Script: `experiments/x5_liao_baseline.py` · CSV: `results/x5_liao_summary.csv` · "ancho / rango": mediana sobre
barras del ancho del intervalo sobre el rango de la solución AC (todos los vértices en IEEE 9 y PEGASE; vértices
extremos y Monte Carlo en IEEE 57)*

- **Las dos derivaciones coinciden.** Los coeficientes de segundo orden de Liao (una resolución por cada par de
  variables), proyectados sobre |V|, son iguales a los de la vía adjunta (una resolución traspuesta por barra).
  El script lo comprueba en cada caso.
- **La Ec. (26) tampoco contiene la solución exacta.** En IEEE 9 con ±10 %, <<X5_STATE_OUT>> de las
  <<X5_STATE_N>> componentes de e y f tienen algún vértice fuera de su intervalo, por hasta <<X5_STATE_SHORT>>
  p.u. Es el tercer orden, en el sistema de prueba del propio artículo.
- **Pasar de e y f a |V| con intervalos no escala.** El intervalo queda a <<X5_LIAO_9>> veces el rango en 9
  barras, <<X5_LIAO_57>> en 57 y <<X5_LIAO_PEG_LO>> a <<X5_LIAO_PEG_HI>> en PEGASE, porque se pierde que e y f
  rotan juntas. La forma directa se mantiene entre <<X5_DIRECT_LO>> y <<X5_DIRECT_HI>>.
- **La vía adjunta no es más barata con pocas fuentes.** Con 8 zonas en PEGASE, la vía de Liao necesita
  <<X5_SOLVES_L8>> resoluciones lineales y la adjunta <<X5_SOLVES_D8>>. La adjunta solo gana cuando hay más
  fuentes que barras a vigilar. Lo que marca la diferencia es proyectar sobre |V| antes de acotar, y eso se puede
  hacer con cualquiera de las dos.

**¿Se reproduce la Tabla 1 del artículo?** No.

{{x5_table1}}

{{x5_variants}}

*CSV: `results/x5_liao_table1.csv`, `results/x5_liao_variants.csv`*

- Los voltajes nominales coinciden con cuatro decimales: es el mismo sistema.
- El artículo habla de 12 variables de perturbación sin listarlas. Con cinco lecturas distintas del conjunto
  incierto, el rango exacto queda a no menos de <<X5_VAR_BEST>> p.u. de su columna de Monte Carlo, que no depende
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

{{x6}}

*Script: `experiments/x6_fixed_point.py` · CSV: `results/x6_fixed_point.csv` · "variación real del Jacobiano":
mayor radio espectral de J⁻¹(J(x* + Cε) − J(x*)) en 64 vértices, sin ninguna cota*

- **Funciona en el alimentador radial.** En IEEE-33 con δ de 5 y 15 % hay punto fijo; la holgura añadida es
  <<X6_ALLOW_33_LO>> a <<X6_ALLOW_33_HI>> % del radio y el encierro queda a <<X6_W_33_LO>>–<<X6_W_33_HI>> veces el
  rango exacto.
- **Donde hay punto fijo, el encierro contiene.** En los <<X6_N_FOUND>> casos con punto fijo (de <<X6_N_CASES>>),
  el número de violaciones en los 256 vértices, sumando el resto del estado, la parte de tercer orden y |V|, es
  <<X6_VIOL>>. La cota de segundo orden sola fallaba en esos mismos casos.
- **En 118 barras solo cierra con δ = 2 %**, y ya suelta: la holgura es <<X6_ALLOW_118_LO>> a
  <<X6_ALLOW_118_HI>> % del radio y el encierro queda a <<X6_W_118_LO>>–<<X6_W_118_HI>> veces el rango. Es el mismo
  orden de caja que usa Zhu (1,5 a 2 % en 39 y 118 barras).
- **En 300 y 2.869 barras no cierra.** El radio de G en PEGASE con δ = 5 % es <<X6_PEG_G>>; en coordenadas
  polares baja a <<X6_PEG_GP>>, que sigue sobre 1.
- **El obstáculo es la cota componente a componente, no la física.** La variación real del Jacobiano en PEGASE
  con δ = 5 %, en polares, tiene radio espectral <<X6_PEG_TP>>: hay margen para un argumento de contracción,
  pero los valores absolutos elemento a elemento lo multiplican por <<X6_PEG_FACTOR>> o más.
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

{{x7}}

*Script: `experiments/x7_block_choices.py` · CSV: `results/x7_block_choices.csv` · necesita `cvxpy`*

- **La relajación semidefinida da prácticamente el rango exacto**: <<X7_SDP_LO>> a <<X7_SDP_HI>> veces en mediana
  (<<X7_SDP_MAX>> en la peor forma). RAA-K queda entre <<X7_RAAK_LO>> y <<X7_RAAK_HI>>, los bloques espectrales
  entre <<X7_SPEC_LO>> y <<X7_SPEC_HI>> y AF2 entre <<X7_AF2_LO>> y <<X7_AF2_HI>>.
- **RAA-K es mucho más rápida**: <<X7_RAAK_MS_LO>> a <<X7_RAAK_MS_HI>> ms por forma contra
  <<X7_SDP_MS_LO>> a <<X7_SDP_MS_HI>> ms. Con 16 símbolos y 2.300 barras, la relajación semidefinida son unos
  cinco minutos; con miles de símbolos no es viable.
- **RAA-K se puede hacer rigurosa con redondeo dirigido** (el grupo ya lo tiene en C++). Una cota semidefinida
  sale de un solver numérico y necesita un certificado dual verificado para ser rigurosa.
- **La cota híbrida no ayuda:** supera a la mejor de sus dos partes solo en <<X7_HYB_WINS>> de los
  <<X7_NCASES>> casos, y por menos de 1 %.

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
"""
