# Trabajos relacionados: primera pasada (Parte 1 de CHECKLIST_AA.md)

Estado: **primera pasada, no cierra la Puerta 1.** Se abrieron unas 35 fuentes. Cuatro se leyeron a texto
completo (Vaccaro 2010, Vaccaro y Cañizares 2017, Schmidt y Schegner 2020 y Liao et al. 2017, este último con el
PDF que entregó la autora); del resto se vio el resumen o los metadatos, porque IEEE, Wiley, Springer y MDPI no
dieron acceso. Las citas de Vaccaro (2010) están contrastadas
contra el PDF. "n.v." significa no verificado: ahí hay que leer el artículo antes de citarlo.

## Las dos preguntas de la Puerta 1

**Q1. ¿Alguien acotó ya el resto cuadrático del flujo rectangular con aritmética afín, o con una forma
cuadrática que la extienda? → En parte sí: Liao et al. (2017) publican la expansión de segundo orden en
rectangulares, evaluada con intervalos.** Con aritmética afín o con bloques no encontré a nadie.

### Liao, Liu, Zhang, Wang y Qin (2017), leído completo

*Interval method for uncertain power flow analysis based on Taylor inclusion function*, IET GTD 11(5):1270-1278.

Qué hace:

- Escribe el flujo en rectangulares "since it allows us to manage quadratic equations" y expande la solución
  x(ξ) en las variables de perturbación hasta segundo orden: flujo nominal (Ec. 18), primeras derivadas con el
  Jacobiano (Ec. 23) y segundas derivadas con las Hessianas constantes (Ec. 25).
- El intervalo sale de evaluar esa expansión con aritmética de intervalos (Ec. 26): x_c + (∂x/∂ξ)[Δξ] +
  ½(∂²x/∂ξ²)[Δξ_iΔξ_j]. Es decir, la forma cuadrática se acota monomio a monomio.
- Afirma: "there is no higher order term than second order in Taylor series expansion of (4)".
- Prueba la exactitud en IEEE 9 y 57 contra Monte Carlo (1.000 y 2.000 muestras) y Krawczyk, y reporta tiempos
  en IEEE 300 y en el sistema polaco de 2.383 nudos (75 a 91 s).

Qué significa para nosotros:

| Aspecto | Liao et al. (2017) | Lo que tenemos |
|---|---|---|
| Idea base | Segundo orden en rectangulares con Jacobiano y Hessiana | **La misma.** H_i de las secciones 17 y 18 es su ∂²x/∂ξ², proyectada sobre \|V_i\|. No es un aporte nuevo |
| Cota del término cuadrático | Intervalos, monomio a monomio, sobre e y f | Equivale a lo que llamamos AF2 sobre la forma ensamblada. Los bloques (teorema k-nomial, RAA-K) son la diferencia: quitan el lado positivo espurio. El lado inferior ya es casi exacto monomio a monomio cuando la forma es semidefinida negativa (sección 18.6 de RESULTS.md) |
| Órdenes superiores | Los da por inexistentes. Es cierto para F(x), que es cuadrática en x, pero no para la solución x(ξ) | **Medidos:** el tercer orden existe y hace que la cota no contenga la solución AC, hasta 0,006 p.u. en PEGASE (sección 18.4). Su intervalo es una aproximación, no un encierro |
| Costo | Un sistema lineal por cada par (i, j) de perturbaciones (Ec. 25): del orden de m² resoluciones | Una resolución por símbolo y una traspuesta por barra de salida; con K zonas, K resoluciones |
| Magnitud de voltaje | Intervalos de e y f; no dice cómo pasa a magnitud y ángulo. Los flujos de línea, con aritmética de intervalos (Ecs. 27 y 28) | Forma cuadrática de \|V_i\| directa, con el término de curvatura; evita recombinar e y f con intervalos |
| Validación | Monte Carlo, hasta 57 nudos | Flujo AC exacto en todos los vértices (hasta 65.536) y 20.000 muestras, 2.869 nudos |
| Garantía que declara | Ninguna formal ("approximate enclosure"); Monte Carlo como referencia | Ninguna todavía: cota de segundo orden con el déficit medido |

Consecuencias:

1. **"Segundo orden exacto en rectangulares" no se puede reclamar.** Hay que citarlo como punto de partida.
2. **Liao es la referencia de comparación natural** (punto 5.1 del checklist), más cercana que Vaccaro.
3. Lo que queda como aporte posible: (a) mostrar con vértices que la expansión de segundo orden no es un
   encierro y controlar el tercer orden; (b) acotar la forma cuadrática con bloques en vez de monomio a
   monomio; (c) el planteamiento adjunto, que baja el costo de m² a m resoluciones; (d) la escala y la
   validación exhaustiva. De estos, (a) con punto fijo es el más fuerte.

Otros cercanos:

| Trabajo | Qué hace | Por qué importa | Qué falta leer |
|---|---|---|---|
| Ding, Bo, Li et al. (2015), IEEE TPWRS 30(1) | Flujo cartesiano como QCQP, relajación lineal y ajuste de cotas (OBBT) | Mismas ecuaciones cuadráticas, acotadas por relajación convexa y no por aritmética afín. Hasta IEEE-118 | Texto completo |
| Schmidt, Schegner (2020), PSCC | Aritmética afín en rectangulares, forma afín por sensibilidades, producto afín estándar | Trata los símbolos nuevos como "bounded linearization errors". Las formas cuadráticas solo las cita para circuitos. Hasta 906 nudos | — (leído completo) |
| Chen, Xiao, Wang (2019), JMPSCE 7(1) | OPF intervalar con aritmética afín y Taylor intervalar | Combina las dos ideas en distribución, hasta 113 nudos | Cómo trata los términos de orden alto |

No apareció ninguna aplicación de Messine y Touhami (forma cuadrática general) ni de Skalna y Hladík
(multiplicación de Chebyshev) a flujo de potencia, ni ninguna comparación entre el producto afín término a
término y la forma cuadrática ensamblada.

**Q2. ¿Está publicado que la forma afín de primer orden no encierra la solución AC? → Sí, de forma cualitativa.**

- Vaccaro, Cañizares y Villacci (2010): "these equations are nonlinear expressions, and hence the obtained affine
  forms are usually an underestimation of the exact result [24]. Thus, to guarantee the inclusion of the solution
  domain, each partial deviation is multiplied by an amplification coefficient [24]." El artículo no da regla
  para elegir el coeficiente.
- El mismo artículo, sobre la validación: "the random, uniformly distributed variation of parameters (with mean
  equal zero) assumed in the Monte Carlo approach tends to underestimate the worst case variations."
- Vaccaro y Cañizares (2017) describen esos métodos de sensibilidad con amplificación como "a rough estimation
  ... using a simple trial-and-error approach" y dicen que ninguno garantiza la inclusión.

Lo que no apareció: una cuantificación del fallo en los vértices de la caja de incertidumbre, ni en sistemas
grandes. Las validaciones leídas usan Monte Carlo uniforme (5.000 muestras), que no visita los vértices.

**Consecuencia para el planteamiento.** El hallazgo "la caja afín no encierra la solución AC" no se puede
presentar como nuevo: es conocido. Tampoco el uso del segundo orden para corregirlo (Liao 2017). Lo que no
encontré publicado es su medición en los vértices y a escala, ni que el segundo orden tampoco alcanza.

**Cómo resuelven los métodos existentes la falta de garantía**

| Vía | Trabajos |
|---|---|
| Amplificar las desviaciones y contraer con programación lineal | Vaccaro 2010; Xu 2016; Schmidt 2020 |
| Optimizar sobre los coeficientes afines | Vaccaro y Cañizares 2017; Vaccaro 2013; Ding 2013; Leng 2020; Zhang 2017 |
| Programas min/max (escenarios optimizados) | Zhang 2018; Cheng 2018; Refs. [10] y [11] |
| Operador de Krawczyk-Moore | Ding 2012; Araujo y da Costa 2018 |
| Relajación lineal y ajuste de cotas | Ding 2015 |
| Alcanzabilidad con sobreaproximación demostrable | Zhou y Zhang 2020, 2021 (hasta 2.000 nudos) |
| Renunciar a la garantía (intervalos de confianza) | Ref. [10] |

## Tabla de trabajos localizados

### Aritmética afín aplicada a flujo de potencia

| Referencia | Método | Producto / no linealidad | Garantía que declara | Mayor sistema | Leído |
|---|---|---|---|---|---|
| Vaccaro, Cañizares, Villacci (2010), IEEE TPWRS 25(2):624-632 | Polares. Centro del flujo nominal, desviaciones por sensibilidades amplificadas, contracción con 2 LP por símbolo | Cada operación no afín agrega un símbolo (Chebyshev) | "reliable estimation of the power flow solution hull"; comparación con Monte Carlo | IEEE 57 | Completo |
| Vaccaro, Cañizares (2017), IEEE TPWRS 32(1):274-288 | Marco unificado PF y OPF como optimización en dos etapas, más PCA para reducir símbolos | Chebyshev con símbolos "endógenos" cuyo radio se minimiza | "guaranteed to include the real solution domains"; admite que no garantiza convergencia | IEEE 118; distribución real de 1.310 nudos con 10 símbolos | Completo |
| Pirnia, Cañizares, Bhattacharya, Vaccaro (2014), IEEE TPWRS 29(6):2775-2783 | Aritmética afín para OPF con generación incierta | n.v. | n.v. | Sistema europeo real de 1.211 nudos | Resumen |
| Vaccaro, Cañizares, Bhattacharya (2013), IEEE TPWRS 28(2):1179-1186 | Aritmética de rangos: dos problemas de optimización con complementariedad | n.v. | "Reliable interval bounds" | n.v. | Resumen |
| Coletta, Vaccaro, Villacci (2019), EPSR 175:105860 | Reformulación con Jacobiano explícito y Newton-Raphson | n.v. | n.v. | n.v. | Resumen |
| Ding, Cui, Gu, Wan (2012), Autom. Electr. Power Syst. 36(13) | Intervalos con Krawczyk-Moore y aritmética afín para estrechar | n.v. | n.v. | IEEE 57 | Resumen |
| Ding, Bo, Guo et al. (2013), Proc. CSEE 33(19) | Aritmética afín no iterativa: LP, NLP y QP | n.v. | n.v. | 57 nudos | Resumen |
| Xu, Gu, Gao et al. (2016), IET GTD 10(15):3910-3918 | Aproximación lineal de la aritmética afín con LP | n.v. | n.v. | n.v. | Resumen |
| Zhang, Chen, Ngan, Yang, Hua (2017), IEEE TPWRS 32(2):1422-1429 | Flujo intervalar mixto rectangular y polar, con optimización | n.v. | n.v. | n.v. | Resumen |
| Araujo, da Costa (2018), IEEE Latin Am. Trans. 16(7) | Inyección de corrientes en rectangulares, Krawczyk en cada iteración | n.v. | n.v. | IEEE 57 | Resumen |
| Adusumilli, Kumar (2018), IET GTD 12(18):4225-4232 | Aritmética afín modificada para flujo de continuación | "reduces the overestimation caused by the AA operations"; mecanismo n.v. | n.v. | Polaco de 2.383 nudos | Resumen |
| Wang, Han, Wu (2015), IEEE TPWRS 30(6):3053-3062 | Aritmética afín compleja en barrido trifásico | n.v. | n.v. | 292 nudos | Resumen |
| Luo, Shi, Ni (2018), IEEE TPWRS 33(1):1113-1115 | Teoría de la evidencia y aritmética afín "extendida"; p-boxes | n.v. | Cotas de probabilidad | n.v. | Resumen |
| Leng, Liu, Ran, Chen, Zhang (2020), IEEE Access 8:60293-60304 | Aritmética afín y optimización no lineal; correlación por paralelogramos | n.v. | n.v. | 118 nudos (distribución) | Resumen |
| Cheng, Cheng, Shi et al. (2018), Energies 11(11):3176 | Aritmética afín con contracción por optimización, y escenarios optimizados | n.v. | n.v. | n.v. | Resumen |
| Chen, Xiao, Wang (2019), JMPSCE 7(1):139-150 | OPF intervalar con aritmética afín y Taylor intervalar | Orden alto n.v. | "strict operational bounds"; detalle n.v. | 113 nudos | Resumen |
| Schmidt, Schegner (2020), PSCC | Estimación de estado trifásica en rectangulares; contracción LP iterada | Producto estándar; símbolos nuevos como error de linealización acotado | "rigorous; enclosure of the true value is guaranteed" | 906 nudos | Completo |

### Flujo por intervalos sin aritmética afín, relajación y alcanzabilidad

| Referencia | Método | Garantía que declara | Mayor sistema | Leído |
|---|---|---|---|---|
| Liao, Liu, Zhang, Wang, Qin (2017), IET GTD 11(5):1270-1278 | Taylor de segundo orden en rectangulares, evaluado con intervalos | Encierro aproximado; Monte Carlo como referencia | Exactitud: IEEE 57. Tiempos: 2.383 nudos | Completo |
| Ding, Bo, Li et al. (2015), IEEE TPWRS 30(1):177-188 | QCQP cartesiano, relajación lineal y OBBT | Envolvente exterior | IEEE 118 | Resumen |
| Zhang, Chen, Shi et al. (2018), IEEE TSG 9(5):5217-5226 | Escenarios optimizados | n.v. | n.v. | Metadatos |
| Zhang, Liu, Zhou et al. (2023), IEEE TSTE 14(1):563-575 — Ref. [10] | Teorema central del límite sobre formas afines | Intervalos de confianza | IEEE-30 (AC), polaco 2.383 (DC) | Completo (PDF de los autores) |
| Wang, Zhang, Li et al. (2024), IEEE TPWRS 39(4):6111-6114 — Ref. [11] | Conjuntos híbridos caja-elipsoide sobre las entradas | Según el conjunto | IEEE-30 | Completo (PDF de los autores) |
| Zhou, Zhang (2020, 2021), IEEE TPWRS 35(4) y 36(3) | Flujo de potencia alcanzable (ReachFlow) | "provably over-approximated enclosure" | 2.000 nudos | Resumen |
| Dhople, Guggilam, Chen (2015), arXiv:1512.06903 | Aproximación lineal en rectangulares | Cota en norma del desajuste de potencia, no del voltaje | Sin casos | Completo (HTML) |
| Bolognani, Zampieri (2016), IEEE TPWRS 31(1):163-172 | Existencia de solución y aproximación lineal en distribución | Cota explícita del error | n.v. | Resumen |

### Herramientas de aritmética afín sin aplicación a flujo de potencia localizada

| Referencia | Qué aporta | Leído |
|---|---|---|
| Messine, Touhami (2006), Reliable Computing 12(3):171-192 | Forma cuadrática general como extensión de la aritmética afín | Metadatos |
| Skalna, Hladík (2017), Numerical Algorithms 76(4):1131-1152 | Multiplicación de Chebyshev de error mínimo para formas reducidas, en tiempo lineal | Resumen |
| Wang, Wang, Wu, Wang (2019), ACM TOMS 45(1) | Aritmética afín polar, aplicada a circuitos | Resumen |
| Bouissou, Goubault, Putot et al. (2016), TACAS | Formas afines probabilísticas con desigualdades de concentración | Resumen |

## Segunda pasada: estado del arte de cada aporte posible

Hecha junto con las pruebas de las secciones 20 a 24 de RESULTS.md. "Leído" indica cuánto: completo (HTML del
preprint), resumen o metadatos. Solo el preprint de Zhu fue contrastado dos veces.

### Aporte 1: encierro por punto fijo sobre el resto

| Trabajo | Qué hace | Qué certifica | Mayor sistema | Leído |
|---|---|---|---|---|
| Zhu (2026), arXiv:2610.26634, *Convex Inner Approximations for Safe AC Dispatch With Device Dynamics* | Rectangulares, x = x⁰ + Mξ + e; parte lineal exacta; resto cuadrático acotado con β(Q) = 1ᵀ\|MᵀQM\|1 + 2ēᵀ\|QM\|1 + ēᵀ\|Q\|ē, "the inverse is applied to the quadratic equations before absolute values are taken" | Brouwer: existencia para toda carga de una caja; unicidad local. Certificado de factibilidad; no entrega intervalos de \|V\| por barra ni compara con el rango real | IEEE 39 y 118, con cajas de 1,5 y 2 %; sintético de 3.200 solo para costo | Completo |
| Lee, Turitsyn, Molzahn, Roald (2021), IEEE TPWRS 36:4953-4966 | Restricción convexa robusta, en polares; residuo con envolventes | Brouwer: existencia dentro de un politopo de estados para toda incertidumbre de un elipsoide; los intervalos de \|V\| en barras PQ son variables de decisión | 118 | Completo |
| Lee, Nguyen, Dvijotham, Turitsyn (2019), arXiv:1803.00818 | Restricción convexa del conjunto factible; polares con admitancia ajustada en fase | Brouwer: condición suficiente de existencia | 118 | Completo |
| Dvijotham, Nguyen, Turitsyn (2018), IEEE L-CSS 2(1):25-30 | Ecuaciones cuadráticas con parámetros afines; resto acotado en norma | Brouwer: existencia a distancia ≤ r de la solución nominal | 18 | Completo |
| Nguyen, Dvijotham, Turitsyn, arXiv:1708.06845 | Politopo auto-mapeado en (θ, log V) | Brouwer: existencia para toda inyección de una caja | 1.354; la holgura reportada cae de 1 (9 barras) a 0,335 (1.354) | Completo |
| Bolognani y Zampieri (2016), IEEE TPWRS 31(1) | Punto fijo sobre el resto del modelo lineal; solo barras PQ | Banach: existencia, unicidad y cota por barra del error del modelo lineal, que ellos mismos llaman conservadora | Banco tipo IEEE 123 | Completo |

Veredicto: **el mecanismo está publicado, incluida la versión rectangular con la inversa aplicada antes de los
valores absolutos.** Todos lo usan para factibilidad y ninguno reporta más de 118 barras con cajas mayores que
2 %. Lo que no encontré: intervalos de |V| por barra tras restar el segundo orden exacto, medidos contra todos
los vértices.

### Aporte 2: ni el primer ni el segundo orden encierran

- Liao et al. (2017) tiene unas 35 citas. En lo accesible (títulos, resúmenes y unos pocos textos completos),
  nadie señala que la solución tiene términos por encima del segundo orden ni valida en vértices; la validación
  habitual es Monte Carlo.
- La única crítica explícita encontrada es de Laib et al. (2018, arXiv:1806.09689), y trata de otra cosa: "the
  obtained bounds in [5] may be conservative or even the set of obtained solutions may be empty".
- La continuación directa del método es la línea de de Sousa y Melo (2018-2023), no leída.

Veredicto: **aparentemente abierto**, con la reserva de que no se pudo leer el texto completo de la mayoría de
los trabajos que citan a Liao.

### Series de orden superior (afecta a los aportes 1 y 4)

- Embebido holomorfo multidimensional (Liu, Wang, Xu, Sun, Bak, 2017, IEEE Access): series multivariables de las
  tensiones hasta orden 14, en 2 a 4 variables.
- Sun, Ding, Qu, Li, Shahidehpour (2021), IEEE TPWRS 36(3): usan esas series para plantear el flujo por
  intervalos como optimización polinómica en caja, con relajación semidefinida. Solo resumen.
- Jabr (2019), IEEE TPWRS 34(6): Taylor de orden 2 y 3 en variable compleja para distribución, sin cota del
  resto. Leído completo en la tercera pasada (más abajo).
- Adjunto de segundo orden con dos resoluciones por dirección: Pacaud et al. (2022), arXiv:2201.00241.

Veredicto: **las series de orden 3 o más existen y ya se usaron para intervalos.** No apareció una cota
rigurosa del truncamiento.

### Aporte 3: acotar una forma cuadrática sobre una caja

- El rango exacto es NP-duro en general. Escalera establecida: Messine y Touhami (término a término), relajación
  semidefinida de Shor, y McCormick/RLT. Adjé, Garoche y Werey (2015) prueban que la semidefinida queda entre el
  rango exacto y la cota de Messine y Touhami.
- Para rango bajo hay algoritmos exactos de costo polinómico (Karystinos y Liavas, 2010); solo fragmento.
- No apareció ninguna descomposición en cuadrados de formas lineales usada en flujo de potencia.

Veredicto: **la matemática está cubierta; el uso en flujo de potencia no apareció.** La sección 23 de RESULTS.md
mide que la relajación semidefinida es más ajustada que RAA-K en estas formas.

### Aporte 5: escoger el peor escenario

- Método de escenarios optimizados (Zhang et al., 2018, IEEE TSG): un programa no lineal por variable, por punto
  interior. Solo resumen.
- Flujo de cargas de frontera (Dimitrovski y Tomsovic, 2004, IEEE TPWRS). Solo resumen.
- Sou, Khong y Sandberg (2026), IEEE TPWRS 41(1):124-141: en redes radiales, un solo flujo determinista da todos
  los voltajes mínimos; hasta 10.476 barras. Leído completo; ver la subsección siguiente.
- Gan, Li, Topcu y Low (arXiv:1311.7170), Lema 1: en redes radiales el modelo lineal acota por arriba los
  voltajes del modelo de flujo por ramas (v ≤ v̂). Es el "resto de un solo signo" observado en IEEE-33.

Veredicto: **en gran parte cubierto, sobre todo en redes radiales.** No apareció un método que use el segundo
orden para elegir el vértice de cada barra en redes malladas con barras PV.

### Encierros en redes radiales: Sou, Khong y Sandberg (2026), leído completo

*Real-time Worst-Case Load Flow for Distribution Systems with Power Injection Uncertainties*, IEEE TPWRS
41(1):124-141, DOI 10.1109/TPWRS.2025.3597048 (versión aceptada, repositorio de KTH). Es el trabajo que más
afecta a la idea de un encierro verificado para alimentadores de distribución.

| Aspecto | Qué hace |
|---|---|
| Modelo | Redes radiales balanceadas, ecuaciones DistFlow, sin elementos en derivación |
| Incertidumbre | La parte incierta de la generación distribuida de cada barra, en un semidisco (Δp ≥ 0, \|Δp + jΔq\| ≤ C). La demanda es determinista. Dicen que se adapta a cajas, sin detallarlo |
| Voltajes mínimos | Exactos con un solo flujo de carga, por monotonía del voltaje respecto de la inyección, que argumentan "for typical values of line parameters and loads" |
| Voltajes y corrientes máximos | Cota inferior con unos pocos flujos de carga y cota superior con una iteración de punto fijo; la diferencia entre ambas certifica el error durante la operación |
| Ajuste reportado (voltaje) | Error relativo certificado medio de hasta 0,37 % y máximo de hasta 1,1 % en las 17 redes de las Tablas VI y VII (el menor es 2·10⁻⁶ %); hasta 2,1 % de máximo en la red real de 533 barras (Tabla II) |
| Escala | Hasta 10.476 barras en 0,53 s |
| Redes malladas | Solo una prueba de concepto en 33 barras, muestreando la frontera del conjunto de inyecciones |

Consecuencias:

1. **En redes radiales con generación distribuida, el encierro de voltajes ya existe**, con cota superior
   garantizada para el máximo, mínimo exacto bajo monotonía, y a una escala que la construcción de este
   repositorio no alcanza (usa matrices densas del tamaño del estado al cuadrado).
2. Lo que ese trabajo no cubre: redes malladas (más allá de la prueba de concepto), elementos en derivación,
   barras con voltaje controlado, y una garantía que no dependa de suponer monotonía.
3. Medido en la sección 25 de RESULTS.md, con su método reimplementado y validado contra su Tabla I: el
   encierro de este repositorio cierra con la incertidumbre grande que ellos usan; su cota superior es más
   ajustada con incertidumbre grande en flujo directo y la de aquí lo es con incertidumbre chica y en flujo
   inverso. En redes radiales su método sigue siendo la referencia por escala y velocidad.
4. Su "error certificado" (cota superior menos cota inferior) lo domina la cota inferior: la cota superior
   sola queda bastante más cerca del máximo alcanzado (sección 25.3).

## Tercera pasada: encierros garantizados fuera del caso radial

Hecha el 15-11-2026 al empezar la opción A, con un agente de búsqueda. Niveles de lectura que reportó el agente:
**TC** = texto completo; **TC-auto** = HTML completo leído por extracción automática; **SR** = solo resumen;
**NA** = no se pudo abrir. De esta pasada revisé yo mismo solo dos textos: el Teorema 1 de Wang et al. (2018),
que está implementado en `cdhim/discs.py` y reproduce el ejemplo numérico del artículo (ρ = 0,0412), y la
Proposición 3.4 de la tesis de Wang, cuya extracción perdió los símbolos; lo demás es lectura del agente.

Pregunta: qué métodos publicados dan cotas exteriores **garantizadas** de |V| sobre un conjunto de inyecciones,
en redes de distribución malladas, con elementos en derivación o con barras de voltaje controlado.

### Familia de punto fijo (barras PQ, topología cualquiera)

| Trabajo | Qué garantiza | Ajuste reportado | Mayor red | Leído |
|---|---|---|---|---|
| Wang, Bernstein, Le Boudec, Paolone (2018), IEEE TSG 9(2):953-962 | Banach: solución única en el disco \|v_i − v̂_i\| ≤ ρ\|w_i\|, con ρ explícito. Admite elementos en derivación | Un solo caso: ρ = 0,0412. No lo comparan con el rango real | IEEE 13 modificado | TC-auto + PDF |
| Bernstein, Wang, Dall'Anese, Le Boudec, Zhao (2018), IEEE TPWRS 33(6):5832-5843 | Multifásico. Unicidad en un disco y cota del error del modelo lineal de punto fijo. Sin cota para el Taylor de primer orden ni para el modelo lineal de magnitudes | Solo errores reales (< 1 % y < 1,4 %); no comparan la cota teórica con el error real | 8.500 nudos | TC-auto + PDF |
| Wang, Le Boudec, Paolone (2019), IEEE TSG 10(2):1349-1362; contenido tomado de la tesis de C. Wang (EPFL 9456, 2019), cap. 3 | "Admisibilidad" de un conjunto de inyecciones (producto de polígonos en el plano P-Q): todo estado respeta los límites de voltaje y corriente. Usa el disco anterior y una cota por barra; si falla, particiona el conjunto. Veredicto unilateral: "admisible" o "no se sabe" | IEEE 37 e IEEE 123 modificados, ambos **radiales**: certifica hasta que el voltaje mínimo queda a 0,017 p.u. del límite | IEEE 123 trifásico, 0,03 s | Artículo NA; tesis TC |
| Wang, Stai, Le Boudec, arXiv:1810.06256 | Admisibilidad con relajaciones semidefinidas de Lasserre de orden 2 | Red **mallada de 4 nudos** e IEEE 13 monofásico | IEEE 13; minutos por caso | TC |
| Nguyen, Dvijotham, Yu, Turitsyn (2019), IEEE TSG 10(5):4827-4837 | Brouwer: existencia con \|V*\|/(1+r) ≤ \|V\| ≤ \|V*\|/(1−r). Desprecia los elementos en derivación | Cubre hasta ~80 % de la región real | 141 nudos, 3,44 ms | TC-auto |
| Cui, Sun, arXiv:1904.08855 | Brouwer y unicidad; cotas explícitas de \|V\|. Generadores como barras de módulo y ángulo fijos, no PV | IEEE 39 en carga base: error de las cotas < 0,1 p.u. | 2.383 barras (transmisión) | TC-auto + PDF |

Dos datos de la tesis de Wang que importan aquí. Reporta un caso (IEEE 37) en que **no** hay monotonía: el voltaje
mínimo se alcanza con potencia activa máxima. Y critica las cotas por relajación convexa de Molzahn y Roald
porque suponen, sin probarlo, que el flujo tiene solución para toda inyección del conjunto.

### Familia de relajación convexa

| Trabajo | Qué garantiza | Ajuste reportado | Mayor red | Leído |
|---|---|---|---|---|
| Molzahn, Roald (2019), HICSS-52 | Mínimo y máximo de \|V\| sobre una caja con relajación QC y ajuste de cotas. Válido **suponiendo** que existe solución de alta tensión para toda inyección | Cualitativo: las cotas inferiores quedan "significantly lower than the baseline voltages" ya con 0 % de variabilidad | 56 barras, **radial** | TC |
| Molzahn, Roald (2018), PSCC | Extremos de cada restricción con relajaciones semidefinida y QC; transmisión con barras PV | Sin dato de ajuste en \|V\| | 39 barras, 10 a 20 min | TC |
| Duan, Jiang, Fang, Liu (2017), IEEE TPWRS 32(1):522-530 | Jerarquía de Lasserre; transmisión con barras PV. Sin prueba de existencia | Orden 2 exacto en 6, 9 y 14 barras | 300 barras solo a orden 1; el orden 2 no pasa de ~40 | TC |

### Aproximados o sin garantía (para no confundirlos con lo anterior)

- **Jabr (2019), IEEE TPWRS 34(6)**: Taylor de orden 2 y 3 desde el punto sin carga, con aritmética circular, en
  redes radiales **y malladas** de 161, 1.464 y 3.146 nudos. Es lo más cercano en espíritu a un modelo de
  segundo orden para distribución, pero **no acota el resto**: "approximate circular voltage interval
  solutions", validadas con 10.000 muestras. TC.
- Zhou y Zhang (2020, 2021), flujo "alcanzable" con zonotopos: afirma un encierro probado, pero el resto de
  Lagrange "can be approximately calculated" y la validación es por Monte Carlo. TC parcial.
- Rodrigues Junior, Melo, Nepomuceno (2022), IJEPES 139: operador de Krawczyk trifásico en redes radiales; sin
  prueba de inclusión ni redondeo dirigido a la vista. TC.
- Buason, Misra, Molzahn (arXiv:2501.13825): lineal a trozos guiado por sensibilidades de segundo orden;
  conservador solo sobre las muestras. TC-auto.
- Raj y Kumar (2020): aritmética afín en distribución débilmente mallada de 33 y 118 barras; el resumen no
  afirma garantía. SR.

### Veredicto de la tercera pasada

- **Malladas con barras PQ:** la garantía teórica existe (discos de Wang et al., admisibilidad de Wang, Le
  Boudec y Paolone, Brouwer en Nguyen et al.). Ninguna fuente abierta reporta, para un alimentador **mallado**,
  el ancho del intervalo certificado frente al rango real sobre un conjunto de inyecciones; la evidencia
  numérica en distribución mallada es una red de 4 nudos.
- **Con elementos en derivación:** cubierto en teoría por la familia de Wang (van dentro de la matriz de
  admitancias); ninguna fuente mide su efecto en el ajuste.
- **Con barras PV:** no cubierto por la familia de punto fijo de distribución (es trabajo futuro en Wang et al.
  y en Nguyen et al.). Solo por relajación convexa en transmisión de hasta 300 barras, suponiendo solvabilidad.
- **Segundo orden con garantía:** no apareció ningún encierro garantizado basado en sensibilidades de segundo
  orden o modelos de Taylor. Jabr (2019) usa órdenes 2 y 3 sin cota del resto.
- La sección 25 de RESULTS.md mide el disco de Wang et al. (2018) en las mismas redes que el encierro de aquí.

Lo que el agente **no** pudo verificar: el artículo de Wang, Le Boudec y Paolone en IEEE TSG (solo su resumen y
la tesis); los textos completos de Ding (2015), Sun (2021), Zhang (2018), Raj y Kumar (2020), Cheng (2023) y
Reihani (2025); y no pudo barrer de forma sistemática 2022-2026. Sou et al. no tiene citas registradas todavía.
Que no haya aparecido un encierro reciente para distribución mallada no prueba que no exista.

## Cuarta pasada: lo que faltaba leer antes de sostener el aporte

Hecha el 15-11-2026 con un agente de búsqueda. Mismos niveles de lectura (TC, SR, NA). Yo revisé las notas del
agente (qué leyó, con qué acceso y qué no pudo verificar), no cada fuente.

### Wang, Le Boudec y Paolone (2019), IEEE TSG 10(2): no se pudo leer completo. **Bloqueado: falta el PDF**

- No existe preprint en arXiv (la consulta a su API no devuelve nada), el registro de EPFL Infoscience no tiene
  archivo y IEEE Xplore es de pago. **NA.**
- Lo que sí se leyó: su resumen; las láminas de Le Boudec en el Simons Institute (marzo de 2018), que traen
  recortes del artículo (su Ec. (6), la tabla de notación y tres figuras) y enuncian sus teoremas 1 y 3;
  arXiv:1810.06256, que cita su Lema 2 y describe el método como heurístico y de tiempo exponencial; y la tesis
  de C. Wang, §3.3 y §3.5. Que la tesis reproduce el artículo es una inferencia (misma fórmula, mismas figuras,
  mismos ejemplos); la tesis no lo dice.
- La cota por nodo del artículo (su Ec. (6), Ec. (3.9) de la tesis) es la que implementa `cdhim/discs.py` en
  lectura monofásica. No apareció una cota de un solo paso más ajustada, ni una que use del conjunto algo más
  que la desviación máxima por entrada.
- Lo que el artículo tiene y `discs.py` no: **subdivisión del conjunto con empalme de dominios**, que es su única
  vía al ajuste; recentrado de la referencia; cota de corrientes; multifase.
- Dos precisiones de rigor sobre el baseline: la unicidad vale en el disco mayor (ρ‡) y la localización en el
  menor (ρ†); y la inyección de referencia debe ser F(v̂), no la especificada, así que el residuo del flujo de
  carga queda fuera de la hipótesis tal como está programado.
- No imprime radios ni cotas con que validar una implementación. Sus ejemplos son los alimentadores IEEE 13, 37
  y 123, todos radiales.
- Medición propia a raíz de esto: la subdivisión con discos frente a la subdivisión con el encierro de aquí,
  sección 26 de RESULTS.md.

### Multifase: Bernstein, Wang, Dall'Anese, Le Boudec y Zhao (2018). TC

- Teorema 2: existencia y unicidad con cargas en estrella y en delta (ξ = ξ^Y + ξ^Δ). Da un radio uniforme; no
  da cota por nodo ni un enunciado sobre conjuntos de inyecciones.
- La ξ trifásica usa valores absolutos por entrada de la matriz 3N × 3N, así que un certificado monofásico (de
  secuencia positiva) no se transfiere: en su ejemplo balanceado da 0,1085 frente a 0,1855.

### Rigor en punto flotante

No apareció ningún trabajo que pruebe, con el redondeo controlado, un encierro del rango de voltaje bajo
inyecciones inciertas en una red de tamaño realista.

- Barboza, Dimuro y Reiser (2004): Newton intervalar con INTLAB en una red de 3 barras. TC.
- Rodrigues Junior, Melo y Nepomuceno (2022): Krawczyk trifásico por inyección de corrientes con INTLAB, IEEE
  13, 33 y 69; validado por Monte Carlo, sin afirmar una inclusión verificada. TC.
- Jabr (2019): aritmética circular con INTLAB sobre una serie truncada sin cota del resto. Es la evaluación
  rigurosa de un modelo que no lo es. TC.
- Pereira, da Costa y Rosa (2012) y Wang y Alvarado (1992): solo fragmentos; no se sabe si controlan el redondeo.
- Modelos de Taylor: nada para el flujo de potencia estacionario.

### Vecinos nuevos, 2023 a 2026

- **Nguyen, Le y Le (2026), preprint en Research Square, sin revisión por pares.** Encierro por punto fijo con un
  radio por barra y subdivisión de la caja, en alimentadores **radiales** de 16 a 141 barras. Declara que usa
  márgenes en doble precisión "en vez de aritmética intervalar con redondeo dirigido". Es el vecino más cercano
  en redes radiales; no cubre mallas. TC.
- Zhu (2026), arXiv:2610.26634: restricción convexa por punto fijo en transmisión mallada con barras PV (IEEE 39,
  118 y una sintética de 3.200 barras). Certifica que se cumplen los límites; no da un encierro ajustado por
  barra. Resumen y lectura rápida.
- Lin, Chen, Li y Zhao (2023), IEEE TPWRS: restricciones lineales a partir de Brouwer y del desarrollo de Taylor
  de segundo orden, IEEE 33 y 136. Son condiciones interiores (suficientes), no un encierro. SR.
- Christianen et al. (2025) y Mavalizadeh y Almassalkhi (2024): aproximaciones interiores. SR y lectura rápida.
- Los 16 trabajos que citan a Jabr (2019) desde 2023 son aproximaciones o revisiones, sin cota del resto según
  sus resúmenes. Sou et al. sigue sin citas indexadas.

**Resultado del barrido:** no apareció ningún trabajo que reporte encierros garantizados por barra en
distribución mallada a pocos puntos porcentuales del rango real. Es un resultado negativo de búsqueda, no una
prueba de ausencia: no hubo acceso a la búsqueda en texto completo de IEEE Xplore ni de Google Scholar.

## Dónde se consultó cada fuente principal

- Vaccaro, Cañizares, Villacci (2010), texto completo: https://dspacemainprd01.lib.uwaterloo.ca/server/api/core/bitstreams/c5d87021-6777-444e-a64d-f5978b782d56/content
- Vaccaro, Cañizares (2017), texto completo: https://dspacemainprd01.lib.uwaterloo.ca/server/api/core/bitstreams/bb1de9a2-08d6-4e28-8faa-5b4aadeae511/content
- Schmidt, Schegner (2020), texto completo: https://pscc-central.epfl.ch/repo/papers/2020/442.pdf
- Liao et al. (2017), texto completo (PDF de la autora); ficha: https://doi.org/10.1049/iet-gtd.2016.1344
- Ding et al. (2015), resumen: https://scholarsmine.mst.edu/ele_comeng_facwork/2670
- Coletta, Vaccaro, Villacci (2019), resumen: https://ic.unicamp.br/~stolfi/EXPORT/bibliography/split-bib/2019/col-vac-vil-19-aa-powflow.bib
- Adusumilli, Kumar (2018), resumen: https://ic.unicamp.br/~stolfi/EXPORT/bibliography/split-bib-5/2018/adu-kum-18-aa-voltst.bib
- Zhou, Zhang (2021), resumen: https://researchconnect.stonybrook.edu/en/publications/reachable-power-flow-theory-to-practice/
- Dhople, Guggilam, Chen (2015): https://ar5iv.labs.arxiv.org/html/1512.06903
- Sou, Khong, Sandberg (2026), texto completo: https://kth.diva-portal.org/smash/get/diva2:1998957/FULLTEXT01
- Zhu (2026), texto completo: https://arxiv.org/html/2610.26634
- Lee et al. (2021): https://ar5iv.labs.arxiv.org/html/2005.04835 · Nguyen et al.: https://ar5iv.labs.arxiv.org/html/1708.06845
- Dvijotham et al. (2018): https://ar5iv.labs.arxiv.org/html/1703.08881 · Bolognani y Zampieri: https://ar5iv.labs.arxiv.org/html/1403.5031
- Adjé et al. (2015): https://ar5iv.labs.arxiv.org/html/1411.5847 · Gan et al.: https://arxiv.org/html/1311.7170v1
- Laib et al. (2018): https://ar5iv.labs.arxiv.org/html/1806.09689 · Pacaud et al. (2022): https://arxiv.org/html/2201.00241
- Tercera pasada. Wang et al. (2018): https://arxiv.org/abs/1602.08372 · Bernstein et al. (2018): https://arxiv.org/abs/1702.03310 ·
  tesis de C. Wang: https://infoscience.epfl.ch/entities/publication/2987d778-4466-47e5-90d7-3dc96197603f ·
  Wang, Stai, Le Boudec: https://arxiv.org/abs/1810.06256 · Nguyen et al. (2019): https://ar5iv.labs.arxiv.org/html/1705.05774 ·
  Cui y Sun: https://arxiv.org/abs/1904.08855 · Molzahn y Roald: https://molzahn.github.io/pubs/molzahn_roald-hicss2019.pdf y
  https://molzahn.github.io/pubs/molzahn_roald-acopf_robust2018.pdf · Duan et al. (2017): https://livrepository.liverpool.ac.uk/3005659/1/main.pdf ·
  Jabr (2019): https://scholarworks.aub.edu.lb/items/1f3ec672-3967-4556-93aa-2ad88b640ce7 · Buason et al.: https://ar5iv.labs.arxiv.org/html/2501.13825
- El resto de las fichas salió de Crossref, OpenAlex, de los repositorios de las universidades de los autores y de la
  bibliografía de aritmética afín de J. Stolfi (ic.unicamp.br).

## Lo que hay que leer completo antes de cerrar la Puerta 1

1. Torelli y Vaccaro (2015), "A second order dynamic power flow model", EPSR 126: Liao lo cita para el uso de
   rectangulares; es de la revista objetivo y hay que ver qué hace con el segundo orden.
2. Coletta, Vaccaro y Villacci (2019): es el antecedente en la misma revista y en rectangulares.
3. Adusumilli y Kumar (2018): cómo modifica el producto y qué garantía declara a 2.383 nudos.
4. Chen, Xiao y Wang (2019): tratamiento del orden alto.
5. Messine y Touhami (2006) y Ninin, Messine y Hansen (2015), para ubicar el k-nomial.
6. Luo, Shi y Ni (2018): qué es la aritmética afín "extendida".
7. La Ref. [12], que sigue sin PDF.

---

## Quinta pasada: la capa probabilística (versión de formas afines, experimentos B4 y B5)

Pregunta: ¿está publicado acotar |V| por barra y cantidades de rama con la forma afín de sensibilidades en el
punto AC, un radio de concentración que no supone normalidad (Bernstein) y una corrección de segundo orden del
error de linearización, validado contra AC a escala de transmisión? Búsqueda del 26-11-2026, no exhaustiva:
nueve búsquedas web y dos artículos leídos completos. Niveles: **TC** texto completo, **R** solo resumen o
fragmentos del buscador.

| Trabajo | Qué hace | Distribución | Linearización | Mayor sistema / validación | Nivel |
|---|---|---|---|---|---|
| Roald y Andersson (2019), IEEE TPWRS 33(3), arXiv 1706.03241 | CC-AC-OPF: el punto previsto con AC completo; el efecto de la incertidumbre linealizado en ese punto. Restricciones de azar sobre **voltajes, corrientes de rama** y generación. Márgenes analíticos o por Monte Carlo, iterando | Gaussiana (Φ⁻¹); dicen que se puede cambiar por reformulaciones distribucionalmente robustas ([29] del artículo) | Primer orden. Escriben que "la imprecisión por la distribución es del mismo orden que la del error de linealización". No corrigen el segundo orden | Polaco 2.383 barras con 941 cargas inciertas; pruebas dentro y fuera de muestra | TC |
| Lubin, Dvorkin y Roald (2019), IEEE TPWRS, arXiv 1803.08754 | CC-AC-OPF como SOCP; restricciones cuadráticas de azar para la potencia aparente | Gaussiana, extensible con DRO | Primer orden en el punto AC | IEEE 118 | TC |
| Halilbasic, Pinson y Chatzivasileiadis (2019), IEEE TPWRS, arXiv 1804.05754 | Relajación SOC del CC-AC-OPF y ajustes por Taylor de primer orden; recuperación de factibilidad AC | — | Primer orden; reconocen el error de linealización como fuente de imprecisión | — | R |
| Summers et al. (2016), NAPS | Regulación de voltaje con restricciones de azar tipo Chebyshev desde media y covarianza estimadas en línea | Sin distribución (solo momentos) | Flujo AC linealizado | Distribución | R (el PDF redirige a otro sitio) |
| Li, Jiang y Mathieu (2019), IEEE TCNS, arXiv 1811.10217 | OPF distribucionalmente robusto con momentos y unimodalidad | DRO, menos conservador que Chebyshev | — | — | R |
| Talkington et al. (2025), arXiv 2510.17798 | Desigualdades de concentración matriciales (Bernstein matricial) para el espectro de la admitancia con **parámetros de red** inciertos; errores de DC y LinDistFlow | Concentración | No es incertidumbre de inyección | Redes IEEE | R |
| Bouissou, Goubault, Putot et al. (2016), TACAS | Formas afines probabilísticas con desigualdades de concentración (ya en la primera pasada) | Concentración | — | Sin aplicación a sistemas de potencia | R |
| Pirnia, Cañizares, Bhattacharya, Vaccaro (2012-2014) | Flujo de potencia y OPF con aritmética afín frente a Monte Carlo: cotas algo más anchas que el rango de MC, sin suponer distribución | Ninguna (peor caso) | Aritmética afín | 1.211 barras | R |
| Nemirovski y Shapiro (2006), SIAM J. Optim. | "Aproximación de Bernstein" de restricciones de azar, convexa | Concentración (cotas de la función generadora de momentos) | — | No es de sistemas de potencia; no encontré quién la aplique a OPF | R |
| Flujo probabilístico por cumulantes (Usaola 2009 y otros) | Momentos de salida desde la linealización y expansión de Gram-Charlier o Cornish-Fisher | Reconstruye la distribución | Primer orden; la linealización se nombra como su principal limitación | — | R |
| PPF acelerado, arXiv 2010.14995 | Las ecuaciones son cuadráticas en coordenadas cartesianas: el desarrollo de segundo orden es exacto ahí | — | Hessiana constante en un modelo reducido | Distribución | R |

### Qué ya está publicado

- **Linealizar en el punto AC y acotar voltajes y corrientes de rama con un cuantil analítico** es la reformulación
  analítica de Roald y Andersson. El "3σ" de este repositorio es esa misma reformulación gaussiana. Ellos ya
  señalan que el error por la distribución y el de linealización son del mismo orden.
- **Reformulaciones sin suponer distribución** con momentos (Chebyshev, Cantelli, DRO, unimodalidad) son
  corrientes en CC-OPF.
- **Concentración sobre formas afines** está publicada en análisis de programas (Bouissou et al., 2016).
- **El peor caso con aritmética afín** frente a Monte Carlo está publicado (Pirnia, Vaccaro).

### Qué no apareció

1. **Corrección de segundo orden del error de linealización en cotas probabilísticas de voltaje a escala de
   transmisión**, usando el rango del término cuadrático (asimétrico y de cola pesada, que un ajuste de media y
   varianza no captura: B4) en vez de iterar con Monte Carlo AC.
2. **Bernstein (soporte acotado más varianza) sobre la forma afín de |V|** como radio sin suponer distribución.
   Es una desigualdad de libro y la idea está en Bouissou; lo nuevo sería la aplicación y su medición frente a
   Chebyshev y Gauss.
3. **Caídas de voltaje por rama desde la forma afín de la diferencia**, y cuánto ganan frente a restar
   intervalos (B5: 4 a 21 veces). Roald y Andersson sí acotan corrientes de rama con la forma lineal, así que la
   idea general de cantidades de rama no es nueva: lo medible es la ganancia y la validación.
4. **Medición de que la caja de peor caso de segundo orden no encierra la solución AC en los vértices** de una red
   de 2.869 barras (B5). Ya estaba anotado en las secciones 17 y 18 de RESULTS.md para el primer orden.

### Lo que falta para cerrar esta pasada

- Leer completo Roald, Oldewurtel, Van Parys y Andersson (2015), con DRO, y la referencia [29] de Roald y
  Andersson. Son la comparación natural para el radio sin distribución: hay que medir Bernstein frente a Cantelli
  (que solo usa la varianza) en los mismos casos.
- Buscar en IEEE Xplore o Scopus "second-order" junto con "chance constraint" y "voltage", y flujo probabilístico
  con Hessiana a escala de transmisión. El buscador web no basta para descartar.
- Leer Summers et al. (2016) completo.

---

## Sexta pasada: estado del arte de la contribución principal (26-11-2026)

**Método.** Búsqueda en OpenAlex (indexa IEEE, Elsevier, IET, MDPI y otros, con resúmenes) por título y resumen,
veinte consultas booleanas ordenadas por relevancia y por citas: segundo orden o Hessiana en flujo
probabilístico; error de linealización en restricciones de azar y en flujo probabilístico; desigualdades de
concentración en sistemas de potencia; radios sin distribución para voltaje; cumulantes, estimación por puntos y
punto de silla con atención a las colas; caos polinomial y procesos gaussianos; aritmética afín probabilística;
caídas de voltaje y corrientes de rama con incertidumbre. 154 trabajos en la primera ronda y unos 170 en la
segunda. Se leyeron los resúmenes de los cercanos; ninguno de los nuevos se leyó completo. **No reemplaza una
búsqueda en IEEE Xplore o Scopus**: OpenAlex no tiene resumen de todos los artículos de Elsevier (por ejemplo, el
de 1985) y su búsqueda por texto es menos fina.

### Los trabajos más cercanos

| Trabajo | Qué hace | Diferencia con la propuesta | Nivel |
|---|---|---|---|
| Brucoli, Torelli, Napoli (1985), "Quadratic probabilistic load flow with linearly modelled dispatch", IJEPES, doi 10.1016/0142-0615(85)90042-0 | Flujo probabilístico con el modelo cuadrático | Antecedente clásico del segundo orden; sin resumen en OpenAlex ni acceso abierto | Ficha |
| "The algorithm of probabilistic load flow retaining nonlinearity" (2002), PowerCon, doi 10.1109/icpst.2002.1047154 | Serie de Taylor con términos no lineales; esperanza y varianza; IEEE 14 | Momentos, no colas | R |
| Sandoval, Gong, Chung (2022), "Three-phase second-order analytic probabilistic load flow with voltage-dependent load", IEEE TPWRS, doi 10.1109/tpwrs.2022.3159392 | Modelo cuadrático exacto por inyección de corriente; **los dos primeros momentos** de voltaje y corriente; distribución | **Es el antecedente directo del segundo orden.** Usa momentos; B4 muestra que corregir con media y varianza no basta al 0,27 % (falla en 15 de 25 celdas de PEGASE). Distribución, una red práctica | R |
| Liu, Sun, Wang, Ju (2018), "Probabilistic power flow analysis using multidimensional holomorphic embedding and generalized cumulants", IEEE TPWRS, doi 10.1109/tpwrs.2018.2846203 | Solución AC analítica no lineal; cumulantes y Gram-Charlier; "mejores colas" | Otra vía a la no linealidad; IEEE 30; la cola sale de Gram-Charlier, que B8 sugiere poco fiable lejos en la cola | R |
| Basirat, Mohammadi, Keihan Asl (2026), "Second order exponential asymptotic expansion for probabilistic load flow analysis", IJEPES, doi 10.1016/j.ijepes.2026.111582 | Expansión de segundo orden, costo independiente del número de variables, CDF sin integrar | Distribuciones, no cotas de cobertura; IEEE 14 y 118 | R |
| Roald y Andersson (2019), IEEE TPWRS (quinta pasada) | Linealización en el punto AC, cuantil gaussiano, márgenes por Monte Carlo iterando | Primer orden; el error de linealización se corrige con Monte Carlo AC, no con el modelo cuadrático | TC |
| "Convex relaxations and approximations of chance-constrained AC-OPF problems" (2019), IEEE TPWRS | Ajustes por Taylor de primer orden | Primer orden | R |
| "Asymptotically tight conic approximations for chance-constrained AC OPF" (2022), EJOR, doi 10.1016/j.ejor.2022.06.020 | Aproximación cónica de restricciones de azar bilaterales | Optimización; sin resumen en OpenAlex | Ficha |
| "Sequential convex approximation approach for CC-ACOPF under arbitrary random distribution" (2026), IEEE TPWRS, doi 10.1109/tpwrs.2026.3670499 | Restricciones de azar con distribuciones arbitrarias, polinomios de Bernstein (no la desigualdad) | Optimización y muestreo; "Bernstein" es otra cosa | R |
| "Nonparametric probabilistic load flow with saddle point approximation" (2017), IEEE TSG | Punto de silla para las colas | Lineal; colas aproximadas, sin garantía | Ficha |
| "Risk assessment of rare events in probabilistic power flow via hybrid multi-surrogate method" (2019), IEEE TSG, doi 10.1109/tsg.2019.2940928 | Caos polinomial más reevaluación AC de las colas | Muestreo con el modelo completo en la cola; costo mayor | R |
| Singh, Moger, Jena (2025), "A modified point estimate-based probabilistic load flow approach for improving tail accuracy", EPSR, doi 10.1016/j.epsr.2025.111606 | Estimación por puntos modificada para las colas | **Hay que leerlo**: es el intento más reciente de arreglar la cola de la estimación por puntos | Ficha |
| "Safety screening for voltage control in active distribution grids via distributionally robust conformal screening" (2026), arXiv | Cobertura por predicción conforme | Basado en muestras; distribución | Ficha |
| "Probabilistic interval power flow ... correlation of distributed wind power" (2024), IJEPES, doi 10.1016/j.ijepes.2024.109827 | Aritmética afín más teoría de la evidencia; correlación con paralelogramo, polígono o elipse | Correlación como restricción geométrica; distribución | R |

### Cómo queda el aporte

- **No se puede reclamar "flujo probabilístico de segundo orden"**: existe desde 1985 y en TPWRS 2022 en forma
  analítica. Hay que citarlo y diferenciarse.
- **Lo que no apareció en ninguna de las dos rondas**:
  1. usar **el rango** del término de segundo orden, no sus momentos, para que una cota linealizada cubra al nivel
     del 0,27 %, con la evidencia de que los momentos no bastan (B4) y de por qué (asimetría y cola);
  2. combinarlo con un radio de concentración que solo supone independencia y soporte acotado, en vez de un cuantil
     gaussiano o una expansión de Cornish-Fisher o Gram-Charlier;
  3. validarlo contra AC en las colas a escala de transmisión (4 sistemas, hasta 2.869 barras, hasta el 95 % del
     límite de cargabilidad), junto con los métodos publicados (cumulantes, estimación por puntos) que no cubren;
  4. corrientes de rama y caídas de voltaje con la misma construcción.
- **Riesgo:** el artículo de TPWRS 2022 también calcula corrientes. Hay que leerlo completo para saber si reporta
  cobertura de cuantiles extremos o solo momentos y densidades.

### Lo que hay que leer completo antes de redactar

1. TPWRS 2022, segundo orden analítico: qué cuantiles reporta y en qué red.
2. EPSR 2025, estimación por puntos para las colas.
3. Brucoli, Torelli y Napoli (1985): ver qué momentos usa (está tras muro de pago).
4. IJEPES 2026, expansión asintótica de segundo orden.
5. TPWRS 2018, holomorphic embedding con cumulantes: si mide el error en la cola.

Acceso: los de 2022 (TPWRS), 2018 (TPWRS) y 1985 (IJEPES) no tienen versión abierta. Los de IJEPES 2026 y EPSR 2025
figuran como abiertos en OpenAlex, pero ScienceDirect y SSRN rechazan la descarga automática (HTTP 403); hay que
bajarlos a mano.

**Lo que sigue haciendo falta de los autores:** una búsqueda en IEEE Xplore o Scopus con acceso a texto completo.
Consultas sugeridas: ("second-order" OR quadratic) AND ("probabilistic power flow" OR "probabilistic load flow")
AND (quantile OR tail OR "confidence interval"); y "linearization error" AND "chance constraint" AND voltage AND
(quadratic OR Hessian).
