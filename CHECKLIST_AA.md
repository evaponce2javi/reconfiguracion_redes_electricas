# Checklist: versión de aritmética afín

Objetivo: pasar de CD-HIM (elipse o caja por rama) a un **flujo de potencia afín de segundo orden**, donde el
resto de la forma afín se acota con los operadores del grupo (RAA, presupuesto diagonal, teorema k-nomial).
Revista objetivo: Electric Power Systems Research.

El proceso tiene tres puertas. Si una no se pasa, se detiene y se revisa el plan antes de seguir (ver "Plan B").

Cada punto dice quién lo hace: **(grupo)** es una decisión o un insumo de los autores, **(código)** es trabajo
en este repositorio, **(texto)** es redacción.

---

## Reglas de trabajo

Siguen las del primer encargo:

- No inventar ni redondear resultados. Lo que no se pueda calcular se reporta como bloqueado, con la causa.
- Semilla fija en todo experimento aleatorio, guardada junto al resultado.
- Cada experimento queda como CSV en `results/`, con un script que lo regenera de principio a fin.
- Si un resultado contradice lo que se esperaba, no se ajusta el código para que calce: se reportan ambas cosas.
- Los números van a `RESULTS.md`; los autores los pasan al manuscrito.

Y se agregan:

- No escribir "garantizado", "verificado" ni "riguroso" mientras no estén hechos el punto fijo y el redondeo
  dirigido (3.4 y 3.6). Hasta entonces: "cota de segundo orden, validada contra AC".
- Toda cota se prueba contra el flujo AC exacto en los vértices extremos, no solo con Monte Carlo. El muestreo
  no alcanza los vértices (sección 17.2 de RESULTS.md).
- Toda mejora se mide contra una referencia de la literatura, no contra una versión débil propia.
- Los trabajos no publicados del grupo no se envían a servicios externos.

---

## Punto de partida (ya hecho, sección 17 de RESULTS.md)

- [x] Modelo de segundo orden en coordenadas rectangulares, sin ensamblar H_i, verificado contra las
      sensibilidades polares y contra el flujo AC (`experiments/x1_aa_second_order.py`).
- [x] Flujo AC exacto en los dos vértices extremos de cada barra, en IEEE-33 y PEGASE-2869.
- [x] Comparación de siete cotas del resto en 32 barras de IEEE-33 y 24 de PEGASE.
- [x] Dependencia por rama: zonotopo de cada pareja y radio de la caída de voltaje.
- [x] Radio de Hoeffding y corrimiento del centro de segundo orden.

Lo que eso deja abierto: la cota buena cuesta N_u² por barra, solo se acota el segundo orden, la cota de
presupuesto diagonal es una lectura simplificada y no hay ninguna comparación con un flujo afín publicado.

---

## Parte 0 — Decisiones del grupo, antes de programar

Siguen abiertas: son del grupo. Para poder correr las Partes 1 y 2 se usó un valor provisional en cada una,
anotado abajo. Ninguno compromete al grupo; si se cambia, se vuelve a correr `experiments/x2_few_sources.py`.

- [ ] **0.1 (grupo) Tipo de garantía.** (a) Encierro verificado: punto fijo más redondeo dirigido. (b) Cota de
      segundo orden validada contra AC en vértices y por muestreo. Recomendación: (b) para este envío y (a) como
      extensión. *Hecho cuando:* está escrito en una frase qué afirma el paper y qué no.
      *Provisional:* (b). La Parte 2 muestra que (b) sola no contiene la solución AC; ver 3.4.
- [ ] **0.2 (grupo) Modelo de incertidumbre.** (a) Miles de símbolos independientes, uno por carga, como ahora.
      (b) Decenas de fuentes correlacionadas: zonas o plantas, cada símbolo mueve un conjunto de barras.
      Recomendación: (b) como caso principal y (a) como prueba de esfuerzo. *Hecho cuando:* están definidos K,
      qué barras mueve cada símbolo y con qué semiancho.
      *Provisional:* (b), con zonas sacadas solo de la topología (`cdhim.spectral_zones`: bisección espectral,
      mismas barras por zona), un símbolo por zona que mueve P y Q de todas sus cargas, y el δ de diseño de
      cada sistema. K = 4, 8, 16 en IEEE-33 y 8, 16, 32, 64 en PEGASE. Las zonas reales las define el grupo.
- [ ] **0.3 (grupo) Barra slack.** Única, como ahora, o distribuida con factores de participación. Con slack
      única el peor caso de PEGASE mueve los ángulos cerca de 1 rad. *Hecho cuando:* está decidido y, si es
      distribuida, de dónde salen los factores.
      *Provisional:* se corrieron las dos. "Participación": la generación sigue al cambio de demanda en partes
      proporcionales al despacho (`cdhim.participation_factors`); la slack solo cubre su parte y las pérdidas.
- [ ] **0.4 (grupo) Salidas del método.** Solo |V| por barra, o también caídas de voltaje y flujos por rama. La
      dependencia entre barras rinde en las cantidades de rama (sección 17.5).
      *Provisional:* solo |V| por barra.
- [ ] **0.5 (grupo) Qué se conserva de CD-HIM.** Propuesta: se retiran la elipse, G y γ_th; se conserva el motor
      de sensibilidades dispersas.
- [ ] **0.6 (grupo) Relación con los trabajos no publicados.** El paper nuevo se apoya en el k-nomial y en RAA.
      Decidir el orden de envío y si habrá un preprint para citar.
- [x] **0.7 (grupo) Operadores.** Entregar la implementación real del presupuesto diagonal y del k-nomial, o
      decidir si se portan a Python o se enlaza el código de IBEX.
      *Resuelto:* la implementación del grupo está en `ibex-lib/raak/affine/ibex_KNomial.h`. Se portó a
      `cdhim/knomial.py` y reproduce sin diferencias los 80 casos de `knomial_golden.txt`. No redondea hacia
      fuera; para rigor hay que enlazar el C++ o emular el redondeo.

---

## Parte 1 — Revisión bibliográfica (Puerta 1: novedad)

**Estado: primera pasada hecha, Puerta 1 sin cerrar y con un solape confirmado.** Está en
[RELATED_WORK.md](RELATED_WORK.md): 35 fuentes abiertas, cuatro leídas completas.

- **Liao et al. (2017) ya publican el segundo orden en rectangulares** (Jacobiano y Hessiana), evaluado con
  intervalos monomio a monomio. Esa idea no se puede reclamar; pasa a ser la referencia de comparación.
- Liao da por inexistentes los órdenes superiores. La sección 18 mide que existen y que por ellos la cota no
  contiene la solución AC.
- Que la forma afín de primer orden subestima también está publicado (Vaccaro 2010), sin cuantificar.
- No apareció nadie que acote la forma cuadrática con bloques ni con aritmética afín.

Lo que queda como aporte posible: control del tercer orden, cota por bloques, costo lineal en vez de cuadrático
en el número de fuentes, y validación exhaustiva a escala. Falta leer completos los artículos de la lista final
de RELATED_WORK.md.

- [ ] **1.1 (texto) Flujo de potencia con aritmética afín.** Vaccaro, Cañizares y Villacci (2010); Vaccaro y
      Cañizares (2017); Coletta, Vaccaro y Villacci (2019, rectangulares); Ding et al. (2012, con Krawczyk);
      aritmética afín modificada y aritmética afín polar; Luo, Shi y Ni (2018). Para cada uno: cómo trata el
      producto de formas afines, qué garantía declara y cuál es el sistema más grande que prueba.
- [ ] **1.2 (texto) Formas cuadráticas en aritmética afín.** Messine y Touhami (forma cuadrática general);
      multiplicación de Chebyshev de formas reducidas; AF2 y rAF2. Ubicar el k-nomial frente a ellos.
- [ ] **1.3 (texto) Flujo de potencia por intervalos y verificado.** Krawczyk, relajaciones con ajuste de cotas,
      modelos de Taylor.
- [ ] **1.4 (texto) Capa probabilística.** Refs. [10] y [11], formas afines probabilísticas con desigualdades de
      concentración. Conseguir la Ref. [12].
- [ ] **1.5 (texto) Dos preguntas que hay que contestar con cita.** ¿Alguien acotó el resto cuadrático del flujo
      rectangular con aritmética afín a esta escala? ¿Está reportado que la caja afín de primer orden no encierra
      la solución AC?

**Puerta 1.** Se pasa si ninguna de las dos respuestas es "sí, igual que aquí". Si lo es, se reposiciona el
aporte antes de programar.

*Entregable:* tabla de trabajos relacionados (método, tipo de garantía, tratamiento del producto, mayor sistema).

---

## Parte 2 — Prueba decisiva con pocas fuentes (Puerta 2: ¿cierra el método?)

**Estado: corrida. La Puerta 2 no se pasa: falla el criterio 1.** Detalle en la sección 18 de RESULTS.md.

- [x] **2.1 (código) Fuentes agregadas.** Extender `build_sources` para que un símbolo mueva varias barras con
      participaciones fijas. Probar K = 8, 16, 32 y 64. *Hecho cuando:* con un símbolo por carga se recuperan
      exactamente los resultados actuales (prueba unitaria).
- [x] **2.2 (código) H_i para todas las barras.** Con K chico, H_i es K × K y se ensambla para todas las barras.
- [x] **2.3 (código) Rango de referencia.** Para K ≤ 16, flujo AC exacto en los 2^K vértices, más búsqueda
      local desde el peor. Reemplaza al "rango alcanzado" de la sección 17, que era solo una cota inferior.
- [x] **2.4 (código) Cotas con los operadores reales del grupo** (0.7): AF2, presupuesto diagonal, k-nomial.
      *Hecho cuando:* cada operador pasa su propia prueba de encierro contra fuerza bruta en formas chicas.
      Hecho con RAA-K portado desde ibex-lib; la prueba de encierro está en `tests/test_knomial.py`.
- [x] **2.5 (código) Medir, por sistema y por K:** exceso de la caja afín en los vértices; ancho de cada cota
      sobre el rango de referencia; ancho frente a la caja 3σ y al rango de Monte Carlo al 99,73 %; tiempo.

**Puerta 2.** Criterios propuestos, a confirmar por el grupo:

1. La cota de segundo orden contiene la solución AC en todos los vértices probados y en todas las muestras,
   en todos los sistemas, en el punto de diseño.
2. Su ancho no supera 1,5 veces el rango de referencia, en mediana.
3. El mejor operador del grupo mejora a AF2 en una cantidad que valga la pena declarar. Referencia de la
   sección 17: 1,17 y 1,51 veces en el ancho del resto, más la eliminación del lado positivo espurio.

*Entregable:* sección nueva en `RESULTS.md` con las tres respuestas. Hecho: sección 18.

**Resultado.**

1. **No.** Una cota de segundo orden ajustada queda corta exactamente por el término de tercer orden. En IEEE-33
   son 6·10⁻⁶ p.u.; en PEGASE con δ = 5 %, hasta 0,005 p.u., y decenas de barras por más de 0,001 p.u. Por
   muestreo casi no falla (0,01 % de las muestras como máximo, contra 24 a 51 % de la caja afín con K = 8).
2. **Sí.** Entre 1,01 y 1,07 veces el rango AC de los vértices, en mediana.
3. **Sí, según el número de fuentes.** Con hasta 8 zonas, RAA-K es la mejor cota (exacta con K = 4 en IEEE-33)
   y supera a los bloques espectrales; con 16 o más pierde ventaja. La ganancia está en el lado superior: el
   inferior, que es el que manda, ya es casi exacto con AF2.

**Consecuencia.** El punto 3.4 pasa a ser obligatorio y la decisión 0.1 hay que tomarla sabiendo esto: o se
controla el tercer orden, o el paper declara una cota de segundo orden con el déficit medido y su ley en δ.

---

## Parte 3 — Método completo

- [ ] **3.1 (código) Pasar el modelo de segundo orden al paquete**, con pruebas: coincidencia con las
      sensibilidades polares; consistencia de H_i con la evaluación término a término; el error después de sumar
      la forma cuadrática se multiplica por 8 cuando δ se duplica (orden δ³).
- [ ] **3.2 (código) Salida de un solo lado.** El resto es negativo, así que la cota superior queda en el radio
      afín y solo se corrige la inferior. *Hecho cuando:* está resuelto qué hacer en las barras donde H_i no es
      semidefinida negativa (en PEGASE hay casos con hasta 11,6 % de parte positiva). Con fuentes por zona
      son la mayoría: H_i es semidefinida negativa solo en el 7 a 42 % de las barras de PEGASE (sección 18.3).
- [ ] **3.3 (código) Solo si se mantiene el caso de miles de símbolos (0.2a): elegir los bloques sin formar H_i.**
      Candidatos: descomposición por diferencias de rama, coordenadas polares con resto acotado, autovectores
      dominantes sin matriz. Es el problema abierto y puede no tener solución barata.
- [ ] **3.4 (código) Tercer orden. Obligatorio tras la Puerta 2.** Hoy queda sin explicar el 3 % del exceso en
      mediana y el 17 % en el percentil 90 (PEGASE), y eso basta para que la cota no contenga. Según 0.1: (a) punto fijo r = −J⁻¹q(Cε + r) con verificación de contracción, o
      (b) margen medido y declarado como límite.
- [ ] **3.5 (código) Cantidades de rama** (si 0.4 las incluye): caída de voltaje, corriente y flujos como formas
      bilineales de las formas afines, acotadas con los mismos operadores.
- [ ] **3.6 (código) Redondeo dirigido** (solo con 0.1a): coeficientes intervalares y residuo riguroso.
- [ ] **3.7 (código) Slack distribuida** (si 0.3 lo decide): en la linealización y en el solver AC por lotes.
      *Avance:* hay un redespacho afín (`participation`), en que la generación sigue a la demanda en partes
      fijas y la slack cubre su parte y las pérdidas. No reparte las pérdidas.
- [ ] **3.8 (grupo) Capa probabilística.** Decidir si entra. Recomendación: fuera, o una subsección corta con el
      radio de Hoeffding; mezclar garantías fue una de las críticas al paper anterior.

---

## Parte 4 — Validación

- [ ] **4.1 (código) Vértices extremos, todas las barras, ambos lados, todos los sistemas.**
- [ ] **4.2 (código) Búsqueda adversarial.** Maximizar la violación con el flujo AC, partiendo del minimizador
      del modelo de segundo orden. *Hecho cuando:* se reporta el peor punto encontrado y si queda dentro.
- [ ] **4.3 (código) Monte Carlo AC**, 20.000 muestras por caso.
- [ ] **4.4 (código) Barrido carga × δ**, igual que C1, hasta el 95 % del límite de cargabilidad. Reportar desde
      qué combinación deja de valer la cota y la tasa de no convergencia.
- [ ] **4.5 (código) Escala.** Los seis sistemas de C7 (33 a 9.241 barras): tiempo, memoria y ancho.
- [ ] **4.6 (código) Sensibilidad a K y a la definición de las zonas.**

*Hecho cuando:* cada afirmación de validez del manuscrito apunta a una fila de un CSV.

---

## Parte 5 — Comparaciones

- [ ] **5.1 (código) El método publicado más cercano, reimplementado: Liao et al. (2017)**, con sus intervalos
      sobre e y f. En segundo lugar, un flujo afín: Vaccaro (2010) o Coletta (2019). Es la comparación que
      decide si hay mejora. *Hecho cuando:* reproduce las Tablas 1 y 2 de Liao (IEEE 9).
- [ ] **5.2 (código) Caja afín de primer orden** (lo que hay hoy).
- [ ] **5.3 (código) Caja 3σ**, llamada por su nombre. Decidir si se implementa la variante AC completa de la
      Ref. [10] en un sistema chico.
- [ ] **5.4 (código) Flujo por intervalos con Krawczyk** (opcional).
- [ ] **5.5 (código) Referencias de verdad:** rango exacto por vértices (K chico) y rango de Monte Carlo.
- [ ] **5.6 (código) Tabla única:** tipo de garantía, ancho sobre el rango de referencia, contención en vértices
      y en muestras, tiempo. Mismas condiciones para todos.

---

## Parte 6 — Manuscrito (Puerta 3: revisión interna)

- [ ] **6.1 (texto) Lista de afirmaciones**, cada una con su CSV y su límite. Sin fila, no entra.
- [ ] **6.2 (texto) Estructura.** Forma afín por sensibilidades dispersas; resto cuadrático exacto; cota con el
      operador elegido; salida de un solo lado; validación; comparaciones; límites.
- [ ] **6.3 (texto) Sección de límites explícita:** orden del resto, slack, modelo de incertidumbre, tamaño de K,
      qué pasa cerca del límite de cargabilidad.
- [ ] **6.4 (texto) Figuras:** exceso de la caja afín en los vértices; ancho de las cotas; escala. Todas
      regenerables con `experiments/figures.py`.
- [ ] **6.5 (texto) Requisitos de la revista:** leer la guía para autores vigente de EPSR (extensión, formato,
      puntos destacados, declaración de datos) antes de dar formato.
- [ ] **6.6 (grupo) Disponibilidad del código** y de los operadores no publicados.

**Puerta 3.** Una persona del grupo que no haya hecho los experimentos revisa la lista de 6.1 contra los CSV,
con las críticas de los revisores de TPWRS en la mano.

---

## Parte 7 — Entregables

- [ ] `RESULTS.md` con una sección por parte (2, 4 y 5), generado desde los CSV.
- [ ] Paquete con el modelo de segundo orden, los operadores y sus pruebas.
- [ ] Tabla de trabajos relacionados (Parte 1).
- [ ] Lista de afirmaciones con su evidencia (6.1).
- [ ] Lista de puntos bloqueados y qué se necesita para destrabarlos.

---

## Parte 8 — Prueba crítica de los cinco aportes posibles

**Estado: hecha una primera vuelta.** Detalle en las secciones 20 a 24 de RESULTS.md y en la segunda pasada de
RELATED_WORK.md. Antes de cada prueba se revisó el código contra la derivación.

- [x] **8.0 Auditoría del modelo de segundo orden.** Pasado al paquete (punto 3.1) con pruebas independientes:
      H_i es la segunda derivada de la solución AC (diferencias finitas) y el error restante escala como δ³.
- [x] **8.1 Aporte 2, las cotas no encierran.** Confirmado con pandapower (diferencia de 10⁻¹⁰ p.u.). El mínimo
      está en un vértice; el máximo, no siempre. La fórmula de Liao falla en su propio IEEE 9. **Se sostiene.**
- [x] **8.2 Aporte 5, peor escenario.** 95 a 98 % de aciertos en PEGASE. En redes radiales ya está resuelto en
      la literatura. **Parcial.**
- [x] **8.3 Aporte 4, formulación.** Las vías de Liao y adjunta dan los mismos coeficientes; la adjunta no ahorra
      con pocas fuentes. La Tabla 1 de Liao no se pudo reproducir (punto 5.1, bloqueado: el artículo no define
      las variables ni el paso a \|V\|). **Parcial.**
- [x] **8.4 Aporte 3, RAA-K.** Mejor que AF2; peor que la relajación semidefinida en ajuste, mejor en costo y
      en que admite redondeo dirigido. **Parcial.**
- [x] **8.5 Aporte 1, punto fijo (punto 3.4a), primera versión.** La construcción ya está publicada (Zhu, 2026).
      Cierra en IEEE-33 hasta δ = 15 % (encierro a 1,002–1,009 veces el rango exacto) y en 118 barras solo con
      δ = 2 %; no cierra en 300 ni en 2.869. **No escala.**

Pendiente de esta parte:

- [ ] **8.6 Punto fijo de tercer orden** (restar también la forma cúbica exacta): debería extender el rango en
      sistemas chicos; no arregla la escala, porque el límite ahí es la matriz G.
- [ ] **8.7 Argumento de contracción sin valores absolutos elemento a elemento** (norma adaptada o diferencias
      por rama, en polares). Es el problema abierto; en PEGASE con δ = 5 % la variación real del Jacobiano tiene
      radio 0,28 y la cota componente a componente da 1,7 o más. Resuelto solo para redes sin barras PV, con
      la forma de inyección de corriente (punto 9.3); con barras PV sigue abierto.
- [ ] **8.8 Redondeo dirigido** de todo lo anterior, para poder decir "verificado" (punto 3.6).
- [ ] **8.9 Leer completos** los trabajos de de Sousa y Melo y Sun et al. (2021). Sou et al. (2026) ya está leído
      (RELATED_WORK.md): cubre el encierro de voltajes en redes radiales con generación distribuida, hasta
      10.476 barras. La opción de "encierro verificado en distribución" queda acotada a lo que ese trabajo no
      cubre: redes malladas, elementos en derivación, barras de voltaje controlado.

**Consecuencia para la decisión 0.1.** El camino (a), encierro verificado, es realista hoy solo para
alimentadores de distribución. A escala de transmisión lo defendible es el camino (b): cota de segundo orden
con el déficit medido.

---

## Parte 9 — Opción A: encierros de voltaje en redes de distribución

Elegida el 15-11-2026. Detalle y cifras en la sección 25 de RESULTS.md; bibliografía en la tercera pasada de
RELATED_WORK.md.

- [x] **9.1 Línea base radial.** Método de Sou, Khong y Sandberg (2026) reimplementado (`cdhim/distflow.py`) y
      validado contra su Tabla I. Sus Tablas VI y VII se reproducen solo en parte: la regla de montaje es ambigua
      y la causa de la diferencia en 33 y 69 barras no está identificada.
- [x] **9.2 Línea base para cualquier topología.** Discos de Wang, Bernstein, Le Boudec y Paolone (2018),
      Teorema 1, y su refinamiento por barra (`cdhim/discs.py`).
- [x] **9.3 Cota del resto que no se degrada con líneas cortas.** Forma de inyección de corriente
      (`InjectionEnclosure`), comprobada pieza por pieza contra el flujo AC. Reemplaza en redes PQ a la cota de la
      sección 22, que es la que fallaba (cubre el punto 8.7 solo para redes sin barras PV).
- [x] **9.4 Comparación en redes radiales.** Hecha (sección 25.3). **Sin aporte ahí:** WCLF es más ajustado con
      incertidumbre grande, da el mínimo exacto y escala.
- [x] **9.5 Comparación en redes malladas y con elementos en derivación.** Hecha (sección 25.4) en cuatro redes
      de hasta 179 barras. **Es donde hay algo medible.** En la red CIGRE de media tensión el encierro se
      degrada con incertidumbre grande y deja de cerrar en el nivel más alto. En CIGRE y Oberrhein el voltaje
      mínimo no está en el punto que supone la monotonía.
- [x] **9.6 Redondeo dirigido** (es el punto 8.8). Hecho: `cdhim/verified.py`, verificación a posteriori con
      existencia y unicidad; las cotas dirigidas se comparan con aritmética racional exacta en las pruebas
      (sección 26.3). Queda una reserva: no usa una biblioteca de intervalos certificada (punto 9.19).
- [x] **9.7 Referencia exacta.** Hecho: cota interior certificada en todos los casos y rango exacto certificado
      por bisección en cajas de pocos símbolos (`cdhim/certified.py`, sección 26.4). No hay rango exacto para los
      semidiscos ni para muchos símbolos.
- [ ] **9.8 Leer el artículo de Wang, Le Boudec y Paolone (2019).** **Bloqueado: falta el PDF** (no hay preprint
      y está tras muro de pago). Por láminas, resumen y tesis, su cota por nodo es la de `cdhim/discs.py` y su vía
      al ajuste es subdividir el conjunto; eso se midió (sección 26.2).
- [x] **9.9 Redes desbalanceadas.** Hecho con supuestos: el encierro se aplica sin cambios a nodos-fase con
      cargas en estrella (`cdhim/threephase.py`, sección 26.7), validado contra `runpp_3ph`. Fuera: cargas en
      delta, transformadores, neutro explícito.
- [x] **9.10 Barras de voltaje controlado.** Hecho en distribución, límite medido en transmisión: forma
      voltaje-corriente (`cdhim/verified_vi.py`, sección 26.6). Verifica con fuentes PV en la red de 33 barras y
      en la de 30; no en las de 118 y 300.
- [x] **9.11 Escala.** Medida (sección 26.5): hasta 2.940 barras con 64 símbolos y hasta 533 con 108. El límite
      es el número de símbolos (memoria m·K²), después las matrices densas de 2m × 2m.
- [x] **9.12 Acotar juntos el primer y el segundo orden.** Hecho: cota conjunta por tangentes de la parte
      cóncava (sección 26.3). Bajar la cota del resto de tercer orden sigue abierto.
- [x] **9.17 Baseline de Krawczyk.** Hecho (`cdhim/krawczyk.py`, sección 26.1). Por inyección de corriente cierra
      en los mismos casos y queda cerca; el encierro es más ajustado en todas las barras y más caro.
- [ ] **9.18 Cota que cierre con barras PV a mayor escala.** En la red de 118 barras la parte lineal contrae y
      lo que falla es el mayorante por sumas de fila de la forma voltaje-corriente: habría que conservar las
      cancelaciones de J⁻¹ΔJ. En la de 300 barras ni la parte lineal contrae.
- [ ] **9.19 Segunda implementación del redondeo** con una biblioteca de intervalos certificada (por ejemplo
      Arb), para no depender de la disciplina de redondeo escrita a mano.
- [ ] **9.20 Revisar si la cota conjunta por tangentes ya existe con otro nombre** (es corriente en optimización
      global) antes de presentarla como propia.
- [ ] **9.21 Residuo con más precisión en la forma voltaje-corriente:** falla en redes PQ con ramas muy cortas
      (141 barras) porque la cota de Y·(I − C·J) crece con las admitancias.
- [ ] **9.23 Bloques de RAA-K dentro de la cota conjunta.** Con la cota conjunta, RAA-K ya no mejora el ancho en
      los casos de la escala (sección 26.5). Falta probar si sus bloques, usados como tangentes en vez de los
      autovectores, dan una cota conjunta igual de buena y más barata o dispersa; sería la forma de que el
      teorema k-nomial del grupo siga dentro del método.
- [ ] **9.22 Detalles de carga de redes** encontrados en la escala: `load_grid` supone la tabla de barras
      ordenada; el lector depende de la codificación del sistema; la URL de MATPOWER está fija en 7.1.
- [ ] **9.13 Decidir el conjunto de incertidumbre del artículo.** Semidiscos con demanda determinista es el caso
      de Sou et al.; demanda incierta por zonas es el caso natural de la aritmética afín (sección 25.5).
- [ ] **9.14 Rehacer las secciones 18 a 23 con RAA-K normalizada.** El operador usa una tolerancia absoluta y con
      las formas del flujo de potencia pierde bloques (sección 25.8): esas secciones lo subestiman. Hay que
      volver a correr `x2` a `x7` pasando `normalise=True`, revisar el texto y el veredicto del aporte 3. Son
      varias horas de cómputo. `x6` hay que regenerarlo de todos modos: `Enclosure.solve` ahora infla el punto
      fijo y comprueba el auto-mapeo, lo que puede mover sus cifras en los últimos dígitos.
- [ ] **9.15 Avisar al grupo** de esa dependencia de escala en `ibex_KNomial.h` (decisión de los autores).
      El 15-11-2026 se consultó en una sesión de trabajo sobre `ibex-lib`. Según esa respuesta (no la verifiqué
      yo), el efecto se reproduce en el binario C++ y se evalúa un parche 0015 con umbrales relativos.
- [ ] **9.16 Si el grupo cambia los umbrales en `ibex_KNomial.h`,** actualizar `cdhim/knomial.py`, regenerar los
      casos dorados y volver a correr `tests/test_knomial.py`; con umbrales relativos `normalise=True` deja de
      hacer falta. Ojo: el umbral de la línea 156 se aplica a un producto de dos diagonales, así que su versión
      relativa lleva la escala al cuadrado.

---

## Orden sugerido

Las estimaciones son gruesas y mías; dependen sobre todo de 0.1, 0.2 y 0.7.

| Paso | Qué | Depende de | Esfuerzo estimado |
|---|---|---|---|
| 1 | Parte 0, decisiones | — | una reunión |
| 2 | Parte 1, bibliografía, en paralelo con el paso 3 | — | 1 a 2 semanas |
| 3 | Parte 2, prueba con pocas fuentes | 0.2, 0.7 | 1 semana |
| — | **Puertas 1 y 2** | pasos 2 y 3 | — |
| 4 | 3.1, 3.2 y 3.4 | Puerta 2 | 1 a 2 semanas |
| 5 | 5.1, referencia de la literatura | Puerta 1 | 1 a 2 semanas |
| 6 | Parte 4, validación | paso 4 | 1 semana |
| 7 | Resto de la Parte 5 | pasos 4 y 5 | 1 semana |
| 8 | Opcionales: 3.3, 3.5, 3.6, 3.7 | según Parte 0 | 1 a 4 semanas cada uno |
| 9 | Parte 6, manuscrito | pasos 6 y 7 | 2 a 3 semanas |

El paso 3 es el que más informa por lo que cuesta: reutiliza casi todo lo que ya existe y contesta si el método
cierra.

---

## Plan B

- **Si falla la Puerta 1** (alguien ya lo publicó): reposicionar como estudio de escala y de reproducibilidad sobre un
  método existente, o llevar la aplicación como sección del paper del k-nomial.
- **Si falla el criterio 1 de la Puerta 2** (la cota no contiene): el tercer orden pesa; pasar a 3.4(a) antes de
  seguir, o restringir el rango de δ y de carga que se declara.
- **Si falla el criterio 2** (cota demasiado ancha): revisar la elección de bloques; si no mejora, el aporte se
  reduce al hallazgo de que la caja afín no encierra la solución AC y a su cuantificación.
- **Si falla el criterio 3** (el operador del grupo no mejora a AF2): el paper de sistemas de potencia se puede
  escribir con AF2 ensamblada, pero deja de ser una aplicación de los operadores del grupo. Decidir si aun así
  interesa.
