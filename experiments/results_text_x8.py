"""Section 25 of RESULTS.md: option A, step 1 (verified voltage enclosures on distribution networks)."""

X8_TEXT = r"""
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
  semidiscos), los nueve flujos de frontera del método de Sou et al. y <<X8_SAMPLES>> inyecciones aleatorias por
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

{{x8_validation}}

*Script: `experiments/x8_feeders.py` · CSV: `results/x8_wclf_validation.csv`*

- **La Tabla I se reproduce** a los dígitos publicados (<<X8_T1_OK>>). Es el único caso con el montaje
  completamente especificado: las cinco fuentes y sus potencias están dibujadas en su Fig. 10.
- **Las Tablas VI y VII se reproducen en <<X8_T67_HIT>> de <<X8_T67_N>> filas** (redes: <<X8_T67_NETS>>). En las
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

{{x8_radial}}

*Script: `experiments/x8_feeders.py` · CSV: `results/x8_feeders.csv` · semilla por fila en el CSV. "Exceso" =
(cota superior − máximo alcanzado) / máximo alcanzado; "holgura" = (mínimo alcanzado − cota inferior) / mínimo
alcanzado; ambos promediados sobre las barras. Los tiempos son de una sola corrida en Python.*

- **El punto fijo cierra en <<X8_RAD_FP>> de los <<X8_RAD_N>> casos radiales**, incluido ρU = 80 %, donde el
  rango de voltaje alcanzado llega a <<X8_WIDTH_MAX>> p.u. La cota de la sección 22 dejaba de cerrar mucho antes (sección 25.6).
- **Ninguna cota queda del lado equivocado de un valor alcanzado**: <<X8_VIOL>> violaciones entre los tres
  métodos en toda la sección.
- **Cota superior.** La del encierro es más ajustada que la de WCLF en <<X8_WIN_ALL>> de <<X8_CNT_ALL>> casos:
  <<X8_WIN_R>> de <<X8_CNT_R>> en flujo inverso y <<X8_WIN_F>> de <<X8_CNT_F>> en flujo directo. Con
  incertidumbre grande (ρU ≥ 50 %) en flujo directo gana WCLF en <<X8_LOSE_LARGE_F>> de <<X8_CNT_LARGE_F>>
  casos; con incertidumbre chica (ρU ≤ 10 %) en flujo inverso gana el encierro en <<X8_WIN_SMALL_R>> de
  <<X8_CNT_SMALL_R>>. El peor caso para el encierro es la red de <<X8_WORST_NET>> barras con
  ρU = <<X8_WORST_RHO>> %: exceso de <<X8_WORST_EXC>> % contra <<X8_WORST_WCLF>> % de WCLF, una razón de
  <<X8_BEST_FOR_WCLF>>. En el otro extremo la razón a favor del encierro llega a <<X8_BEST_FOR_US>>.
- **El "error certificado" de WCLF lo domina su cota inferior.** Esa cota sale de nueve flujos con la misma
  dirección para todas las fuentes; la búsqueda local encuentra máximos mayores, y esa diferencia explica entre
  <<X8_GAP_SHARE_LO>> y <<X8_GAP_SHARE_HI>> % del error que el método certifica. Su cota superior sola es bastante
  mejor que lo que su propia métrica sugiere.
- **Cota inferior.** WCLF da el mínimo con un solo flujo en (Δp, Δq) = (0, −C), suponiendo monotonía. En estos
  cuatro alimentadores la suposición se cumple: su valor coincide con el mínimo alcanzado en todas las barras
  (diferencia máxima <<X8_VMIN_GAP>> p.u.). El encierro da una cota inferior que no depende de esa suposición, con
  holgura de hasta <<X8_LOW_SMALL_MAX>> % para ρU ≤ 10 % y de hasta <<X8_LOW_RAD_MAX>> % en el peor caso. En otras
  redes la suposición no se cumple (sección 25.4).
- **Discos de Wang et al.** Aplican en todos los casos salvo <<X8_DISC_FAIL>>. El disco uniforme es entre
  <<X8_DISCU_LO>> y <<X8_DISCU_HI>> veces el rango alcanzado y el refinado por barra entre <<X8_DISC_LO>> y
  <<X8_DISC_HI>> (todas las redes de la sección).
- **Costo.** El encierro tarda entre <<X8_MS_ENC_LO>> ms y <<X8_MS_ENC_HI_S>> s (todas las redes de la sección) y
  la reimplementación de WCLF entre <<X8_MS_WCLF_LO>> y <<X8_MS_WCLF_HI>> ms. Sus autores reportan 3 a 15 ms en estas redes y 0,53 s en 10.476
  barras; el encierro usa matrices densas y aquí no pasa de <<X8_OB_BUS>> barras.

**Lectura.** En redes radiales con este tipo de incertidumbre, WCLF sigue siendo la referencia: es más ajustado
cuando la incertidumbre es grande, da el mínimo exacto y escala a diez mil barras. Donde el encierro gana
(incertidumbre chica, flujo inverso) las diferencias son de milésimas de punto porcentual de voltaje, sin
consecuencia práctica. **Aquí no hay aporte.**

### 25.4 Redes malladas y con elementos en derivación: donde WCLF no aplica

{{x8_general}}

*Script: `experiments/x8_feeders.py` · CSV: `results/x8_feeders.csv`*

- **El encierro cierra en <<X8_GEN_FP>> de los <<X8_GEN_N>> casos.** No cierra en: <<X8_GEN_NOFP>>. Donde
  cierra queda entre <<X8_W_GEN_LO>> y <<X8_W_GEN_HI>> veces el rango alcanzado (mediana por caso;
  <<X8_WMAX_GEN>> en la peor barra).
- **Se degrada cuando el rango de voltaje es grande.** Con ρU ≤ 20 % el encierro queda a lo más a
  <<X8_W_GEN_SMALL_HI>> veces el rango en las cuatro redes. El peor caso es la red CIGRE, donde la regla de
  montaje pone toda la potencia en dos fuentes: con ρU = 50 % el rango alcanzado llega a <<X8_CIGRE_WIDTH>> p.u.,
  el encierro queda a <<X8_CIGRE_W>> veces ese rango, y con ρU = 80 % ni el encierro ni los discos cierran.
- **La cota de segundo orden sola no alcanza:** sin la holgura del punto fijo, la cota inferior deja fuera algún
  valor alcanzado en <<X8_2ND_LOW>> de los <<X8_N_OK>> casos de la sección (por hasta <<X8_2ND_LOW_MAX>> p.u.), y
  la superior en <<X8_2ND_UP>>. Es el mismo hallazgo de las secciones 18 y 20, ahora en distribución.
- **La garantía publicada para estas redes es bastante más ancha.** Los discos de Wang et al. contienen todo lo
  alcanzado, pero en estas cuatro redes su ancho es de <<X8_DISC_GEN_LO>> a <<X8_DISC_GEN_HI>> veces el rango
  (refinados por barra) y de <<X8_DISCU_GEN_LO>> a <<X8_DISCU_GEN_HI>> veces sin refinar.
- **Tamaño.** La red Oberrhein tiene <<X8_OB_BUS>> barras y <<X8_OB_K>> símbolos; el encierro tarda hasta
  <<X8_OB_MS>> s. En un perfil aparte (sin CSV) casi todo ese tiempo se va en la cota RAA-K de las formas
  cuadráticas.

**La monotonía no se cumple en todas las redes.** El mínimo de voltaje de cada barra, ¿está en (0, −C)?

{{x8_monotone}}

*Script: `experiments/x8_feeders.py monotone` · CSV: `results/x8_monotonicity.csv`. Se compara el voltaje en
(0, −C) con el menor valor que encuentra la búsqueda local; uno menor lo alcanza una inyección factible.*

- En <<X8_MONO_OK_NETS>> el mínimo está en (0, −C) en todas las barras.
- En <<X8_MONO_NETS>> no: en <<X8_MONO_CASES>> de los <<X8_MONO_N>> casos hay barras con un mínimo menor, hasta
  <<X8_MONO_BUSES>> barras a la vez y por hasta <<X8_MONO_GAP>> p.u. (<<X8_MONO_GAP_PCT>> % del rango de esa
  barra). Son las dos redes con transformadores y cables con capacidad. La diferencia es chica, pero un método
  que tome el mínimo de un solo flujo deja de ser exacto ahí. La tesis de C. Wang reporta lo mismo en el IEEE 37.

### 25.5 Cajas de demanda: las dos cotas del resto

El montaje de la sección 22 (demanda incierta ±δ en 8 zonas, caja simétrica), ahora con las dos cotas del resto
y contra los 256 vértices.

{{x8_boxes}}

*Script: `experiments/x8_feeders.py` · CSV: `results/x8_load_boxes.csv`*

- La cota de la sección 22 cierra en <<X8_BOX_OLD>> de <<X8_BOX_N>> casos; la nueva, en <<X8_BOX_NEW>>, hasta
  δ = <<X8_BOX_DMAX>> %, con el encierro entre <<X8_BOX_W_LO>> y <<X8_BOX_W_HI>> veces el rango exacto en
  vértices y <<X8_BOX_VIOL>> barras con algún vértice fuera.

### 25.6 Qué cambió en el método

**Por qué fallaba la cota de la sección 22.** Escribía el resto en coordenadas rectangulares y acotaba el término
cuadrático con valores absolutos de J⁻¹ΔJ. Cuando dos barras están unidas por una rama muy corta, esa cota
trata sus voltajes como independientes y la corriente entre ellas se dispara. Era un defecto de la cota, no del
problema: deja de cerrar aunque la parte lineal contraiga con holgura (radio de G desde <<X8_RB_RADIUS_LO>> en
los casos en que falla).

{{x8_remainder}}

*Script: `experiments/x8_feeders.py remainder` · CSV: `results/x8_remainder_bounds.csv`. Mismos casos que las
tablas de 25.3 y 25.4.*

- La cota de la sección 22 cierra en <<X8_RB_OLD>> de <<X8_RB_N>> casos y la nueva en <<X8_RB_NEW>>. La de la
  sección 22 no cierra en ningún caso de: <<X8_RB_NEVER>>.
- Donde cierran las dos, la cota del tercer orden de la sección 22 es entre <<X8_RB_RATIO_LO>> y
  <<X8_RB_RATIO_HI>> veces la nueva.

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
   <<X8_W_GEN_LO>>–<<X8_W_GEN_HI>> veces el rango alcanzado, donde WCLF no aplica y la garantía publicada es de
   <<X8_DISC_GEN_LO>> a <<X8_DISC_GEN_HI>> veces el rango. En lo que se pudo leer, nadie reporta esa comparación en una red de distribución mallada. La ventaja
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
   - escala: todo es denso y se probó hasta <<X8_OB_BUS>> barras.
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

{{x8_scale}}

*Script: `experiments/x8_feeders.py` · CSV: `results/x8_raak_scale.csv`*

- **Normalizada, la cota es la misma fracción de AF2 para todo δ**, como dice la teoría: entre <<X8_SC_NRM_LO>> y
  <<X8_SC_NRM_HI>>. Sin normalizar va de <<X8_SC_RAW_LO>> a <<X8_SC_RAW_HI>>, y con δ chico no mejora nada a AF2.
- **En los δ de diseño**, que son los que usan las secciones 18 a 23, la cota sin normalizar es
  <<X8_SC_DES_RAW_LO>> a <<X8_SC_DES_RAW_HI>> veces AF2 y la normalizada <<X8_SC_DES_NRM_LO>> a
  <<X8_SC_DES_NRM_HI>>. Frente al rango exacto en vértices, con 8 zonas: <<X8_SC_VR_RAW_LO>> a
  <<X8_SC_VR_RAW_HI>> sin normalizar y <<X8_SC_VR_NRM_LO>> a <<X8_SC_VR_NRM_HI>> normalizada. Hay formas que
  se estrechan hasta <<X8_SC_NARROW_MAX>> %.
- **No es un problema de validez.** Las dos cotas son válidas; la normalizada contiene los valores en los vértices
  en todas las formas con 8 zonas (<<X8_SC_SHORT>> excepciones). La búsqueda es voraz y no es monótona: en
  <<X8_SC_WIDER>> de las <<X8_SC_FORMS>> formas la cota normalizada sale más ancha que la otra.
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
"""
