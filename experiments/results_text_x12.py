"""Section 26 of RESULTS.md: the checks left open by section 25 (X9, X10, X12, X13, X14)."""

X12_TEXT = r"""
---

## 26. Las comprobaciones que faltaban para sostener la opción A

La sección 25.7 listó lo que faltaba. Esta sección lo hace, o dice por qué no se pudo.

| Pendiente de 25.7 | Estado | Resultado en una línea | Dónde |
|---|---|---|---|
| Comparar con Krawczyk | Hecho | Un Krawczyk bien armado cierra en los mismos casos y queda cerca; el encierro es más ajustado en todas las barras, por un orden de exactitud, y más caro | 26.1 |
| Leer Wang, Le Boudec y Paolone (2019) | **Bloqueado: falta el PDF** | Por fuentes indirectas, su cota es la ya implementada y su vía al ajuste es subdividir; eso se midió | 26.2 |
| Redondeo dirigido | Hecho | El encierro está verificado en punto flotante, con unicidad, sin perder ancho | 26.3 |
| Referencia exacta o certificada | Hecho | Cota interior certificada en todos los casos; rango exacto certificado en cajas de pocos símbolos | 26.4 |
| Escala | Hecho, con un límite medido | Verifica hasta <<X13_BUS_MAX>> barras con <<X13_K_AT_BUS_MAX>> símbolos y hasta <<X13_BUS_AT_K_MAX>> con <<X13_K_MAX>>; lo que limita es el número de símbolos | 26.5 |
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

{{x9_variants}}

*Script: `experiments/x9_krawczyk.py` · CSV: `results/x9_krawczyk.csv`. Los 80 casos de la sección 25.*

- **Por desbalance de potencia no cierra en ningún caso** (<<X9_POWER>> de 320 corridas). Trata como independientes
  los voltajes de barras vecinas; es el mismo defecto de la cota de la sección 22.
- **Por inyección de corriente cierra:** en <<X9_CUR_SLOPE>> de <<X9_N>> casos con pendientes y en <<X9_CUR_DER>>
  con derivadas. Con pendientes falla en los mismos dos casos que el encierro.
- **Ninguna violación** en ninguna variante (<<X9_VIOL>>).
- Con demanda incierta en 8 zonas (`results/x9_krawczyk_load_boxes.csv`) cierra en <<X9_LB_CLOSES>> de
  <<X9_LB_N>> casos y mide entre <<X9_LB_LO>> y <<X9_LB_HI>> veces el rango exacto en vértices.

**Frente al encierro verificado**, con la mejor variante (corriente, pendientes, ejes rotados) y contra la misma
referencia certificada de 26.4:

{{x12_krawczyk}}

*Script: `experiments/x12_verified.py verified` · CSV: `results/x12_verified.csv`.*

- **Krawczyk no es más ajustado que el encierro en ninguna cota de ninguna barra** (<<X12_KR_TIGHTER>>).
- **Su intervalo es entre <<X12_KR_OVER_LO>> y <<X12_KR_OVER_HI>> veces el del encierro** (mediana por caso). En
  exceso sobre el rango exacto, que es lo que distingue a dos métodos ajustados, el de Krawczyk es entre
  <<X12_KR_EXC_LO>> y <<X12_KR_EXC_HI>> veces el del encierro.
- **El encierro es más caro:** <<X12_VER_MS>> ms de mediana frente a <<X12_KR_MS>> ms, y la diferencia crece con el
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

{{x12_subdivision}}

*Script: `experiments/x12_verified.py subdivision` · CSV: `results/x12_subdivision.csv`. A lo más 3.000 encierros
por caso y método. Krawczyk y el disco no controlan el redondeo: sus columnas valen en aritmética exacta.*

- **Encierros necesarios para llegar a 1 %:**
  - el verificado: <<X12_SUB_V>>, es decir, ya llega sin subdividir;
  - Krawczyk: <<X12_SUB_K>> en los casos en que llega; en <<X12_SUB_K_FAIL>> de <<X12_SUB_N>> no llega dentro del
    presupuesto y queda a <<X12_SUB_K_END>> %;
  - el disco: <<X12_SUB_D_LO>> a <<X12_SUB_D_HI>> en los casos en que llega; en <<X12_SUB_D_FAIL>> de <<X12_SUB_N>>
    no llega.
- La razón es el orden del exceso en el tamaño de la hoja: primero para el disco, segundo para Krawczyk y tercero
  para el encierro.
- Sin subdividir: el disco mide entre <<X12_SUB_D_ROOT_LO>> y <<X12_SUB_D_ROOT_HI>> veces el rango exacto,
  Krawczyk entre <<X12_SUB_K_ROOT_LO>> y <<X12_SUB_K_ROOT_HI>>, y el encierro entre <<X12_SUB_V_ROOT_LO>> y
  <<X12_SUB_V_ROOT_HI>>.
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

{{x12_verified}}

*Script: `experiments/x12_verified.py verified` · CSV: `results/x12_verified.csv`. Los 80 casos de la sección 25
y los 10 de la red trifásica de 26.7.*

- **Cierra en <<X12_VER>> de <<X12_N>> casos de la sección 25 y prueba unicidad en <<X12_UNI>>.** El encierro de
  esa sección cerraba en <<X12_X8>>. <<X12_FAIL>>
- **El rigor no cuesta ancho.** Con la cota separada, el encierro verificado mide entre <<X12_SEPW_LO>> y
  <<X12_SEPW_HI>> veces el de la sección 25: sale igual o más angosto porque la cota del resto es otra (el radio
  del tubo es entre <<X12_TUBE_LO>> y <<X12_TUBE_HI>> veces la cota del resto de esa sección).
- **La cota conjunta lo estrecha más.** En las redes de la sección 25 la cota superior sale de ella en
  <<X12_JOINT_UP>> y la inferior en <<X12_JOINT_LO>>; en la red trifásica, las dos. El ancho sobre el rango exacto
  baja de <<X12_SEP_LO>>–<<X12_SEP_HI>> a <<X12_WC_LO>>–<<X12_WC_HI>> (mediana por caso).
- **El redondeo no pesa.** El radio del tubo, redondeo incluido, es a lo más <<X12_TUBE_W>> % del mayor ancho
  alcanzado del caso; el piso numérico se mide aparte en 26.5.
- **Costo.** La verificación tarda entre <<X12_MS_LO>> ms y <<X12_MS_HI>> s por caso, en una corrida a solas. En
  los casos lentos casi todo es la cota RAA-K de las formas cuadráticas, que resulta prescindible (26.5).

### 26.4 Referencia certificada

Hasta aquí el ancho se comparaba con valores alcanzados por búsqueda. Ahora hay dos referencias con garantía
(`cdhim/certified.py`).

**Cota interior certificada.** El mejor punto encontrado para cada barra se verifica: se comprueba en aritmética
racional que pertenece al conjunto y se encierra su solución del flujo de potencia con el mismo verificador,
sin incertidumbre. El extremo inferior de ese encierro es una cota inferior demostrada del máximo exacto, y lo
análogo para el mínimo. Entonces el cociente entre el ancho del encierro y el ancho entre esas dos cotas es
una **cota superior demostrada** de cuánto sobreestima el encierro.

- Puntos que no se pudieron certificar: <<X12_PTS_BAD>>. El valor certificado difiere del alcanzado en a lo más
  <<X12_CERT_GAP>> p.u.
- **Sobreestimación demostrada:** el encierro mide a lo más entre <<X12_WC_LO>> y <<X12_WC_HI>> veces el rango
  exacto (mediana por caso; la peor barra, <<X12_WC_MAX>>). Con ρU de hasta 20 %, a lo más <<X12_WC_20>>.
- **Dos demostraciones independientes no se contradicen:** en <<X12_CONTRA>> barras el encierro deja fuera un
  punto certificado.

**Rango exacto.** En cajas con pocos símbolos se puede acotar el rango exacto por los dos lados: por dentro con
puntos certificados y por fuera con la unión de los encierros verificados de una partición de la caja, que se
refina por bisección donde sobra. El exceso de cada encierro es de tercer orden en el tamaño de la caja, así que
bastan pocas cajas.

{{x12_exact}}

*Script: `experiments/x12_verified.py exact` · CSV: `results/x12_exact_range.csv`. La bisección se detiene cuando
ninguna hoja excede los valores alcanzados en más de 10⁻⁴ del ancho, o tras 3.000 encierros.*

- En los <<X12_EX_N>> casos el rango exacto queda determinado con un error relativo de a lo más
  <<X12_EX_GAP>> %, con entre <<X12_EX_BOX_LO>> y <<X12_EX_BOX_HI>> encierros y a lo más <<X12_EX_SEC>> s.
- **Demanda incierta en 8 zonas:** el encierro mide entre <<X12_EX_LOAD_LO>> y <<X12_EX_LOAD_HI>> veces el rango
  exacto.
- **Cajas envolventes de los semidiscos:** entre <<X12_EX_BB_LO>> y <<X12_EX_BB_HI>> veces.
- Límite: la bisección necesita cajas, así que no sirve tal cual para los semidiscos, y su costo crece con el
  número de símbolos. Para las redes grandes queda solo la cota interior certificada.

### 26.5 Escala

**Qué se midió** (`experiments/x13_scale.py`). El encierro verificado en redes PQ de distribución más grandes:
los casos de MATPOWER 7.1 de 118, 136 y 533 barras con sus enlaces cerrados, la de 141 barras, Oberrhein y la red
de baja tensión de Schutterwald de pandapower. Tres diseños de fuentes: una cada 10 barras, como en la sección
25 (K crece con la red), y 8 o 32 fuentes (K fijo), para separar el efecto del tamaño de la red del de K.
Flujos F y R con ρU de 10 y 50 %. La cota de los módulos se calcula de tres maneras: con la forma cuadrática
acotada término a término, con bloques de hasta 3 símbolos y con RAA-K completa.

{{x13}}

*Script: `experiments/x13_scale.py` · CSV: `results/x13_scale.csv`. Tiempos de una corrida a solas. El rango
alcanzado es el de flujos de frontera y muestras, sin búsqueda local, para que el tiempo no se lo lleve la
referencia.*

- **Verifica, con unicidad, en <<X13_VER>> de los <<X13_OK>> casos corridos:** hasta <<X13_BUS_MAX>> barras con
  <<X13_K_AT_BUS_MAX>> símbolos y hasta <<X13_BUS_AT_K_MAX>> barras con <<X13_K_MAX>>. Ninguna violación
  (<<X13_VIOL>>).
- **Lo que limita es el número de símbolos, no el de barras.** El modelo de segundo orden guarda un tensor de
  barras × K × K. Con una fuente cada 10 barras en la red de <<X13_BUS_MAX>> barras son <<X13_BLOCKED_K>> símbolos
  y un solo tensor ocupa <<X13_BLOCKED_GB>> GB: no se intentó (<<X13_BLOCKED_N>> casos).
- **Costo en la red más grande:** verificar el tubo toma de <<X13_TUBE_BIG_LO>> a <<X13_TUBE_BIG_HI>> s, por las
  matrices densas de 2m × 2m, y la memoria pico es de <<X13_RSS_MAX>> GB.
- **RAA-K es lo más caro y, con la cota conjunta, aquí no aporta nada.** La versión completa llega a
  <<X13_FULL_MAX>> s donde la cota término a término toma <<X13_TERM_AT_FULL_MAX>> s y la de bloques de hasta 3,
  <<X13_K3_AT_FULL_MAX>> s. El ancho que entregan, sobre el de la cota término a término: <<X13_FULL_GAIN>> la
  completa y <<X13_K3_GAIN>> la de bloques de hasta 3, en los <<X13_OK>> casos. Una explicación coherente con
  esto, que no se midió por separado: RAA-K ayudaba en la cota superior, que ahora sale de la cota conjunta, y
  para una forma cóncava su cota inferior coincide con la de término a término.
- **Sin RAA-K el encierro es mucho más barato.** En Oberrhein con <<X13_OB_K>> símbolos, modelo, tubo y módulos
  suman <<X13_OB_TERM>> s con la cota término a término y <<X13_OB_FULL>> s con RAA-K completa.
- **Piso numérico.** Sin incertidumbre, el radio del tubo es a lo más <<X13_FLOOR>> p.u., que es a lo más
  <<X13_FLOOR_SHARE>> % del mayor radio del caso. El redondeo no limita el tamaño de la red.
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

**Distribución con fuentes de voltaje controlado.** Las <<X14_D_PV>> barras con generación de la regla de la
sección 25 pasan a ser PV, con su potencia activa incierta (<<X14_D_K>> símbolos):

{{x14_distribution}}

*Script: `experiments/x14_pv.py` · CSV: `results/x14_pv_distribution.csv`.*

- **Verifica, con unicidad, en <<X14_D_VER>> de <<X14_D_N>> casos**, sin violaciones (<<X14_D_VIOL>>).
- **Es ajustado con incertidumbre moderada y deja de serlo con mucha.** Con ρU de hasta 20 % mide entre
  <<X14_D_W20_LO>> y <<X14_D_W20_HI>> veces el rango exacto en vértices; con 50 y 80 %, entre <<X14_D_WBIG_LO>> y
  <<X14_D_WBIG_HI>>, y hasta <<X14_D_WORST>> veces en la peor barra. Es peor que con barras PQ a igual ρU. No se
  aisló cuánto se debe a que el rango de las barras vecinas a una PV es muy chico y cuánto a la cota.

**Transmisión**, con demanda incierta en 8 zonas:

{{x14_transmission}}

*Script: `experiments/x14_pv.py` · CSV: `results/x14_pv_transmission.csv`. «Radio de la parte lineal»: radio
espectral de la suma de |J⁻¹ΔJ| de la sección 22; por encima de 1, ningún punto fijo con mayorante por entradas
puede cerrar.*

- **Verifica en <<X14_T_NETS_OK>>** (<<X14_T_VER>> de <<X14_T_N>> casos), entre <<X14_T_W_LO>> y <<X14_T_W_HI>>
  veces el rango exacto en vértices. **No verifica en <<X14_T_NETS_BAD>>, ni con δ = 2 %.**
- **En la red de 118 barras el límite es de esta cota, no del problema.** La parte lineal contrae (radio de
  <<X14_T_G118_LO>> a <<X14_T_G118_HI>>), pero el mayorante de la forma voltaje-corriente, que suma normas por bloque
  sobre todas las barras, tiene radio de <<X14_T_M118_LO>> para arriba. Una cota que conserve las cancelaciones de
  la parte lineal podría cerrar ahí; no se implementó.
- **En la de 300 barras ni la parte lineal contrae** (radio desde <<X14_T_G300_LO>> con δ = 2 %). Con esta familia
  de cotas no hay encierro verificado a esa escala.

**Las dos formas en las redes PQ de la sección 25:**

{{x14_forms}}

*Script: `experiments/x14_pv.py` · CSV: `results/x14_forms.csv`.*

- La forma de corriente verifica en <<X14_F_CUR>> de <<X14_F_N>> casos y la de voltaje-corriente en
  <<X14_F_VI>>. Esta falla donde la otra no en: <<X14_F_ONLY_CUR>>. Al revés, en <<X14_F_ONLY_VI>> casos.
- En la red de 141 barras la causa está aislada: el bloque que acota Y·(I − C·J) crece con las admitancias de las
  ramas muy cortas, y sin esos términos de redondeo el radio del mayorante baja a la sexta parte. Haría falta
  evaluar ese residuo con más precisión.
- Donde verifican las dos, el radio del tubo es el mismo (razón de <<X14_F_RATIO_LO>> a <<X14_F_RATIO_HI>>).
- Consecuencia: en redes PQ se usa la forma de corriente; la de voltaje-corriente queda para cuando hay barras PV.

### 26.7 Redes desbalanceadas

**Observación.** En el dominio de fases, con el neutro a potencial de tierra en cada barra y cargas de potencia
constante entre fase y neutro, cada nodo-fase cumple (YV)ₙ + cₙ = conj(Sₙ/Vₙ): la misma forma que una barra
monofásica. El encierro se aplica entonces sin cambios, con la matriz de admitancias de fases (bloques 3 × 3 con
acoplamiento mutuo) y la fuente en 1∠0°, 1∠−120° y 1∠120°. Se verificó en vez de darla por buena
(`cdhim/threephase.py`, `tests/test_threephase.py`).

**Red.** El lado de baja tensión del IEEE European LV Test Feeder, tal como viene en pandapower: 906 barras a
416 V y 55 cargas monofásicas en estrella. Se usa completa (<<X10_FULL_NODES>> nodos-fase PQ) y reducida de forma
exacta a <<X10_BUSES>> barras (<<X10_NODES>> nodos-fase). Los alimentadores IEEE 13, 37 y 123 se descartaron:
traen reguladores, transformadores y cargas en delta o de corriente constante.

{{x10_validation}}

*Script: `experiments/x10_three_phase.py validation` · CSV: `results/x10_model_validation.csv`.*

- **El modelo coincide con `runpp_3ph` de pandapower** a <<X10_VAL_DIFF>> p.u., y a <<X10_VAL_RES>> p.u. si se
  descuenta el residuo de convergencia de pandapower.
- **Es desbalanceada de verdad:** hasta <<X10_SPREAD>> p.u. de diferencia de módulo entre fases, secuencia negativa
  de hasta <<X10_NEG>> % y secuencia cero de hasta <<X10_ZERO>> %.
- Dejar fuera el transformador Dyn cambia los módulos en a lo más <<X10_TRAFO>> p.u.

**Encierro sin redondeo dirigido** (el de la sección 25), <<X10_SRC>> inversores monofásicos repartidos
<<X10_PHASES>> entre las fases, <<X10_K>> símbolos:

{{x10_cases}}

*Script: `experiments/x10_three_phase.py cases` · CSV: `results/x10_three_phase.csv`.*

- **Cierra en <<X10_FP>> de <<X10_N>> casos** y queda entre <<X10_W_LO>> y <<X10_W_HI>> veces el rango alcanzado;
  los discos por nodo, entre <<X10_D_LO>> y <<X10_D_HI>>. Ninguna violación (<<X10_VIOL>>).
- **Con redondeo dirigido y cota conjunta** (`results/x12_verified.csv`): verifica en <<X12_3P_VER>> de
  <<X12_3P_N>> casos, con unicidad en <<X12_3P_UNI>>, y mide a lo más entre <<X12_3P_WC_LO>> y <<X12_3P_WC_HI>>
  veces el rango exacto. Krawczyk, entre <<X12_3P_KR_LO>> y <<X12_3P_KR_HI>>.
- **La monotonía no se cumple nunca:** en los <<X10_MONO>> casos el mínimo de todos los nodos-fase está fuera de
  (0, −C), a una distancia de <<X10_MONO_GAP>> % del rango. Una fuente en una fase mueve
  las otras dos en sentido contrario, así que un mínimo tomado de un solo flujo no sirve en trifásico.
- **Tamaño completo** (`results/x10_full_model.csv`): con <<X10_FULL_NODES>> nodos-fase cierra en <<X10_FULL_FP>>
  de <<X10_FULL_N>> casos en <<X10_FULL_S_LO>> a <<X10_FULL_S_HI>> s, con las mismas cotas que la red reducida
  (diferencia máxima <<X10_FULL_DIFF>> p.u.) y <<X10_FULL_VIOL>> violaciones en las muestras.

**Hasta dónde llega** (capacidad total multiplicada, flujo inverso, ρU = 80 %):

{{x10_limits}}

*Script: `experiments/x10_three_phase.py limits` · CSV: `results/x10_limits.csv`.*

- El encierro y los discos dejan de cerrar en el mismo escalón (capacidad <<X10_LIM_FIRST_FAIL>> veces la carga),
  aunque el flujo AC sigue teniendo solución en todas las muestras. Con <<X10_LIM_LAST>> veces, los voltajes ya
  están fuera de cualquier rango de operación.

**Qué queda fuera:** cargas en delta, transformadores, la impedancia de cortocircuito de la fuente y el neutro
explícito. Las líneas no transpuestas y las fases faltantes están probadas solo en una red sintética, sin
contraste con un solver externo.

### 26.8 Lectura crítica: qué cambia

1. **El aporte resistió la comparación que podía tumbarlo, y quedó mejor delimitado.** Frente a Krawczyk por
   inyección de corriente, que es el mejor baseline intervalar que se pudo armar y más fuerte que lo publicado,
   el encierro es más ajustado en todas las barras: su exceso sobre el rango exacto es entre <<X12_KR_EXC_LO>> y
   <<X12_KR_EXC_HI>> veces menor. Lo que ya no se puede decir es que la garantía publicada queda a 1,24–2,17 veces
   el rango: esa es la de los discos, y Krawczyk queda a <<X12_KR_WC_LO>>–<<X12_KR_WC_HI>>.
2. **Qué es el método, dicho sin adorno:** Krawczyk un orden más arriba. Modelo de segundo orden, resto de
   tercer orden acotado por un punto fijo en forma de inyección de corriente, y una cota conjunta de los
   términos de primer y segundo orden. Un revisor puede leerlo como incremental; el argumento es que ese orden
   de más cambia el exceso en un factor de <<X12_KR_EXC_LO>> a <<X12_KR_EXC_HI>> y la subdivisión necesaria de
   decenas de cajas a una.
3. **Ahora es un encierro verificado de verdad:** existencia y unicidad con el redondeo controlado, y la
   sobreestimación demostrada contra puntos certificados, no medida contra una búsqueda. En la revisión de
   literatura no apareció nada así a tamaño realista; es un resultado negativo de búsqueda, no una prueba.
4. **Alcance comprobado:** redes PQ radiales, malladas y con elementos en derivación; trifásicas desbalanceadas
   con cargas en estrella; fuentes de voltaje controlado en distribución, con menos ajuste; hasta
   <<X13_BUS_MAX>> barras con <<X13_K_AT_BUS_MAX>> símbolos.
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
7. **RAA-K deja de aportar en este encierro.** Con la cota conjunta, en los <<X13_OK>> casos de la escala el ancho
   es el mismo con RAA-K completa que acotando la forma término a término (razón <<X13_FULL_GAIN>>), a un costo de
   hasta <<X13_FULL_MAX>> s. Quitarlo abarata el encierro (en Oberrhein, de <<X13_OB_FULL>> a <<X13_OB_TERM>> s) y
   acerca su costo al de Krawczyk. Para un artículo cuyo interés es la aritmética afín del grupo, eso hay que
   decirlo tal cual: hoy el resultado no se sostiene en el operador. Lo que sí queda de esa línea es la estructura: la cota
   conjunta es una descomposición en bloques (uᵀε)², como la del teorema k-nomial, acotados junto con el término
   lineal en vez de uno por uno. Usar ahí los bloques de RAA-K en lugar de autovectores no se probó (punto 9.23).
"""
