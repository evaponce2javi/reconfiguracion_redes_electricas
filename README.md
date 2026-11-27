# CD-HIM

Correlation-Driven Hybrid Interval Method: cotas de voltaje bajo incertidumbre de inyección, con un modelo
AC linealizado (sensibilidades dispersas) y selección Elipse/Caja por rama.

- [METODO.md](METODO.md): el encierro verificado de voltajes, explicado: objetivos, proceso, diagrama y resultados.
- [CAMBIOS.md](CAMBIOS.md): todo lo que se cambió desde el código original, con la razón de cada cambio, y lo que falta.
- [AUDIT.md](AUDIT.md): auditoría del código original (`finalFinalCode.py`) frente al manuscrito rechazado.
- [RESULTS.md](RESULTS.md): resultados de todos los experimentos, los pedidos por los revisores y los de la versión de aritmética afín.
- [CHECKLIST_AA.md](CHECKLIST_AA.md): pasos para la versión orientada a aritmética afín, con el estado de cada uno.
- [RELATED_WORK.md](RELATED_WORK.md): revisión bibliográfica, con el nivel de lectura de cada fuente.

## Instalación

```bash
python -m venv .venv
.venv\Scripts\activate            # Windows;  source .venv/bin/activate en Linux/macOS
pip install -r requirements.txt
```

Probado con Python 3.12.10, pandapower 3.5.6, NumPy 2.4.6 y SciPy 1.18.1.

## Uso

```bash
python -m cdhim ieee33                       # resumen del método en el punto de diseño
python -m cdhim pegase2869 --radius coordinate --csv cotas.csv
python -m cdhim ieee33 --legacy              # reproduce las cifras del manuscrito enviado
python -m pytest                             # 326 pruebas unitarias
python experiments/run_all.py                # regenera todos los CSV de results/ (2 a 2,5 h; C1 es lo más lento)
python experiments/run_all.py b1 c5          # solo algunos experimentos
python experiments/run_all.py x1             # exploratorio: resto de segundo orden y aritmética afín (unos 20 min)
python experiments/run_all.py x2 x3          # Puerta 2 del checklist: pocas fuentes por zona (unos 40 min); x3 usa los CSV de x2
python experiments/run_all.py x4 x5 x6 x7    # prueba crítica de los cinco aportes posibles (usan los CSV de x2; x7 necesita cvxpy)
python experiments/run_all.py x8             # opción A: encierros en redes de distribución frente a lo publicado (descarga 4 casos de MATPOWER 7.1)
python experiments/run_all.py x9 x10         # baseline de Krawczyk y red trifásica desbalanceada (secciones 26.1 y 26.7)
python experiments/run_all.py x12            # encierro verificado, referencia certificada y subdivisión (unas 2 h; usa los módulos de x9 y x10)
python experiments/run_all.py x13 x14        # escala (más de 1 h, a solas para los tiempos; descarga 3 casos de MATPOWER) y barras PV
python experiments/run_all.py b4 b5          # corrección de segundo orden (B4: punto de diseño, barrido y tiempos, unas 2,5 h) y formas afines (B5: punto de diseño)
python experiments/b5_affine_forms.py sweep  # B5: barrido carga x delta (unas 2 h); también: surrogate, radii
python experiments/b6_correlated.py          # B6: fuentes correlacionadas como símbolos compartidos (unos 30 min)
python experiments/b7_generality.py          # B7: cobertura simultánea e IEEE-118 / PEGASE-1354 (unas 2 h)
python experiments/b8_ppf_baselines.py       # B8: cumulantes y estimación por puntos (unos 40 min)
python experiments/b9_branch_currents.py     # B9: corrientes de rama (unos 40 min)
python experiments/b10_marginals.py          # B10: distribuciones no uniformes (unos 30 min)
python experiments/figures_paper.py compute  # datos de la Fig. 1 del paper (forma del término cuadrático)
python experiments/figures_paper.py          # figuras del paper en paper/figures/
python experiments/build_report.py           # reconstruye RESULTS.md a partir de los CSV
python experiments/figures.py                # figuras en results/figures/
```

Todos los parámetros (δ, γ_th, δ_stab, nivel de confianza, tamaños de muestra, semilla) están en
[config.yaml](config.yaml). Cada CSV de `results/` lleva al lado un `.meta.json` con la semilla, el tamaño
de muestra, el script que lo generó y las versiones de las librerías.

## Estructura

| Ruta | Contenido |
|---|---|
| `cdhim/geometry.py` | Cajas (covarianza, 3σ, afín), elipse, G, κ, radios, Mahalanobis. Funciones puras. |
| `cdhim/core.py` | Algoritmo 1: selección por rama e intervalos de salida. Definiciones de CRI. |
| `cdhim/powerflow.py` | Jacobiano disperso, Newton-Raphson y solver AC por lotes. |
| `cdhim/study.py` | Punto de operación → LU disperso (SuperLU) → coeficientes afines. |
| `cdhim/uncertainty.py` | Distribuciones de los símbolos de ruido, Σ_ε, cópula gaussiana. |
| `cdhim/validation.py` | Monte Carlo: violación por barra, por rama y a nivel sistema, contra el modelo lineal y contra AC. |
| `cdhim/network.py` | Carga de casos (pandapower, PGLib y alimentadores de MATPOWER 7.1, a los que aplica la conversión de unidades escrita en sus archivos). |
| `cdhim/zones.py` | Zonas para fuentes agregadas (bisección espectral de la topología). |
| `cdhim/second_order.py` | Modelo de segundo orden de \|V\| (forma cuadrática exacta por barra) y cotas de referencia. |
| `cdhim/enclosure.py` | Encierro por punto fijo (sin redondeo dirigido): cota del resto en rectangulares (sección 22) y en forma de inyección de corriente (sección 25), y `voltage_enclosure`. |
| `cdhim/halfdisk.py` | Incertidumbre en semidiscos (generación con inversor): estudio sobre las cajas que los cubren, primer orden exacto sobre los semidiscos, búsqueda local de extremos. |
| `cdhim/distflow.py` | Línea base para redes radiales: reimplementación de Sou, Khong y Sandberg (2026), con su regla de montaje. |
| `cdhim/discs.py` | Línea base para cualquier topología con barras PQ: discos de Wang, Bernstein, Le Boudec y Paolone (2018). |
| `cdhim/verified.py` | Encierro **verificado** en punto flotante para redes PQ: tubo alrededor del modelo de segundo orden (existencia y unicidad), cota conjunta de primer y segundo orden y cotas de \|V\| con redondeo dirigido (sección 26.3). |
| `cdhim/verified_vi.py` | La misma verificación en forma voltaje-corriente, que admite barras PV (sección 26.6). |
| `cdhim/certified.py` | Referencias certificadas: puntos interiores verificados y rango exacto por bisección de la caja (sección 26.4). |
| `cdhim/krawczyk.py` | Línea base intervalar: operador de Krawczyk paramétrico, por desbalance de potencia o por inyección de corriente (sección 26.1). |
| `cdhim/threephase.py` | Redes trifásicas desbalanceadas en el dominio de fases, con cargas en estrella (sección 26.7). |
| `cdhim/correction.py` | Corrección de segundo orden de las cotas CD-HIM: media exacta, varianza y rango muestreados del término cuadrático, q en los vértices, elipse exterior a E + rango (sección 27). |
| `cdhim/affine_forms.py` | Salidas desde las formas afines sin elipse: barras y caídas por rama, radios afín y de Bernstein, modelo de segundo orden por objetivo sin formar H, búsqueda del peor vértice (sección 28). |
| `cdhim/knomial.py` | RAA-K: reescritura k-nomial por presupuesto diagonal y su cota. Port de `ibex-lib/raak/affine/ibex_KNomial.h`, validado contra `tests/data/knomial_golden.txt`. |
| `experiments/` | Un script por experimento (B1 a C8; B4 y B5 son posteriores a la revisión: secciones 27 a 29), `build_report.py` (arma RESULTS.md) y `figures.py`. Los scripts `x1_` a `x14_` son de la versión de aritmética afín (secciones 17 a 26 de RESULTS.md): no los pidieron los revisores. `report_x12.py` arma las tablas de la sección 26. |
| `tests/` | Pruebas unitarias. |
| `finalFinalCode.py` | Código original, sin modificar. Se conserva como referencia de la auditoría. |

## Modo `legacy`

El código original evalúa el Jacobiano en arranque plano y con un error de signo en ∂P/∂θ (ver AUDIT.md).
El paquete nuevo usa por defecto el Jacobiano exacto en el punto AC convergido. Para reproducir las cifras del
manuscrito enviado: `linearization.point: flat` y `linearization.jacobian: legacy_signbug` en `config.yaml`,
o la opción `--legacy`. La prueba `tests/test_powerflow.py::test_legacy_mode_reproduces_the_submitted_figures`
verifica que ese modo da 16/31 ramas, CRI 86,70 % y la rama 5–6 de la Fig. 4.

## Datos

`pglib_opf_case2869_pegase.m` se descarga de [PGLib-OPF](https://github.com/power-grid-lib/pglib-opf) a `data/`
en el primer uso (misma fuente que el script original). Los demás casos vienen con pandapower.
