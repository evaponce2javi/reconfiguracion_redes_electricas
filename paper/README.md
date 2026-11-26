# Paper de formas afines: esqueleto

Revista objetivo: **International Journal of Electrical Power & Energy Systems (IJEPES)**, Elsevier.

| Archivo | Contenido |
|---|---|
| `main.tex` | Esqueleto con la clase `elsarticle`: título, resumen y highlights en borrador; la proposición de cobertura con su prueba; el plan de cada sección en rojo (`\plan`, `\nota`) y la fuente de cada cifra en azul (`\src`). |
| `refs.bib` | Bibliografía del kit de CD-HIM más las referencias nuevas (DOI verificados en OpenAlex). Las marcadas `TODO: verify` hay que revisarlas. |
| `figures/` | Vacía. Las figuras se generan con un script por escribir (`experiments/figures_paper.py`), desde los CSV. |

Compilar: `cd paper && latexmk -pdf main.tex` (en el contenedor: `apt-get install texlive-latex-extra texlive-publishers
latexmk`; `elsarticle` viene en `texlive-publishers`), o en Overleaf, que la trae incluida. Los archivos intermedios y
el PDF no se versionan.

## Reglas

- Ninguna cifra a mano: cada número sale de `RESULTS.md` (secciones 27 a 29) o de `results/*.csv`. La lista de
  afirmaciones con su respaldo y su límite es la **sección 29 de RESULTS.md**.
- No se menciona CD-HIM ni la elipse: es un trabajo nuevo. De `caja-elipse.tex` se reutiliza la introducción
  (primer párrafo), la formulación linealizada y las sensibilidades dispersas (§II-A, §II-B, §III-A).
- El modelo de segundo orden no es aporte: se cita a Brucoli et al. (1985) y Sandoval et al. (2022). El aporte es
  usar su **rango** con un radio de concentración, la **garantía** y la **validación de cobertura** a escala.

## Lo que hay que verificar en la guía de IJEPES

No se pudo leer en línea: ScienceDirect rechaza la descarga automática. Por el artículo de Basirat et al. (IJEPES
2026), la revista es **de acceso abierto (CC BY)**. **APC: 0 para este envío** (confirmado por los autores).
- **Extensión:** el artículo de referencia tiene unas 10.000 palabras y 13 páginas.
- **Resumen y highlights:** límite de palabras del resumen (el de referencia tiene unas 190), número de highlights y caracteres por viñeta (habitual en Elsevier: 3 a 5, hasta 85 caracteres).
- **Revisión:** simple o doble ciego. Si es doble ciego, va una página de título aparte.
- **Declaraciones:** CRediT, conflictos de interés, disponibilidad de datos y uso de IA generativa.

## Pendientes antes de redactar

1. Elegir el título (tres opciones en `main.tex`).
2. Script de figuras: histograma del término cuadrático frente a la normal; mapa carga × δ; ancho de las corrientes.
   Antes hay que generar el CSV del diagnóstico de asimetría y curtosis, que hoy no tiene uno.
3. Decidir la declaración de uso de IA y el repositorio de datos.
4. Revisar las referencias marcadas `TODO`.
