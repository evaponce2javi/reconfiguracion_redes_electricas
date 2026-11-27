# Paper de formas afines: esqueleto

Revista objetivo: **Electric Power Systems Research (EPSR)**, Elsevier, por la vía de **suscripción (sin APC)**. EPSR es
híbrida: el acceso abierto es opcional y cobra (unos US$ 2.500 a 3.000); por suscripción no hay cargo. IJEPES se
descartó porque desde 2024 es de acceso abierto completo con APC de unos US$ 3.500 a 4.190. Alternativa sin APC:
Journal of Modern Power Systems and Clean Energy (acceso abierto sin cargo según DOAJ; verificar en su sitio).

| Archivo | Contenido |
|---|---|
| `main.tex` | Esqueleto con la clase `elsarticle`: título, resumen y highlights en borrador; la proposición de cobertura con su prueba; el plan de cada sección en rojo (`\plan`, `\nota`) y la fuente de cada cifra en azul (`\src`). |
| `refs.bib` | Bibliografía del kit de CD-HIM más las referencias nuevas (DOI verificados en OpenAlex). Las marcadas `TODO: verify` hay que revisarlas. |
| `cover_letter.md` | Borrador de la carta al editor (en inglés). |
| `tables/` | Tablas 1 a 3 en LaTeX, generadas por `experiments/tables_paper.py` desde los CSV. |
| `figures/` | Fig. 1 a 3 en PDF, generadas por `experiments/figures_paper.py` desde los CSV de `results/` (la Fig. 1 primero con `compute`). |

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

## Lo que hay que verificar en la guía de EPSR

No se pudo leer en línea: ScienceDirect rechaza la descarga automática. Lo que se sabe por el artículo de referencia
de la revista (Singh et al. 2025, EPSR 245, 111606):

- Publicado por suscripción ("All rights are reserved"): esa es la vía sin APC. **Al enviar, no marcar la opción de
  acceso abierto.**
- 5 palabras clave, resumen de unas 160 palabras (el nuestro tiene unas 220), unas 7.700 palabras en total.
- Lleva CRediT, declaración de conflictos de interés y disponibilidad de datos.
- Referencias numeradas (compatibles con `elsarticle-num`).

Por verificar en la guía: límite del resumen, highlights (número y caracteres), extensión máxima, revisión simple o
doble ciego, declaración de uso de IA generativa, y la política de preprint (si se puede subir a arXiv o al
repositorio de la PUCV).

## Comprobación de cifras

`python experiments/check_paper_numbers.py` recalcula desde los CSV las 45 cifras citadas en el texto y verifica que
`main.tex` diga exactamente eso. Correrlo después de editar el texto o de regenerar resultados.

## Pendientes antes de enviar

1. Revisión del texto completo por el grupo; la introducción es la que más pesa para el editor.
2. Elegir el título (tres opciones comentadas en `main.tex`).
3. CRediT (roles de cada autor), y ajustar los borradores de las declaraciones de uso de IA y de datos.
4. Publicar el repositorio con un DOI y poner el enlace en la declaración de datos y en la carta.
5. Verificar la guía de EPSR (sección anterior) y quitar las marcas en rojo y azul (`\plan`, `\nota`, `\src`).
6. Al enviar, no marcar la opción de acceso abierto.
