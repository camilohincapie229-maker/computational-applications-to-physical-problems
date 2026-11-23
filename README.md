# Aplicaciones computacionales a problemas físicos / Computational Applications to Physical Problems

## Español

Soluciones de cinco casos de estudio de la asignatura Aplicaciones Computacionales a Problemas Físicos de la Universidad Industrial de Santander (UIS).

| Caso | Tema | Documentos | Código |
| --- | --- | --- | --- |
| 1 | Intercepción de un objeto en movimiento | [Reporte](Caso_1/reporte.pdf) · [Guía](Caso_1/guia.pdf) | [Sin gravedad](Caso_1/modelo_computacional.py) · [Con gravedad](Caso_1/persecucion.py) |
| 2 | Modos normales de una cuerda no homogénea | [Reporte](Caso_2/reporte.pdf) · [Guía](Caso_2/guia.pdf) | [Script](Caso_2/cuerda_no_homogenea.py) · [Notebook](Caso_2/Caso_2.ipynb) |
| 3 | Propiedades ópticas y electrónicas de la plata | [Reporte](Caso_3/reporte.pdf) · [Guía](Caso_3/guia.pdf) | [Notebook](Caso_3/Caso_3.ipynb) |
| 4 | Oscilaciones forzadas de una cuerda fija-libre | [Reporte](Caso_4/reporte.pdf) · [Guía](Caso_4/guia.pdf) | [Notebook](Caso_4/Caso_4.ipynb) |
| 5 | Temperatura estacionaria en una placa bidimensional | [Reporte](Final/reporte.pdf) · [Guía](Final/guia.pdf) | [Notebook](Final/Caso_5.ipynb) |

Desde la raíz, con Python 3.11 o posterior, instala las dependencias con `python -m pip install -r requirements.txt`. Abre los notebooks con `python -m notebook` y ejecuta todas las celdas en orden. El caso 3 requiere conexión a internet para descargar los datos experimentales.

Ejecuta los scripts con `python Caso_1/modelo_computacional.py`, `python Caso_1/persecucion.py` o `python Caso_2/cuerda_no_homogenea.py -N 6 -L 1 --mu1 1 --mu2 2`. Para usar `--simulate` en el caso 2, instala las dependencias gráficas con `python -m pip install -r Caso_2/requirements.txt`.

## English

Solutions to five case studies from the Computational Applications to Physical Problems course at Universidad Industrial de Santander (UIS).

| Case | Topic | Documents | Code |
| --- | --- | --- | --- |
| 1 | Interception of a moving object | [Report](Caso_1/reporte.pdf) · [Guide](Caso_1/guia.pdf) | [Without gravity](Caso_1/modelo_computacional.py) · [With gravity](Caso_1/persecucion.py) |
| 2 | Normal modes of a nonuniform string | [Report](Caso_2/reporte.pdf) · [Guide](Caso_2/guia.pdf) | [Script](Caso_2/cuerda_no_homogenea.py) · [Notebook](Caso_2/Caso_2.ipynb) |
| 3 | Optical and electronic properties of silver | [Report](Caso_3/reporte.pdf) · [Guide](Caso_3/guia.pdf) | [Notebook](Caso_3/Caso_3.ipynb) |
| 4 | Forced oscillations of a fixed-free string | [Report](Caso_4/reporte.pdf) · [Guide](Caso_4/guia.pdf) | [Notebook](Caso_4/Caso_4.ipynb) |
| 5 | Steady-state temperature in a two-dimensional plate | [Report](Final/reporte.pdf) · [Guide](Final/guia.pdf) | [Notebook](Final/Caso_5.ipynb) |

From the repository root, use Python 3.11 or later and install the dependencies with `python -m pip install -r requirements.txt`. Open the notebooks with `python -m notebook` and run all cells in order. Case 3 requires an internet connection to download the experimental data.

Run the scripts with `python Caso_1/modelo_computacional.py`, `python Caso_1/persecucion.py`, or `python Caso_2/cuerda_no_homogenea.py -N 6 -L 1 --mu1 1 --mu2 2`. To use `--simulate` in case 2, install the graphical dependencies with `python -m pip install -r Caso_2/requirements.txt`.
