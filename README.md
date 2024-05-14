
## Acerca de

Directorio donde se ha llevado a cabo la parte métodos y obtención de resultados del TFM.

Consiste en 5 scripts que descargan alphamissense (0), procesan las variantes de acuerdo al criterio ACMG (1), transforman las variantes a formato tabular aceptado por VEP (2) para llevar a cabo la anotación con los predictores de interés (3), seguido de un procesado de los resultados obtenidos con vep (4). Por ultimo se anota (merge) los resultados descargados de AMS para las variantes de interés (5).

El archivo resultante se utilizará en `analysis/sarcomeric_validation.ipynb` para hacer un último procesado, analisis exploratorio y obtener resultados de curvas PR y ROC.
Posteriormente se utilizará R en `R_figures/figures.Rmd` para obtener las gráficas finales.

## Contenidos

- _in/ --> archivos de entrada para los scripts, principalmente los raw obtenidos de la web de AlphaMissense (/raw/...) y las variantes en los genes sarcoméricos (`sarcomeric_filtered_hg38.csv`).
- analysis/ --> python notebooks donde se hace un análisis exploratorio y se obtienen las curvas correspondientes. También contiene la carpeta R_figures donde se obtienen los gráficos que se utilizarán en la memoria del TFM
- VEP_format/ --> archivos para dar forma y resultados de la anotación de vep



