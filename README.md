# KineticBox 1D

KineticBox 1D es una simulación paralela escrita en C y OpenMP de un gas unidimensional confinado, con interacciones aleatorias en los límites del recipiente y análisis estadístico de las distribuciones de posición y momento lineal.

El programa está diseñado para simular grandes cantidades de partículas y estudiar la evolución de sus distribuciones espaciales y de momento a medida que interactúan con los límites del sistema.

## Modelo

Cada partícula está caracterizada por:

- una posición x;
- un momento lineal p.

Entre interacciones con los límites, las partículas evolucionan mediante movimiento libre:

```tex
x ← x + p · Δt / m
```

Cuando una partícula alcanza uno de los límites del recipiente, su trayectoria se refleja y se introducen perturbaciones aleatorias.

## Incertidumbre en la posición del límite

La posición después de una interacción con la pared incluye una perturbación aleatoria cuya amplitud está controlada por el parámetro:

```tex
sigmaL
```

Este parámetro puede configurarse desde el archivo `datos.in`.

## Variación del momento

Durante una interacción con un límite también se introduce una modificación aleatoria del momento.

La amplitud utilizada actualmente por la simulación es:

```tex
DeltaE = alfa · ((p - pmin) · (pmax - p))²
```

y el nuevo valor del momento se obtiene mediante una perturbación aleatoria basada en DeltaE.

Los parámetros principales definidos actualmente en el código son:

```tex
pmin = 3.0e-26
pmax = 3.0e-23
alfa = 1.5e42
```

## Análisis estadístico

Durante la simulación se generan histogramas correspondientes a las distribuciones de:

- posición x;
- momento p.

El programa calcula distintas medidas basadas en χ² para comparar las distribuciones obtenidas con las distribuciones esperadas y analizar propiedades como uniformidad y simetría.

Los resultados se escriben en archivos `.dat`.

## Requisitos

Para compilar KineticBox 1D se necesita:

- GCC con soporte para OpenMP;
- GNU Make;
- una plataforma compatible con C99.

## Compilación

Clona el repositorio:

```bash
git clone https://github.com/ignabelitzky/kinetic-box-1d.git
```

Ingresa al directorio:

```bash
cd kinetic-box-1d
```

Compila el programa:

```bash
make
```

Esto genera el ejecutable:

```bash
main
```

El proyecto se compila utilizando optimizaciones de GCC y soporte para OpenMP.

## Configuración

La simulación utiliza el archivo:

```bash
datos.in
```

para definir sus principales parámetros.

Por ejemplo:

```bash
Número de partículas (N_PART): 2097152
(Número de bins - 1) / 2 (BINS): 500
Delta t en seg (DT): 2.215491295991E-7
Masa de las partículas en kg (M): 6.646473667973E-27
Nº de hilos para la corrida (N_THREADS): 32
```

También pueden configurarse:

- la cantidad de pasos de evolución;
- el archivo desde el cual retomar la simulación;
- el archivo utilizado para guardar el estado;
- la incertidumbre `sigmaL` asociada a los límites del recipiente.

El valor predeterminado de N_PART es:

```tex
2097152 = 2²¹ partículas
```

## Ejecución

Una vez compilado el programa y configurado `datos.in`, ejecuta:

```bash
./main
```

Durante la ejecución se muestra información como:

- número de pasos realizados;
- energía total del sistema;
- parámetros utilizados por la simulación.

## Archivos generados

La simulación puede generar diferentes tipos de archivos.

### .dat

Contienen histogramas de las distribuciones de posición y momento junto con los valores estadísticos calculados durante la simulación.

Ejemplo:
```bash
X1000000.dat
```

### .dmp

Contienen el estado de la simulación y permiten guardar o retomar una ejecución.

Para la configuración predeterminada, estos archivos contienen información correspondiente a 2²¹ partículas.

### Limpieza

Para eliminar el ejecutable y los archivos objeto generados durante la compilación:

```bash
make clean
```

## Estructura del proyecto

```bash
kinetic-box-1d/
├── include/
│   ├── constants.h
│   └── utils.h
├── src/
│   ├── g1D-sp-075d.c
│   └── utils.c
├── datos.in
├── Makefile
└── README.md
```