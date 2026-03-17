# OptKnock-GECKO

![Status](https://img.shields.io/badge/Status-In%20Development-yellow)
![MATLAB](https://img.shields.io/badge/MATLAB-R2021a%2B-blue)
![Gurobi](https://img.shields.io/badge/Solver-Gurobi-green)

Este repositorio contiene una implementación avanzada de **OptKnock Dual** adaptada para modelos con restricciones enzimáticas (**GECKO 3.0**). El software permite identificar intervenciones genéticas (knockouts de enzimas) que obligan a la célula a producir compuestos de interés industrial mientras maximiza su crecimiento.

---

# Tabla de Contenidos
- [Requisitos](#-requisitos)
- [Estructura del Proyecto](#-estructura-del-proyecto)
- [Instalación](#-instalación)
- [Guía de Uso](#-guía-de-uso)
- [Lógica del Algoritmo](#-lógica-del-algoritmo)

---

## Requisitos

Para utilizar este programa, es necesario contar con el siguiente stack tecnológico:

### 1. Software y Solvers
* **MATLAB**: Versión R2021a o superior.
* **Gurobi Optimizer**: Motor de optimización MILP. Es fundamental para resolver el problema bi-nivel.
    * Ejecutar `gurobi_setup` en MATLAB tras instalarlo.

### 2. Toolboxes de Biología de Sistemas
* **COBRA Toolbox**: Manejo de redes metabólicas a escala genómica.
* **GECKO Toolbox v3.0**: Necesario para la integración de datos enzimáticos ($k_{cat}$).

### 3. Modelo Metabólico
* **ecYeastGEM**: Modelo enzimáticamente restringido de *S. cerevisiae*. Debe estar en formato `.yml`.

---

##  Estructura del Proyecto

El repositorio se organiza de la siguiente manera:

* `geckoOptKnockDualFull.m`: Función principal que construye la matriz de optimización MILP utilizando la restricción de **Dualidad Fuerte**.
* `GECKONOCKfull.m`: Script para cargar el modelo, cantidad de knockouts y producto a maximizar. Resuelve el problema llamando al script anterior.
* `ecYeastGEM.yml`: (No incluido por peso) Archivo del modelo que debe posicionarse en la raíz.

---

##  Instalación

1.  **Clonar el repositorio:**
    ```bash
    git clone [https://github.com/tu-usuario/CICV-Beauchef.git](https://github.com/tu-usuario/CICV-Beauchef.git)
    cd CICV-Beauchef/Proyecto_3_GeckoKnock
    ```
2.  **Configurar Rutas en MATLAB:**
    Asegurarse de que las carpetas de COBRA y GECKO estén en el *Path* de MATLAB:
    ```matlab
    addpath(genpath('C:/ruta/a/cobratoolbox'))
    addpath(genpath('C:/ruta/a/GECKO'))
    ```

---

##  Guía de Uso

Para ejecutar un diseño de producción, sigue estos pasos:

1.  Abre `GECKONOCKfull.m`.
2.  Configura el **Target** (ej. `r_2056` para Succinato) y el número de **Knockouts**.
3.  Ejecuta el script. El programa realizará un barrido de biomasa automático.

### Ejemplo de salida en consola:
> `Biomasa Min: 0.01 -> Producción: 0.8523`  
> `Biomasa Min: 0.10 -> Producción: 0.4211`

---

##  Lógica del Algoritmo

El programa resuelve un problema de optimización bi-nivel. Mientras el usuario busca maximizar el producto, la célula busca maximizar su biomasa. 



Utilizamos la **Dualidad Fuerte** para convertir este problema en un solo nivel:

$$ \text{max } z = c^T v + \epsilon \cdot \text{biomass} $$

Sujeto a:
* Balance de masa: $S \cdot v = 0$
* Optimalidad Dual: $S^T w + \lambda - \mu = c_{biom}$
* Restricción de presupuesto: $\sum (1 - y_i) = K$

---


