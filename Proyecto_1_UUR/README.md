# Algoritmo de Tasas de Consumo No Convencionales (UOF)

Este proyecto implementa un enfoque alterno para estudiar complejas cuyo objetivo no es maximizar su biomasa, utilizando el modelo de **Tasas de Consumo No Convencionales (UOF)**. El método identifica conexiones relevantes entre metabolitos aprovechando la sensibilidad de los precios duales.

---

## 1. Problema a resolver

Cuando trabajamos con células complejas cuyo objetivo no es maximizar su biomasa las herramientas como el **Análisis de Balance de Flujos (FBA)** pueden no ser exactas. Este algoritmo propone una solución en dos etapas:

### Etapa 1: Essential Nutrient Minimization (ENM)
El algoritmo recibe un modelo metabólico y plantea el siguiente **problema lineal** para cada nutriente:

$$
\begin{aligned}
    & \text{Minimizar} & &| v_{\text{nutriente}}| \\
    & \text{sujeto a}  & & S \cdot v = 0, \\
    &                  & & v_{\text{biomasa}}=\mu_{\text{obs}}\\
    &                  & & v_{\text{proteina}}=p_{\text{obs}}\\
    &                  & & l_j \le v_j \le u_j, \quad \forall j \in\text{Reacciones}\setminus \text{Nutrientes}\\
    &                  & & v_i \in \mathbb{R} \quad \forall i \text{ Nutriente}
\end{aligned}
$$

Si el óptimo es distinto de cero, el nutriente se clasifica como **esencial** y se registra su consumo mínimo.

### Etapa 2: Uptakes Objective Functions (UOF)
Utiliza la lista de esenciales de la etapa anterior para resolver el siguiente problema para cada nutriente **no esencial**:

$$
\begin{aligned}
    & \text{Minimizar} & &| v_{\text{nutriente}}| \\
    & \text{sujeto a}  & & S \cdot v = 0, \\
    &                  & & v_{\text{biomasa}}=\mu_{\text{obs}}\\
    &                  & & v_{\text{proteina}}=p_{\text{obs}}\\
    &                  & & l_j \le v_j \le u_j, \quad \forall j \in\text{Reacciones}\setminus (\text{Esenciales}\cup \{v_{\text{nutrientes}}\})\\
    &                  & & v_i = v_i^{\text{ENM}} \quad \forall i \text{ Esenciales}\\
    &                  & & v_{\text{nutriente}}\in \mathbb{R}
\end{aligned}
$$

Finalmente, el algoritmo retorna los óptimos y el **precio dual** $\lambda$ de cada esencial respecto a cada no esencial.


### ¿Qué es el Precio Dual?

Imaginemos que tenemos un problema de optimización, por ejemplo, queremos maximizar un beneficio:

$$\max p(x) \quad \text{sujeto a} \quad g(x) \leq 0$$

Si tomamos una restricción $g_i(x) \leq 0$ (un límite de presupuesto o de un recurso) y la **relajamos** permitiendo un pequeño margen $u$ (es decir, $g_i(x) \leq u$), el resultado final de nuestro problema cambiará. 

Si llamamos $p(u)$ al nuevo valor óptimo, existe una relación matemática fundamental:

$$\frac{dp(0)}{du} = \lambda$$



### Como se interpreta

El **Precio Dual ($\lambda$)** mide la **sensibilidad**. Nos indica cuánto cambiaría nuestro resultado final si tuviéramos "un poquito más" de un recurso que nos está limitando.


**En nuestro caso:**
Se interpreta como el cambio en el valor óptimo del **nutriente no esencial** frente a una pequeña variación en la cantidad disponible de un **nutriente esencial**.





## 2. Instalación

A continuación se muetra como instalar y utilizar el programa que aplica lo explicado anteriormente. Por motivos ilustrativos haremos la muestra con el modelo  `iCHOv1_DG44.xml` de las células CHO.

1.  **Solver:** Instala tu solver de preferencia (por ejemplo, **Gurobi**).
2.  **Modelo:** Descarga el modelo metabólico (ej. `iCHOv1_DG44.xml`) y guárdalo en una ubicación conocida.
3.  **Script:** Descarga el archivo `Algoritmo.py` disponible en el GitHub e impórtalo en tu entorno de Python.
4.  **Configuración:** En el script, debemos decirle al programa cuál es el modelo con el que estamos trabajando, para esto debemos decirle donde está guardado el mismo en nuetra computadora. Ve al archivo, click derecho y "Copiar como ruta" ("Copy as path") y deberias tener algo como "C:\Users\tu_usuario\Desktop\iCHOv1_DG44.xml".  Define la ruta de tu modelo anteponiendo una `r`, debería quedar algo como:


`ruta_modelo = r"C:\Users\tu_usuario\Desktop\iCHOv1_DG44.xml"`

Con esto el codigo ya esta funcionando con el modelo que elegiste.


## 3. Cómo usarlo

Recordemos que el algoritmo consta de dos etapas. Cada etapa es una función, ENM o UOF, además hay una tercera función llamado  `mapa_calor`. Nosotros solamente llamaremos a esta última.


Para la ejecución de la función, es necesario suministrar cuatro argumentos principales:

1. **`modelo_base`:** Corresponde al modelo metabólico a escala genómica cargado previamente (por ejemplo, el modelo `iCHOv1\_DG44`).


3. **`proteinas`**: Esta sección define las reacciones de demanda asociadas a la síntesis de proteínas específicas y su flujo metabólico. Para implementarlo, utilizamos un diccionario, una estructura de datos que vincula cada identificador de reacción con su respectiva cota de consumo (tasa de producción). De esta manera, el diccionario permite fijar el flujo deseado para cada reacción de la siguiente forma:


  `tasas_sintesis = {
    "ID_REACCION_1": tasa_fijada_1,
    "ID_REACCION_2": tasa_fijada_2,
    # ...
    "ID_REACCION_N": tasa_fijada_n
  }`


  
  En el ejemplo del modelo `iCHOv1\_DG44` nos queda el siguiente diccionario:

 `proteinas ={ 'igg_hc_1':2.6060000000000005e-05 , 'igg_lc_1':2.606e-05, 'igg_formation':1.303e-05, 'DM_igg_g':1.303e-05}`


  4.**Lista de Nuetrientes (`lista_nutrientes` ó `medio`)**:  Define las condiciones del medio extracelular. Esta entrada es flexible y permite entregar tanto una lista de identificadores de reacciones como un diccionario con límites específicos de consumo. Es fundamental especificar este parámetro debido a que los modelos suelen ser "generalistas" y permiten el consumo de una vasta gama de nutrientes. Mediante esta entrada, se acota la disponibilidad de componentes para simular un escenario biológico realista; la distinción técnica entre el uso de una lista o un diccionario se profundizará más adelante.



  
1. **`modelo_base`**: Corresponde al modelo metabólico a escala genómica cargado previamente (por ejemplo, el modelo `iCHOv1_DG44`).

3. **`proteinas`**: Esta sección define las reacciones de demanda asociadas a la síntesis de proteínas específicas y su flujo metabólico. Para implementarlo, utilizamos un diccionario, una estructura de datos que vincula cada identificador de reacción con su respectiva cota de consumo (tasa de producción). De esta manera, el diccionario permite fijar el flujo deseado para cada reacción de la siguiente forma:

   ```python
   tasas_sintesis = {
       "ID_REACCION_1": tasa_fijada_1,
       "ID_REACCION_2": tasa_fijada_2,
       # ...
       "ID_REACCION_N": tasa_fijada_n
   }

```

En el ejemplo del modelo `iCHOv1_DG44` nos queda el siguiente diccionario:

```python
proteinas = {
    'igg_hc_1': 2.6060000000000005e-05, 
    'igg_lc_1': 2.606e-05, 
    'igg_formation': 1.303e-05, 
    'DM_igg_g': 1.303e-05
}

```

4. **Lista de Nutrientes (`lista_nutrientes` ó `medio`)**: Define las condiciones del medio extracelular. Esta entrada es flexible y permite entregar tanto una lista de identificadores de reacciones como un diccionario con límites específicos de consumo. Es fundamental especificar este parámetro debido a que los modelos suelen ser "generalistas" y permiten el consumo de una vasta gama de nutrientes. Mediante esta entrada, se acota la disponibilidad de componentes para simular un escenario biológico realista; la distinción técnica entre el uso de una lista o un diccionario se profundizará más adelante.

```



