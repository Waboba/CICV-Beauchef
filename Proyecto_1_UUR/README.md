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

Finalmente, el algoritmo retorna los óptimos y el **precio dual** de cada esencial respecto a cada no esencial.


### ¿Qué es el Precio Dual?

Imaginemos que tenemos un problema de optimización, por ejemplo, queremos maximizar un beneficio:

$$\max p(x) \quad \text{sujeto a} \quad g(x) \leq 0$$

Si tomamos una restricción $g_i(x) \leq 0$ (un límite de presupuesto o de un recurso) y la **relajamos** permitiendo un pequeño margen $u$ (es decir, $g_i(x) \leq u$), el resultado final de nuestro problema cambiará. 

Si llamamos $p^*(u)$ al nuevo valor óptimo, existe una relación matemática fundamental:

$$\frac{dp^*(0)}{du} = \lambda^*$$

---

### En palabras sencillas (sin tecnicismos)

El **Precio Dual ($\lambda^*$)** mide la **sensibilidad**. Nos indica cuánto mejoraría nuestro resultado final si tuviéramos "un poquito más" de un recurso que nos está limitando.

**En tu caso práctico:**
Se interpreta como el cambio en el valor óptimo del **nutriente no esencial** frente a una pequeña variación en la cantidad disponible de un **nutriente esencial**.

* **Si el Precio Dual es alto:** Significa que ese nutriente esencial es un "cuello de botella". Si consiguiéramos un poco más, nuestro objetivo final mejoraría mucho.
* **Si el Precio Dual es cero:** Significa que ese nutriente nos sobra; aunque tuviéramos más, el resultado final no cambiaría porque no es lo que nos está frenando.


---

## 2. Instalación

1.  **Solver:** Instala tu solver de preferencia (por ejemplo, **Gurobi**).
2.  **Modelo:** Descarga el modelo metabólico (ej. `iCHOv1_DG44.xml`) y guárdalo en una ubicación conocida.
3.  **Script:** Descarga el archivo `Algoritmo.py` e impórtalo en tu entorno de Python.
4.  **Configuración:** En el script, define la ruta de tu modelo anteponiendo una `r` para evitar errores de caracteres:

```python
ruta_modelo = r"C:\Users\tu_usuario\Desktop\iCHOv1_DG44.xml"
