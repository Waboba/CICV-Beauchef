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

---

## 2. Instalación

1.  **Solver:** Instala tu solver de preferencia (por ejemplo, **Gurobi**).
2.  **Modelo:** Descarga el modelo metabólico (ej. `iCHOv1_DG44.xml`) y guárdalo en una ubicación conocida.
3.  **Script:** Descarga el archivo `Algoritmo.py` e impórtalo en tu entorno de Python.
4.  **Configuración:** En el script, define la ruta de tu modelo anteponiendo una `r` para evitar errores de caracteres:

```python
ruta_modelo = r"C:\Users\tu_usuario\Desktop\iCHOv1_DG44.xml"
