# CICV-Beauchef
<img src="[TU_URL_DE_IMAGEN_AQUI](https://pbs.twimg.com/profile_images/1509939244382490624/HpHnMwQC.jpg)" alt="CICV" width="500">
![CICV]( | width=50) ![FCFM](https://ingenieria.uchile.cl/.resources/portal-ingenieria/images/fcfm.png | width=50)
Repositorio de proyectos de colaboración CICV-Beauchef en 2026. Los proyectos están separados por carpeta.
Los proyectos de 1 a 3 son en base a modelos metabólicos de balance de flujo y cómo adaptar sus limitaciones:

## Análisis de Balance de Flujos (Flux Balance Analysis - FBA)

### 1.- Deducción
1.  **Balance de masa:** Se define por la ley:\\
    $$\frac{dx}{dt} = S \cdot v$$\\
    Donde $x$ es la concentración de metabolitos, $v$ es el vector de flujos metabólicos y S una matriz $m\times n$ llamada "Matriz estequiométrica". $m$ es el número de metabolitos y $n$ es el número de reacciones.
2.  **Estado estacionario:** Se asume que las concentraciones internas no cambian en el tiempo ($\frac{dx}{dt} = 0$), lo que simplifica el sistema a:
    $$S \cdot v = 0$$
3.  **Optimización:** Como suele haber más reacciones que metabolitos, el sistema está subdeterminado. Se utiliza **Programación Lineal (LP)** para encontrar una solución que maximice o minimice una función objetivo ($Z$), generalmente la biomasa:
    $$\max Z = c^T \cdot v$$

---
## 2.- Supuestos Principales

Para que el modelo sea funcional, el FBA se apoya en tres pilares fundamentales:

* **Estado Estacionario (Steady State):** Se asume que el sistema está en equilibrio dinámico; es decir, la velocidad a la que se produce un metabolito es igual a la velocidad a la que se consume.
* **Optimality (Optimización):** Se parte de la premisa de que los organismos biológicos han evolucionado para maximizar su eficiencia (como el crecimiento celular o la producción de ATP).
* **Limitaciones de Capacidad:** Aunque no conocemos las cinéticas, el modelo impone restricciones de capacidad a los flujos ($v_{min} \leq v \leq v_{max}$), basándose en la disponibilidad de nutrientes o la termodinámica.

---

## 3.- Ventajas y Limitaciones

| Ventajas | Limitaciones |
| :--- | :--- |
| No requiere parámetros cinéticos ($K_m, V_{max}$). | No predice concentraciones de metabolitos. |
| Computacionalmente eficiente para redes grandes. | No considera la regulación genética o alostérica. |
| Ideal para ingeniería metabólica. | El supuesto de optimización no siempre es real. |
