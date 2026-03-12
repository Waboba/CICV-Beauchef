# CICV-Beauchef
Repositorio de proyectos de colaboración CICV-Beauchef en 2026. Los proyectos están separados por carpeta.
Los proyectos de 1 a 3 son en base a modelos metabólicos de balance de flujo:

## Análisis de Balance de Flujos (Flux Balance Analysis - FBA)
1.  **Balance de masa:** Se define por la ecuación diferencial:
    $$\frac{dx}{dt} = S \cdot v$$
    Donde $x$ es la concentración de metabolitos y $v$ es el vector de flujos metabólicos.
2.  **Estado estacionario:** Se asume que las concentraciones internas no cambian en el tiempo ($\frac{dx}{dt} = 0$), lo que simplifica el sistema a:
    $$S \cdot v = 0$$
3.  **Optimización:** Como suele haber más reacciones que metabolitos, el sistema está subdeterminado. Se utiliza **Programación Lineal (LP)** para encontrar una solución que maximice o minimice una función objetivo ($Z$), generalmente la biomasa:
    $$\max Z = c^T \cdot v$$

---
