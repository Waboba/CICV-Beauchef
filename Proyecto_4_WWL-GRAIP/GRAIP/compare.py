import numpy as np

# Cargar matrices
dist_julia = np.loadtxt("dist_julia.csv", delimiter=',')
dist_python = np.loadtxt("dist_python.csv", delimiter=',')

# Verificar dimensiones
if dist_julia.shape != dist_python.shape:
    print(f"Error: dimensiones diferentes {dist_julia.shape} vs {dist_python.shape}")
    exit(1)

# Comparar con tolerancia
if np.allclose(dist_julia, dist_python, rtol=1e-5, atol=1e-8):
    print("¡Las matrices son iguales! Traducción correcta.")
else:
    # Mostrar estadísticas de diferencia
    diff = np.abs(dist_julia - dist_python)
    print(f"Diferencias encontradas:")
    print(f"  Máxima diferencia absoluta: {np.max(diff)}")
    print(f"  Media de diferencia absoluta: {np.mean(diff)}")
    print(f"  Desviación estándar de la diferencia: {np.std(diff)}")
    
    # Opcional: guardar diferencias para inspección
    np.savetxt("diff.csv", diff, delimiter=',')
    print("Matriz de diferencias guardada en diff.csv")