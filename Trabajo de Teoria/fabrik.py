import numpy as np

def fabrik(joints, target, lengths, tolerance=0.01, max_iter=10):
    """
    Implementación del algoritmo FABRIK.
    joints: Lista de coordenadas [x, y] de las articulaciones.
    target: Coordenada [x, y] del objetivo.
    lengths: Distancias rígidas entre articulaciones.
    """
    joints = np.array(joints, dtype=float)
    target = np.array(target, dtype=float)
    base_origin = joints[0].copy()
    
    # 1. Comprobar si el objetivo está fuera de alcance
    dist_to_target = np.linalg.norm(joints[0] - target)
    if dist_to_target > sum(lengths):
        print("Objetivo fuera de alcance: Alineando articulaciones.")
        for i in range(len(joints) - 1):
            r = np.linalg.norm(target - joints[i])
            l = lengths[i] / r
            joints[i+1] = (1 - l) * joints[i] + l * target
        return joints

    # 2. Proceso Iterativo (Dentro del alcance)
    for _ in range(max_iter):
        # Condición de parada: ¿Estamos lo suficientemente cerca?
        if np.linalg.norm(joints[-1] - target) < tolerance:
            break

        # FASE BACKWARD: Del extremo a la base
        joints[-1] = target
        for i in range(len(joints) - 2, -1, -1):
            r = np.linalg.norm(joints[i+1] - joints[i])
            l = lengths[i] / r
            joints[i] = (1 - l) * joints[i+1] + l * joints[i]

        # FASE FORWARD: De la base al extremo
        joints[0] = base_origin
        for i in range(len(joints) - 1):
            r = np.linalg.norm(joints[i+1] - joints[i])
            l = lengths[i] / r
            joints[i+1] = (1 - l) * joints[i] + l * joints[i+1]

    return joints

# --- Ejemplo de uso ---
puntos_iniciales = [[0, 0], [0, 2], [0, 4], [0, 6]] # 3 segmentos de longitud 2
longitudes = [2, 2, 2]
objetivo = [4, 3]

resultado = fabrik(puntos_iniciales, objetivo, longitudes)
print("Nuevas posiciones de las articulaciones:\n", resultado)