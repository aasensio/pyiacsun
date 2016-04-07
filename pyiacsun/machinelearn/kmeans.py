import numpy as np


def distancia(media, matriz):
    """Calcula la distancia euclideana desde la media hasta
    cada vector en la matriz
    """
    res = []
    for i in matriz:
        res.append(np.sum((media - i)**2.))
    return res


def kmeans(datos_fila, nGrupos, umbral=None):
    """Calcula los vectores clase de la muestra para nGrupos

    Args:
        datos_fila (array [vectores x lambdas]): matriz de datos
        nGrupos (int): Numero de grupos para la clasificacion
        umbral (None, optional): Error relativo

    Returns:
        dict: Los diccionarios con los grupos y las medias
    """
    if umbral is None:
        umbral = 1e-6

    media = np.zeros((nGrupos, datos_fila.shape[1]))
    media[0:nGrupos, :] = datos_fila[0:nGrupos, :]
    difrel = 1.
    contador = 0
    # Creamos nGrupos grupos:
    grupos = {}
    for i in range(nGrupos):
        grupos['grupo' + str(i)] = []
    # Creamos nGrupos medias:
    medias = {}
    for i in range(nGrupos):
        medias['media' + str(i)] = []
    # Creo la matriz de distancias para todos los vectores
    distancias = np.zeros((nGrupos, datos_fila.shape[0]))

    # Comienza las iteraciones hasta que no lleguen a la convergencia
    while difrel > umbral:

        # Reseteamos los grupos:
        for i in range(nGrupos):
            grupos['grupo' + str(i)] = []
        # Calculamos sus distancias
        for grupo in range(nGrupos):
            distancias[grupo, :] = distancia(media[grupo, :], datos_fila)
        # Inserto cada vector con distancia minima al vector medio
        for i in range(len(datos_fila)):
            grupos['grupo' + str(np.argmin(distancias[:, i]))].append(datos_fila[i, :])
        # Calculo el nuevo vector medio
        for grupo in range(nGrupos):
            medias['media' + str(grupo)] = np.mean(np.array(grupos['grupo' + str(grupo)]), axis=0)
        # Calculo las diferencias relativas con la iteracion anterior
        difrel = 0.
        for grupo in range(nGrupos):
            difrel += np.max(np.abs(medias['media'+str(grupo)]-media[grupo, :])/media[grupo, :])
        # Inserto la nueva media:
        for grupo in range(nGrupos):
            media[grupo, :] = medias['media' + str(grupo)]
        contador += 1

    return [medias, grupos]
