#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Robótica Computacional - 
# Grado en Ingeniería Informática (Cuarto)
# Práctica: Resolución de la cinemática inversa mediante CCD
#           (Cyclic Coordinate Descent).

import sys
from math import *
import numpy as np
import matplotlib.pyplot as plt
import colorsys as cs

# ******************************************************************************
# Declaración de funciones

def muestra_origenes(O,final=0):
  # Muestra los orígenes de coordenadas para cada articulación
  print('Origenes de coordenadas:')
  for i in range(len(O)):
    print('(O'+str(i)+')0\t= '+str([round(j,3) for j in O[i]]))
  if final:
    print('E.Final = '+str([round(j,3) for j in final]))

def muestra_robot(O,obj):
  # Muestra el robot graficamente
  plt.figure()
  plt.xlim(-L,L)
  plt.ylim(-L,L)
  T = [np.array(o).T.tolist() for o in O]
  for i in range(len(T)):
    plt.plot(T[i][0], T[i][1], '-o', color=cs.hsv_to_rgb(i/float(len(T)),1,1))
  plt.plot(obj[0], obj[1], '*')
  plt.pause(0.0001)
  plt.show()
  plt.close()

def matriz_T(d,th,a,al):
  # Calcula la matriz de transformación usando parámetros DH
  return [[cos(th), -sin(th)*cos(al),  sin(th)*sin(al), a*cos(th)]
         ,[sin(th),  cos(th)*cos(al), -sin(al)*cos(th), a*sin(th)]
         ,[      0,          sin(al),          cos(al),         d]
         ,[      0,                0,                0,         1]
         ]

def cin_dir(th,a):
  # Sea 'th' el vector de thetas y 'a' el vector de longitudes
  T = np.identity(4)
  o = [[0,0]]
  for i in range(len(th)):
    T = np.dot(T,matriz_T(0,th[i],a[i],0))
    tmp=np.dot(T,[0,0,0,1])
    o.append([tmp[0],tmp[1]])
  return o

# ******************************************************************************
# Cálculo de la cinemática inversa de forma iterativa por el método CCD

# Valores articulares iniciales y longitudes de los eslabones
th=[0.,0.,0.]
a =[5.,5.,5.]
L = sum(a) # Variable para representación gráfica
EPSILON = .01

# Introducción del punto objetivo
if len(sys.argv) != 3:
  sys.exit("python " + sys.argv[0] + " x y")
objetivo=[float(i) for i in sys.argv[1:]]

O=cin_dir(th,a)
print ("- Posicion inicial:")
muestra_origenes(O)

dist = float("inf")
prev = 0.
iteracion = 1

# Bucle principal de CCD
while (dist > EPSILON and abs(prev-dist) > EPSILON/100.):
  prev = dist
  O_historial = [cin_dir(th,a)]
  
  # IMPLEMENTACIÓN CCD: De la articulación más distal a la más proximal
  for i in range(len(th)-1, -1, -1):
    # 1. Calcular posición actual con cinemática directa
    O_actual = cin_dir(th, a)
    pos_articulacion = np.array(O_actual[i])
    pos_extremo = np.array(O_actual[-1])
    pos_objetivo = np.array(objetivo)

    # 2. Vectores desde la articulación al extremo y al objetivo
    v_extremo = pos_extremo - pos_articulacion
    v_objetivo = pos_objetivo - pos_articulacion

    # 3. Cálculo del ángulo de corrección mediante atan2
    alfa_extremo = atan2(v_extremo[1], v_extremo[0])
    alfa_objetivo = atan2(v_objetivo[1], v_objetivo[0])
    
    # 4. Actualización del ángulo
    th[i] = th[i] + (alfa_objetivo - alfa_extremo)
    
    # 5. NORMALIZACIÓN: Elemento clave de la implementación para mantener ángulos en [-pi, pi]
    th[i] = (th[i] + pi) % (2 * pi) - pi

    # Guardar estado para la visualización
    O_historial.append(cin_dir(th,a))

  # Actualizar distancia y mostrar progreso
  dist = np.linalg.norm(np.subtract(objetivo, O_historial[-1][-1]))
  print ("\n- Iteracion " + str(iteracion) + ':')
  muestra_origenes(O_historial[-1])
  muestra_robot(O_historial, objetivo)
  print ("Distancia al objetivo = " + str(round(dist,5)))
  iteracion += 1

# Resultados finales
if dist <= EPSILON:
  print ("\n" + str(iteracion-1) + " iteraciones para converger.")
else:
  print ("\nNo hay convergencia tras " + str(iteracion-1) + " iteraciones.")

print ("- Umbral de convergencia epsilon: " + str(EPSILON))
print ("- Distancia al objetivo:          " + str(round(dist,5)))
print ("- Valores finales de las articulaciones:")
for i in range(len(th)):
  print ("  theta" + str(i+1) + " = " + str(round(th[i],3)))
for i in range(len(th)):
  print ("  L" + str(i+1) + "     = " + str(round(a[i],3)))