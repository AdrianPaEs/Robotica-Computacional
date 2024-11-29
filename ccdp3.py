#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Universidad de La Laguna
# Escuela Superior de Ingeniería y Tecnología
# Grado en Ingeniería Informática
# Asignatura: Robótica Computacional
# Curso: 4º
# Práctica 2: Cinemática Inversa CCD prismática limites
# Autor: Adrián Padrón Espinosa (alu01013232@ull.edu.es)

# Ejecutar python ccdp3.py 10 10

import sys
from math import *
import numpy as np
import matplotlib.pyplot as plt
import colorsys as cs

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
  
  #input("Pulse una tecla para continuar...")
  plt.close()

def matriz_T(d,th,a,al):
   
  return [[cos(th), -sin(th)*cos(al),  sin(th)*sin(al), a*cos(th)]
         ,[sin(th),  cos(th)*cos(al), -sin(al)*cos(th), a*sin(th)]
         ,[      0,          sin(al),          cos(al),         d]
         ,[      0,                0,                0,         1]
         ]

def cin_dir(th,a):
  #Sea 'th' el vector de thetas
  #Sea 'a'  el vector de longitudes
  T = np.identity(4)
  o = [[0,0]]
  for i in range(len(th)):
    T = np.dot(T,matriz_T(0,th[i],a[i],0))
    tmp=np.dot(T,[0,0,0,1])
    o.append([tmp[0],tmp[1]])
  return o

# ******************************************************************************
# Cálculo de la cinemática inversa de forma iterativa por el método CCD

# valores articulares arbitrarios para la cinemática directa inicial
th=[0.,0.,0.]
a =[5.,5.,5.]
L = sum(a) # variable para representación gráfica
EPSILON = .01

tipo_articulacion = [0, 0, 0] # Tipo de articulación - (0) rotacional, (1) prismática
limite_sup = [45, 90, 90]     # Límite superior de las articulaciones, por encima eje X
limite_inf = [-45, -90, -90]    # Límite inferior de las articulaciones, por debajo eje X

plt.ion() # Modo Interactivo

# introducción del punto para la cinemática inversa
if len(sys.argv) != 3:
  sys.exit("python " + sys.argv[0] + " x y")
objetivo=[float(i) for i in sys.argv[1:]]
O=cin_dir(th,a)

# Calculamos la posicion inicial
print ("- Posicion inicial:")
muestra_origenes(O)

dist = float("inf")
prev = 0.
iteracion = 1
while (dist > EPSILON and abs(prev-dist) > EPSILON/100.):
  prev = dist           
  O=[cin_dir(th,a)]
  num_articulaciones = len(th)
  for i in range(num_articulaciones):
    j = num_articulaciones - i - 1  # Articulación actual
    if [tipo_articulacion[j] == 0 for j in range(num_articulaciones)][i]: # Articulacion rotatoria
      
      v1 = np.subtract(objetivo, O[i][j]) # Vector entre el punto final y el objetivo
      v2 = np.subtract(O[i][-1], O[i][j]) # Vector entre el origen de la articulación y el punto final
      
      alpha1 = atan2(v1[1], v1[0]) # Ángulo entre el eje X y el vector v1
      alpha2 = atan2(v2[1], v2[0]) # Ángulo entre el eje X y el vector v2
      
      th[j] += alpha1 - alpha2 # Actualizamos theta
      
      # while th[j] > pi: th[j] -= 2 * pi
      # while th[j] < -pi: th[j] += 2 * pi
      
      if (th[j] > np.radians(limite_sup[j])):
        th[j] = np.radians(limite_sup[j]) # Comprobamos los límites
      elif (th[j] < np.radians(limite_inf[j])):
        th[j] = np.radians(limite_inf[j]) # Comprobamos los límites

    elif [tipo_articulacion[j] == 1 for j in range(num_articulaciones)][i]: # Articulacion prismatica
      print("Articulacion prismatica")
      w = np.sum(th[:j + 1]) # Sumatorio de 'th' de 0 hasta 'actual' (incluido)
      d = np.dot([np.cos(w), np.sin(w)], np.subtract(objetivo, O[i][-1]))  
      a[j] += d   # Actualizamos 'a'
      
      # Comprobamos los límites
      if (a[j] > limite_sup[j]):
        a[j] = limite_sup[j]
      elif (a[j] < limite_inf[j]):
        a[j] = limite_inf[j]    
    else:
      print("Tipo de articulacion no reconocido, tiene que ser rotatoria (0) o prismatica (1)")
      exit(1)
      
    O.append(cin_dir(th,a))

  dist = np.linalg.norm(np.subtract(objetivo,O[-1][-1]))
  print ("\n- Iteracion " + str(iteracion) + ':')
  muestra_origenes(O[-1])
  muestra_robot(O,objetivo)
  print ("Distancia al objetivo = " + str(round(dist,5)))
  iteracion+=1
  O[0]=O[-1]

if dist <= EPSILON:
  print ("\n" + str(iteracion) + " iteraciones para converger.")
else:
  print ("\nNo hay convergencia tras " + str(iteracion) + " iteraciones.")
print ("- Umbral de convergencia epsilon: " + str(EPSILON))
print ("- Distancia al objetivo:          " + str(round(dist,5)))
print ("- Valores finales de las articulaciones:")
for i in range(len(th)):
  print ("  theta" + str(i+1) + " = " + str(round(th[i],3)))
for i in range(len(th)):
  print ("  L" + str(i+1) + "     = " + str(round(a[i],3)))
