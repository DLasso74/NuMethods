# -*- coding: utf-8 -*-
"""
Proyecto Métodos Numéricos - WaveSolve
Soluciona la ecuación de onda para una malla de entrada con elementos triangulares.

@authors: Melissa Acosta - Alejandra Garzón - Darío Lasso
"""
#%% Librerías
from __future__ import division
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import numpy as np
import sympy
import meshio
#%% Funciones y setup
# Jacobiano 
def Jacob(dS, coord_el):
    return sympy.simplify(dS * coord_el)
# Jacobiano inv
def Jacob_inv(dS, coord_el):
    jac = Jacob(dS, coord_el)
    return sympy.Matrix([[jac[1, 1], -jac[0, 1]], [-jac[1, 0], jac[0, 0]]])/jac.det()
# Cuadratura Gaussiana
W1 = 1/6
W2 = 1/6
W3 = 1/6
Xi1 = 0.5
Xi2 = 0
Xi3 = 0.5
Nu1 = 0
Nu2 = 0.5
Nu3 = 0.5
GQPoints = sympy.Matrix([[Xi1, Nu1],[Xi2, Nu2],[Xi3, Nu3]])
# Funciones de Forma
S1 =  sympy.Matrix([0.5,0.5,0])
S2 =  sympy.Matrix([0.5,0,0.5])
S3 =  sympy.Matrix([0,0.5,0.5])
S = sympy.Matrix([[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5]]) # Matriz de Coordenadas nodales
dS = sympy.Matrix([[-1, 1, 0],[-1, 0, 1]])
# Velocidad de Onda
prof = float(input("Introduzca la profundidad del lago en metros: "))
c = np.sqrt(9.78*prof) # En m/s

#%% Preprocesamiento
# Malla
points, cells, point_data, cell_data, field_data = meshio.read("input.msh") # Se requiere una malla con triángulos uniformes
Xp = points[:, 0] # Puntos en X
Yp = points[:, 1] # Puntos en Y
Zp = points[:, 2] # Puntos en Z
connections =  cells["triangle"] # Conexiones de los nodos
dt = 0.1 # Time step
dy = Yp[connections[1,1]] - Yp[connections[1,0]] # Space step basandose en la distancia entre nodos de los triángulos

#%% Procesamiento
TNodes = np.size(points,0)
TElements = np.size(connections,0)
Kg = sympy.zeros(TNodes)
Mg = sympy.zeros(TNodes)
Hg = sympy.zeros(TNodes)
Fg = sympy.zeros(TNodes,1)
bg = sympy.zeros(TNodes,1)
for n in range(TNodes):
    Node0 = connections[n,0] # Nodo 1 de la celda
    Node1 = connections[n,1] # Nodo 2 de la celda
    Node2 = connections[n,2] # Nodo 3 de la celda
    x1 = Xp[Node0]
    x2 = Xp[Node1]
    x3 = Xp[Node2]
    y1 = Yp[Node0]
    y2 = Yp[Node1]
    y3 = Yp[Node2]
    coord_el = sympy.Matrix([[x1,y1],[x2,y2],[x3,y3]])
    # Rigidez local
    K = 0.5*sympy.det(Jacob(dS, coord_el)) * sympy.Matrix([[2, -1, -1],[-1, 2, -1],[-1, -1, 2]])
    # K Global
    Kg[Node0,Node0] =+ K[0,0] 
    Kg[Node0,Node1] =+ K[0,1] 
    Kg[Node0,Node2] =+ K[0,2] 
    Kg[Node1,Node0] =+ K[1,0] 
    Kg[Node1,Node1] =+ K[1,1] 
    Kg[Node1,Node2] =+ K[1,2] 
    Kg[Node2,Node0] =+ K[2,0] 
    Kg[Node2,Node1] =+ K[2,1] 
    Kg[Node2,Node2] =+ K[2,2] 
    # H Local
    H = sympy.det(Jacob(dS, coord_el))/24 * sympy.Matrix([[2, 1, 1],[1, 2, 1],[1, 1, 2]])
    # H Global
    Hg[Node0,Node0] =+ H[0,0] 
    Hg[Node0,Node1] =+ H[0,1] 
    Hg[Node0,Node2] =+ H[0,2] 
    Hg[Node1,Node0] =+ H[1,0] 
    Hg[Node1,Node1] =+ H[1,1] 
    Hg[Node1,Node2] =+ H[1,2] 
    Hg[Node2,Node0] =+ H[2,0] 
    Hg[Node2,Node1] =+ H[2,1] 
    Hg[Node2,Node2] =+ H[2,2]
    # M local
    M = 1/(c**2) * H
    # M Global
    Mg[Node0,Node0] =+ M[0,0] 
    Mg[Node0,Node1] =+ M[0,1] 
    Mg[Node0,Node2] =+ M[0,2] 
    Mg[Node1,Node0] =+ M[1,0] 
    Mg[Node1,Node1] =+ M[1,1] 
    Mg[Node1,Node2] =+ M[1,2] 
    Mg[Node2,Node0] =+ M[2,0] 
    Mg[Node2,Node1] =+ M[2,1] 
    Mg[Node2,Node2] =+ M[2,2]

B = 2*dt**2/dy**2 *Hg   
A = 2*Mg - dt**2 *Kg - B
Iteraciones = TNodes

U = sympy.zeros(Iteraciones, TNodes)
longitud = np.size(connections,0)
y = np.linspace(0, longitud, TNodes)
t = np.linspace(0, dt*Iteraciones, Iteraciones)
U[1, 25] = 1 # Condicion inicial

for cont_t in range(2, Iteraciones): # Solución para cada tiempo mayor a 2*dr
    for cont_y in range(1, TNodes - 1):
        U[cont_t - 1, cont_y] = A * U[cont_t, cont_y] + B * (U[cont_t, cont_y + 1] + U[cont_t, cont_y - 1]) - Mg * U[cont_y + 1, cont_y]
#%% Postprocesamiento

# Exportación
tri_mesh = {
        'points': points,
        'cells': {'triangle': cells}}
point_data = {'data': data}

meshio.write(
    "Datos.vtk",
    tri_mesh["points"],
    tri_mesh["cells"],
    point_data=point_data)

# Gráfico
figure = plt.figure() # Grafico de la superficie
ax = figure.gca(projection='3d')
ax.plot_trisurf(x, y, data, cmap=plt.cm.viridis)
