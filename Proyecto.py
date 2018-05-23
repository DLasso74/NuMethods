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
S = sympy.Matrix([[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5]])
dS = sympy.Matrix([[-1, 1, 0],[-1, 0, 1]])
# Velocidad de Onda
prof = float(input("Introduzca la profundidad del lago en metros: "))
c = np.sqrt(9.78*prof) # En m/s

x, y, t = sympy.symbols("x y t")
#%% Preprocesamiento
# Malla
points, cells, point_data, cell_data, field_data = meshio.read("input.msh")
Xp = points[:, 0] # Puntos en X
Yp = points[:, 1] # Puntos en Y
Zp = points[:, 2] # Puntos en Z
connections =  cells["triangle"] # Conexiones de los nodos

#%% Procesamiento
TNodes = np.size(points,0)
Kg = sympy.zeros(TNodes)
Mg = sympy.zeros(TNodes)
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
    # Matriz B
    B = Jacob_inv(dS, coord_el) * dS
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
    # M local
    M = sympy.det(Jacob(dS, coord_el))/24 * sympy.Matrix([[2, 1, 1],[1, 2, 1],[1, 1, 2]])
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
scale = 2
plot_args = {'rstride': 1, 'cstride': 1, 'cmap':"viridis",
             'linewidth': 0.1, 'antialiased': True, 'edgecolor': '#1e1e1e',
             'shade': True, 'alpha': 1.0, 'vmin': 0, 'vmax':scale}
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ani = animation.FuncAnimation(fig, wave_iter, range(nframes), blit=False,
                              fargs=(ax, X, Y, Z, Z0, dt, ntime_anim, L, scale,
                                     plot_args))
plt.show()