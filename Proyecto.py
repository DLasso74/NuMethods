# -*- coding: utf-8 -*-
"""
Proyecto Métodos Numéricos


@authors: Melissa Acosta - Alejandra Garzón - Darío Lasso
"""
#%% Librerías
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import meshio
import math

#%% Malla
points, cells, point_data, cell_data, field_data = meshio.read("input.msh")
x = points[:, 0] # Puntos en X
y = points[:, 1] # Puntos en Y
z = points[:, 2] # Puntos en Z
connections =  cells["quad"] # Conexiones de los nodos
#%% Cálculos

#%% Gráfico
