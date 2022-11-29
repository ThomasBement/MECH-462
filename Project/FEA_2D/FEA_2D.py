# ---------------------------------------- #
# FEA_2D [Python File]
# Written By: Thomas Bement
# Created On: 2022-11-14
# ---------------------------------------- #

"""
IMPORTS
"""
import os

import numpy as np
import matplotlib.pyplot as plt

"""
GLOBAL CONSTANTS
"""
save_path = './DATA/STORAGE'
img_path = './IMG'

# Derivative mapping matrix
T = np.array([
    [1,  0,  0,  0],
    [0,  0,  0,  1],
    [0,  1,  1,  0]])  

# Gauss point matrix
g = np.array([
    [-1, -1],
    [ 1, -1],
    [-1,  1],
    [ 1,  1]], dtype='d')
g /= np.sqrt(3)


"""
FUNCTIONS
"""
def gen_mesh(len_x, len_y, elem_x, elem_y):
    # Number of nodes
    node_x = elem_x+1
    node_y = elem_y+1
    # Node position axies
    node_pos_x = np.linspace(0, len_x, node_x)
    node_pos_y = np.linspace(0, len_y, node_y)
    # Assign array variables and initialize to nan for safe storage
    x_vec_n = np.nan*np.ones((node_x, node_y, 2))
    x_vec_e = np.nan*np.ones((elem_x, elem_y, 2))
    for ey in range(elem_y):
        for ex in range(elem_x):
            # Calculate nodal positions
            x1 = node_pos_x[ex];    x2 = node_pos_x[ex+1];  x3 = x1;                x4 = x2
            y1 = node_pos_y[ey];    y2 = y1;                y3 = node_pos_y[ey+1];  y4 = y3
            # Calculate element positions
            x_e = 0.25*(x1+x2+x3+x4)
            y_e = 0.25*(y1+y2+y3+y4)
            # Write node vector
            x_vec_n[ex,ey,:] = [x1, y1];        x_vec_n[ex+1,ey,:] = [x2, y2]
            x_vec_n[ex,ey+1,:] = [x3, y3];      x_vec_n[ex+1,ey+1,:] = [x4, y4]
            # Write element vector
            x_vec_e[ex,ey,:] = [x_e, y_e]
    return x_vec_n, x_vec_e

def plot_deformed(node_pos, deformed_pos):
    plt.scatter(node_pos[:,:,0], node_pos[:,:,1])
    plt.show()
    plt.close()

def elasticity_matrix(E, V, condition='strain'):
    if (condition == 'strain'):
        # Plane strain elastic matrix
        Em  = np.array([
                [(1-V),     V,      0],
                [V,         (1-V),  0],
                [0,         0,      (1-2*V)/2]])
        Em *= E/((1+V)*(1-2*V))
    return Em

def assemble_k(elem_x, elem_y, x_vec_n, Em):
    N = (elem_x+1)*(elem_y+1)               # Total number of nodes
    K = np.nan*np.ones((2*N,2*N))           # Stifness Matrix
    for ey in range(elem_y):                #
        for ex in range(elem_x):            #
            i = ey*(nx+1)+ex                # Index of bottom left node
            j = i+(nx+1)                    # Index of top left node
            [x1, y1] = x_vec_n[ex,ey,:];  [x2, y2] = x_vec_n[ex+1,ey,:]
            [x3, y3] = x_vec_n[ex,ey+1,:];  [x4, y4] = x_vec_n[ex+1,ey+1,:]
            Kq = np.zeros((8,8))            #
            for gg in range(4):             # Integral aproximation using Gauss points
                xx = g[gg,0]                # Gauss point for X
                yy = g[gg,1]                # Gauss point for Y
                Nxx = np.array([-0.25*(1-yy), 0.25*(1-yy), -0.25*(1+yy), 0.25*(1+yy)])
                Nyy = np.array([-0.25*(1-xx), -0.25*(1+xx), 0.25*(1-xx), 0.25*(1+xx)])
                J = np.array([ # Jacobian Matrix
                    [Nxx@[x1,x2,x3,x4], Nxx@[y1,y2,y3,y4]],
                    [Nyy@[x1,x2,x3,x4], Nyy@[y1,y2,y3,y4]]])
                dJ = np.linalg.det(J)
                iJ = np.linalg.inv(J)
                A = np.zeros((4,4)); A[0:2, 0:2] = iJ; A[2:4, 2:4] = iJ
                A = T@A
                G = np.array([ # G Matrix
                    [Nxx[0], 0,  Nxx[1], 0,  Nxx[2], 0,  Nxx[3], 0],
                    [Nyy[0], 0,  Nyy[1], 0,  Nyy[2], 0,  Nyy[3], 0],
                    [0,  Nxx[0], 0,  Nxx[1], 0,  Nxx[2], 0,  Nxx[3]],
                    [0,  Nyy[0], 0,  Nyy[1], 0,  Nyy[2], 0,  Nyy[3]]])
                B = A@G # B Matrix
                Kq += dJ*B.T@(Em@B) # ELement Stiffness Matrix
                # Bg[]:,:,gg] = B         # ????
            # EDOFS - DOF indicies for elemnt q
            edofs = np.array([2*i,2*i+1,2*(i+1),2*(i+1)+1,2*j,2*j+1,2*(j+1),2*(j+1)+1])
            K[edofs.reshape(-1, 1), edofs.reshape(1, -1)] += Kq
    return K

"""
CONSTANTS
"""
#-Problem Definition--------#---------------------------------------# [units]
a = 7.5e-3                  # Punch radius                          # [m]
delta = -1e-5               # Punch displacement                    # [m]
R = 20*a                    # Sample radius                         # [m]
H = R                       # Sample height                         # [m]
E_s = 70e9                  # Elastic modulus                       # [Pa]
V_s = 0.3                   # Poissons ratio                        # [N.a.]
#-Plot parameters-----------#---------------------------------------# [units]
nms = 3                     # Node marker size                      #
scale = 1e3                 # Deformed scale                        #
#-Mesh----------------------#---------------------------------------# [units]
nx = 50                     # Number of X elements                  # [#]
ny = 50                     # Number of Y elements                  # [#]

"""
MAIN
"""
node_pos, elem_pos = gen_mesh(R, H, nx, ny)
E_M = elasticity_matrix(E_s, V_s, 'strain')
K_M = assemble_k(nx, ny, node_pos, E_M)
print(K_M)