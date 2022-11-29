"""
IMPORTS
"""
import numpy as np
import matplotlib.pyplot as plt
"""
FUNCTIONS
"""
# Theoretical pressure distribution for punch problem
def pressure(r, a, delta, E, V):
    G = E/(2*(1+V))
    ans = np.zeros_like(r)
    for i in range(len(r)):
        if (r[i] > a):
            ans[i] = np.nan
        else:
            ans[i] = (2*G*delta)/(np.pi*(1-V)*np.sqrt(a**2-r[i]**2))
    return ans
# Theoretical displacement for punch problem
def disp_z(r, a, delta, E, V):
    G = E/(2*(1+V))
    ans = np.zeros_like(r)
    for i in range(len(r)):
        if (r[i] > a):
            ans[i] = (2*delta/np.pi)*np.arcsin(a/r[i])
        else:   
            ans[i] = delta
    return ans
# Function to find index of nearest value to a search value in a given array
def find_nearest(array, value):
    diff = abs(array[0]-value)
    idx = 0
    for i in range(len(array)):
        if (abs(array[i]-value) < diff):
            diff = abs(array[i]-value)
            idx = i
    return idx
# Color map plot to show deformed and undefomed plots of the sample area with a variable of interest
def color_map(x, y, data, data_label, zoom, show=False, undeformed=True):
    tag = 'DEF'
    if undeformed:
        y, x = np.meshgrid(x,y)
        tag = 'UNDEF'
    plt.contourf(x, y, data, 50)
    cbar = plt.colorbar()
    cbar.set_label(data_label, rotation=270, labelpad=20)
    plt.xlabel('Radial Position [m]')
    plt.ylabel('Axial Position [m]')
    plt.tight_layout()
    plt.axis('scaled')
    if not undeformed:
        for i in range(np.shape(x)[0]):
            plt.plot(x[i,:], y[i,:], linewidth=0.5, color='#1a1a1a')
            plt.plot(x[:,i], y[:,i], linewidth=0.5, color='#1a1a1a')
    plt.xlim([-0.1*max(x[-1,:]), 1.1*max(x[-1,:])])
    plt.ylim([-0.1*max(y[:,-1]), 1.1*max(y[:,-1])])
    plt.savefig('%s/CMAP_%s_%s' %(img_path, data_label, tag), bbox_inches='tight')
    plt.xlim(zoom[0])
    plt.ylim(zoom[1])
    plt.savefig('%s/CMAP_ZOOM_%s_%s' %(img_path, data_label, tag), bbox_inches='tight')
    if show:
        plt.show()
    plt.close()
"""
CONSTANTS MAIN
"""
#-Problem Definition----------------#-----------------------------------# [units]
a = 7.5e-3                          # Punch Radius                      # [m]
delta = -1e-5                       # Punch Displacement                # [m]
R = 20*a                            # Sample Radius                     # [m]
H = R                               # Sample Height                     # [m]
E_s = 70e9                          # Elastic Modulus                   # [Pa]
V_s = 0.3                           # Poissons Ratio                    # [N.a.]
#-Plot Parameters-------------------#-----------------------------------# [units]
nms = 3                             # Node Marker Size 
scale = 1.5e3                       # Deformed Scale
img_path = './IMG/TEMP'             # Save path for images
#-Material Definition---------------#-----------------------------------# [units]
E = E_s                             # Elastic Modulus                   # [Pa]
v = V_s                             # Poisson's Ratio                   # [N.a.]
#-Structure Definition--------------#-----------------------------------# [units]
Lx = R                              # Length in X                       # [m]
Ly = H                              # Length in Y                       # [m]
#-Mesh------------------------------#-----------------------------------# [units]
nx = 50                             # Number of X Elements              # [#]
ny = 50                             # Number of Y Elements              # [#]
m = nx*ny                           # Total Number of Elements          # [#]
Dx = Lx/nx                          # Size of the Element in X          # [m]
Dy = Ly/ny                          # Size of the Element in Y          # [m]
#-Nodes-----------------------------#-----------------------------------# [units]
N = (nx+1)*(ny+1)                   # Total Number of Nodes             # [#]
ix = np.arange(nx + 1)              #
iy = np.arange(ny + 1)              #
#-Node Coordinates------------------#-----------------------------------# [units]
x = ix*Dx                           # Node X Positions                  # [m]
y = iy*Dy                           # Node Y Positions                  # [m]
xn, yn = np.meshgrid(x,y)           # Node Mesh with X and Y            # [m]
#-Stiffness Matrix------------------#-----------------------------------# [units] 
K = np.zeros((2*N,2*N))             # Stifness Matrix
#-Plane Strain Elasticity Matrix----#-----------------------------------# [units]
Em  = np.array([                    # Elastic Matrix
    [(1-v),   v,      0],           #
    [v,      (1-v),   0],           #
    [0,      0,      (1-2*v)/2]])   #
Em *= E/((1+v)*(1-2*v))             #
T = np.array([                      # Derivative Mapping Matrix
    [1,  0,  0,  0],                #
    [0,  0,  0,  1],                #
    [0,  1,  1,  0]])               #
#-Gauss points----------------------#-----------------------------------# [units]
g = np.array([                      # Gauss point matrix
    [-1, -1],                       #
    [ 1, -1],                       #
    [-1,  1],                       #
    [ 1,  1]], dtype='d')           #
g /= np.sqrt(3)                     #
#-----------------------------------#-----------------------------------#
"""
MAIN
"""
#-----------------------------------------------------------------------#
# K Matrix Assembly                                                     #
#-----------------------------------------------------------------------#
# If broken switch order of loops
for ey in range(ny):
    for ex in range(nx):
        i = ey*(nx+1)+ex        # Index of bottom left node
        j = i+(nx+1)            # Index of top left node
        x1 = xn[ey,ex]; x2 = xn[ey,ex+1]; x3 = x1; x4 = x2
        y1 = yn[ey,ex]; y2 = y1; y3 = yn[ey+1,ex]; y4 = y3
        Kq = np.zeros((8,8))
        # Integral aproximation using Gauss points
        for gg in range(4):     # gg == gauss point index
            xx = g[gg,0]        # Gauss point for X
            yy = g[gg,1]        # Gauss point for Y
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
        # EDOFS - DOF indicies for elemnt q
        edofs = np.array([2*i,2*i+1,2*(i+1),2*(i+1)+1,2*j,2*j+1,2*(j+1),2*(j+1)+1])
        K[edofs.reshape(-1, 1), edofs.reshape(1, -1)] += Kq
#-----------------------------------------------------------------------#
# Boundary Conditions                                                   #
#-----------------------------------------------------------------------#
u = np.nan*np.ones(2*N)
for i in range(N):
    i_x = 2*i
    i_y = 2*i+1
    row = int(np.floor(i/(nx+1)))
    col = i%(ny+1)
    # Fixed Bottom Row
    if (row == 0/Dy):
        u[i_x] = 0
        u[i_y] = 0
    # Symmetry
    elif (col == 0/Dy) and not (row == H/Dy):
        u[i_x] = 0
    # Perscribed Displacements for Punch
    elif (row == H/Dy) and (col <= a/Dx):
        u[i_x] = 0
        u[i_y] = delta
#-----------------------------------------------------------------------#
# Loads                                                                 #
#-----------------------------------------------------------------------#
f = np.zeros(2*N)
#-----------------------------------------------------------------------#
# DOF Indicies                                                          #
#-----------------------------------------------------------------------#
fixed_dofs = np.array([i for i in range(len(u)) if np.isnan(u[i])!=np.isnan(np.nan)])
free_dofs = np.array([i for i in range(len(u)) if np.isnan(u[i])==np.isnan(np.nan)])
#-----------------------------------------------------------------------#
# Solution                                                              #
#-----------------------------------------------------------------------#
# Perscribed displacement solution notation for displacement
u_1 = u[free_dofs]
u_2 = u[fixed_dofs]
# Perscribed displacement solution notation for forces
f_1 = f[free_dofs]
f_2 = f[fixed_dofs]
# Perscribed displacement solution notation for stiffness matrix
K_11 = K[free_dofs.reshape(-1, 1), free_dofs.reshape(1, -1)]
K_12 = K[free_dofs.reshape(-1, 1), fixed_dofs.reshape(1, -1)]
K_21 = K[fixed_dofs.reshape(-1, 1), free_dofs.reshape(1, -1)]
K_22 = K[fixed_dofs.reshape(-1, 1), fixed_dofs.reshape(1, -1)]
# Perscribed displacement solution
u[free_dofs] = np.linalg.solve(K_11, (f_1-K_12@u_2))
#-----------------------------------------------------------------------#
# Plotting Displacements                                                #
#-----------------------------------------------------------------------#
# Python Displacements
dx = u[0::2]
dy = u[1::2]
# Abaqus Data
x_a = np.loadtxt('./Abaqus_Validation/Abaqus_X_Data.csv', delimiter=',')
dy_a = np.loadtxt('./Abaqus_Validation/Abaqus_Displacement_Data.csv', delimiter=',')
sy_a = np.loadtxt('./Abaqus_Validation/Abaqus_Stress_Data.csv', delimiter=',')
dy_a_s = np.loadtxt('./Abaqus_Validation/Abaqus_Displacement_Slip_Data.csv', delimiter=',')
sy_a_s = np.loadtxt('./Abaqus_Validation/Abaqus_Stress_Slip_Data.csv', delimiter=',')
# Theoretical Displacements
dy_t = disp_z(x, a, delta, E, v)
# Node position vectors
x_vec = np.zeros_like(u)
for i in range(N):
    i_x = 2*i
    i_y = 2*i+1
    row = int(np.floor(i/(nx+1)))
    col = i%(ny+1)
    x_vec[i_x] = col*Dx
    x_vec[i_y] = row*Dy
x_x = x_vec[0::2]
x_y = x_vec[1::2]
# Displacement plotting
plt.plot([], [], ' ', label='Scale Factor: %.2E' %(scale))
plt.plot(x_x[-(nx+1):], x_y[-(ny+1):]+dy_t*scale, color='r', alpha=0.5, label='Theoretical deformation')
plt.plot(x_a, x_y[-(ny+1)]+dy_a*scale, color='b', alpha=0.5, label='Abaqus deformation')
plt.scatter(x_x, x_y, s=nms, label='Undeformed')
plt.scatter(x_x+dx*scale, x_y+dy*scale, s=nms, label='Deformed Python')
plt.xlabel('Radial Position [m]')
plt.ylabel('Axial Position [m]')
plt.legend(loc='upper left', bbox_to_anchor=(1.04, 1))
plt.axis('scaled')
plt.savefig('%s/Displacements_Plot' %(img_path), bbox_inches='tight')
plt.show()
plt.close()
#-----------------------------------------------------------------------#
# Plotting Stress                                                       #
#-----------------------------------------------------------------------#
strain = np.nan*np.ones((nx,ny,3))
stress = np.nan*np.ones_like((strain))
stress_vm = np.nan*np.ones((nx,ny))
Uv = np.nan*np.ones((nx,ny))
ex_vec = np.nan*np.ones((nx,ny,2))
deformations = np.nan*np.ones((nx,ny,3))
ex_vec_def = np.nan*np.ones((nx,ny,2))
for ey in range(ny):
    for ex in range(nx):
        i = ey*(nx+1)+ex        # Index of bottom left node
        j = i+(nx+1)            # Index of top left node
        i_x = 2*i; i_y = 2*i+1
        j_x = 2*j; j_y = 2*j+1
        # Position Values
        x1 = x_vec[i_x]; x2 = x_vec[i_x+2]; x3 = x_vec[j_x]; x4 = x_vec[j_x+2]
        y1 = x_vec[i_y]; y2 = x_vec[i_y+2]; y3 = x_vec[j_y]; y4 = x_vec[j_y+2]
        ex_vec[ex,ey,0] = 0.25*(x1+x2+x3+x4)
        ex_vec[ex,ey,1] = 0.25*(y1+y2+y3+y4)
        # Displacement Values
        u1 = u[i_x]; u2 = u[i_x+2]; u3 = u[j_x]; u4 = u[j_x+2]
        v1 = u[i_y]; v2 = u[i_y+2]; v3 = u[j_y]; v4 = u[j_y+2]
        # Deformations
        deformations[ex,ey,0] = 0.25*(u1+u2+u3+u4)
        deformations[ex,ey,1] = 0.25*(v1+v2+v3+v4)
        deformations[ex,ey,2] = np.sqrt(deformations[ex,ey,0]**2+deformations[ex,ey,1]**2)
        # Deformed positions
        ex_vec_def[ex,ey,0] = ex_vec[ex,ey,0]+deformations[ex,ey,0]*scale
        ex_vec_def[ex,ey,1] = ex_vec[ex,ey,1]+deformations[ex,ey,1]*scale
        # Strain Values
        strain[ex,ey,0] = 0.5*((u2-u1)/(x2-x1)+(u4-u3)/(x4-x3))                                  # In x-x
        strain[ex,ey,1] = 0.5*((v3-v1)/(y3-y1)+(v4-v2)/(y4-y2))                                  # In y-y
        strain[ex,ey,2] = 0.5*((u3-u1)/(y3-y1)+(u4-u2)/(y4-y2)+(v2-v1)/(x2-x1)+(v4-v3)/(x4-x3))  # In x-y
        # Stress Values
        stress[ex,ey,:] = Em@strain[ex,ey,:]
        stress_vm[ex,ey] = np.sqrt(stress[ex,ey,0]**2+stress[ex,ey,1]**2-(stress[ex,ey,0]*stress[ex,ey,1])+3*stress[ex,ey,2]**2)
        #Strain Energy Values
        Uv[ex,ey] = 0.5*strain[ex,ey,:]@stress[ex,ey,:]
plt.plot(ex_vec[:,-1,0], stress[:,-1,1], label='Y-Y Stress Python')
plt.plot(ex_vec[:,-1,0], pressure(ex_vec[:,-1,0], a, delta, E, v), label='Y-Y Pressure Theory')
plt.plot(x_a, sy_a, color='b', alpha=0.5, label='Y-Y Stress Abaqus')
plt.xlabel('Radial Position [m]')
plt.ylabel('Stress Magnitude [Pa]')
plt.legend(loc='upper left', bbox_to_anchor=(1.04, 1))
plt.savefig('%s/Stress_Plot' %(img_path), bbox_inches='tight')
plt.show()
plt.close()
#-----------------------------------------------------------------------#
# Plotting Deformed Values                                              #
#-----------------------------------------------------------------------#
# Deformation color maps
color_map(ex_vec_def[:,:,0], ex_vec_def[:,:,1], deformations[:,:,0], 'X Deformation [m]', [[0, R*0.2], [H*0.8, H]], undeformed=False)
color_map(ex_vec_def[:,:,0], ex_vec_def[:,:,1], deformations[:,:,1], 'Y Deformation [m]', [[0, R*0.2], [H*0.8, H]], undeformed=False)
color_map(ex_vec_def[:,:,0], ex_vec_def[:,:,1], deformations[:,:,2], 'Total Deformation [m]', [[0, R*0.2], [H*0.8, H]], undeformed=False)
# Strain color maps
color_map(ex_vec_def[:,:,0], ex_vec_def[:,:,1], strain[:,:,0], 'X-X Strain [Pa]', [[0, R*0.2], [H*0.8, H]], undeformed=False)
color_map(ex_vec_def[:,:,0], ex_vec_def[:,:,1], strain[:,:,1], 'Y-Y Strain [Pa]', [[0, R*0.2], [H*0.8, H]], undeformed=False)
color_map(ex_vec_def[:,:,0], ex_vec_def[:,:,1], strain[:,:,2], 'X-Y Strain [Pa]', [[0, R*0.2], [H*0.8, H]], undeformed=False)
# Stress color maps
color_map(ex_vec_def[:,:,0], ex_vec_def[:,:,1], stress[:,:,0], 'X-X Stress [Pa]', [[0, R*0.2], [H*0.8, H]], undeformed=False)
color_map(ex_vec_def[:,:,0], ex_vec_def[:,:,1], stress[:,:,1], 'Y-Y Stress [Pa]', [[0, R*0.2], [H*0.8, H]], undeformed=False)
color_map(ex_vec_def[:,:,0], ex_vec_def[:,:,1], stress[:,:,2], 'X-Y Stress [Pa]', [[0, R*0.2], [H*0.8, H]], undeformed=False)
color_map(ex_vec_def[:,:,0], ex_vec_def[:,:,1], stress_vm, 'Von Mises Stress [Pa]', [[0, R*0.2], [H*0.8, H]], undeformed=False)
#-----------------------------------------------------------------------#
# Plotting Undeformed Values                                            #
#-----------------------------------------------------------------------#
# Deformation color maps
color_map(ex_vec[:,-1,0], ex_vec[-1,:,1], deformations[:,:,0], 'X Deformation [m]', [[0, R*0.2], [H*0.8, H]])
color_map(ex_vec[:,-1,0], ex_vec[-1,:,1], deformations[:,:,1], 'Y Deformation [m]', [[0, R*0.2], [H*0.8, H]])
color_map(ex_vec[:,-1,0], ex_vec[-1,:,1], deformations[:,:,2], 'Total Deformation [m]', [[0, R*0.2], [H*0.8, H]])
# Strain color maps
color_map(ex_vec[:,-1,0], ex_vec[-1,:,1], strain[:,:,0], 'X-X Strain [Pa]', [[0, R*0.2], [H*0.8, H]])
color_map(ex_vec[:,-1,0], ex_vec[-1,:,1], strain[:,:,1], 'Y-Y Strain [Pa]', [[0, R*0.2], [H*0.8, H]])
color_map(ex_vec[:,-1,0], ex_vec[-1,:,1], strain[:,:,2], 'X-Y Strain [Pa]', [[0, R*0.2], [H*0.8, H]])
# Stress color maps
color_map(ex_vec[:,-1,0], ex_vec[-1,:,1], stress[:,:,0], 'X-X Stress [Pa]', [[0, R*0.2], [H*0.8, H]])
color_map(ex_vec[:,-1,0], ex_vec[-1,:,1], stress[:,:,1], 'Y-Y Stress [Pa]', [[0, R*0.2], [H*0.8, H]])
color_map(ex_vec[:,-1,0], ex_vec[-1,:,1], stress[:,:,2], 'X-Y Stress [Pa]', [[0, R*0.2], [H*0.8, H]])
color_map(ex_vec[:,-1,0], ex_vec[-1,:,1], stress_vm, 'Von Mises Stress [Pa]', [[0, R*0.2], [H*0.8, H]])
#-----------------------------------------------------------------------#
# Validation Error Output                                               #
#-----------------------------------------------------------------------#
# Model Displacement Comparison
error = np.nan*np.ones_like(dy_t)
error_a = np.nan*np.ones_like(dy_t)
idx_start = int(H/Dy)*(nx+1)
for i in range(len(error)):
    idx_y = i+idx_start
    error[i] = abs((dy[idx_y] - dy_t[i])/dy_t[i]*100)
    dy_a_current = dy_a[find_nearest(x_a, x_x[-(nx+1):][i])]
    error_a[i] = abs((dy[idx_y] - dy_a_current)/dy_a_current*100)
print('DISPLACEMENTS:')
print('Theoretical Error:', max(error))
print('Abaqus Error:', max(error_a))
# Stress Error
x_pos = ex_vec[:,-1,0]
stress_YY = stress[:,-1,1]
stress_YY_t = pressure(x_pos, a, delta, E, v)
error = np.nan*np.ones_like(stress_YY)
error_a = np.nan*np.ones_like(error)
for i in range(len(stress_YY)):
    if (x_pos[i] <= a):
        error[i] = abs((stress_YY_t[i] - stress_YY[i])/stress_YY_t[i]*100)
        sy_a_current = sy_a_s[find_nearest(x_a, x_pos[i])]
        error_a[i] = abs((sy_a_current - stress_YY[i])/sy_a_current*100)
print('STRESS:')
print('Theoretical Error:', max(error))
print('Abaqus Error:', max(error_a))    