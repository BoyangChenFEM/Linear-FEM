# ====================================
# Finite Element Programme
# for 1D linear axial tapered bar
# Author: Boyang CHEN, TU Delft
# Time:   Sep 2019
# ====================================
import numpy as np
import numpy.linalg as la

# Define geometrical parameters
L  = 300   # mm
w1 = 45    # mm
w2 = 22    # mm
t  = 3     # mm

# Define material parameters
E = 70000 # MPa

# Define load parameters
P = 3000  # N 

# Define the mesh
nelems = 4
nnodes = nelems+1
node_coords = np.linspace(0,L,nnodes)
print('Nodal coordinates (mm):')
print(node_coords)

# Define the element connectivity matrix
elem_cnc = np.zeros([nelems,2]).astype(int)
for i in range(nelems):
    elem_cnc[i,:] = np.array([i,i+1]).astype(int)
print('Connectivity matrix: ')
print(elem_cnc+1)  

# Define the boundary conditions (bcd)
# dof: degree of freedom
bcd_dof   = 1 # dof number with imposed bcd
bcd_value = 0 # value of the imposed bcd

# Define the loads
load_dof   = nnodes # dof number with imposed load
load_value = P      # value of the imposed load

# Initialize components of the global matrix equation K a = f
# since only axial dof is considered, no. of dofs = no. of nodes
K = np.zeros([nnodes,nnodes])
f = np.zeros([nnodes,1])

# Calculate element stiffness matrix K_e and assemble to K using elem_cnc
for e in range(nelems):
    # calculate the length of the bar element
    L_e = node_coords[elem_cnc[e,1]]-node_coords[elem_cnc[e,0]]
    # calculate the cross-sectional area of the bar element using mid 
    # coordinate of the bar
    x_middle = (node_coords[elem_cnc[e,1]]+node_coords[elem_cnc[e,0]])/2
    A_cross  = (w1 + (w2-w1)*x_middle/L)*t
    # calculate the local stiffness matrix K_e
    K_e = E*A_cross/L_e*np.matrix([[1,-1],[-1,1]])
    # assemble K_e to K
    K[np.ix_(elem_cnc[e,:],elem_cnc[e,:])] += K_e[:,:]
    
# cross-out rows with imposed displacement; here set them all to zero
# store this row for postprocessing later before crossing it out
Krow_bcd = np.zeros([1,nnodes])
Krow_bcd[:] = K[bcd_dof-1,:] # storage
K[bcd_dof-1,:] = 0        # cross-out row in K

# if imposed displacement is not zero, then move U_D * column(bcd_node) to the RHS and set the column wrt U_D to zero
if (bcd_value != 0):
    f[:] -= bcd_value * K[:,bcd_dof-1].reshape((-1,1))

# cross-out column of K wrt node with imposed displacement
K[:,bcd_dof-1] = 0

# Apply the loads on the force vector
f[load_dof-1] = load_value

# set diagonal term of K wrt imposed node to a dummy stiffness to avoid singularity
K[bcd_dof-1,bcd_dof-1] = 1 # dummy non-zero stiffness
# Solve for a
a = la.solve(K,f)

# update a to include the imposed displacements
a[bcd_dof-1] = bcd_value

print('Axial displacements (mm):')
print(a)

# obtain reaction force at node with imposed displacement
#RF = np.matmul(Krow_bcd,a)
RF = Krow_bcd@a
#print('Reaction force:')
#print(RF)

# derive strain
epsilon = np.zeros([nelems,1])
du = np.zeros([nelems,1])
Le = np.zeros([nelems,1])
du[:]=a[elem_cnc[:,1]]-a[elem_cnc[:,0]]
Le[:]=(node_coords[elem_cnc[:,1]]-node_coords[elem_cnc[:,0]]).reshape((-1,1))
epsilon = du/Le
print('Axial strains:')
print(epsilon)

# derive stress
sigma = E*epsilon
print('Axial stresses (MPa) :')
print(sigma)
