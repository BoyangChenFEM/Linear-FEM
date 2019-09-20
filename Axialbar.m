% ====================================
% Finite Element Programme
% for 1D linear axial tapered bar
% Author: Boyang CHEN, TU Delft
% Time:   Sep 2019
% ====================================

clc;
clear;
close all;

% Define geometrical parameters
L  = 300;   % mm
w1 = 45;    % mm
w2 = 22;    % mm
t  = 3;     % mm

% Define material parameters
E = 70000; % MPa

% Define load parameters
P = 3000;  % N 

% Define the mesh
nelems = 4;
nnodes = nelems+1;
node_coords = linspace(0,L,nnodes);
disp('Nodal coordinates (mm):')
disp(node_coords)

% Define the element connectivity matrix
elem_cnc    = zeros(nelems,2);
for i = 1: nelems
    elem_cnc(i,:) = [i,i+1];
end
disp('Connectivity matrix: ')
disp(elem_cnc)

% Define the boundary conditions (bcd); 
% dof: degree of freedom
bcd_dof   = 1; % dof number with imposed bcd
bcd_value = 0; % value of the imposed bcd

% Define the loads
load_dof   = nnodes; % dof number with imposed load 
load_value = P;      % value of the imposed load


% Initialize components of global matrix equation K a = f
% since only axial dof is considered, no. of dofs = no. of nodes
K = zeros(nnodes,nnodes);
f = zeros(nnodes,1);


% Calculate element stiffness matrix K_e and assemble to K using elem_cnc
for e = 1: nelems
    % calculate the length of the bar element
    L_e = node_coords(elem_cnc(e,2))-node_coords(elem_cnc(e,1));
    % calculate the cross-sectional area of the bar element using mid
    % coordinate of the bar
    x_middle = (node_coords(elem_cnc(e,1))+node_coords(elem_cnc(e,2)))/2;
    A_cross  = (w1 + (w2-w1)*x_middle/L)*t;
    % calculate the local stiffness matrix K_e
    K_e = E*A_cross/L_e*[1,-1;-1,1];
    % assemble K_e to K
    K(elem_cnc(e,:),elem_cnc(e,:))=K(elem_cnc(e,:),elem_cnc(e,:))+K_e(:,:);
end 

% cross-out rows with imposed displacement; here set them all to zero
% store this row for postprocessing later before crossing it out
Krow_bcd = K(bcd_dof,:); % storage
K(bcd_dof,:) = 0; % cross-out row in K
% f(bcd_dof)   = 0; % cross-out row in f; it's already zero

% if imposed displacement is not zero, then move U_D * column(bcd_node) to the RHS and set the column wrt U_D to zero
if (bcd_value ~= 0)
    f(:) = f(:) - bcd_value * K(:,bcd_dof);
end
% cross-out column of K wrt node with imposed displacement
K(:,bcd_dof) = 0;

% Apply the loads on the force vector
f(load_dof) = load_value;

% set diagonal term of K wrt imposed node to a dummy stiffness to avoid singularity
K(bcd_dof,bcd_dof) = 1; % dummy non-zero stiffness

% Solve for a:
a = K\f;

% update a to include the imposed displacements
a(bcd_dof) = bcd_value;

% obtain reaction force at node with imposed displacement
RF = Krow_bcd*a;

% Postprocessing
% Primary output
disp('Axial displacements (mm):')
disp(a')

% Derive strain
epsilon = zeros(nelems, 1);
epsilon(:)=(a(elem_cnc(:,2))-a(elem_cnc(:,1)))...
    ./(node_coords(elem_cnc(:,2))-node_coords(elem_cnc(:,1))).';
% for e = 1 : nelems
%     epsilon(e) = (a(elem_cnc(e,2))-a(elem_cnc(e,1)))...
%     /(node_coords(elem_cnc(e,2))-node_coords(elem_cnc(e,1)));
% end
disp('Axial strain in the elements:')
disp(epsilon')

% Derive stress
sigma = E * epsilon;
disp('Axial stress in the elements (MPa):')
disp(sigma')



