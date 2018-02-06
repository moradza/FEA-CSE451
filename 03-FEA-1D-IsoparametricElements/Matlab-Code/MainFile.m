
clc; clear all;

%% ME 471 / AE 420 / 
% Programming Assignment #2: 1D elastostatic FE code
% includes: 
% - isoparametric linear and quadratic elements
% - numerical integration using Ng = 1,2,3,4 gauss points
% - DLOAD function for user-defined body load b(x)
% - UMAT function for user-defined elasticity modulus E(x)


% Variables definition

% nnode_ele, the (integer) number of nodes per element
% node_dof, the (integer) number of degrees of freedom per node
% edof, the (integer) number of degrees of freedom per element


filename = 'input.dat'; 

