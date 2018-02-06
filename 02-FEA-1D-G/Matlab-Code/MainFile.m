
clc; clear all;

%% ME 471 / AE 420 / 
% Programming Assignment #2: 1D elastostatic FE code
% Prof. Mariana Silva

% Variables definition

% nnode_ele, the (integer) number of nodes per element
% node_dof, the (integer) number of degrees of freedom per node
% edof, the (integer) number of degrees of freedom per element

NNODE_ELE = 2; %Number of nodes per element
DOF_NODE = 1; %Number of degrees of freedom per node
EDOF = NNODE_ELE*DOF_NODE; %number of degrees of freedom per element

filename = 'input.dat'; 

