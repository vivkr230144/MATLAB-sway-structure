% Name: Alicia Jiayun LAW
% CID: 01105518
% Coursework: Modelling of a Sway Structure

%%
function [ sway ] = Radius (type, radius, E, NODES, ELEMENTS, P, ...
    fdofs_restrained, fdofs_free, fdof_topload, load_fdofs, ...
    dofs_restrained, dofs_free, dof_topload, load_dofs, combimin, combimax)

% This function is used by fzero to determine the radius required for the
% lateral displacement of the top left node to be = 25mm

% Inputs;
% 1st row: inputs used for both truss and frames
% 2nd row: inputs used for frames
% 3rd row: inputs used for trusses

% Update to new derived parameters
A = pi*radius^2;           % mm2  - cross-sectional area of each member
EA = E*A;                  % N    - axial rigidity of all members
I = (1/4)*pi*radius^4 ;    % mm4  - second moment of area for a circle
EI = E*I ;                 % Nmm2 - flexural rigidity of all members

% Loading, elements, nodes and dofs are the same.

% Find the radius that gives 25mm lateral sway at the point of
% application of the top-most load
% (i) for FRAMES
if type == 1
    
    A = FRAMES ('5x3 frames question 5', NODES, ELEMENTS, fdofs_restrained, fdofs_free);
    A = A.assemble (EA,EI,P,load_fdofs); % obtain K and Fa
    A = A.solve (); % obtain D and F
    
    sway = A.D (fdof_topload) - 25 ;
    
% (ii) for TRUSS MIN SWAY
elseif type == 2
    
    B = TRUSS('5x3 min configuration', NODES, ELEMENTS, dofs_restrained, dofs_free);
    B = B.AddCrossMembers(combimin);
    B = B.assemble (EA, P, load_dofs);
    B = B.solve ();
    
    sway = B.U (dof_topload) - 25 ;
    
% (iii) for TRUSS MAX SWAY
elseif type == 3
    
    C = TRUSS('5x3 max configuration', NODES, ELEMENTS, dofs_restrained, dofs_free);
    C = C.AddCrossMembers(combimax);
    C = C.assemble (EA, P, load_dofs);
    C = C.solve ();
    
    sway = C.U (dof_topload) - 25 ;
    
else
    
    sway = NaN;
    
end


end

