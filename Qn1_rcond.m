clear
clc

% CI2-221: Computational Methods II
% Coursework - Modeling of a Sway Structure
% Name: Alicia Jiayun LAW

%%

% Starting Values
alpha = mod(12,5) + 3 ;
beta = mod(15,5) + 3 ;
N = min(alpha,beta); % Number of bays
M = max(alpha,beta); % Number of Storeys
n_nodes = (N+1)*(M+1); % Number of Nodes

elements = 35; % Number of Elements w/o cross members

%%

% Parameters - all in N and mm for simplicity
L = 1000; % mm - length of horizontal and vertical members
radius = 25; % mm - radius of each member
E = 2e5; % N/mm2 - Youngs Modulus
P = 100e3; % N - nodal loading

% Derived Parameters
A = pi*radius^2; % mm2 - cross-sectional area of each member
EA = E*A; % N - axial rigidity of all members

%% SPECIFYING NODES AND ELEMENTS

% Initialize
NODES.coords = zeros (n_nodes, 2);
NODES.dofs = zeros (n_nodes, 2);

for i = 1 : M+1     % loop over rows
    for j = 1 : N+1 % loop over columns
        
        node_number = 6*(j-1)+i; % Current node number being specified
        
        % (i)   Specifying nodal x-y coordinates
        NODES.coords (node_number,:) = [L*(j-1) , L*(i-1)];
        
        % (ii)  Specifying nodal dofs
        NODES.dofs (node_number,:) = [node_number*2 - 1 , node_number*2];
        
    end
end

% (iii) Specifying Element-Nodal connectivity
% Initialize
ELEMENTS.nodes = zeros (elements,2);
ELEMENTS.nodes (1,:) = [1 2]; % The nodes for the first element
for e = 2: elements
    
    % (for most vertical bars, with some exceptions:)
    % The first node of my current element is the last node of my previous element
    % The second node of my current element is +1 of my first node
    ELEMENTS.nodes (e,:) = [ELEMENTS.nodes(e-1,2),ELEMENTS.nodes(e-1,2)+1];
    
    % For vertical elements
    if e < 20 && mod(e-1,5) == 0 % Elements in a New Column
        ELEMENTS.nodes (e,:) = [1+6*(e-1)/5,2+6*(e-1)/5];
        
    % For horizontal elements
    elseif e > 20 && e <= 35
        
        if mod(e-21,3)==0 % Elements in a New Row
            ELEMENTS.nodes (e,:) = [2+(e-21)/3,2+(e-21)/3+6];
        else
            ELEMENTS.nodes (e,:) = [ELEMENTS.nodes(e-1,2),ELEMENTS.nodes(e-1,2)+6];
        end
        
    else
        
    end
end

% Manually input and omit cross member's element-nodal connectivity with 
% comments where necessary:

% For one cross member for each storey
ELEMENTS.nodes (end+1,:) = [1 8];
ELEMENTS.nodes (end+1,:) = [2 9];
ELEMENTS.nodes (end+1,:) = [3 10];
ELEMENTS.nodes (end+1,:) = [4 11];
ELEMENTS.nodes (end+1,:) = [5 12];

% 5 cross members but not one in each storey
% ELEMENTS.nodes (end+1,:) = [7 14];

% Degrees of freedom & other parameters
dofs_free = [3:12, 15:24, 27:36, 39:48];
dofs_restrained = [1,2,13,14,25,26,37,38]; % known nodal x-y dofs due to fixed ends at nodes 1, 7, 13 and 19

%% Assemble the stiffness matrix, K, and force vector, fa
load_dofs = [3, 5, 7, 9, 11];

obj = TRUSS('5x3', NODES, ELEMENTS, dofs_restrained, dofs_free);
obj = obj.assemble (EA,P,load_dofs);
obj = obj.solve ();

% Extract KFF for assessment
KFF = obj.K(dofs_free, dofs_free);

%% Assess the stability of the truss
% A. Computationally, by checking rcond value of KFF
rcond_KFF = rcond(KFF) ;

fprintf('\nAssessing truss stability computationally:\nrcond(KFF) = %g\n',rcond_KFF);

%% Results:

% When placed one in every floor
% No. of Additional Cross Members and the Corresponding rcond value
% 1 : rcond(KFF) = 2.38962e-34
% 2 : rcond(KFF) = 2.07753e-34
% 3 : rcond(KFF) = 2.07753e-34
% 4 : rcond(KFF) = 2.07753e-34
% 5 : rcond(KFF) = 0.000163204

% When 5 cross members are placed but not one in each storey:
% rcond(KFF) = 2.0933e-34

% Given that eps is 1e-15, the structure is only stable (ie. greater than
% eps) when 5 cross members are used and for one in each floor.

