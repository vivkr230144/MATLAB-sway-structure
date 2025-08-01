clear
clc

% CI2-221: Computational Methods II
% Coursework - Modeling of a Sway Structure
% Name: Alicia Jiayun LAW

%% 
% Assigned Structure
alpha = mod(12,5) + 3 ;
beta = mod(15,5) + 3 ;

N = min(alpha,beta); % Number of bays
M = max(alpha,beta); % Number of Storeys
n_nodes = (N+1)*(M+1); % Number of Nodes
elements = 35; % Number of Elements w/o cross members

%% Question One
% Parameters - all in N and mm for simplicity

L = 1000; % mm - length of horizontal and vertical members
radius = 25; % mm - radius of each member
E = 2e5; % N/mm2 - Youngs Modulus
P = 100e3; % N - nodal loading

% Derived Parameters
A = pi*radius^2; % mm2 - cross-sectional area of each member
EA = E*A; % N - axial rigidity of all members


% SPECIFYING NODES AND ELEMENTS
% Initialize
NODES.coords = zeros (n_nodes, 2);
NODES.dofs = zeros (n_nodes, 2);

for i = 1 : M+1     % loop over rows
    for j = 1 : N+1 % loop over columns
        
        node_number = 6*(j-1)+i; % Current node number being specified
        
        % (i)   Specifying nodal x-y coordinates
        NODES.coords (node_number,:) = [L*(j-1) , L*(i-1)];
        
        % (ii)  Specifying nodal dofs
        % Only 2 dof per node for trusses
        NODES.dofs (node_number,:) = [node_number*2 - 1 , node_number*2];
        
    end
end

% (iii) Specifying Element-Nodal connectivity
% Initialize
ELEMENTS.nodes = zeros (elements,2);
ELEMENTS.nodes (1,:) = [1 2]; % The nodes for the first element

for e = 2: elements
    
    % (for most vertical bars, except those starting from a new bay:)
    % The first node of my current element is the last node of my previous element
    % The second node of my current element is +1 of my first node
    ELEMENTS.nodes (e,:) = [ELEMENTS.nodes(e-1,2),ELEMENTS.nodes(e-1,2)+1];
    
    % For vertical elements
    if e < 20 && mod(e-1,5) == 0 % Elements in a New Bay
        
        ELEMENTS.nodes (e,:) = [1+6*(e-1)/5,2+6*(e-1)/5];
        
    % For horizontal elements
    elseif e > 20 && e <= 35
        
        if mod(e-21,3)==0 % Elements in a New Row
            
            ELEMENTS.nodes (e,:) = [2+(e-21)/3,2+(e-21)/3+6];
            
        else % all other horizontal elements
            
            ELEMENTS.nodes (e,:) = [ELEMENTS.nodes(e-1,2),ELEMENTS.nodes(e-1,2)+6];
            
        end
        
    else
        
    end
end

% Degrees of freedom & other parameters for trusses
dofs_free = [3:12, 15:24, 27:36, 39:48];
dofs_restrained = [1,2,13,14,25,26,37,38]; 

% Assemble the stiffness matrix, K, and force vector, fa for trusses
load_dofs = [3, 5, 7, 9, 11];

obj = TRUSS('5x3', NODES, ELEMENTS, dofs_restrained, dofs_free);
obj = obj.assemble (EA,P,load_dofs);
obj = obj.solve ();

% Extract KFF for assessment
KFF = obj.K(dofs_free, dofs_free);

% Assessment of the stability of the truss can be found in code 'Qn1_rcond'

%% Questions Two and Three

% Permutate every combination using nested loops
% Initial values 
maxdisp_trussmax = 0;
maxdisp_trussmin = 1000000;

% N=3 and each nested loop has 6 possible variations.
% - for N = 1 to 3, the cross members are inclined to the right
% - for N = 4 to 6, the cross members are inclined to the left

for M1 = 1 : 2*N % Floor 1
    for M2 = 1 : 2*N % Floor 2
        for M3 = 1 : 2*N % Floor 3
            for M4 = 1 : 2*N % Floor 4
                for M5 = 1 : 2*N % Floor 5

                    obj2 = TRUSS('5x3 cross', NODES, ELEMENTS, dofs_restrained, dofs_free);
                    
                    % Obtain the element-nodal connectivity for cross members
                    if M1 <= N
                        nodescross(1,1) = 6*(M1-1)+1; nodescross(1,2) = nodescross(1,1) + 7;
                    elseif M1 > N
                        nodescross(1,1) = 6*(M1-3)-4; nodescross(1,2)= nodescross(1,1) + 5;
                    end

                    
                    if M2 <= N
                        nodescross(2,1) = 2+6*(M2-1);  nodescross(2,2) = nodescross(2,1) + 7;
                    elseif M2 > N
                        nodescross(2,1) = 6*(M2-3)-3;  nodescross(2,2)= nodescross(2,1) + 5;
                    end

                    
                    if M3 <= N
                        nodescross(3,1) = 3+6*(M3-1); nodescross(3,2) = nodescross(3,1) + 7;
                    elseif M3 > N
                        nodescross(3,1) = 6*(M3-3)-2; nodescross(3,2)= nodescross(3,1) + 5;
                    end

                    
                    if M4 <= N
                        nodescross(4,1) = 4+6*(M4-1); nodescross(4,2) = nodescross(4,1) + 7;
                    elseif M4 > N
                        nodescross(4,1) = 6*(M4-3)-1; nodescross(4,2) = nodescross(4,1) + 5;
                    end

                    
                    if M5 <= N
                        nodescross(5,1) = 5+6*(M5-1); nodescross(5,2) = nodescross(5,1) + 7;
                    elseif M5 > N
                        nodescross(5,1) = 6*(M5-3); nodescross(5,2)= nodescross(5,1) + 5;
                    end

                    % Add on the cross member configuration to the
                    % element-nodal connectivity matrix
                    obj2 = obj2.AddCrossMembers(nodescross);
                    
                    % obtain K and Fa
                    obj2 = obj2.assemble (EA,P,load_dofs); 
                    
                    % obtain U and F
                    obj2 = obj2.solve (); 
                    
                    % obtain the maximum lateral displacement for this combination
                    obj2 = obj2.FindMaxLatDisp (); 

                    % Loop over all permutations to determine the

                    if obj2.maxlatdisp > maxdisp_trussmax
                        
                        % (i) configuration for max lateral displacement
                        maxdisp_trussmax = obj2.maxlatdisp;
                        combimax = nodescross;
                        
                    end

                    if obj2.maxlatdisp < maxdisp_trussmin
                        
                        % (ii) configuration for min lateral displacement
                        maxdisp_trussmin = obj2.maxlatdisp;
                        combimin = nodescross;
                        
                    end

                end
            end
        end
    end
end

%% Question Two
% Plot the configuration that provides minimum lateral sway
amp = 1;
obj3 = TRUSS('5x3 min configuration', NODES, ELEMENTS, dofs_restrained, dofs_free);
obj3 = obj3.AddCrossMembers(combimin);
obj3 = obj3.assemble (EA, P, load_dofs);
obj3 = obj3.solve ();
obj3 = obj3.newcoord (amp);
obj3 = obj3.plotting ();

% Check the reaction forces (for report)
reactions_minsway = obj3.F(dofs_restrained);

%% Question Three
% Plot the configuration that provides maximum lateral sway
amp = 1;
obj4 = TRUSS('5x3 max configuration', NODES, ELEMENTS, dofs_restrained, dofs_free);
obj4 = obj4.AddCrossMembers(combimax);
obj4 = obj4.assemble (EA, P, load_dofs);
obj4 = obj4.solve ();
obj4 = obj4.newcoord (amp);
obj4 = obj4.plotting ();

% Check the reaction forces (for report)
reactions_maxsway = obj4.F(dofs_restrained);

%% Question Four

% Additional Parameters for Frames
I = (1/4)*pi*radius^4 ; % mm4 - second moment of area for a circle
EI = E*I ; % Nmm2 - flexural rigidity of all members

% Specifying new nodal dofs for frames
NODES.fdofs = zeros (n_nodes, 3); % Initialize

for i = 1 : M+1     % loop over rows
    for j = 1 : N+1 % loop over columns
        
        node_number = 6*(j-1)+i; % Current node number being specified
        % Now there're 3 dof for frames
        NODES.fdofs (node_number,:) = [node_number*3 - 2 , node_number*3 - 1, node_number*3];
   
    end
end

% Specifying new nodal load dof for frames
load_fdofs = [4, 7, 10, 13, 16];

% Specifying new free and restrained degrees of freedom & other parameters
% for frames
fdofs_free = [4:18, 22:36, 40:54, 58:72];
fdofs_restrained = [1,2,3,19,20,21,37,38,39,55,56,57]; % known nodal x-y dofs due to fixed ends at nodes 1, 7, 13 and 19

% Plot the deformed shape for the frame
amp = 1;
obj5 = FRAMES ('5x3 frames', NODES, ELEMENTS, fdofs_restrained, fdofs_free);
obj5 = obj5.assemble (EA,EI,P,load_fdofs); % obtain K and Fa
obj5 = obj5.solve (); % obtain D and F
obj5 = obj5.FindMaxLatDisp (); % obtain the maximum lateral displacement for this combination
obj5 = obj5.newcoord (amp);
obj5 = obj5.plotting (amp);

% Check the reaction force (for report)
reactions_frames = obj5.F(fdofs_restrained);

% Check maximum lateral displacement
maxdisp_frame = obj5.maxlatdisp;

%% Question Five
fun = @Radius ;

% Point of Interest
% For frames:
fdof_topload = 16 ; 
% For trusses:
dof_topload = 11 ; 

% Find the radius where lateral sway = 25mm
% (i) For Frames: 
x0_frame = 25 ; 
r25_frame = fzero(@(radius) fun (1, radius, E, NODES, ELEMENTS, P, fdofs_restrained, fdofs_free, fdof_topload, load_fdofs, ...
                            dofs_restrained, dofs_free, dof_topload, load_dofs, combimin, combimax),x0_frame);     

% (ii) For Most Efficient Truss (Min Sway)
% Current lateral sway at the point
x0_trussmin = 25 ; 
r25_trussmin = fzero(@(radius) fun (2, radius, E, NODES, ELEMENTS, P,fdofs_restrained, fdofs_free, fdof_topload, load_fdofs, ...
                     dofs_restrained, dofs_free, dof_topload, load_dofs, combimin, combimax),x0_trussmin);     

% (iii) For Least Efficient Truss (Max Sway)
x0_trussmax  = 25 ;
r25_trussmax = fzero(@(radius) fun (3, radius, E, NODES, ELEMENTS, P, fdofs_restrained, fdofs_free, fdof_topload, load_fdofs, ...
                            dofs_restrained, dofs_free, dof_topload, load_dofs, combimin, combimax),x0_trussmax ); 

% Find the volume of steel
Lcross = sqrt(2)*1e3; % mm
vol_frame = (pi*r25_frame^2*L)*35 ; 
vol_trussmin = (pi*r25_trussmin^2*L)*35 + (pi*r25_trussmin^2*Lcross)*5 ;
vol_trussmax = (pi*r25_trussmax^2*L)*35 + (pi*r25_trussmax^2*Lcross)*5 ;

%% Question Six

% (i) For Frames:
% Derived arameters
A_frame = pi*r25_frame^2;
EA_frame = E*A_frame; 
I_frame =(1/4)*pi*r25_frame^4 ;
EI_frame = E*I_frame ;   

% Call for a new class frame
amp = 1;
obj6 = FRAMES ('5x3 FRAME Q6', NODES, ELEMENTS, fdofs_restrained, fdofs_free);
obj6 = obj6.assemble (EA_frame,EI_frame,P,load_fdofs); % obtain K and Fa
obj6 = obj6.solve (); % obtain D and F
obj6 = obj6.newcoord (amp); % obtain new coordinates
obj6 = obj6.axial (EA_frame);

memberfail_frame = []; % stores a vector of members that exceed the critical load
MoSmin_frame = 100000; 
for i = 1:elements

    Lcrit = L/2;
    Pcrit = EI_frame*pi^2/Lcrit^2 ;
    
    % Check if it's less than the critical load
    if -obj6.Faxial(i) > Pcrit % Because we are only interested with members in compression
        
        memberfail_frame(end+1)=i; % if members have failed, store the member
        
    else
        
    end
    
    % Find the minimum margin of safety for compression members
    if obj6.Faxial(i) < 0 % identify compression member
        MoS_frame = Pcrit/abs(obj6.Faxial(i)) - 1; % calculate its MoS
        if MoS_frame < MoSmin_frame
            MoSmin_frame = MoS_frame ; % the smallest MoS of the structure will be the MoS that defines the structure
        end
    end
    
end

% (ii) For Most Efficient Truss (Min Sway)
% Derived arameters
A_trussmin = pi*r25_trussmin^2;
EA_trussmin = E*A_trussmin; 
I_trussmin =(1/4)*pi*r25_trussmin^4 ;
EI_trussmin = E*I_trussmin ;   

% Call for a new class truss
amp = 1;
obj7 = TRUSS ('5x3 TRUSS MIN Q6', NODES, ELEMENTS, dofs_restrained, dofs_free);
obj7 = obj7.AddCrossMembers(combimin);
obj7 = obj7.assemble (EA_trussmin,P,load_dofs); % obtain K and Fa
obj7 = obj7.solve (); % obtain D and F
obj7 = obj7.newcoord (amp); % obtain new coordinates
obj7 = obj7.axial (EA_trussmin);

memberfail_trussmin = []; %stores a vector of elements that exceeded the critical load
MoSmin_trussmin = 100000; 
for i = 1:40
    
    if i > 35 % cross members 
        Lcrit = Lcross ;
    else % vertical and horizontal members
        Lcrit = L ;
    end
    
    Pcrit = EI_trussmin*pi^2/Lcrit^2 ;
    
    if -obj7.Faxial(i) > Pcrit % Because we are only interested with members in compression
        memberfail_trussmin(end+1) = i; % if members have failed, store the member
    else
    end
    
    % Find the minimum margin of safety for compression members
    if obj7.Faxial(i) < 0 % identify compression member
        MoS_trussmin = Pcrit/abs(obj7.Faxial(i)) - 1 ; % calculate its MoS
        if MoS_trussmin < MoSmin_trussmin
            MoSmin_trussmin = MoS_trussmin ; % the smallest MoS of the structure will be the MoS that defines the structure
        end
    end
    
end

% (iii) For Least Efficient Truss (Max Sway)
% Derived arameters
A_trussmax = pi*r25_trussmax^2;
EA_trussmax = E*A_trussmax; 
I_trussmax =(1/4)*pi*r25_trussmax^4 ;
EI_trussmax = E*I_trussmax ;   

% Call for a new class truss
amp = 1;
obj8 = TRUSS ('5x3 TRUSS MAX Q6', NODES, ELEMENTS, dofs_restrained, dofs_free);
obj8 = obj8.AddCrossMembers(combimax);
obj8 = obj8.assemble (EA_trussmax,P,load_dofs); % obtain K and Fa
obj8 = obj8.solve (); % obtain D and F
obj8 = obj8.newcoord (amp); % obtain new coordinates
obj8 = obj8.axial (EA_trussmax);

% Check if it's less than the critical load
memberfail_trussmax = []; %stores a vector of elements that exceeded the critical load
MoSmin_trussmax = 100000; 
for i = 1:40
    
    if i > 35 % cross members 
        Lcrit = Lcross ;
    else % vertical and horizontal members
        Lcrit = L ;
    end
    
    Pcrit = EI_trussmax*pi^2/Lcrit^2 ;
    
    if -obj8.Faxial(i) > Pcrit % Because we are only interested with members in compression
        memberfail_trussmax(end+1) = i; % if members have failed, store the member
    else
        
    end
    
    % Find the minimum margin of safety for compression members
    if obj8.Faxial(i) < 0 % identify compression member
        MoS_trussmax = Pcrit/abs(obj8.Faxial(i)) - 1; % calculate its MoS
        if MoS_trussmax < MoSmin_trussmax
            MoSmin_trussmax = MoS_trussmax ; % the smallest MoS of the structure will be the MoS that defines the structure
        end
    end
    
end

fprintf('The margin of safety is: frames = %g, truss min sway = %g, truss max sway = %g',MoSmin_frame,MoSmin_trussmin,MoSmin_trussmax)

% Plot the elements that failed (for report)
obj7 = obj7.buckleplot(memberfail_trussmin);
