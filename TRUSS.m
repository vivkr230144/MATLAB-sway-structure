% CI2-221: Computational Methods II
% Alicia Jiayun LAW
% 01105518

%%
classdef TRUSS
    
    properties

        % Initial input
        NAME              % String 
        NODES             % Struct 
        ELEMENTS          % Struct 
        dofs_free         % Vector 
        dofs_restrained   % Vector        

        % New properties to be defined later
        K
        Fa
        U
        F
        maxlatdisp
        Faxial
        
    end
    
    methods
%%        
        % Create an instance of the TRUSS class
        function obj = TRUSS(NAME, NODES, ELEMENTS, dofs_restrained, dofs_free)
            
            obj.NAME = NAME;
            obj.NODES  = NODES;
            obj.ELEMENTS = ELEMENTS;
            obj.dofs_restrained = dofs_restrained;
            obj.dofs_free = dofs_free;
            
            % 'sub properties' of NODES that we will define later 
            obj.NODES.amp_coords = [];
            obj.NODES.new_coords = [];
            
        end
        
        %% %%%%%%%%%%%%%% PART ONE: PRELIMINARY STEPS %%%%%%%%%%%%%%%%%% %%
        % Adding cross members (ELEMENTS.nodescross) 
        % to the basic structure (ELEMENTS.nodesbasic)
        
        function obj = AddCrossMembers (obj,nodescross)
            
            obj.ELEMENTS.nodes = [obj.ELEMENTS.nodes; nodescross];
            
        end
        
        %% %%%%%%%%%%%%%%%%%%%% PART TWO: ASSEMBLE %%%%%%%%%%%%%%%%%%%%% %%
        % (i)  The Stiffness Matrix, K, and
        % (ii) The Force Vector, fa
        
        function obj = assemble (obj,EA,P,load_dofs)
            % P is the nodal force
            % load_dofs are the dofs that experience direct loading
            
            nodes = size(obj.NODES.coords,1);      % no. of nodes
            elements = size(obj.ELEMENTS.nodes,1); % no. of elements
            
            % Construct
            obj.K = zeros(2*nodes);    % Initialize square matrix the size of the no. of dof
            obj.Fa = zeros(2*nodes,1); % Initialize vector with the length of the no. of dof
            
            for E = 1:elements % loop through all elements & build stiffness matrix
                
                n1 = obj.ELEMENTS.nodes(E,1); n2 = obj.ELEMENTS.nodes(E,2); % identify element node numbers
                x1 = obj.NODES.coords(n1,1); y1 = obj.NODES.coords(n1,2); % element node 1 - x,y coordinates
                x2 = obj.NODES.coords(n2,1); y2 = obj.NODES.coords(n2,2); % element node 2 - x,y coordinates
                dof11 = obj.NODES.dofs(n1,1); dof12 = obj.NODES.dofs(n1,2); % element node 1 - dofs
                dof21 = obj.NODES.dofs(n2,1); dof22 = obj.NODES.dofs(n2,2); % element node 2 - dofs
                alpha = atan2(y2-y1,x2-x1); % angle of inclination relative to the POSITIVE direction of the x axis
                c = cos(alpha); c2 = c*c; s = sin(alpha); s2 = s*s; cs = c*s; % angle parameters
                Le = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ); % element length
                ke = EA/Le; % element axial stiffness
                
                % Updating global stiffness matrix [K] coefficients
                % Note that for each element you still have to 'think locally'
                % Row 1 - element dof11
                obj.K(dof11,dof11) = obj.K(dof11,dof11) + ke*c2; % Col 1 - element dof11
                obj.K(dof11,dof12) = obj.K(dof11,dof12) + ke*cs; % Col 2 - element dof12
                obj.K(dof11,dof21) = obj.K(dof11,dof21) - ke*c2; % Col 3 - element dof21
                obj.K(dof11,dof22) = obj.K(dof11,dof22) - ke*cs; % Col 4 - element dof22
                
                % Row 2 - element dof12
                obj.K(dof12,dof11) = obj.K(dof12,dof11) + ke*cs; % Col 1 - element dof11
                obj.K(dof12,dof12) = obj.K(dof12,dof12) + ke*s2; % Col 2 - element dof12
                obj.K(dof12,dof21) = obj.K(dof12,dof21) - ke*cs; % Col 3 - element dof21
                obj.K(dof12,dof22) = obj.K(dof12,dof22) - ke*s2; % Col 4 - element dof22
                
                % Row 3 - element dof21
                obj.K(dof21,dof11) = obj.K(dof21,dof11) - ke*c2; % Col 1 - element dof11
                obj.K(dof21,dof12) = obj.K(dof21,dof12) - ke*cs; % Col 2 - element dof12
                obj.K(dof21,dof21) = obj.K(dof21,dof21) + ke*c2; % Col 3 - element dof21
                obj.K(dof21,dof22) = obj.K(dof21,dof22) + ke*cs; % Col 4 - element dof22
                
                % Row 4 - element dof22
                obj.K(dof22,dof11) = obj.K(dof22,dof11) - ke*cs; % Col 1 - element dof11
                obj.K(dof22,dof12) = obj.K(dof22,dof12) - ke*s2; % Col 2 - element dof12
                obj.K(dof22,dof21) = obj.K(dof22,dof21) + ke*cs; % Col 3 - element dof21
                obj.K(dof22,dof22) = obj.K(dof22,dof22) + ke*s2; % Col 4 - element dof22
                
                % There are no applied member forces
                
            end
            
            %Updating the global applied nodal force vector {fa} coefficients
            
            obj.Fa(load_dofs) = P;
            
        end
        
        %% %%%%%%%%%%%%%%%%%%%%  PART THREE: SOLVE  %%%%%%%%%%%%%%%%%%%% %%
        % (i) For the nodal displacements, U, and
        % (ii) the nodal reactions, F
        
        function obj = solve (obj)
            
            nodes = size(obj.NODES.coords,1); % no. of nodes
            
            % Specification of submatrices
            KRR = obj.K(obj.dofs_restrained , obj.dofs_restrained);
            KRF = obj.K(obj.dofs_restrained , obj.dofs_free);
            KFR = obj.K(obj.dofs_free , obj.dofs_restrained);
            KFF = obj.K(obj.dofs_free , obj.dofs_free);
            fF = obj.Fa(obj.dofs_free);
            
            % Solution for the unknown nodal dofs
            uR = zeros(size(obj.dofs_restrained)); % BC - zero displacement on restrained nodes
            uF = KFF\(fF - KFR*uR'); % 1st matrix equation
            obj.U = zeros(2*nodes,1); % full nodal dof vector
            obj.U(obj.dofs_free) = uF'; obj.U(obj.dofs_restrained) = uR';
            
            % Solution for the unknown reactions
            fR = KRF*uF + KRR*uR'; % 2nd matrix equation
            obj.F = zeros(2*nodes,1); % full nodal force vector
            obj.F(obj.dofs_free) = fF'; obj.F(obj.dofs_restrained) = fR';
            
        end
        
        %% %%%%%%%%%%%%  PART FOUR: UPDATE NEW COORDINATES  %%%%%%%%%%%% %%
        
        function obj = newcoord (obj,amp)
            % amp = amplification factor for plotting purposes only
            
            % Initialize
            % (i) actual coordinates
            obj.NODES.new_coords = zeros(size(obj.NODES.coords)); 
            
            % (ii) amplified coordinates
            obj.NODES.amp_coords = zeros(size(obj.NODES.coords)); 
            
            for I = 1:size(obj.NODES.coords,1) % x-coordinate
                for J = 1:size(obj.NODES.coords,2) % y-coordinate
                    
                    % (i) Actual Co-ordinates
                    obj.NODES.new_coords(I,J) = obj.NODES.coords(I,J) + obj.U(obj.NODES.dofs(I,J)); 
                    
                    % (ii) Amplified Co-ordinates
                    obj.NODES.amp_coords(I,J) = obj.NODES.coords(I,J) + obj.U(obj.NODES.dofs(I,J))*amp; 
                    
                end
            end
            
        end
        
        %% %%%%%%%%%%%%  PART FIVE: PLOT NEW COORDINATES  %%%%%%%%%%%%%% %%
        function obj = plotting (obj)
            
            % Plotting
            figure; hold all; grid on;
            
            % Axis Limits
            x_min = min(obj.NODES.coords(:,1)); x_max = max(obj.NODES.coords(:,1)); x_range = x_max - x_min; x_ext = 0.1*x_range;
            y_min = min(obj.NODES.coords(:,2)); y_max = max(obj.NODES.coords(:,2)); y_range = y_max - y_min; y_ext = 0.2*y_range;
            axis([x_min-x_ext  x_max+x_ext  y_min-y_ext  y_max+y_ext])
            
            elements = size(obj.ELEMENTS.nodes,1); % no. of elements
            
            for E = 1:elements
                n1 = obj.ELEMENTS.nodes(E,1); n2 = obj.ELEMENTS.nodes(E,2); % identify element node numbers
                
                % Plotting original structure
                x1 = obj.NODES.coords(n1,1); y1 = obj.NODES.coords(n1,2); % element node 1 - x,y original coordinates
                x2 = obj.NODES.coords(n2,1); y2 = obj.NODES.coords(n2,2); % element node 2 - x,y original coordinates
                Le = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ); % element length
                plot([x1,x2],[y1,y2],'k','Linewidth',3);
                
                % Check on changes in member lengths and plotting amplified deformed structure
                x1_new = obj.NODES.new_coords(n1,1); y1_new = obj.NODES.new_coords(n1,2); % element node 1 - x,y actual deformed coordinates
                x2_new = obj.NODES.new_coords(n2,1); y2_new = obj.NODES.new_coords(n2,2); % element node 2 - x,y actual deformed coordinates
                x1_amp = obj.NODES.amp_coords(n1,1); y1_amp = obj.NODES.amp_coords(n1,2); % element node 1 - x,y amplified deformed coordinates
                x2_amp = obj.NODES.amp_coords(n2,1); y2_amp = obj.NODES.amp_coords(n2,2); % element node 2 - x,y amplified deformed coordinates
                Le_new = sqrt( (x2_new - x1_new)^2 + (y2_new - y1_new)^2 ); % actual new element length
                if Le_new <= Le; % element length has decreased - member in compression
                    plot([x1_amp,x2_amp],[y1_amp,y2_amp],'b','Linewidth',3); % blue colour
                elseif Le_new > Le; % element length as increased - member in tension
                    plot([x1_amp,x2_amp],[y1_amp,y2_amp],'r','Linewidth',3); % red colour
                end
                
                % Plotting nodes last
                % Original Nodes
                plot(x1,y1,'ko','Markersize',7,'MarkerFaceColor','w');
                plot(x2,y2,'ko','Markersize',7,'MarkerFaceColor','w');
                
                % New Nodes
                plot(x1_amp,y1_amp,'ko','Markersize',7,'MarkerFaceColor','y');
                plot(x2_amp,y2_amp,'ko','Markersize',7,'MarkerFaceColor','y');
                xlabel('x-coordinate [mm]')
                ylabel('y-coordinate [mm]')
                
            end
            
        end
        
        %% %%%%  PART SIX: OBTAIN THE MAXIMUM LATERAL DISPLACEMENT  %%%% %%
        
        function obj = FindMaxLatDisp (obj)
        
            obj.maxlatdisp = max (obj.U(1:2:end-1));
        
        end
        
        %% %%%%%%%%  PART SEVEN: OBTAIN THE MEMBER AXIAL FORCE  %%%%%%%% %%
        
        function obj = axial (obj,EA)
            
            elements = size(obj.ELEMENTS.nodes,1);  % no. of elements
            obj.Faxial = size(obj.ELEMENTS.nodes,1); % initialize
            
            for E = 1:elements
                
                n1 = obj.ELEMENTS.nodes(E,1); n2 = obj.ELEMENTS.nodes(E,2); % identify element node numbers
                
                % Original structure
                x1 = obj.NODES.coords(n1,1); y1 = obj.NODES.coords(n1,2); % element node 1 - x,y original coordinates
                x2 = obj.NODES.coords(n2,1); y2 = obj.NODES.coords(n2,2); % element node 2 - x,y original coordinates
                dof11 = obj.NODES.dofs(n1,1); dof12 = obj.NODES.dofs(n1,2); % element node 1 - dofs
                dof21 = obj.NODES.dofs(n2,1); dof22 = obj.NODES.dofs(n2,2); % element node 2 - dofs
                alpha = atan2(y2-y1,x2-x1); % angle of inclination relative to the POSITIVE direction of the x axis
                c = cos(alpha); s = sin(alpha); % angle parameters
                Le = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ); % element length
                
                % Obtain displacement
                % for first node
                U11 = obj.U(dof11); U12 = obj.U(dof12); 
                % for second node
                U21 = obj.U(dof21); U22 = obj.U(dof22);
                
                % Transformation matrix
                T = [c s 0 0; 0 0 c s];
                
                % Transfrom gloabal dof to local dofs
                localdof = T*[U11;U12;U21;U22];
                
                % Obtain local displacement
                u1 = localdof(1);
                u2 = localdof(2);
                
                % Obtain epsilon (strain)
                eps = (u2-u1)/Le;
                
                % Obtain axial force
                obj.Faxial (E) = EA*eps ;
            end
            
        end
        
        %% %%%%%%%%%%  PART EIGHT: Plot Elements that buckle  %%%%%%%%%% %%
        
        function obj = buckleplot (obj,memberfail)
            
            % Plotting
            figure; hold all; grid on;
            
            % Axis Limits
            x_min = min(obj.NODES.coords(:,1)); x_max = max(obj.NODES.coords(:,1)); x_range = x_max - x_min; x_ext = 0.1*x_range;
            y_min = min(obj.NODES.coords(:,2)); y_max = max(obj.NODES.coords(:,2)); y_range = y_max - y_min; y_ext = 0.2*y_range;
            axis([x_min-x_ext  x_max+x_ext  y_min-y_ext  y_max+y_ext])
            
            elements = size(obj.ELEMENTS.nodes,1); % no. of elements
            
            for E = 1:elements
                n1 = obj.ELEMENTS.nodes(E,1); n2 = obj.ELEMENTS.nodes(E,2); % identify element node numbers
                
                % Plotting original structure
                x1 = obj.NODES.coords(n1,1); y1 = obj.NODES.coords(n1,2); % element node 1 - x,y original coordinates
                x2 = obj.NODES.coords(n2,1); y2 = obj.NODES.coords(n2,2); % element node 2 - x,y original coordinates
                Le = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ); % element length
                
                plot([x1,x2],[y1,y2],'k','Linewidth',3);
                
                % Plot buckling elements, override plot
                for i = 1:length(memberfail)
                    
                    % for members that fail, use red
                    if E == memberfail(i)
                        plot([x1,x2],[y1,y2],'r','Linewidth',3);
                        break
                    end
                    
                end
                
                % Plotting original nodes
                plot(x1,y1,'ko','Markersize',7,'MarkerFaceColor','w');
                plot(x2,y2,'ko','Markersize',7,'MarkerFaceColor','w');
                
            end
            
            xlabel('x coordinate [mm]','FontSize',16);
            ylabel('y coordinate [mm]','FontSize',16);
            title(['Buckling Members'],'FontSize',16);
            
        end
        
    end
end
