% CI2-221: Computational Methods II
% Alicia Jiayun LAW
% 01105518

%%
classdef FRAMES
    
    
    properties
        
        % Initial input
        NAME              % String 
        NODES             % Struct 
        ELEMENTS          % Struct 
        fdofs_free        % Vector 
        fdofs_restrained  % Vector
        
        % New properties to be defined later
        K   
        fa  
        D   
        F
        maxlatdisp
        Faxial
        
    end
    
    methods
        
        % Create an instance of the TRUSS class
        function obj = FRAMES(NAME, NODES, ELEMENTS, fdofs_restrained, fdofs_free)
            
            obj.NAME = NAME;
            obj.NODES  = NODES;
            obj.ELEMENTS = ELEMENTS;
            obj.fdofs_restrained = fdofs_restrained;
            obj.fdofs_free = fdofs_free;
            
            % 'sub properties' of NODES that we will define later 
            obj.NODES.amp_coords = [];
            obj.NODES.new_coords = [];
            
        end
        
        %% %%%%%%%%%%%%%%%%%%%% PART ONE: ASSEMBLE %%%%%%%%%%%%%%%%%%%%% %%
        % (i)  The Stiffness Matrix, K, and
        % (ii) The Force Vector, fa
        
        function obj = assemble (obj,EA,EI,P,load_fdofs)
            % P is the nodal force
            % load_dofs are the dofs that experience direct loading
            
            nodes = size(obj.NODES.coords,1);      % no. of nodes
            elements = size(obj.ELEMENTS.nodes,1); % no. of elements
            
            % Construct
            obj.K = zeros(3*nodes);    % Initialize square matrix the size of the no. of dof
            obj.fa = zeros(3*nodes,1); % Initialize vector with the length of the no. of dof
            
            for E = 1:elements % loop through all elements & build stiffness matrix
                
                n1 = obj.ELEMENTS.nodes(E,1); n2 = obj.ELEMENTS.nodes(E,2); % identify element node numbers
                x1 = obj.NODES.coords(n1,1); y1 = obj.NODES.coords(n1,2); % element node 1 - x,y coordinates
                x2 = obj.NODES.coords(n2,1); y2 = obj.NODES.coords(n2,2); % element node 2 - x,y coordinates
                dof11 = obj.NODES.fdofs(n1,1); dof12 = obj.NODES.fdofs(n1,2); dof13 = obj.NODES.fdofs(n1,3); % element node 1 - dofs
                dof21 = obj.NODES.fdofs(n2,1); dof22 = obj.NODES.fdofs(n2,2); dof23 = obj.NODES.fdofs(n2,3); % element node 2 - dofs
                alpha = atan2(y2-y1,x2-x1); % angle of inclination of element to POSITIVE direction of the horizontal x axis
                c = cos(alpha); c2 = c*c; s = sin(alpha); s2 = s*s; cs = c*s; % angle parameters
                Le = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ); % element length
                ka = EA/Le; % generic stretching stiffness term
                kb = EI/(Le*Le*Le); % generic bending stiffness term
                
                % Updating global stiffness matrix [K] coefficients
                % Note that for each element you still have to 'think locally'
                % Row 1 - element dof11 - u1 of element
                obj.K(dof11,dof11) = obj.K(dof11,dof11) + ka*c2 + 12*kb*s2;      % Col 1 - element dof11 - u1 of element
                obj.K(dof11,dof12) = obj.K(dof11,dof12) + ka*cs - 12*kb*cs;      % Col 2 - element dof12 - v1 of element
                obj.K(dof11,dof13) = obj.K(dof11,dof13) - 6*Le*kb*s;             % Col 3 - element dof13 - th1 of element
                obj.K(dof11,dof21) = obj.K(dof11,dof21) - ka*c2 - 12*kb*s2;      % Col 4 - element dof21 - u2 of element
                obj.K(dof11,dof22) = obj.K(dof11,dof22) - ka*cs + 12*kb*cs;      % Col 5 - element dof22 - v2 of element
                obj.K(dof11,dof23) = obj.K(dof11,dof23) - 6*Le*kb*s;             % Col 6 - element dof23 - th2 of element
                
                % Row 2 - element dof12 - v1 of element
                obj.K(dof12,dof11) = obj.K(dof12,dof11) + ka*cs - 12*kb*cs;      % Col 1 - element dof11 - u1 of element
                obj.K(dof12,dof12) = obj.K(dof12,dof12) + ka*s2 + 12*kb*c2;      % Col 2 - element dof12 - v1 of element
                obj.K(dof12,dof13) = obj.K(dof12,dof13) + 6*Le*kb*c;             % Col 3 - element dof13 - th1 of element
                obj.K(dof12,dof21) = obj.K(dof12,dof21) - ka*cs + 12*kb*cs;      % Col 4 - element dof21 - u2 of element
                obj.K(dof12,dof22) = obj.K(dof12,dof22) - ka*s2 - 12*kb*c2;      % Col 5 - element dof22 - v2 of element
                obj.K(dof12,dof23) = obj.K(dof12,dof23) + 6*Le*kb*c;             % Col 6 - element dof23 - th2 of element
                
                % Row 3 - element dof13 - th1 of element
                obj.K(dof13,dof11) = obj.K(dof13,dof11) - 6*Le*kb*s;             % Col 1 - element dof11 - u1 of element
                obj.K(dof13,dof12) = obj.K(dof13,dof12) + 6*Le*kb*c;             % Col 2 - element dof12 - v1 of element
                obj.K(dof13,dof13) = obj.K(dof13,dof13) + 4*Le*Le*kb;            % Col 3 - element dof13 - th1 of element
                obj.K(dof13,dof21) = obj.K(dof13,dof21) + 6*Le*kb*s;             % Col 4 - element dof21 - u2 of element
                obj.K(dof13,dof22) = obj.K(dof13,dof22) - 6*Le*kb*c;             % Col 5 - element dof22 - v2 of element
                obj.K(dof13,dof23) = obj.K(dof13,dof23) + 2*Le*Le*kb;            % Col 6 - element dof23 - th2 of element
                
                % Row 4 - element dof21 - u2 of element
                obj.K(dof21,dof11) = obj.K(dof21,dof11) - ka*c2 - 12*kb*s2;      % Col 1 - element dof11 - u1 of element
                obj.K(dof21,dof12) = obj.K(dof21,dof12) - ka*cs + 12*kb*cs;      % Col 2 - element dof12 - v1 of element
                obj.K(dof21,dof13) = obj.K(dof21,dof13) + 6*Le*kb*s;             % Col 3 - element dof13 - th1 of element
                obj.K(dof21,dof21) = obj.K(dof21,dof21) + ka*c2 + 12*kb*s2;      % Col 4 - element dof21 - u2 of element
                obj.K(dof21,dof22) = obj.K(dof21,dof22) + ka*cs - 12*kb*cs;      % Col 5 - element dof22 - v2 of element
                obj.K(dof21,dof23) = obj.K(dof21,dof23) + 6*Le*kb*s;             % Col 6 - element dof23 - th2 of element
                
                % Row 5 - element dof22 - v2 of element
                obj.K(dof22,dof11) = obj.K(dof22,dof11) - ka*cs + 12*kb*cs;      % Col 1 - element dof11 - u1 of element
                obj.K(dof22,dof12) = obj.K(dof22,dof12) - ka*s2 - 12*kb*c2;      % Col 2 - element dof12 - v1 of element
                obj.K(dof22,dof13) = obj.K(dof22,dof13) - 6*Le*kb*c;             % Col 3 - element dof13 - th1 of element
                obj.K(dof22,dof21) = obj.K(dof22,dof21) + ka*cs - 12*kb*cs;      % Col 4 - element dof21 - u2 of element
                obj.K(dof22,dof22) = obj.K(dof22,dof22) + ka*s2 + 12*kb*c2;      % Col 5 - element dof22 - v2 of element
                obj.K(dof22,dof23) = obj.K(dof22,dof23) - 6*Le*kb*c;             % Col 6 - element dof23 - th2 of element
                
                % Row 6 - element dof23 - th2 of element
                obj.K(dof23,dof11) = obj.K(dof23,dof11) - 6*Le*kb*s;             % Col 1 - element dof11 - u1 of element
                obj.K(dof23,dof12) = obj.K(dof23,dof12) + 6*Le*kb*c;             % Col 2 - element dof12 - v1 of element
                obj.K(dof23,dof13) = obj.K(dof23,dof13) + 2*Le*Le*kb;            % Col 3 - element dof13 - th1 of element
                obj.K(dof23,dof21) = obj.K(dof23,dof21) + 6*Le*kb*s;             % Col 4 - element dof21 - u2 of element
                obj.K(dof23,dof22) = obj.K(dof23,dof22) - 6*Le*kb*c;             % Col 5 - element dof22 - v2 of element
                obj.K(dof23,dof23) = obj.K(dof23,dof23) + 4*Le*Le*kb;            % Col 6 - element dof23 - th2 of element
                
                % There are no applied member forces
                
            end
            
            %Updating the global applied nodal force vector {fa} coefficients
            
            obj.fa(load_fdofs) = P;
            
        end
        
        %% %%%%%%%%%%%%%%%%%%%%  PART TWO: SOLVE  %%%%%%%%%%%%%%%%%%%%%% %%
        % (i) For the nodal displacements, U, and
        % (ii) the nodal reactions, F
        
        function obj = solve (obj)
            
            nodes = size(obj.NODES.coords,1); % no. of nodes
            
            % Specification of submatrices
            KRR = obj.K(obj.fdofs_restrained , obj.fdofs_restrained);
            KRF = obj.K(obj.fdofs_restrained , obj.fdofs_free);
            KFR = obj.K(obj.fdofs_free , obj.fdofs_restrained);
            KFF = obj.K(obj.fdofs_free , obj.fdofs_free);
            faR = obj.fa(obj.fdofs_restrained);
            faF = obj.fa(obj.fdofs_free);
            frF = zeros(size(faF));
            
            % Solution for the unknown nodal dofs
            dR = zeros(size(obj.fdofs_restrained))'; % BC - zero displacement on restrained nodes
            dF = KFF\(frF + faF - KFR*dR); % 1st matrix equation
            obj.D = zeros(size(obj.K,1),1); obj.D(obj.fdofs_free) = dF'; obj.D(obj.fdofs_restrained) = dR'; % full nodal dof vector
            
            % Solution for the unknown reactions
            frR = KRF*dF + KRR*dR - faR; % 2nd matrix equation
            obj.F = zeros(size(obj.K,1),1); obj.F(obj.fdofs_free) = frF' + faF'; obj.F(obj.fdofs_restrained) = frR' + faR'; % full reconstituted nodal force vector
            
        end
        
        %% %%%%%%%%%%%%  PART THREE: UPDATE NEW COORDINATES %%%%%%%%%%%% %%
        
        function obj = newcoord (obj,amp)
            % amp = amplification factor for plotting purposes only
            
            for I = 1:size(obj.NODES.coords,1)
                
                % Amplified Coordinates
                % Nodal x-coordinate gets updated with the axial deformation u dof
                obj.NODES.amp_coords(I,1) = obj.NODES.coords(I,1) + obj.D(3*I-2)*amp;
                % Nodal y-coordinate gets updated with the transverse deformation v dof
                obj.NODES.amp_coords(I,2) = obj.NODES.coords(I,2) + obj.D(3*I-1)*amp;
                
                % New Original Coordinates
                % Nodal x-coordinate gets updated with the axial deformation u dof
                obj.NODES.new_coords(I,1) = obj.NODES.coords(I,1) + obj.D(3*I-2);
                % Nodal y-coordinate gets updated with the transverse deformation v dof
                obj.NODES.new_coords(I,2) = obj.NODES.coords(I,2) + obj.D(3*I-1);
                
            end
        end
        
        %% %%%%%%%%%%%%%  PART FOUR: PLOT NEW COORDINATES  %%%%%%%%%%%%% %%
        function obj = plotting (obj,amp)
            
            figure; hold all; grid on;
            elements = size(obj.ELEMENTS.nodes,1); % no. of elements
            
            for E = 1:elements
                n1 = obj.ELEMENTS.nodes(E,1); n2 = obj.ELEMENTS.nodes(E,2); % identify element node numbers
                
                % Plotting original structure
                x1 = obj.NODES.coords(n1,1); y1 = obj.NODES.coords(n1,2); % element node 1 - x,y original coordinates
                x2 = obj.NODES.coords(n2,1); y2 = obj.NODES.coords(n2,2); % element node 2 - x,y original coordinates
                Le = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ); % element length
                alpha = atan2(y2 - y1,x2 - x1); % angle of inclination of element to the POSITIVE direction of the horizontal x axis
                plot([x1,x2],[y1,y2],'k','Linewidth',3);
                
                % Plotting the deformed structure
                x1_amp = obj.NODES.amp_coords(n1,1); y1_amp = obj.NODES.amp_coords(n1,2); % element node 1 - x,y amplified deformed coordinates
                x2_amp = obj.NODES.amp_coords(n2,1); y2_amp = obj.NODES.amp_coords(n2,2); % element node 2 - x,y amplified deformed coordinates
                % This is a little harder this time around, because the shape functions
                % are no longer linear but cubic. Recall from the lectures that:
                % the transverse deflection is v(Xi) = N1(Xi)v1 + N2(Xi)th1 + N3(Xi)v2 + N4(Xi)th2
                % the axial deflection would be u(Xi) = N5(Xi)u1 + N6(Xi)u2 where N5 and N6 are linear in Xi
                Xi = [-1:0.1:+1]; % The intrinsic Xi coordinate is always between -1 and +1
                
                % The following are the shape functions and 2nd derivatives evaluated at these Xi
                N1 = (1/4)*(1-Xi).*(1-Xi).*(2+Xi);  ddN1ddXi = (3/2)*Xi;
                N2 = (Le/8)*(1-Xi).*(1-Xi).*(1+Xi); ddN2ddXi = (1/4*(1+Xi))*Le-(1/2)*Le*(1-Xi);
                N3 = (1/4)*(1+Xi).*(1+Xi).*(2-Xi);  ddN3ddXi = -(3/2)*Xi;
                N4 = (Le/8)*(1+Xi).*(1+Xi).*(Xi-1); ddN4ddXi = (1/4)*Le*(-1+Xi)+(1/2*(1+Xi))*Le;
                N5 = (1-Xi)/2;
                N6 = (1+Xi)/2;
                
                % The following are the computed global nodal dofs
                ux1 = obj.D(obj.NODES.fdofs(n1,1)); uy1 = obj.D(obj.NODES.fdofs(n1,2)); thz1 = obj.D(obj.NODES.fdofs(n1,3)); % element node 1 - global dofs
                ux2 = obj.D(obj.NODES.fdofs(n2,1)); uy2 = obj.D(obj.NODES.fdofs(n2,2)); thz2 = obj.D(obj.NODES.fdofs(n2,3)); % element node 2 - global dofs
                % Transformation from global to local nodal dofs
                c = cos(alpha); s = sin(alpha);
                u1 = c*ux1 + s*uy1; v1 = -s*ux1 + c*uy1; th1 = thz1; % element node 1 - local dofs
                u2 = c*ux2 + s*uy2; v2 = -s*ux2 + c*uy2; th2 = thz2; % element node 1 - local dofs
                
                % The following are the v(Xi) transverse and u(Xi) axial variations within the element
                vXi = N1*v1 + N2*th1 + N3*v2 + N4*th2;
                uXi = N5*u1 + N6*u2;
                % The following is the curvature d2v/dXi2
                vppXi = ddN1ddXi*v1 + ddN2ddXi*th1 + ddN3ddXi*v2 + ddN4ddXi*th2;
                
                % The following is a transformation of the local coordinate Xi into a
                % local element coordinate equal to 0 at its node 1 and Le_new at its node 2.
                % This is necessary because the structure is now 2D. Call this coord 's'
                s = 0.5*Le*(1+Xi);
                x = x1 + s*cos(alpha) - vXi*sin(alpha)*amp + uXi*cos(alpha)*amp;
                y = y1 + s*sin(alpha) + vXi*cos(alpha)*amp + uXi*sin(alpha)*amp;
                
                % The tricky part - although we have just *computed* v in terms of Xi,
                % we will actually plot those same values in terms of the global x coordinate!
                for I = 2:length(x);
                    % The following checks whether the element length has increased
                    % (TENSION - red) or decreased (COMPRESSION - blue)
                    tol = 1e-6;
                    if u2 - u1 > tol; colourA = 'ro-'; else colourA = 'bo-'; end
                    if abs(u2 - u1) < tol; colourA = 'ko-'; end
                    
                    % The aim of computing the 2nd derivative of v(Xi) is to obtain
                    % the curvature. The bending moment is directly proportional to the
                    % curvature, so that if we know that the 2nd derivative is
                    % positive, the bending moment is also positive (SAGGING) and if
                    % the 2nd derivative is negative, the bending moment is also
                    % negative (HOGGING) - relative to how we define nodes 1 & 2
                    if vppXi(I) >= 0; colourB = 'r'; else colourB = 'b'; end
                    
                    % This ensures that the icon is plotted red for sagging and blue for hogging.
                    % Not only that, the line connecting the icons is red for axial tension and blue
                    % for axial compression. Neat!
                    plot([x(I) x(I-1)],[y(I) y(I-1)],colourA,'Linewidth',3,'Markersize',12,'MarkerFaceColor',colourB);
                end
                % This is much easier than having to change each N from Xi to x and
                % computing v in terms of x
                
                % Plotting nodes last!
                plot(x1,y1,'ko','Markersize',7,'MarkerFaceColor','w');
                plot(x2,y2,'ko','Markersize',7,'MarkerFaceColor','w');
                plot(x1_amp,y1_amp,'ko','Markersize',5,'MarkerFaceColor','y');
                plot(x2_amp,y2_amp,'ko','Markersize',5,'MarkerFaceColor','y');
                
            end
            
            xlabel('x coordinate','FontSize',16);
            ylabel('y coordinate','FontSize',16);
            title(['2D frame under horizontal point loads: deformed shape amplification factor = ',num2str(amp)],'FontSize',16);
            
        end
        %% %%%%% PART FIVE: OBTAIN THE MAXIMUM LATERAL DISPLACEMENT %%%% %%
        
        function obj = FindMaxLatDisp (obj)
            
            obj.maxlatdisp = max (obj.D(1:3:end-2));
            
        end
        
        %% %%%%%%%%%  PART SIX: OBTAIN THE MEMBER AXIAL FORCE %%%%%%%%%% %%
        
        function obj = axial (obj,EA)
            
            elements = size(obj.ELEMENTS.nodes,1);   % no. of elements
            obj.Faxial = size(obj.ELEMENTS.nodes,1); % initialize
            
            for E = 1:elements
                n1 = obj.ELEMENTS.nodes(E,1); n2 = obj.ELEMENTS.nodes(E,2); % identify element node numbers
                
                % Original structure
                x1 = obj.NODES.coords(n1,1); y1 = obj.NODES.coords(n1,2); % element node 1 - x,y original coordinates
                x2 = obj.NODES.coords(n2,1); y2 = obj.NODES.coords(n2,2); % element node 2 - x,y original coordinates
                dof11 = obj.NODES.fdofs(n1,1); dof12 = obj.NODES.fdofs(n1,2); dof13 = obj.NODES.fdofs(n1,3); % element node 1 - dofs
                dof21 = obj.NODES.fdofs(n2,1); dof22 = obj.NODES.fdofs(n2,2); dof23 = obj.NODES.fdofs(n2,3); % element node 2 - dofs
                alpha = atan2(y2-y1,x2-x1); % angle of inclination of element to POSITIVE direction of the horizontal x axis
                c = cos(alpha); s = sin(alpha); % angle parameters
                Le = sqrt( (x2 - x1)^2 + (y2 - y1)^2 ); % element length
                
                % Obtain displacements
                % for first node
                D11 = obj.D(dof11); D12 = obj.D(dof12); D13 = obj.D(dof13);
                % for second node
                D21 = obj.D(dof21); D22 = obj.D(dof22); D23 = obj.D(dof23);
                
                % Transformation matrix
                T = [c s 0 0 0 0 ; 
                    -s c 0 0 0 0 ; 
                     0 0 1 0 0 0 ;
                     0 0 0 c s 0 ; 
                     0 0 0 c s 0 ;
                     0 0 0 -s c 0; 
                     0 0 0 0 0 1]; 
                 
                % Transfrom gloabal dof to local dofs
                localdof = T*[D11;D12;D13;D21;D22;D23];
                
                % Obtain local displacement
                d1 = localdof(1);
                d2 = localdof(4);
                
                % Obtain epsilon (strain)
                eps = (d2 - d1)/Le;
                
                % Obtain axial force
                obj.Faxial(E) = EA*eps ; 
                
            end
        end
    end
end