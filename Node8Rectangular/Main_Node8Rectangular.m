clear; close all; clc

%% Material Data
E = 1e7; nu = 0.3;
D = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];  % Plane Stress

%% Geometric Parameter
nx = 1; ny = 1; 
Lx = 2; Ly = 1;
X_origin = 0; Y_origin = 0;

[Coordinate, Nodes] = generate_Q8_mesh(Lx,Ly,nx,ny);

GaussNumber = 3;

%% Assemble Total K
nElem = size(Nodes,1);
PerNodeDof = 2;
TotalNodes = max(Nodes(:));
TotalDofs = PerNodeDof*TotalNodes;
K = zeros(TotalDofs,TotalDofs);
F = zeros(TotalDofs,1);
U = zeros(TotalDofs,1);

%% Gauss Points
switch GaussNumber
    case 2
        gaussPts = [-1/sqrt(3), 1/sqrt(3)];   
        weights = [1, 1];                     
    case 3
        gaussPts = [-sqrt(3/5), 0, sqrt(3/5)];  
        weights = [5/9, 8/9, 5/9];              
end

%% Loop over elements
for i = 1:nElem
    ke = zeros(16,16);
    Now_Coord = Coordinate(Nodes(i,:),:);
    Now_Dofs = TransNode2Dof(Nodes(i,:),PerNodeDof);
    
    for j = 1:GaussNumber
        for k = 1:GaussNumber
            xi  = gaussPts(j);
            eta = gaussPts(k);
            w   = weights(j)*weights(k);

            [N, dN_dxi, dN_deta] = q8_shape(xi, eta);

            J = [dN_dxi; dN_deta]*Now_Coord;
            detJ = det(J);
            invJ = inv(J);
            dN_dx = invJ * [dN_dxi; dN_deta];

            % Assemble B matrix
            B = zeros(3,16);
            for m = 1:8
                B(1,2*m-1) = dN_dx(1,m);
                B(2,2*m)   = dN_dx(2,m);
                B(3,2*m-1) = dN_dx(2,m);
                B(3,2*m)   = dN_dx(1,m);
            end

            ke = ke + B' * D * B * detJ * w;
        end
    end
    
    % Assemble global K
    K(Now_Dofs,Now_Dofs) = K(Now_Dofs,Now_Dofs) + ke;
end

ZerosColIdx = [];
for i = 1:size(K,1)
    if K(i,i) == 0
        ZerosColIdx = [ZerosColIdx;i];
    end
end

%% Boundary Conditions
FixedNode = find(abs(Coordinate(:,1)) < 1e-10);
FixedDofs = TransNode2Dof(FixedNode,PerNodeDof);
ActiveDofs = setdiff(1:TotalDofs, [FixedDofs,ZerosColIdx']);

FNode = find(abs(Coordinate(:,1)-1) < 1e-10);
FdDofs = TransNode2Dof(FNode,PerNodeDof);
F(FdDofs(2:2:end)) = -1;  % 只给竖向力

%% Solve
U(ActiveDofs) = K(ActiveDofs,ActiveDofs)\F(ActiveDofs);

%% Post-process
ux = U(1:2:end);
uy = U(2:2:end);
u_magnitude = sqrt(ux.^2 + uy.^2);

%% Write VTK File
writeVTK_Q8('results.vtk', Nodes, Coordinate, u_magnitude, 'disp_magnitude');

PlotDisplacement(Coordinate,Nodes,u_magnitude)
