clear;close all;clc
%% Material Data
E = 1e7; nu = 0.3;
D = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];  % Plane Stress Modal

%% 自己导入结构算
% [Coordinate,Nodes] = Obj2Matlab('MBB.obj');
% Coordinate = Coordinate(:,1:2);

%% Geometric Parameter
nx = 50; ny = 50; 
Lx = 2; Ly = 1;
dx = Lx / nx; dy = Ly / ny;
GaussNmuber = 3;
%% 矩形
nodeID = @(i,j) (j-1)*(nx+1)+i;
Coordinate = zeros((nx+1)*(ny+1), 2);
for j = 1:(ny+1)
    for i = 1:(nx+1)
        id = nodeID(i,j);
        Coordinate(id,:) = [(i-1)*dx, (j-1)*dy];
    end
end

%% 单元连接生成
Nodes = zeros(nx*ny, 4);
elem = 0;
for j = 1:ny
    for i = 1:nx
        elem = elem + 1;
        n1 = nodeID(i,j);
        n2 = nodeID(i+1,j);
        n3 = nodeID(i+1,j+1);
        n4 = nodeID(i,j+1);
        Nodes(elem,:) = [n1 n2 n3 n4];
    end
end

%% Assemble Total K
nElem = size(Nodes,1);
PerNodeDof = 2;
TotalDofs = PerNodeDof*size(Coordinate,1);
K = zeros(TotalDofs,TotalDofs);
F = zeros(TotalDofs,1);
U = zeros(TotalDofs,1);

%% Gauss Point
switch GaussNmuber
    case 2
        gaussPts = [-1/sqrt(3), 1/sqrt(3)];   % 积分域的坐标
        weights = [1, 1];                     % 对应的权重
    case 3
        gaussPts = [-0.775 0 0.775];          % 积分域的坐标
        weights = [5/9 8/9 5/9];              % 对应的权重
end

gaussPts = [-1/sqrt(3), 1/sqrt(3)];   % 积分域的坐标
weights = [1, 1];                     % 对应的权重

for i = 1:nElem
    ke = zeros(8, 8);
    Now_Coord = Coordinate(Nodes(i,:),:);
    Now_Dofs = TransNode2Dof(Nodes(i,:),PerNodeDof);
    % Assume Gauss Points 2*2
    for j = 1:length(gaussPts)
        for k = 1:length(gaussPts)
            s = gaussPts(j);
            t  = gaussPts(k);
            w = weights(j)*weights(k);
            
            % Generate B Matrix dN/dx
            dN_dX = 1/4*[-(1+t) -(1-t) 1-t 1+t;1-s -(1-s) -(1+s) 1+s];  % 形函数关于积分域求导
            J = dN_dX*Now_Coord;
            detJ = det(J);
            invJ = inv(J);
            dN_dx = invJ*dN_dX;   % 形状函数关于全局坐标系求导

            B = zeros(3, 8);
            for m = 1:4
                B(1,2*m-1) = dN_dx(1,m);
                B(2,2*m)   = dN_dx(2,m);
                B(3,2*m-1) = dN_dx(2,m);
                B(3,2*m)   = dN_dx(1,m);
            end
            ke = ke + B' * D * B * detJ  * w;

        end
    end
    K(Now_Dofs,Now_Dofs) = K(Now_Dofs,Now_Dofs) + ke;
end

% 4----3-----6
% |    |     |
% 1----2-----5
FixedNode = find(abs(Coordinate(:,1)) < 1e-10);
FixedDofs = TransNode2Dof(FixedNode,PerNodeDof);
ActiveDofs = setdiff((1:TotalDofs),FixedDofs);
FNode = find(abs(Coordinate(:,1) - Lx) < 1e-10 & abs(Coordinate(:,2) - Ly) < 1e-10);
FdDofs = TransNode2Dof(FNode,PerNodeDof);

F((FdDofs(2:2:length(FdDofs))),1) = -1;
U(ActiveDofs) = K(ActiveDofs,ActiveDofs)\F(ActiveDofs);


%% === 位移云图 ===
ux = U(1:2:end);
uy = U(2:2:end);
u_magnitude = sqrt(ux.^2 + uy.^2);  % 每个节点的总位移
umax = max(u_magnitude)
% Write VTK File
% writeVTK('results.vtk', Nodes, Coordinate, u_magnitude, 'disp_magnitude');

PlotDisplacement(Coordinate,Nodes,u_magnitude)

%% === 可视化约束 ===
hold on;
plot(Coordinate(FixedNode,1), Coordinate(FixedNode,2), ...
     'ko', 'MarkerFaceColor','k', 'MarkerSize',6);  % 红色实心圆表示固定
%% === 可视化载荷 ===
hold on;
quiver(Coordinate(FNode,1), Coordinate(FNode,2), ...
       zeros(size(FNode)), -ones(size(FNode))*0.5, ...
       0, 'r', 'LineWidth', 1.5);  % 红色向下箭头