function [Coordinate, Nodes] = generate_Q9_mesh(L, W, m, n)

    nx = 2*m + 1;
    ny = 2*n + 1;
    totalNodes = nx * ny; 

    Coordinate = zeros(totalNodes, 2);
    
    dx = L / (2*m);
    dy = W / (2*n);
    
    % 生成节点坐标（按行优先：从下到上，从左到右）
    nodeIdx = 1;
    for iy = 0:ny-1
        y = iy * dy;
        for ix = 0:nx-1
            x = ix * dx;
            Coordinate(nodeIdx, :) = [x, y];
            nodeIdx = nodeIdx + 1;
        end
    end
    
    % 初始化单元节点矩阵
    totalElements = m * n;
    Nodes = zeros(totalElements, 9);
    
    % 生成每个单元的节点编号
    elemIdx = 1;
    for iy = 0:n-1       % 单元行索引 (y方向)
        for ix = 0:m-1   % 单元列索引 (x方向)
            % 计算全局节点索引（基于节点网格）
            baseRow = 2*iy;      % 单元底部行在节点网格中的行索引
            baseCol = 2*ix;      % 单元左侧列在节点网格中的列索引
            
            % 计算8个节点的全局索引（MATLAB索引从1开始）
            bottomLeft  = baseRow * nx + baseCol + 1;
            bottomRight = baseRow * nx + baseCol + 3;
            topRight    = (baseRow+2) * nx + baseCol + 3;
            topLeft     = (baseRow+2) * nx + baseCol + 1;
            bottomMid   = baseRow * nx + baseCol + 2;
            rightMid    = (baseRow+1) * nx + baseCol + 3;
            topMid      = (baseRow+2) * nx + baseCol + 2;
            leftMid     = (baseRow+1) * nx + baseCol + 1;
            MidQ9       = (baseRow+1) * nx + baseCol + 2;
            
            % 按Q8单元顺序存储节点编号
            Nodes(elemIdx, :) = [bottomLeft, bottomRight, topRight, topLeft, ...
                                 bottomMid, rightMid, topMid, leftMid MidQ9];
            elemIdx = elemIdx + 1;
        end
    end
end