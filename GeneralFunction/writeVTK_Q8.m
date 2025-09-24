function writeVTK_Q8(filename, Nodes, Coordinate, scalarField, fieldName)
    % filename: 输出文件名
    % Nodes: Q8 单元拓扑 (nElem x 8)
    % Coordinate: 节点坐标 (nNode x 2)
    % scalarField: 节点标量 (nNode x 1)
    % fieldName: 标量场名称

    nPoints = size(Coordinate,1);
    nElems = size(Nodes,1);

    triIdx = [
        1 5 4;
        4 7 5;
        7 5 2;
        7 3 2;
        ];

    nTriPerElem = size(triIdx,1);
    totalTri = nElems * nTriPerElem;

    % 打开文件
    fid = fopen(filename,'w');
    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'Q8 mesh as triangles\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');

    % 输出节点坐标
    fprintf(fid, 'POINTS %d float\n', nPoints);
    fprintf(fid, '%f %f %f\n', [Coordinate, zeros(nPoints,1)]');

    % 输出三角形单元
    fprintf(fid, 'CELLS %d %d\n', totalTri, totalTri*4);
    for e = 1:nElems
        elemNodes = Nodes(e,:);
        for t = 1:nTriPerElem
            fprintf(fid, '3 %d %d %d\n', elemNodes(triIdx(t,:)) - 1); % VTK 从 0 开始
        end
    end

    % 输出单元类型
    fprintf(fid, 'CELL_TYPES %d\n', totalTri);
    fprintf(fid, '%d\n', repmat(5, totalTri,1)); % VTK_TRIANGLE = 5

    % 输出节点标量数据
    fprintf(fid, 'POINT_DATA %d\n', nPoints);
    fprintf(fid, 'SCALARS %s float 1\n', fieldName);
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '%f\n', scalarField);

    fclose(fid);
end
