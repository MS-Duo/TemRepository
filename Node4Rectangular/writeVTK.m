function writeVTK(filename, Nodes, Coordinate, scalarField, fieldName)
    % filename: 输出文件名（含 .vtk）
    % Nodes: 元素拓扑连接，四边形面
    % Coordinate: 顶点坐标 (n × 2)
    % scalarField: 节点上的标量（n × 1）
    % fieldName: 这个场的名称（如 'disp'）

    nPoints = size(Coordinate,1);
    nElems = size(Nodes,1);
    fid = fopen(filename,'w');

    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, 'FEM results\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');

    % 输出坐标（加一列 0，变成3D）
    fprintf(fid, 'POINTS %d float\n', nPoints);
    fprintf(fid, '%f %f %f\n', [Coordinate, zeros(nPoints,1)]');

    % 输出单元（四边形）
    % 每行: 5 n0 n1 n2 n3
    totalSize = 5 * nElems;
    fprintf(fid, 'CELLS %d %d\n', nElems, totalSize);
    fprintf(fid, '4 %d %d %d %d\n', (Nodes-1)');

    % 输出单元类型（9 对应 VTK_QUAD）
    fprintf(fid, 'CELL_TYPES %d\n', nElems);
    fprintf(fid, '%d\n', repmat(9, nElems, 1));

    % 输出节点数据
    fprintf(fid, 'POINT_DATA %d\n', nPoints);
    fprintf(fid, 'SCALARS %s float 1\n', fieldName);
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '%f\n', scalarField);

    fclose(fid);
end
