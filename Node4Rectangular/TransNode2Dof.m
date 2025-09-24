function DofsIdx = TransNode2Dof(Nodes,PerNodeDof)
if size(Nodes,1)>size(Nodes,2)
    Nodes = Nodes';
end
DofsIdx = kron(Nodes,ones(1,PerNodeDof)*PerNodeDof);
DofsIdx = DofsIdx + kron(ones(1,length(Nodes)),-PerNodeDof+1:0);
end