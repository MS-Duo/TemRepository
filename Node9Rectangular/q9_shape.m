function [N, dN_dxi, dN_deta] = q9_shape(xi, eta)

% shape functions
N1 = 0.25*xi*eta*(xi-1)*(eta-1);
N2 = 0.25*xi*eta*(xi+1)*(eta-1);
N3 = 0.25*xi*eta*(xi+1)*(eta+1);
N4 = 0.25*xi*eta*(xi-1)*(eta+1);

N5 = -0.5*eta*(xi+1)*(xi-1)*(eta-1);
N6 = -0.5*xi*(xi+1)*(eta+1)*(eta-1);
N7 = -0.5*eta*(xi+1)*(xi-1)*(eta+1);
N8 = -0.5*xi*(xi-1)*(eta+1)*(eta-1);
N9 = (xi+1)*(xi-1)*(eta+1)*(eta-1);
N = [N1 N2 N3 N4 N5 N6 N7 N8 N9];

dN_dxi = [
  0.25*eta*(2*xi-1)*(eta-1);
  0.25*eta*(2*xi+1)*(eta-1);
  0.25*eta*(2*xi+1)*(eta+1);
  0.25*eta*(2*xi-1)*(eta+1);
  -xi*eta*(eta-1);
  -0.5*(2*xi+1)*(eta^2-1);
  -xi*eta*(eta+1);
  -0.5*(2*xi-1)*(eta^2-1);
  2*xi*(eta^2-1)
];

dN_deta = [
  0.25*xi*(xi-1)*(2*eta-1);
  0.25*xi*(xi+1)*(2*eta-1);
  0.25*xi*(xi+1)*(2*eta+1);
  0.25*xi*(xi-1)*(2*eta+1);
  -0.5*(xi^2-1)*(2*eta-1);
  -xi*(xi+1)*eta;
  -0.5*(xi^2-1)*(2*eta+1);
  -xi*(xi-1)*eta;
  2*eta*(xi^2-1)
];
end
