function [N, dN_dxi, dN_deta] = q8_shape(xi, eta)

% shape functions
N1 = -0.25*(1-xi)*(1-eta)*(1+xi+eta);
N2 = -0.25*(1+xi)*(1-eta)*(1-xi+eta);
N3 = -0.25*(1+xi)*(1+eta)*(1-xi-eta);
N4 = -0.25*(1-xi)*(1+eta)*(1+xi-eta);

N5 = 0.5*(1-xi^2)*(1-eta);
N6 = 0.5*(1+xi)*(1-eta^2);
N7 = 0.5*(1-xi^2)*(1+eta);
N8 = 0.5*(1-xi)*(1-eta^2);

N = [N1 N2 N3 N4 N5 N6 N7 N8];

% dN/dxi
dN1_dxi =  0.25*(1-eta).*(2*xi + eta);
dN2_dxi =  0.25*(1-eta).*(2*xi - eta);
dN3_dxi =  0.25*(1+eta).*(2*xi + eta);
dN4_dxi =  0.25*(1+eta).*(2*xi - eta);
dN5_dxi = -xi.*(1-eta);
dN6_dxi =  0.50*(1 - eta.^2);
dN7_dxi = -xi.*(1+eta);
dN8_dxi = -0.50*(1 - eta.^2);

dN_dxi = [dN1_dxi dN2_dxi dN3_dxi dN4_dxi dN5_dxi dN6_dxi dN7_dxi dN8_dxi];

% dN/deta
dN1_deta = 0.25*(1-xi).*(xi + 2*eta);
dN2_deta = 0.25*(1+xi).*(-xi + 2*eta);
dN3_deta = 0.25*(1+xi).*(xi + 2*eta);
dN4_deta = 0.25*(1-xi).*(2*eta - xi);
dN5_deta = -0.50*(1 - xi.^2);
dN6_deta = -(1+xi).*eta;
dN7_deta =  0.50*(1 - xi.^2);
dN8_deta = -(1 - xi).*eta;

dN_deta = [dN1_deta dN2_deta dN3_deta dN4_deta dN5_deta dN6_deta dN7_deta dN8_deta];
end
