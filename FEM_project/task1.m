% Problem Parameters
a = 0.6;   % Thermal conductivity of water
c = 10;    % Conductivity coefficient for boundary
u0 = 300;  % Surrounding temperature (K)
f = 100;   % Heat source

% Create PDE Model and Mesh
model = createpde();
geometryFromEdges(model, @circleg); % Circular domain (radius = 1)
generateMesh(model, 'Hmax', 0.2, 'GeometricOrder','linear');   % Generate mesh
[p, e, t] = meshToPet(model.Mesh); % Extract mesh points, edges, and elements


% a(phi, u) = L(phi)
% a(phi,u) = A + B, A = int_D, B = int_randentillD
% Assemble System Matrices and Vectors
A = IntMatrix(p, t, a);         % Stiffness matrix
B = BdyMatrix(p, e, c);         % Boundary stiffness matrix
F = IntVector(p, t, f);         % Internal load vector
G = BdyVector(p, e, c, u0);     % Boundary load vector

% Solve the Linear System
K = A + B;
% K = K + 1e-12 * speye(size(K));
% Global stiffness matrix
RHS = F + G;                    % Global load vector
U = K \ RHS;                    % Solve for nodal temperatures

% Visualize Results
pdeplot(model, 'XYData', U, 'Mesh', 'off'); % Plot solution
title('Heat Distribution');
xlabel('x'); ylabel('y');
