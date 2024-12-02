% Problem Parameters
a = 0.6;           % Thermal conductivity
u0 = 300;          % Surrounding temperature (K)
f = @(x, y) 100 * sin(6 * y); % Source function

% Define Circle with pdecirc
geometry = @circleg; % Define circular geometry

% Mesh Refinement
meshSize = 0.2; % Fixed mesh size

% Range of c values
c_values = logspace(-2, 2, 100); % Log-spaced c values from 0.01 to 100
condition_numbers = zeros(size(c_values)); % To store condition numbers

% Loop Over c Values
for i = 1:length(c_values)
    c = c_values(i); % Current value of c

    % Create PDE Model
    model = createpde();
    geometryFromEdges(model, geometry);
    generateMesh(model, 'Hmax', meshSize, 'GeometricOrder', 'linear'); % Mesh refinement
    [p, e, t] = meshToPet(model.Mesh); % Extract mesh points, edges, and elements

    % Assemble System Matrices and Vectors
    A = IntMatrix(p, t, a);         % Stiffness matrix
    B = BdyMatrix(p, e, c);         % Boundary stiffness matrix
    F = IntVectorQuad(p, t, f);     % Internal load vector using f(x, y)
    G = BdyVector(p, e, c, u0);     % Boundary load vector

    % Global System
    K = A + B;                      % Global stiffness matrix

    % Compute Condition Number
    condition_numbers(i) = cond(K); % Condition number of matrix K
end

% Log-Plot of Condition Number vs c
figure;
loglog(c_values, condition_numbers, 'b-', 'LineWidth', 1.5);
title('Condition Number of A vs c');
xlabel('c (Boundary Conductivity)');
ylabel('Condition Number');
grid on;

% Special Cases: c -> 0 and c -> Inf
% Case 1: c -> 0 (Weak boundary condition)
c = 1e-6; % Small c
model = createpde();
geometryFromEdges(model, geometry);
generateMesh(model, 'Hmax', meshSize, 'GeometricOrder', 'linear');
[p, e, t] = meshToPet(model.Mesh);
A = IntMatrix(p, t, a);
B = BdyMatrix(p, e, c);
F = IntVectorQuad(p, t, f);
G = BdyVector(p, e, c, u0);
K = A + B;
RHS = F + G;
U_weak = K \ RHS;

% Case 2: c -> Inf (Strong boundary condition)
c = 1e6; % Large c
model = createpde();
geometryFromEdges(model, geometry);
generateMesh(model, 'Hmax', meshSize, 'GeometricOrder', 'linear');
[p, e, t] = meshToPet(model.Mesh);
A = IntMatrix(p, t, a);
B = BdyMatrix(p, e, c);
F = IntVectorQuad(p, t, f);
G = BdyVector(p, e, c, u0);
K = A + B;
RHS = F + G;
U_strong = K \ RHS;

% Plot Results
figure;
subplot(1, 2, 1);
pdeplot(model, 'XYData', U_weak, 'Mesh', 'off');
title('Solution for c -> 0');
xlabel('x'); ylabel('y');
colormap hot;
axis equal;

subplot(1, 2, 2);
pdeplot(model, 'XYData', U_strong, 'Mesh', 'off');
title('Solution for c -> Inf');
xlabel('x'); ylabel('y');
colormap hot;
axis equal;
