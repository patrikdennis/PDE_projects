% Problem Parameters
a = 0.6; % Thermal conductivity
u0 = 300; % Boundary temperature
f_handle = @(x, y) 100 * sin(6 * y); % Source function

% Geometry and Mesh
geometry = @circleg; % Circular geometry
model = createpde();
geometryFromEdges(model, geometry);
generateMesh(model, 'Hmax', 0.2, 'GeometricOrder', 'linear'); % Mesh
[p, e, t] = meshToPet(model.Mesh); % Extract mesh points, edges, and elements

% Range of c values
c_values = logspace(-5, 5, 50); % Logarithmic scale for c
condition_numbers = zeros(size(c_values)); % Store condition numbers

% Loop over c values
for i = 1:length(c_values)
    c = c_values(i);

    % Assemble System Matrices
    A = IntMatrix(p, t, a); % Stiffness matrix
    B = BdyMatrix(p, e, c); % Boundary stiffness matrix
    K = A + B; % Global stiffness matrix

    % Compute condition number
    condition_numbers(i) = cond(K);
end

% Plot Condition Number vs c
figure;
loglog(c_values, condition_numbers, 'b-', 'LineWidth', 1.5);
xlabel('Boundary Conductivity (c)');
ylabel('Condition Number of A');
title('Condition Number vs c');
grid on;

% Solve for u with c ≈ 0 and c → ∞
c_approx_zero = 1e-5; % c ≈ 0
c_large = 1e5; % c → ∞

% Case: c ≈ 0
B_zero = BdyMatrix(p, e, c_approx_zero);
K_zero = A + B_zero;
F = IntVectorQuad(p, t, f_handle); % Internal load vector
G_zero = BdyVector(p, e, c_approx_zero, u0); % Boundary load vector
RHS_zero = F + G_zero;
K_zero = K_zero + 1e-12 * speye(size(K_zero)); % Regular




% Compute RHS for c ≈ 0
u_zero = K_zero \ RHS_zero; % Solve the linear system for u (c ≈ 0)

% Case: c → ∞
B_large = BdyMatrix(p, e, c_large);
K_large = A + B_large;
G_large = BdyVector(p, e, c_large, u0);
RHS_large = F + G_large;
K_large = K_large + 1e-12 * speye(size(K_large)); % Regularization
u_large = K_large \ RHS_large; % Solve for u (c → ∞)

% Plot solutions
figure;
subplot(1, 2, 1);
pdeplot(model, 'XYData', u_zero, 'Mesh', 'off');
title('Solution for c ≈ 0');
xlabel('x'); ylabel('y');
colormap hot;

subplot(1, 2, 2);
pdeplot(model, 'XYData', u_large, 'Mesh', 'off');
title('Solution for c → ∞');
xlabel('x'); ylabel('y');
colormap hot;

