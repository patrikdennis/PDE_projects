%% deluppgift (a)
% Problem Parameters
% Problem Parameters
a = 0.6;      % Thermal conductivity
c = 9.6;       % Boundary conductivity
u0 = 300;     % Surrounding temperature (K)
f = 48;       % Heat source
f = @(x, y) 320*a*(x.^2+y.^2);
% Define Circle with pdecirc
geometry = @circleg; % Define circular geometry

% Mesh Refinement Levels
meshSizes = [0.4, 0.2, 0.1, 0.05]; % Different levels of refinement
errors = zeros(size(meshSizes)); % To store maximum absolute errors

% Loop Over Mesh Sizes
for i = 1:length(meshSizes)
    % Create PDE Model
    model = createpde();
    geometryFromEdges(model, geometry);
    generateMesh(model, 'Hmax', meshSizes(i), 'GeometricOrder', 'linear'); % Refine the mesh
    [p, e, t] = meshToPet(model.Mesh); % Extract mesh points, edges, and elements

    % Assemble System Matrices and Vectors
    A = IntMatrix(p, t, a);         % Stiffness matrix
    B = BdyMatrix(p, e, c);         % Boundary stiffness matrix
    F = IntVectorQuad(p, t, f); % Internal load vector using constant f
    G = BdyVector(p, e, c, u0);     % Boundary load vector

    % Solve the Linear System
    K = A + B;                      % Global stiffness matrix
    RHS = F + G;                    % Global load vector
    U = K \ RHS;                    % Solve for nodal temperatures

    % Analytical Solution
    U_analytical = 325 - 20 * (p(1, :).^2 + p(2, :).^2).^2; % Exact solution
    U_analytical = U_analytical'; % Ensure correct orientation

    % Compute Maximum Absolute Error
    errors(i) = max(abs(U - U_analytical)); % Ensure correct dimensions

    % Plot Numerical and Analytical Solutions for the First Mesh
    if i == 2 % Visualize for the second mesh refinement
        figure;
        plot(1:length(U), U, 'b-', 'DisplayName', 'Numerical Solution'); % Numerical Solution at Nodes
        hold on;
        plot(1:length(U_analytical), U_analytical, 'r-', 'DisplayName', 'Analytical Solution'); % Analytical Solution at Nodes
        legend show;
        title('Comparison of Numerical and Analytical Solutions');
        xlabel('Node Index');
        ylabel('Temperature');
    end
end


% Plot Convergence (Log-Log Plot)
figure;
loglog(meshSizes, errors, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
title('Convergence Analysis: Error vs Mesh Size');
xlabel('Mesh Size (h)');
ylabel('Maximum Absolute Error');
grid on;

set(gca, 'XDir', 'reverse');


% Estimate Convergence Rate
p = polyfit(log(meshSizes), log(errors), 1); % Linear fit in log-log scale
convergenceRate = p(1); % Slope of the line
disp(['Estimated Convergence Rate (N): ', num2str(convergenceRate)]);


% Analytical Solution Plot
figure;
subplot(1, 2, 1); % Create a 1x2 grid of subplots, select the first subplot
pdeplot(model, 'XYData', U_analytical, 'Mesh', 'off');
title('Heat Distribution Analytical Solution');
xlabel('x'); ylabel('y');
colormap hot;
axis equal; % Ensure the axes are squared
axis tight; % Adjust the axis limits to fit the data tightly

% Numerical Solution Plot
subplot(1, 2, 2); % Select the second subplot
pdeplot(model, "XYData", U, 'Mesh', 'off');
title('Heat Distribution Numerical Solution');
xlabel('x'); ylabel('y');
colormap hot;
axis equal; % Ensure the axes are squared
axis tight; % Adjust the axis limits to fit the data tightly


%% deluppgift (b)
%% deluppgift (b)
% Beer-Lambert source function and binary search for f0
mu = log(2); % Attenuation constant
f0_min = 0;  % Lower bound for f0
f0_max = 100; % Upper bound for f0
tol = 1e-3;  % Convergence tolerance for f0

% Define Circle with pdecirc
geometry = @circleg; % Define circular geometry

% Binary search loop for f0
while abs(f0_max - f0_min) > tol
    f0 = (f0_min + f0_max) / 2; % Midpoint of current range

    % Generate new mesh for each iteration
    model = createpde();
    geometryFromEdges(model, geometry);
    generateMesh(model, 'Hmax', 0.05, 'GeometricOrder', 'linear'); % Refine the mesh
    [p, e, t] = meshToPet(model.Mesh); % Extract mesh points, edges, and elements

    % Assemble System Matrices and Vectors
    A = IntMatrix(p, t, a);         % Stiffness matrix
    B = BdyMatrix(p, e, c);         % Boundary stiffness matrix
    K = A + B;                      % Global stiffness matrix
    F = IntVectorQuad(p, t, @(x, y) f0 * exp(-mu * sqrt(max(0, 1 - y.^2 - x)))); % Internal load vector
    G = BdyVector(p, e, c, u0);     % Boundary load vector

    % Solve the system
    RHS = F + G;
    K = K + 1e-12 * speye(size(K)); % Regularize to avoid numerical issues
    U = K \ RHS;

    % Check maximum temperature
    if max(U) > 320
        f0_max = f0; % Too hot, decrease f0
    else
        f0_min = f0; % Too cold, increase f0
    end
end

% Final f0 value
f0 = (f0_min + f0_max) / 2;
disp(['f0 for max temperature 320: ', num2str(f0)]);

% Plot the temperature distribution
figure;
pdeplot(model, 'XYData', U, 'Mesh', 'off');
title('Temperature Distribution with Beer-Lambert Source');
xlabel('x');
ylabel('y');
colorbar;
colormap hot;


%% redundant slow code does the same as the above 

% Problem Parameters
a = 0.6;          % Thermal conductivity
c = 10;           % Boundary conductivity
u0 = 300;         % Surrounding temperature (K)
mu = log(2);      % Attenuation coefficient

% Define Circle with pdecirc
geometry = @circleg; % Define circular geometry

% Mesh Refinement Levels
meshSizes = [0.4, 0.2, 0.1, 0.05]; % Different levels of refinement
errors = zeros(size(meshSizes)); % To store maximum absolute errors

% Find f0 for Maximum Temperature 320
f0 = 1; % Initial guess for f0
tolerance = 1e-3; % Convergence tolerance for f0
maxTemp = 0; % Initialize max temperature
while abs(maxTemp - 320) > tolerance
    % Define Source Function
    f = @(x, y) f0 * exp(-mu * (sqrt(1 - y.^2) - x));
    
    % Loop Over Mesh Sizes
    for i = 1:length(meshSizes)
        % Create PDE Model
        model = createpde();
        geometryFromEdges(model, geometry);
        generateMesh(model, 'Hmax', meshSizes(i), 'GeometricOrder', 'linear'); % Refine the mesh
        [p, e, t] = meshToPet(model.Mesh); % Extract mesh points, edges, and elements

        % Assemble System Matrices and Vectors
        A = IntMatrix(p, t, a);         % Stiffness matrix
        B = BdyMatrix(p, e, c);         % Boundary stiffness matrix
        F = IntVectorQuad(p, t, f); % Internal load vector using source function
        G = BdyVector(p, e, c, u0);     % Boundary load vector

        % Solve the Linear System
        K = A + B;                      % Global stiffness matrix
        RHS = F + G;                    % Global load vector
        U = K \ RHS;                    % Solve for nodal temperatures

        % Update Max Temperature
        maxTemp = max(U); % Maximum temperature in the domain

        % Break loop if mesh refinement is sufficient
        if i == length(meshSizes)
            break;
        end
    end
    
    % Adjust f0 based on max temperature
    if maxTemp < 320
        f0 = f0 * 1.01; % Increase f0
    else
        f0 = f0 * 0.9; % Decrease f0
    end
end

% Display final f0 value
disp(['Optimal f0: ', num2str(f0)]);

% Plot Heat Distribution for Final Solution
figure;
subplot(1, 2, 1); % Create a 1x2 grid of subplots, select the first subplot
pdeplot(model, 'XYData', U, 'Mesh', 'off');
title('Heat Distribution Numerical Solution');
xlabel('x'); ylabel('y');
colormap hot;
axis equal; % Ensure the axes are squared
axis tight; % Adjust the axis limits to fit the data tightly

% Analytical Solution Plot
U_analytical = 325 - 20 * (p(1, :).^2 + p(2, :).^2).^2; % Exact solution
subplot(1, 2, 2); % Select the second subplot
pdeplot(model, "XYData", U_analytical, 'Mesh', 'off');
title('Heat Distribution Analytical Solution');
xlabel('x'); ylabel('y');
colormap hot;
axis equal; % Ensure the axes are squared
axis tight; % Adjust the axis limits to fit the data tightly

