% Problem Parameters
a = 0.6;      % Thermal conductivity
c = 4.8;      % Boundary conductivity
u0 = 300;     % Surrounding temperature (K)
f = 48;       % Heat source

% Define Circle with pdecirc
R = 1; % Radius of the circle
geometry = @circleg; % Define circular geometry

% Mesh Refinement Levels
meshSizes = [0.4, 0.2, 0.1, 0.05]; % Different levels of refinement
errors = zeros(size(meshSizes)); % To store maximum absolute errors

% Loop Over Mesh Sizes
for i = 1:length(meshSizes)
    % Create PDE Model
    model = createpde();
    geometryFromEdges(model, geometry);
    generateMesh(model, 'Hmax', meshSizes(i) , 'GeometricOrder','linear'); % Refine the mesh
    [p, e, t] = meshToPet(model.Mesh); % Extract mesh points, edges, and elements

    % Assemble System Matrices and Vectors
    A = IntMatrix(p, t, a);         % Stiffness matrix
    B = BdyMatrix(p, e, c);         % Boundary stiffness matrix
    F = IntVector(p, t, f);         % Internal load vector
    G = BdyVector(p, e, c, u0);     % Boundary load vector

    % Solve the Linear System
    K = A + B;                      % Global stiffness matrix
    RHS = F + G;                    % Global load vector
    K = K + 1e-12 * speye(size(K)); % Regularize to prevent numerical issues
    U = K \ RHS;                    % Solve for nodal temperatures

    % Analytical Solution
    % Ensure U_analytical is computed at the same node locations as U
    U_analytical = 325 - 20 * (p(1, :).^2 + p(2, :).^2);
    U_analytical = U_analytical';
    % Compute Maximum Absolute Error
    % Ensure U and U_analytical have the same length
    if length(U) == length(U_analytical)
        errors(i) = max(abs(U - U_analytical));
    else
        error('Mismatch in the length of numerical and analytical solution vectors');
    end

    % Plot Numerical and Analytical Solutions for the First Mesh
    if i == 2
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


% figuren ovan make it make sense

%%

% Plot Convergence (Log-Log Plot)
figure;
plot(meshSizes, errors, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
title('Convergence Analysis: Error vs Mesh Size');
xlabel('Mesh Size (h)');
ylabel('Maximum Absolute Error');
grid on;

% Reverse the x-axis direction
set(gca, 'XDir', 'reverse');

%%

% Estimate Convergence Rate
p = polyfit(log(meshSizes), log(errors), 1); % Linear fit in log-log scale
convergenceRate = p(1); % Slope of the line
disp(['Estimated Convergence Rate (N): ', num2str(convergenceRate)]);

pdeplot(model, 'XYData', U_analytical, 'Mesh', 'off'); % Plot solution
title('Heat Distribution Analytical Solution');
xlabel('x'); ylabel('y');
colormap hot;

%%

pdeplot(model, "XYData",U, 'Mesh', 'off');
title('Heat Distribution Numerical Solution');
xlabel('x'); ylabel('y');
colormap hot; 

% figuren ovan make it make sense


figure; % Create a new figure

% Analytical Solution Plot
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


%% error plot
% Create a tiled layout for larger subplots
% Create a figure and set its size
figure;
set(gcf, 'Position', [100, 100, 1200, 400]); % [x, y, width, height]

% Analytical Solution Plot
subplot(1, 3, 1); % Create a 1x3 grid of subplots, select the first subplot
pdeplot(model, 'XYData', U_analytical, 'Mesh', 'off');
title('Heat Distribution Analytical Solution');
xlabel('x'); ylabel('y');
colormap hot;
axis equal; % Ensure the axes are squared
%axis tight; % Adjust the axis limits to fit the data tightly

% Numerical Solution Plot
subplot(1, 3, 2); % Select the second subplot
pdeplot(model, "XYData", U, 'Mesh', 'off');
title('Heat Distribution Numerical Solution');
xlabel('x'); ylabel('y');
colormap hot;
axis equal; % Ensure the axes are squared
%axis tight; % Adjust the axis limits to fit the data tightly

% Error Plot
% Compute the pointwise absolute error
error_distribution = abs(U - U_analytical);

subplot(1, 3, 3); % Select the third subplot
pdeplot(model, 'XYData', error_distribution, 'Mesh', 'off');
title('Error Distribution');
xlabel('x'); ylabel('y');
colormap hot;
axis equal; % Ensure the axes are squared
%axis tight; % Adjust the axis limits to fit the data tightly
