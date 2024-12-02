% Problem Parameters
c = 10;         % Boundary conductivity
u0 = 250;       % Outside temperature (K)
a_ice = 2.2;    % Thermal conductivity for ice
a_water = 0.6;  % Thermal conductivity for water
tol = 1e-3;     % Convergence tolerance for fixed-point iteration
tol_split = 1e-3; % Tolerance for 50/50 ice-water split
max_iter = 100; % Maximum number of iterations for fixed-point

% Geometry and Mesh
model = createpde();
geometryFromEdges(model, @circleg); % Circular domain
generateMesh(model, 'Hmax', 0.05, 'GeometricOrder', 'linear'); % Initial mesh
[p, e, t] = meshToPet(model.Mesh); % Extract mesh points, edges, and elements

% Define Bounds for Heat Source f
% Lower bound for heat source
f_min = 0;

% Upper bound for heat source
f_max = 10000;  

% Initialize Variables
%Initial guess for heat source
f = (f_min + f_max) / 2; 
a = a_water * ones(size(p, 2), 1); % Initial thermal conductivity

% Binary Search for Heat Source f
while abs(f_max - f_min) > tol_split
    % Initialize Fixed-Point Iteration for Updating a(u)
    for iter = 1:max_iter
        % Assemble System Matrices
        A = IntMatrix2(p, t, a);         % Stiffness matrix
        B = BdyMatrix(p, e, c);          % Boundary stiffness matrix
        G = BdyVector(p, e, c, u0);      % Boundary load vector
        F = IntVector(p, t, f);          % Internal load vector

        % Solve Linear System
        K = A + B;
        RHS = F + G;
        U = K \ RHS;

        % Update Thermal Conductivity
        a_new = a_water * ones(size(p, 2), 1); % Default to water
        a_new(U < 273) = a_ice; % Update nodes with ice
        if norm(a_new - a) < tol % Convergence check for fixed-point
            break;
        end
        a = a_new; % Update a
    end

    % Calculate Ice Fraction Based on Area
    ice_area = 0; % Reset ice area
    total_area = 0; % Reset total area

    for K = 1:size(t, 2)
        % Get the nodes of the current triangle
        nodes = t(1:3, K);
        coords = p(:, nodes);

        % Compute the area of the triangle
        x = coords(1, :);
        y = coords(2, :);
        area = 0.5 * abs(x(1)*(y(2)-y(3)) + x(2)*(y(3)-y(1)) + x(3)*(y(1)-y(2)));
        total_area = total_area + area;

        % Interpolate temperatures at the triangle's centroid
        U_centroid = mean(U(nodes)); % Average temperature at the triangle's centroid

        % Compute the fraction of the triangle that is ice
        if U_centroid < 273

            % Entire triangle is ice --> set 1 
            ice_fraction_triangle = 1; 

        elseif U_centroid > 273

            %Entire triangle is water
            ice_fraction_triangle = 0; 

        else

            % Linear interpolation ice fraction for partially ice-covered triangles
            ice_fraction_triangle = (273 - min(U(nodes))) / (max(U(nodes)) - min(U(nodes)));

        end

        % Add the ice contribution of this triangle
        ice_area = ice_area + ice_fraction_triangle * area;

    end

    ice_fraction = ice_area / total_area;

    % Adjust Heat Source f Based on Ice Fraction
    if abs(ice_fraction - 0.5) < tol_split
        disp(['Found f = ', num2str(f), ' for 50/50 ice/water split.']);
        break;

    elseif ice_fraction > 0.5
        f_min = f; % Too much ice, increase f
    else
        f_max = f; % Too little ice, decrease f
    end

    % Update f for the next iteration
    f = (f_min + f_max) / 2;
end

% Final Output
disp(['Final value of f: ', num2str(f)]);
disp(['Final ice fraction: ', num2str(ice_fraction)]);

% Visualize Final Temperature Distribution
figure;
pdeplot(model, 'XYData', U, 'Mesh', 'off');
title('Final Temperature Distribution');
xlabel('x'); ylabel('y');
colormap hot;
colorbar;
