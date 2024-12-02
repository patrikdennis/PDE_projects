function A = IntMatrixQuadrature(p, t, a)
    % IntMatrixQuadrature: Assembles the stiffness matrix using quadrature
    % Inputs:
    %   p - Node coordinates
    %   t - Element connectivity
    %   a - Thermal conductivity at nodes (vector, size: nNodes x 1)
    % Output:
    %   A - Stiffness matrix

    nNodes = size(p, 2); % Total number of nodes
    A = sparse(nNodes, nNodes); % Initialize stiffness matrix (sparse for efficiency)

    for K = 1:size(t, 2) % Loop over elements
        nodes = t(1:3, K); % Nodes of the current element
        coords = p(:, nodes); % Coordinates of the nodes
        x1 = coords(1, 1); y1 = coords(2, 1);
        x2 = coords(1, 2); y2 = coords(2, 2);
        x3 = coords(1, 3); y3 = coords(2, 3);
        area = 0.5 * abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)); % Element area

        % Compute gradients of shape functions (constant per element)
        gradN = [y2-y3, y3-y1, y1-y2; x3-x2, x1-x3, x2-x1] / (2*area);

        % Quadrature point: center of the element
        x_quad = mean(coords(1, :)); % x-coordinate of the centroid
        y_quad = mean(coords(2, :)); % y-coordinate of the centroid

        % Interpolate a(u) at the quadrature point using shape functions
        shape_func = [1/3; 1/3; 1/3]; % Shape functions evaluated at centroid
        a_quad = shape_func' * a(nodes); % Interpolated value of a(u) at the centroid

        % Local stiffness matrix
        localA = a_quad * (gradN' * gradN) * area;

        % Assemble into global stiffness matrix
        A(nodes, nodes) = A(nodes, nodes) + localA;
    end
end
