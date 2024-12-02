function A = IntMatrix2(p, t, a)
    % IntMatrix: Assembles the stiffness matrix for variable a(u)
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

        % Average a(u) for the element (can be refined to quadrature-based integration)
        a_elem = mean(a(nodes)); % Average of a(u) at the nodes of the element

        % Local stiffness matrix
        localA = a_elem.* (gradN' * gradN) * area;

        % Assemble into global stiffness matrix
        A(nodes, nodes) = A(nodes, nodes) + localA;
    end
end
