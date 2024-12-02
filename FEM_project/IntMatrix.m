function A = IntMatrix(p, t, a)
    % Assemble the stiffness matrix A
    % p: Node coordinates (2 x number of nodes)
    % t: Element connectivity (3 x number of elements)
    % a: Thermal conductivity
    % A: Assembled stiffness matrix

    % phi = a_i + b_i * x_i + c_i * y_i 
    
    np = size(p, 2); % Number of nodes
    A = zeros(np, np); % Initialize sparse matrix

    for K = 1:size(t, 2) % Loop over elements
        nodes = t(1:3, K); % Nodes of the current element
        coords = p(:, nodes); % Coordinates of the nodes

        % Compute the area of the triangular element using a handwritten formula
        x1 = coords(1, 1); y1 = coords(2, 1); % Vertex 1
        x2 = coords(1, 2); y2 = coords(2, 2); % Vertex 2
        x3 = coords(1, 3); y3 = coords(2, 3); % Vertex 3
        area = 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));

        % Compute the gradients of the shape functions (Element Gradient)
        gradN = [y2 - y3, y3 - y1, y1 - y2;
                 x3 - x2, x1 - x3, x2 - x1] / (2 * area);

        % Contribution to the stiffness matrix for this element
        A(nodes, nodes) = A(nodes, nodes) + a * (gradN' * gradN) * area;
    end
end

