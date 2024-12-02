function F = IntVectorQuad(p, t, f_handle)
    % Assemble the load vector F using the quadrature rule (D.2)
    % p: Node coordinates (2 x number of nodes)
    % t: Element connectivity (3 x number of elements)
    % f_handle: Function handle for the source term f(x, y)

    np = size(p, 2); % Number of global nodes
    F = zeros(np, 1); % Initialize global source vector

    for K = 1:size(t, 2) % Loop over elements
        % Nodes of the current triangle
        nodes = t(1:3, K); % Indices of the triangle's nodes
        coords = p(:, nodes); % Coordinates of the triangle vertices

        % Triangle vertices
        x1 = coords(1, 1); y1 = coords(2, 1);
        x2 = coords(1, 2); y2 = coords(2, 2);
        x3 = coords(1, 3); y3 = coords(2, 3);

        % Compute the area of the triangle
        area = 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));

        % Quadrature points (midpoints of edges)
        q1 = [(x1 + x2) / 2; (y1 + y2) / 2];
        q2 = [(x2 + x3) / 2; (y2 + y3) / 2];
        q3 = [(x3 + x1) / 2; (y3 + y1) / 2];

        % Evaluate f at quadrature points
        f1 = f_handle(q1(1), q1(2));
        f2 = f_handle(q2(1), q2(2));
        f3 = f_handle(q3(1), q3(2));

        % Quadrature contribution to the load vector
        localF = (area / 3) * [f1; f2; f3];

        % Assemble contributions into the global load vector
        for i = 1:3
            F(nodes(i)) = F(nodes(i)) + localF(i);
        end
    end
end
