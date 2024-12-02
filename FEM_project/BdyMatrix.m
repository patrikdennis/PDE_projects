function B = BdyMatrix(p, e, c)
    % Assemble boundary stiffness matrix B (Boundary matrix)
    np = size(p, 2); % Number of nodes
    B = zeros(np, np); 

    for k = 1:size(e, 2) % Loop over boundary edges
        nodes = e(1:2, k); % Nodes of the edge
        coords = p(:, nodes); % Coordinates of the edge
        edgeLength = norm(diff(coords, 1, 2)); % Length of the edge

        % Contribution to the boundary stiffness matrix
        B(nodes, nodes) = B(nodes, nodes) + c * edgeLength / 6 * [2, 1; 1, 2];
    end
end
