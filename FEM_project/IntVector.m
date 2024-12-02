function F = IntVector(p, t, f)
    % Assemble internal load vector F
    np = size(p, 2); % Number of nodes
    F = zeros(np, 1); % Initialize load vector

    for K = 1:size(t, 2) % Loop over elements
        nodes = t(1:3, K); % Nodes of the current element
        coords = p(:, nodes); % Coordinates of the nodes
        area = polyarea(coords(1, :), coords(2, :)); % Element area
       
        % Contribution to the load vector for this element
        F(nodes) = F(nodes) + f * area / 3; % Uniform load
    end
end
