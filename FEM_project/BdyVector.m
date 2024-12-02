% function G = BdyVector(p, e, c, u0)
%     % Assemble boundary load vector G
%     np = size(p, 2); % Number of nodes
%     G = zeros(np, 1); % Initialize load vector
% 
%     for k = 1:size(e, 2) % Loop over boundary edges
%         nodes = e(1:2, k); % Nodes of the edge
%         coords = p(:, nodes); % Coordinates of the edge
%         edgeLength = norm(diff(coords, 1, 2)); % Length of the edge
% 
%         % Contribution to the boundary load vector
%         G(nodes) = G(nodes) + c * u0 * edgeLength / 2 * [1; 1];
%     end
% end

% Boundary Vector
function G = BdyVector(p, e, c, u0)
    % Assemble boundary load vector G
    np = size(p, 2); % Number of nodes
    G = zeros(np, 1); % Initialize boundary load vector

    for k = 1:size(e, 2) % Loop over boundary edges
        nodes = e(1:2, k); % Nodes of the edge
        coords = p(:, nodes); % Coordinates of the edge
        edgeLength = norm(diff(coords, 1, 2)); % Length of the edge

        % Contribution to the boundary load vector
        G(nodes) = G(nodes) + c * u0 * edgeLength / 2 * [1; 1];
    end
end