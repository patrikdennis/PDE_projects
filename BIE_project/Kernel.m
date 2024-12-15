function val = Kernel(i, j, y1, y2, nu1, nu2, rvec, rprimevec, rbisvec, k, N)
    % Kernel function (based on D.2)
    if i == j
        % Diagonal correction using curvature formula (Example D.2.2)
        val = (rprimevec(i)^2 - rvec(i) * rbisvec(i) / 2 + rvec(i)^2 / 2) / (2 * pi * (rprimevec(i)^2 + rvec(i)^2)^1.5);
    else
        % Off-diagonal kernel
        dx = y1(j) - y1(i);
        dy = y2(j) - y2(i);
        rdist = sqrt(dx^2 + dy^2); % Distance between points
        %grad_phi = besselh(1, k * rdist) * (k * [dx; dy] / rdist); % Gradient of Green's function
        grad_phi = (1i * k / 4) * (1 / rdist) * besselh(1, 1, k * rdist) * [dx; dy];
        val = grad_phi(1) * nu1(j) + grad_phi(2) * nu2(j); % Contribution
    end
end
