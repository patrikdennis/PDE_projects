% Wavenumber
k = 1;

% Number of discretization points on the boundary
N = 1000;

% Number of field points for plotting solution
M = 100; 

%% Boundary discretization

% Parameterization variable
tvec = linspace(-pi + N, pi, N);

% Radial distance
rvec = 3 + cos(4*tvec + pi);

% First derivative w.r.t. t
rprimevec = -4*sin(4*tvec + pi);

% Second derivative w.r.t. t
rbisvec = -16*cos(4*tvec + pi);

% Boundary coordinates

% x-coordinates of boundary points
y1 = rvec .* cos(tvec); 

% y-coordinates of boundary points
y2 = rvec .* sin(tvec); 

% Derivatives w.r.t. t

% dx/dt
x_t = rprimevec .* cos(tvec) - rvec .* sin(tvec); 

% dy/dt
y_t = rprimevec .* sin(tvec) + rvec .* cos(tvec); 

% Arclength element
dsdt = sqrt(rprimevec.^2 + rvec.^2); 

% Outward unit normal vectors
nu1 = rvec .* cos(tvec) + rprimevec .* sin(tvec);
nu2 = rvec .* sin(tvec) - rprimevec .* cos(tvec);
normalizer = sqrt(rvec.^2 + rprimevec.^2);

% Normalize normal vector
nu1 = nu1 ./ normalizer; 
nu2 = nu2 ./ normalizer;

% Assemble the BIE matrix
Kmat = zeros(N, N);
for i = 1:N
    for j = 1:N
        Kmat(i, j) = Kernel(i, j, y1, y2, nu1, nu2, rvec, rprimevec, rbisvec, k, N);
    end
end

% Visualization
figure;
imagesc(tvec, tvec, real(Kmat.'));
axis xy;
colorbar;
title('Boundary Integral Kernel Matrix');

figure;
imagesc(tvec,tvec, imag(Kmat.'));
axis xy;
colorbar;


%% Compute Neumann data (Example C.1)

% Example point for Neumann data computation
p = [0; 0]; 

% Initialize Neumann data
gvec = zeros(N, 1); 

for i = 1:N
    xi = [y1(i); y2(i)];

    % Distance from boundary point to p
    diffp = xi - p; 
    rp = norm(diffp);
    
    % Green's function
    H1p = besselh(1, 1, k * rp); 

    % Gradient at point p
    gradv1 = (1i * k / 4) * (diffp / rp) * H1p; 

    % Neumann data
    gvec(i) = gradv1.' * [nu1(i); nu2(i)]; 
end

% Plot Neumann data
figure;
plot(tvec, real(gvec));
title('Neumann Data');
xlabel('t (boundary parameter)');
ylabel('g (Neumann data)');
grid on;


%% Kernel function (based on D.2)
function val = Kernel(i, j, y1, y2, nu1, nu2, rvec, rprimevec, rbisvec, k, N)
    if i == j
        % Diagonal correction using curvature formula (Example D.2.2)
        val = (rprimevec(i)^2 - rvec(i) * rbisvec(i) / 2 + rvec(i)^2 / 2) / ...
              (2 * pi * (rprimevec(i)^2 + rvec(i)^2)^1.5);
    else
        % Off-diagonal kernel
        dx = y1(j) - y1(i);
        dy = y2(j) - y2(i);

        % Distance between points
        rdist = sqrt(dx^2 + dy^2); 

        % Gradient of Green's function
        %grad_phi = besselh(1, 1, k * rdist) * (-k * [dx; dy] / rdist); 
        grad_phi = (1i * k / 4) * (1 / rdist) * besselh(1, 1, k * rdist) * [dx; dy];

        % Contribution
        val = grad_phi(1) * nu1(j) + grad_phi(2) * nu2(j); 
    end
end

