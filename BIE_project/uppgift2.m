% Wavenumber
k = 1;

% Number of discretization points on the boundary
N = 1000; 
M = 100; % Number of field points in D^-
p = [0; 0]; % Fixed point in D^+

%% Boundary discretization

% Parameterization variable
tvec = linspace(-pi, pi, N);

% Radial distance
rvec = 3 + cos(4 * tvec + pi);

% First derivative r'
rprimevec = -4 * sin(4 * tvec + pi);

% Second derivative r''
rbisvec = -16 * cos(4 * tvec + pi); 

% x-coordinates of boundary
y1 = rvec .* cos(tvec);

% y-coordinates of boundary
y2 = rvec .* sin(tvec); 

% Arc length element
dsdt = sqrt(rprimevec.^2 + rvec.^2); 

% Outward normal vectors
nu1 = rvec .* cos(tvec) + rprimevec .* sin(tvec);
nu2 = rvec .* sin(tvec) - rprimevec .* cos(tvec);
normalizer = sqrt(rvec.^2 + rprimevec.^2);
nu1 = nu1 ./ normalizer;
nu2 = nu2 ./ normalizer;

%% Solve BIE for density h

% Exact Neumann data
gvec = zeros(N, 1);
for i = 1:N
    xi = [y1(i); y2(i)];
    diffp = xi - p;
    rp = norm(diffp);
    H1p = besselh(1,1, k * rp);
    gradv1 = -k * (diffp / rp) * H1p;
    gvec(i) = gradv1.' * [nu1(i); nu2(i)];
end

% build matrix Kmat
Kmat = zeros(N, N);
for i = 1:N
    for j = 1:N
        Kmat(i, j) = Kernel(i, j, y1, y2, nu1, nu2, rvec, rprimevec, rbisvec, k, N);
    end
end

% solve (1/2 I + (2pi/N)K^T*dsdt)h = g
A = 0.5 * eye(N) + (2*pi / N) * Kmat.' * diag(dsdt);
h = A \ gvec;

%% Compute numerical solution in D^-
xfield = linspace(-4, 4, M); % Extend grid for exterior domain
yfield = linspace(-4, 4, M);
ufield_numerical = zeros(M, M); % Initialize solution in D^-

for ix1 = 1:M
    for ix2 = 1:M
        x1 = xfield(ix1); % Field point x-coordinate
        x2 = yfield(ix2); % Field point y-coordinate

        % Check if point (x1, x2) is inside the boundary (D^+)
        % Mask the interior points based on inpolygon
        % Only process exterior points
        if ~inpolygon(x1, x2, y1, y2) 

            % Compute Green's function contributions
            % Initialize kernel contributions
            phi_vec = zeros(N, 1); 

            for j = 1:N

                % Boundary point
                yj = [y1(j); y2(j)];

                % Distance between (x1, x2) and yj
                rdist = sqrt((x1 - yj(1))^2 + (x2 - yj(2))^2); 

                % Green's function
                phi_vec(j) = -(1i / 4) * besselh(0, k * rdist);             
            end

            % Use the trapezoidal rule to compute u(x) in the exterior
            ufield_numerical(ix1, ix2) = sum(phi_vec .* (h .* dsdt.')) * (2 * pi / N); % Trapezoidal integration
        else
            % Mask interior points
            ufield_numerical(ix1, ix2) = NaN;
        end
    end
end

%% Exact solution

% Initialize ufield for the exact solution
xfield = linspace(-4,4,M);
yfield = linspace(-4,4,M);
ufield_exact = zeros(M, M);

% Compute exact solution in D^-
for ix1 = 1:M
    for ix2 = 1:M
        x1 = xfield(ix1); % Field point x-coordinate
        x2 = yfield(ix2); % Field point y-coordinate

        % Check if the point is inside the boundary
        if ~inpolygon(x1, x2, y1, y2) % Only process exterior points
            rdist = norm([x1; x2] - p); % Distance from field point to fixed point p
            ufield_exact(ix1, ix2) =  besselh(0, k * rdist); % Exact solution
        else
            ufield_exact(ix1, ix2) = NaN; % Mask interior points as NaN
        end
    end
end


%% Plots

figure;

% Plot real part of numerical solution
subplot(2, 2, 1);
imagesc(xfield, yfield, real(ufield_numerical.'));
axis xy;
colormap turbo;
pbaspect([1 1 1]);
colorbar;
title('Real Part of Numerical Solution in D^-');
xlabel('x'); ylabel('y');

% Plot imaginary part of numerical solution
subplot(2, 2, 2);
imagesc(xfield, yfield, imag(ufield_numerical.'));
axis xy;
colormap turbo;
pbaspect([1 1 1]);
colorbar;
title('Imaginary Part of Numerical Solution in D^-');
xlabel('x'); ylabel('y');

% Plot real part of exact solution
subplot(2, 2, 3);
imagesc(xfield, yfield, real(ufield_exact.'));
axis xy;
colormap turbo;
pbaspect([1 1 1]);
colorbar;
title('Real Part of Exact Solution in D^-');
xlabel('x'); ylabel('y');

% Plot imaginary part of exact solution
subplot(2, 2, 4);
imagesc(xfield, yfield, imag(ufield_exact.'));
axis xy;
colormap turbo;
pbaspect([1 1 1]);
colorbar;
title('Imaginary Part of Exact Solution in D^-');
xlabel('x'); ylabel('y');

%% Plot absolute error

% Compute the absolute error across the domain
error_abs = abs(ufield_numerical - ufield_exact);

% Compute the log10 of the absolute error
log_error = log10(error_abs);

% Mask invalid values (e.g., NaN from interior points)
log_error(isnan(log_error)) = -Inf; % Set NaN to -Inf for better visualization in log scale

% Plot the log10 of the absolute error
figure;
imagesc(xfield, yfield, log_error.');
axis xy;
colormap turbo;
pbaspect([1 1 1]);
colorbar;
caxis([-15, 0]); % Adjust color range for better visualization (adjust as needed)
title('Log10 of Absolute Error');
xlabel('x');
ylabel('y');

%% error convergence study

% Parameters
k = 1; % Wavenumber
M = 100; % Number of field points in D^-
p = [0; 0]; % Fixed point in D^+

N_values = [200, 400, 800, 1600]; % Different values of N
errors = zeros(size(N_values)); % Store average errors

for idx = 1:length(N_values)
    % Current N
    N = N_values(idx);
    
    % Boundary discretization
    tvec = linspace(-pi, pi, N);
    rvec = 3 + cos(4 * tvec + pi);
    rprimevec = -4 * sin(4 * tvec + pi);
    rbisvec = -16 * cos(4 * tvec + pi);
    y1 = rvec .* cos(tvec);
    y2 = rvec .* sin(tvec);
    dsdt = sqrt(rprimevec.^2 + rvec.^2);
    nu1 = rvec .* cos(tvec) + rprimevec .* sin(tvec);
    nu2 = rvec .* sin(tvec) - rprimevec .* cos(tvec);
    normalizer = sqrt(rvec.^2 + rprimevec.^2);
    nu1 = nu1 ./ normalizer;
    nu2 = nu2 ./ normalizer;

    % Exact Neumann data
    gvec = zeros(N, 1);
    for i = 1:N
        xi = [y1(i); y2(i)];
        diffp = xi - p;
        rp = norm(diffp);
        H1p = besselh(1,1, k * rp);
        gradv1 = -k * (diffp / rp) * H1p;
        gvec(i) = gradv1.' * [nu1(i); nu2(i)];
    end

    % Solve BIE for density h
    Kmat = zeros(N, N);
    for i = 1:N
        for j = 1:N
            Kmat(i, j) = Kernel(i, j, y1, y2, nu1, nu2, rvec, rprimevec, rbisvec, k, N);
        end
    end
    A = 0.5 * eye(N) + (2*pi / N) * Kmat.' * diag(dsdt);
    h = A \ gvec;

    % Compute numerical solution in D^-
    xfield = linspace(-4, 4, M);
    yfield = linspace(-4, 4, M);
    ufield_numerical = zeros(M, M);
    for ix1 = 1:M
        for ix2 = 1:M
            x1 = xfield(ix1);
            x2 = yfield(ix2);
            if ~inpolygon(x1, x2, y1, y2)
                phi_vec = zeros(N, 1);
                for j = 1:N
                    yj = [y1(j); y2(j)];
                    rdist = sqrt((x1 - yj(1))^2 + (x2 - yj(2))^2);
                    phi_vec(j) = -(1i / 4) * besselh(0, k * rdist);
                end
                ufield_numerical(ix1, ix2) = sum(phi_vec .* (h .* dsdt.')) * (2 * pi / N);
            else
                ufield_numerical(ix1, ix2) = NaN;
            end
        end
    end

    % Compute exact solution in D^-
    ufield_exact = zeros(M, M);
    for ix1 = 1:M
        for ix2 = 1:M
            x1 = xfield(ix1);
            x2 = yfield(ix2);
            if ~inpolygon(x1, x2, y1, y2)
                rdist = norm([x1; x2] - p);
                ufield_exact(ix1, ix2) = besselh(0, k * rdist);
            else
                ufield_exact(ix1, ix2) = NaN;
            end
        end
    end

    % Compute average absolute error in D^- (away from boundary)
    error_abs = abs(ufield_numerical - ufield_exact);
    error_abs(isnan(error_abs)) = 0; % Ignore NaN values (masked points)
    errors(idx) = mean(error_abs(:)); % Average error
end

%% Log-log plot of errors vs. N
figure;
loglog(N_values, errors, '-o', 'LineWidth', 2);
xlabel('Number of Boundary Points (N)');
ylabel('Average Absolute Error');
title('Convergence of Error with Increasing N');
grid on;

%% Fit log-log data to find convergence rate k
coeffs = polyfit(log10(N_values), log10(errors), 1);
k = -coeffs(1); % Convergence rate (slope of the line)
fprintf('Estimated convergence rate: O(1/N^%.2f)\n', k);


%% Error for real and imaginary seperately
% Parameters
k = 1; % Wavenumber
M = 100; % Number of field points in D^-
p = [0; 0]; % Fixed point in D^+

N_values = [200, 400, 800, 1600]; % Different values of N
errors_real = zeros(size(N_values)); % Store real part average errors
errors_imag = zeros(size(N_values)); % Store imaginary part average errors

for idx = 1:length(N_values)
    % Current N
    N = N_values(idx);
    
    %  Boundary discretization
    tvec = linspace(-pi, pi, N);
    rvec = 3 + cos(4 * tvec + pi);
    rprimevec = -4 * sin(4 * tvec + pi);
    rbisvec = -16 * cos(4 * tvec + pi);
    y1 = rvec .* cos(tvec);
    y2 = rvec .* sin(tvec);
    dsdt = sqrt(rprimevec.^2 + rvec.^2);
    nu1 = rvec .* cos(tvec) + rprimevec .* sin(tvec);
    nu2 = rvec .* sin(tvec) - rprimevec .* cos(tvec);
    normalizer = sqrt(rvec.^2 + rprimevec.^2);
    nu1 = nu1 ./ normalizer;
    nu2 = nu2 ./ normalizer;

    % Exact Neumann data
    gvec = zeros(N, 1);
    for i = 1:N
        xi = [y1(i); y2(i)];
        diffp = xi - p;
        rp = norm(diffp);
        H1p = besselh(1,1, k * rp);
        gradv1 = -k * (diffp / rp) * H1p;
        gvec(i) = gradv1.' * [nu1(i); nu2(i)];
    end

    % Solve BIE for density h
    Kmat = zeros(N, N);
    for i = 1:N
        for j = 1:N
            Kmat(i, j) = Kernel(i, j, y1, y2, nu1, nu2, rvec, rprimevec, rbisvec, k, N);
        end
    end
    A = 0.5 * eye(N) + (2*pi / N) * Kmat.' * diag(dsdt);
    h = A \ gvec;

    % Compute numerical solution in D^-
    xfield = linspace(-4, 4, M);
    yfield = linspace(-4, 4, M);
    ufield_numerical = zeros(M, M);
    for ix1 = 1:M
        for ix2 = 1:M
            x1 = xfield(ix1);
            x2 = yfield(ix2);
            if ~inpolygon(x1, x2, y1, y2)
                phi_vec = zeros(N, 1);
                for j = 1:N
                    yj = [y1(j); y2(j)];
                    rdist = sqrt((x1 - yj(1))^2 + (x2 - yj(2))^2);
                    phi_vec(j) = -(1i / 4) * besselh(0, k * rdist);
                end
                ufield_numerical(ix1, ix2) = sum(phi_vec .* (h .* dsdt.')) * (2 * pi / N);
            else
                ufield_numerical(ix1, ix2) = NaN;
            end
        end
    end

    % Compute exact solution in D^-
    ufield_exact = zeros(M, M);
    for ix1 = 1:M
        for ix2 = 1:M
            x1 = xfield(ix1);
            x2 = yfield(ix2);
            if ~inpolygon(x1, x2, y1, y2)
                rdist = norm([x1; x2] - p);
                ufield_exact(ix1, ix2) = besselh(0, k * rdist);
            else
                ufield_exact(ix1, ix2) = NaN;
            end
        end
    end

    % Compute real and imaginary part errors
    error_real = abs(real(ufield_numerical) - real(ufield_exact));
    error_imag = abs(imag(ufield_numerical) - imag(ufield_exact));
    
    % Mask NaN values (interior points)
    error_real(isnan(error_real)) = 0;
    error_imag(isnan(error_imag)) = 0;

    % Average errors in D^- (excluding boundary points)
    errors_real(idx) = mean(error_real(:));
    errors_imag(idx) = mean(error_imag(:));
end

%% Log-log plots for real and imaginary part errors
figure;
loglog(N_values, errors_real, '-o', 'LineWidth', 2);
hold on;
loglog(N_values, errors_imag, '-s', 'LineWidth', 2);
xlabel('Number of Boundary Points (N)');
ylabel('Average Absolute Error');
legend('Real Part Error', 'Imaginary Part Error');
title('Convergence of Real and Imaginary Part Errors');
grid on;

%% Fit log-log data to find convergence rates
coeffs_real = polyfit(log10(N_values), log10(errors_real), 1);
coeffs_imag = polyfit(log10(N_values), log10(errors_imag), 1);
k_real = -coeffs_real(1); % Convergence rate for real part
k_imag = -coeffs_imag(1); % Convergence rate for imaginary part
fprintf('Estimated convergence rate for real part: O(1/N^%.2f)\n', k_real);
fprintf('Estimated convergence rate for imaginary part: O(1/N^%.2f)\n', k_imag);

