% Number of boundary points
N = 1000;

% Far-field directions
%theta_values = linspace(0, 2*pi, 360); 
theta_values = linspace(-pi + 2*pi/N,pi,N);

% Chosen k-values to verify the optic theorem in R^2
k_values = [0.00000001, 0.0001, 0.05, 0.5, 1, 5];           

% Domain discretization of smooth X-shape
tvec = linspace(-pi + 2*pi/N, pi, N);
%theta_values = linspace(-pi,pi, 360);
rvec = 3 + cos(4 * tvec + pi);
rprimevec = -4 * sin(4 * tvec + pi);
rbisvec = -16 * cos(4 * tvec + pi);
y1 = rvec .* cos(tvec);
y2 = rvec .* sin(tvec);
dsdt = sqrt(rprimevec.^2 + rvec.^2);
nu1 = (rvec .* cos(tvec) + rprimevec .* sin(tvec)) ./ dsdt;
nu2 = (rvec .* sin(tvec) - rprimevec .* cos(tvec)) ./ dsdt;

%%

% small enough k 
k_small = 1e-10;

c_const = (1 + 1i)/(4*sqrt(pi));

theta_values = linspace(-pi + 2*pi/N,pi, N);

area_dplus = compute_area(y1, y2, dsdt, nu1, nu2 , rvec, rprimevec, rbisvec, k_small, c_const, theta_values,N);

fprintf('Area of D+: %g \n', area_dplus);

%%

% Far-field directions
%theta_values = linspace(0, 2*pi, 360); 
%theta_values = linspace(-pi,pi,360);

% Chosen k-values to verify the optic theorem in R^2
k_values = [0.00000001, 0.0001, 0.05, 0.5, 1, 5];           

% Construct the Laplace K-matrix (k=0)
Kmat0 = zeros(N, N);
for i =  1:N
    yi = [y1(i); y2(i)];
    for j = 1:N
        if i == j
            % Principal value integral handled by the -1/2 on the diagonal
            Kmat0(i,j) = (rprimevec(i)^2 - rvec(i) * rbisvec(i) / 2 + rvec(i)^2 / 2) / (2 * pi * (rprimevec(i)^2 + rvec(i)^2)^1.5);;
            %Kmat0(i,j) = 0;
        else
            yj = [y1(j); y2(j)];
            dx = yi - yj;
            rdist2 = dx(1)^2 + dx(2)^2;
            dotprod = dx(1)*nu1(j) + dx(2)*nu2(j);
            % Laplace double-layer kernel
            %dotprod = dx(1)*nu1(i) + dx(2)*nu2(i);
            Kmat0(i,j) = (1/(2*pi)) * (dotprod / rdist2);
           %Kmat0(i,j) = (1/(2*pi)) * (dotprod / rdist2);
        end
    end
end

% Solve (-1/2 I + K)*h0 = nu1
% 0.5 I --> shows good convergence
% - 0.5 I --> shows bad convergence
% which one is it? Aren't we still solving exterior Neumman --> 0.5I?
A0 = 0.5*eye(N) + (2*pi/N)*Kmat0.'*diag(dsdt);
h0 = A0 \ nu1(:);

% Compute b1, b2
b1 = sum(y1(:).*h0.*dsdt(:))*(2*pi/N);
b2 = sum(y2(:).*h0.*dsdt(:))*(2*pi/N);



%%
k_values = [0.00000001, 0.0001, 0.05, 0.5, 1, 5];
% Compute |D+|
% For the given X-shape, the area is:
% |D+| = (1/2) * ∫_0^{2π} r(t)^2 dt
D_plus_area = 0.5 * sum((rvec.^2).* (2*pi/N)); 

% Choose theta_values for plotting a_0(theta)
%theta_values = linspace(0, 2*pi, 360);

% Compute a_0(theta)
a0_theta = b1*cos(theta_values) + b2*sin(theta_values) - D_plus_area;

c_const = (1 + 1i)/(4*sqrt(pi));

figure; hold on;
title('Convergence of a(\theta)/(c k^{3/2}) to a_0(\theta)');
xlabel('\theta'); ylabel('Amplitude');
grid on;

% Plot a_0(theta) as a black thick line for reference
plot(theta_values, a0_theta, 'k', 'LineWidth',1.5, 'DisplayName','a_0(\theta)');


% Loop over different wavenumbers k
for k = k_values

    % Incident plane wave from the left
    gvec = zeros(N, 1);

    for i = 1:N
        % If wave is in x-direction
        plane_wave = exp(1i * k * y1(i)); 
        gvec(i) = -1i * k * nu1(i) * plane_wave;
    end

    % Build the matrix for the integral equation
    Kmat = zeros(N, N);
    for i = 1:N
        for j = 1:N
            Kmat(i, j) = Kernel(i, j, y1, y2, nu1, nu2, rvec, rprimevec, rbisvec, k, N); 
        end
    end

    % Solve the BIE for h
    A = 0.5 * eye(N) + (2*pi/N) * Kmat.' * diag(dsdt);
    h = A \ gvec;

    % Compute a(theta)
    % ensure column vectors
    h = h(:);     
    y1 = y1(:);
    y2 = y2(:);
    dsdt = dsdt(:);

    a_theta = zeros(size(theta_values));
    for idx = 1:length(theta_values)
        theta = theta_values(idx);
        x_hat = [cos(theta); sin(theta)];

        integrand = h .* exp(-1i * k * (x_hat(1)*y1 + x_hat(2)*y2)) .* dsdt;
        integral_val = sum(integrand) * (2*pi/N);
        a_theta(idx) = (-(1+1i)/(4*sqrt(pi*k))) * integral_val;
    end

    ratio = a_theta./(c_const*(k^(3/2)));
    %ratio_real = -real(ratio) - 60 ;
    ratio_real = real(ratio);

    plot(theta_values, ratio_real, 'LineWidth', 1.5, 'DisplayName', sprintf('k=%.3g', k));

end
legend();

% we have computed a(theta) correctly. now we only have to compute
% a_0(theta)

L2_errors = zeros(length(k_values), 1);
L_inf_errors = zeros(length(k_values), 1);

% Loop over different wavenumbers k
for kk = 1:length(k_values)
    k = k_values(kk);

    
    % Incident plane wave from the left
    gvec = zeros(N, 1);

    for i = 1:N
        % If wave is in x-direction
        plane_wave = exp(1i * k * y1(i)); 
        gvec(i) = -1i * k * nu1(i) * plane_wave;
    end

    % Build the matrix for the integral equation
    Kmat = zeros(N, N);
    for i = 1:N
        for j = 1:N
            Kmat(i, j) = Kernel(i, j, y1, y2, nu1, nu2, rvec, rprimevec, rbisvec, k, N); 
        end
    end

    % Solve the BIE for h
    A = 0.5 * eye(N) + (2*pi/N) * Kmat.' * diag(dsdt);
    h = A \ gvec;

    % Compute a(theta)
    % ensure column vectors
    h = h(:);     
    y1 = y1(:);
    y2 = y2(:);
    dsdt = dsdt(:);

    a_theta = zeros(size(theta_values));
    for idx = 1:length(theta_values)
        theta = theta_values(idx);
        x_hat = [cos(theta); sin(theta)];

        integrand = h .* exp(-1i * k * (x_hat(1)*y1 + x_hat(2)*y2)) .* dsdt;
        integral_val = sum(integrand) * (2*pi/N);
        a_theta(idx) = (-(1+1i)/(4*sqrt(pi*k))) * integral_val;
    end

    ratio = a_theta./(c_const*(k^(3/2)));
    %ratio_real = -real(ratio) - 60 ;
    ratio_real = real(ratio);


    % Compute errors
    L2_errors(kk) = sqrt(sum((ratio_real - a0_theta).^2) * (2 * pi / length(theta_values)));
    L_inf_errors(kk) = max(abs(ratio_real - a0_theta));

    % Plot as before
    plot(theta_values, ratio_real, 'LineWidth', 1.5, 'DisplayName', sprintf('k=%.3g', k));
end

% Display errors
fprintf('k-values\t L2 Error\t L∞ Error\n');
for kk = 1:length(k_values)
    fprintf('%.1e\t %.4e\t %.4e\n', k_values(kk), L2_errors(kk), L_inf_errors(kk));
end

%% Plot to show the convergence of L2 and Linf

% Generate a dense set of k-values (logarithmic spacing)
k_values = logspace(-10, 0, 10);  % From 1e-10 to 1 with 50 points

% Initialize arrays for errors
L2_errors = zeros(length(k_values), 1);
L_inf_errors = zeros(length(k_values), 1);
MAE_errors = zeros(length(k_values), 1);

% Loop over different wavenumbers k
for kk = 1:length(k_values)
    k = k_values(kk);

    % Incident plane wave from the left
    gvec = zeros(N,1);
    for i = 1:N
        plane_wave = exp(1i * k * y1(i)); 
        gvec(i) = -1i * k * nu1(i) * plane_wave;
    end

    % Build the matrix for the integral equation
    Kmat = zeros(N, N);
    for i = 1:N
        for j = 1:N
            Kmat(i, j) = Kernel(i, j, y1, y2, nu1, nu2, rvec, rprimevec, rbisvec, k, N);
        end
    end

    % Solve the BIE for h
    A = 0.5 * eye(N) + (2 * pi / N) * Kmat.' * diag(dsdt);
    h = A \ gvec;

    % Compute a(theta)
    h = h(:); y1 = y1(:); y2 = y2(:); dsdt = dsdt(:);
    a_theta = zeros(size(theta_values));

    for idx = 1:length(theta_values)
        theta = theta_values(idx);
        x_hat = [cos(theta); sin(theta)];
        integrand = h .* exp(-1i * k * (x_hat(1)*y1 + x_hat(2)*y2)) .* dsdt;
        integral_val = sum(integrand) * (2 * pi / N);
        a_theta(idx) = (-(1+1i)/(4*sqrt(pi*k))) * integral_val;
    end

    % Compute ratio
    ratio = a_theta ./ (c_const * (k^(3/2)));
    ratio_real = real(ratio);

    % Compute errors
    L2_errors(kk) = sqrt(sum((ratio_real - a0_theta).^2) * (2 * pi / length(theta_values)));
    L_inf_errors(kk) = max(abs(ratio_real - a0_theta));
    MAE_errors(kk) = mean(abs(ratio - a0_theta));
end


figure;
loglog(k_values, MAE_errors, '-o', 'LineWidth', 1.5, 'DisplayName', 'MAE');
grid on;
xlabel('$k$', 'Interpreter', 'latex');
ylabel('Mean Absolute Error', 'Interpreter', 'latex');
title('Convergence of Mean Absolute Error', 'Interpreter', 'latex');
legend();

% Plot L2 and L_inf errors vs k on log-log scale
figure;
loglog(k_values, L2_errors, '-o', 'LineWidth', 1.5, 'DisplayName', 'L_2 Fel');
hold on;
loglog(k_values, L_inf_errors, '-s', 'LineWidth', 1.5, 'DisplayName', 'L_\infty Fel');
grid on;
set(gca, 'XDir', 'reverse');
xlabel('$k$', Interpreter='latex');
ylabel('Error');
title('$L_2$ och $L_\infty$ mot $k$', Interpreter='latex');
legend show;
