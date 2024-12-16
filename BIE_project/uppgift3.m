% Number of boundary points
N = 4000;

% Far-field directions
theta_values = linspace(0, 2*pi, 360); 

% Chosen k-values to verify the optic theorem in R^2
k_values = [0.5, 5];           

% Domain discretization of smooth X-shape
tvec = linspace(-pi, pi, N);
rvec = 3 + cos(4 * tvec + pi);
rprimevec = -4 * sin(4 * tvec + pi);
rbisvec = -16 * cos(4 * tvec + pi);
y1 = rvec .* cos(tvec);
y2 = rvec .* sin(tvec);
dsdt = sqrt(rprimevec.^2 + rvec.^2);
nu1 = (rvec .* cos(tvec) + rprimevec .* sin(tvec)) ./ dsdt;
nu2 = (rvec .* sin(tvec) - rprimevec .* cos(tvec)) ./ dsdt;

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

    % Plot Real and Imaginary Parts of a(theta)
    figure;
    subplot(3,1,1);
    plot(theta_values, real(a_theta), 'LineWidth', 2);
    xlabel('$\theta$', 'Interpreter', 'latex'); 
    ylabel('$Re(a(\theta))$', 'Interpreter','latex');
    title(['Real Part of $a(\theta)$, k=', num2str(k)], 'Interpreter','latex');
    grid on;
    
    subplot(3,1,2);
    plot(theta_values, imag(a_theta), 'LineWidth', 2);
    xlabel('$\theta$',Interpreter='latex'); 
    ylabel('$Im(a(\theta))$', Interpreter='latex');
    title(['Imaginary Part of $a(\theta), k=$', num2str(k)], Interpreter="latex");
    grid on;
    
    % Plot |a|^2
    subplot(3,1,3);
    plot(theta_values, abs(a_theta).^2, 'LineWidth', 2);
    xlabel('$\theta$', Interpreter='latex'); 
    ylabel('$|a(\theta)|^2$', Interpreter = 'latex');
    title(['$|a(\theta)|^2, k=$', num2str(k)], Interpreter = 'latex');
    grid on;


    % Verify (6.10)
    dtheta = theta_values(2) - theta_values(1); 
    lhs = sum(abs(a_theta).^2) * dtheta;
    a_0 = a_theta(1);  % assuming theta_values(1) = 0
    rhs = imag(2*(1 - 1i)*sqrt(pi/k)*a_0);

    % Display results
    fprintf('\nVerifying (6.10) for k = %g:\n', k);
    fprintf('LHS = ∫_0^{2π} |a(θ)|^2 dθ = %.6e\n', lhs);
    fprintf('RHS = Im(2(1-i)*√(π/k)*a(0)) = %.6e\n', rhs);
    fprintf('Difference = %.6e\n', abs(lhs - rhs));
end
