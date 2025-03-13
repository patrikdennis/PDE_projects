function area_dplus = compute_area(y1, y2, dsdt, nu1, nu2 , rvec, rprimevec, rbisvec, k_small, c_const, theta_values,N)
    % Compute a(theta) for k_small using the Helmholtz BIE solution
    % You need to implement or have a function `compute_a_theta` that:
    %  - solves Helmholtz BIE for density h
    %  - computes a(theta)

    
    % Incident plane wave from the left
    gvec = zeros(N, 1);
    
    for i = 1:N
        % If wave is in x-direction
        plane_wave = exp(1i * k_small * y1(i)); 
        gvec(i) = -1i * k_small * nu1(i) * plane_wave;
    end
    
    % Build the matrix for the integral equation
    Kmat = zeros(N, N);
    for i = 1:N
        for j = 1:N
            Kmat(i, j) = Kernel(i, j, y1, y2, nu1, nu2, rvec, rprimevec, rbisvec, k_small, N); 
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
    
    a_theta_small_k = zeros(size(theta_values));
    for idx = 1:length(theta_values)
        theta = theta_values(idx);
        x_hat = [cos(theta); sin(theta)];
    
        integrand = h .* exp(-1i * k_small * (x_hat(1)*y1 + x_hat(2)*y2)) .* dsdt;
        integral_val = sum(integrand) * (2*pi/N);
        a_theta_small_k(idx) = (-(1+1i)/(4*sqrt(pi*k_small))) * integral_val;
    end
    
    % Integrate a(theta) over -π to π.
    dtheta = theta_values(2)-theta_values(1);
    I =  sum(a_theta_small_k)*dtheta; 
    
    % |D+| = -(1/(2π))*(1/(c*k_small^(3/2)))*I
    area_dplus = -(1/(2*pi))*(I/(c_const*(k_small^(3/2))));
    area_dplus = real(area_dplus);

end

