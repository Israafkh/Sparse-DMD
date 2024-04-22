% Applying sparse dmd on v velocity component of a heated cylinder flow
v_velocity_comp = ncread("boussinesq.nc","v");
nx = 150;
ny = 450;
velocity = zeros(nx*ny,150);
for i=1:150
    num = i*10;
    v = v_velocity_comp(:,:,num);
    velocity(:,i) = reshape(v,nx*ny,1);
   
end
X = velocity(:, 1:end-1);
Y = velocity(:, 2:end);
% Initialize vectors
gamma_vector = logspace(log10(10),log10(70),60);
nz_values = zeros(1, length(gamma_vector));
performance_values = zeros(1, length(gamma_vector));
% Iterate over gamma
for i = 1:length(gamma_vector)
    gamma = gamma_vector(i);
    [nz,~,~,~,~,performance,~,~,~,~] = Sparse_DMD(X, Y, gamma);
    nz_values(i) = nz;
    performance_values(i) = performance;
    
end

% Plotting non-zero elements vs. gamma
figure;
semilogx(gamma_vector, nz_values, 'o-', 'LineWidth', 2);
xlabel('Gamma');
ylabel('Non-zero elements');
title('Non-zero elements vs. Gamma');
grid on;

% Plotting performance vs. gamma
figure;
semilogx(gamma_vector, performance_values, 'o-', 'LineWidth', 2);
xlabel('Gamma');
ylabel('Performance');
title('Performance vs. Gamma');
grid on;

% Plotting performance vs. non-zero elements
figure
plot(nz_values, performance_values, 'o-', 'LineWidth', 2);
xlabel('Non-zero elements');
ylabel('Performance');
title('Performance vs. Non-zero elements');
grid on;
 
% Choosing a specific gamma to plot dmd spectrum and modes
[nz,a_opt,phi,lambda,lambda_chosen,performance,nonzero_indices,V_r,S_r,W_r] = Sparse_DMD(X, Y, 700);

% Plot DMD spectrum
figure
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'k--') % plot unit circle
hold on, grid on

% Plot lambda with circles
scatter(real(lambda),imag(lambda),'ok')

% Plot lambda_chosen with crosses
scatter(real(lambda_chosen),imag(lambda_chosen),'xk')

axis([-1.1 1.1 -1.1 1.1]);
legend({'Unit Circle', 'Lambda', 'Lambda\_chosen'}, 'Location', 'best');
% Plot DMD modes
phi1 = phi(:,nonzero_indices);
figure;
for i = 1:1:nz
    residual = norm(Y*V_r/S_r*W_r(:,nonzero_indices(i))-lambda_chosen(i)*phi1(:,i),2);
    
    % Plot real part of phi
    figure;
    imagesc(reshape(real(phi1(:,i)), nx, ny));
    title(['Real Part of \phi - Residual: ' num2str(residual)]);
    colorbar; % Add colorbar for reference
    text(5, 5, ['Residual: ' num2str(residual)], 'Color', 'white', 'FontSize', 12); % Add residual as text
    
    % Plot imaginary part of phi
    figure;
    imagesc(reshape(imag(phi1(:,i)), nx, ny));
    title(['Imaginary Part of \phi - Residual: ' num2str(residual)]);
    colorbar; % Add colorbar for reference
    text(5, 5, ['Residual: ' num2str(residual)], 'Color', 'white', 'FontSize', 12); % Add residual as text
    
    pause(0.5); % Add a delay between plotting to see the images
end

