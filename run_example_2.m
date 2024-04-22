% Applying sparse dmd on vorticity data of a flow around cylinder
clear all; close all; clc;
load /Users/israafakih/Downloads/DATA/FLUIDS/CYLINDER_ALL.mat;
X = VORTALL(:, 1:end-1);
Y = VORTALL(:, 2:end);
%initializing vectors
gamma_vector = logspace(log10(0.005),log10(50),15);
nz_values = zeros(1, length(gamma_vector));
performance_values = zeros(1, length(gamma_vector));
%Iterating over gamma
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
figure;
plot(nz_values, performance_values, 'o-', 'LineWidth', 2);
xlabel('Non-zero elements');
ylabel('Performance');
title('Performance vs. Non-zero elements');
grid on;
% Choosing a specific gamma to plot the spectrum and modes
[nz,a_opt,phi,lambda,lambda_chosen,performance,nonzero_indices,V_r,S_r,W_r] = Sparse_DMD(X, Y,log10(20) );

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
%% Plot DMD modes
phi = phi(:,nonzero_indices);
for i=1:3:nz
    residual = norm(Y*V_r/S_r*W_r(:,nonzero_indices(i))-lambda_chosen(i)*phi(:,i),2);
    plotCylinder(reshape(real(phi(:,i)),nx,ny),nx,ny);
    plotCylinder(reshape(imag(phi(:,i)),nx,ny),nx,ny);
    %disp(residual);
end


