%Applying sparse DMD on synthetic example
%Generate a random matrix A
n=1000;
m=400;
A=randn(n,n);
A=A/norm(A,2);
temp=randn(n,1);
%Construct the snapshot matrices
X=zeros(n,m);
Y=zeros(n,m);
X(:,1)=temp;
for i=2:m+1
    temp=A*temp;
    if i~= m+1
        X(:,i)=temp;
    end
    Y(:,i-1)=temp;
end
%Initialize variables
gamma_vector = logspace(log10(0.5),log10(2),40);
nz_values = zeros(1, length(gamma_vector));
performance_values = zeros(1, length(gamma_vector));
%Iterate over gamma values 
for i = 1:length(gamma_vector)
    gamma = gamma_vector(i);
    [nz,~,~,~,~,performance,~,~,~,~] = Sparse_DMD(X, Y, gamma);
    nz_values(i) = nz;
    performance_values(i) = performance;
end
%Plotting nonzero values vs gamma
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

%Choosing a certain gamma to plot the DMD spectrum
% Plot DMD spectrum
[nz,~,~,lambda,lambda_chosen,~,~,~,~,~] = Sparse_DMD(X, Y, 2);
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
