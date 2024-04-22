function [nz,a_opt,phi,lambda,lambda_chosen,performance,nonzero_indices,V_r,S_r,W_r] = Sparse_DMD(X, Y, gamma,rho, max_iter, eps_abs, eps_rel)
    
    % Check if any optional argument is not provided and assign default values
    if nargin < 7
        eps_rel = 1.e-4;
        if nargin < 6
            eps_abs = 1.e-6;
            if nargin < 5
                max_iter = 10000;
                if nargin < 4
                    rho = 1;
                end
            end
        end
    end
    
    % Getting the SVD of the matrix X
    [U, S, V] = svd(X, 'econ');
    r=0;
    sz=size(Y,1);
    for i=1:length(S)

       if S(i,i)>sz*S(1,1)*eps
            r=r+1;
       end
    end
    % truncate to rank-r
    U_r = U(:, 1:r); 
    S_r = S(1:r, 1:r);
    V_r = V(:, 1:r);
    % Forming the DMD matrix
    F_DMD = (U_r)'*Y*V_r;
    for i=1:r
    F_DMD(:,i) = F_DMD(:,i)/S_r(i,i);
    end
    %Calculating eigenvales and eigenvectors of F_DMD
    [W_r, D] = eig(F_DMD);

    lambda = diag(D);
    
    % Calculating DMD modes
    phi =  Y * V_r / S_r *W_r;
    
    % Calculating the vandermonde matrix determined by lambda
    V_mode = zeros(r,size(V,1));
    for i = 1:size(X,2)
    
    V_mode(:,i) =lambda.^(i-1);
    
    end

    % Calculating the components P, q and s
    component_1 = (W_r)'*W_r;
    component_2 = conj(V_mode*(V_mode)');
    P = component_1 .* component_2;
    q = conj(diag(V_mode*V_r*(S_r)'*W_r));
    s = trace((S_r*(V_r)')'*(S_r*(V_r)'));

    % we start now by the iterations for the sparsity
    %use ADMM for the augmented Lagrangian minimization
    n = length(q);
    I = eye(n);
    %Initialize the variables beta and lambda
    beta = zeros(n,1);
    lam = zeros(n,1);
    %Compute the iterates 
    for step = 1:max_iter
    alpha = (P+(rho/2)*I)\(q+(rho/2)*(beta - (1/rho)*lam));
    kappa = (gamma/rho)*ones(n,1);
    v = alpha + (1/rho)*lam;
    beta_new = ((1-kappa./abs(v)).*v).*(abs(v)>kappa);
    lam = lam + rho*(alpha-beta_new);

    % computing epsilon_prime and epsilon_dual
    eps_prim = sqrt(n) * eps_abs + eps_rel * max([norm(alpha),norm(beta_new)]);
    eps_dual = sqrt(n) * eps_abs + eps_rel * norm(lam);

    % computing the two residuals res_prim and res_dual
    res_prim = norm(alpha-beta_new);
    res_dual = rho*norm(beta_new-beta);
    
    % Checking the stopping criteria
        if (res_prim<eps_prim) && (res_dual<eps_dual)
            
            break;
        end
        beta = beta_new;
    end
    % Recording the results of the first step
    nonzero_indices = find(abs(beta) >1.e-12 );
    nz = length(nonzero_indices);
    a_sp = beta;
    
    % now after having the sparsity we need to find the amplitudes after polishing 
    % step1 is to compute the constraint matrix E
    ind_zero = find( abs(beta) < 1.e-12);
    m = length(ind_zero);
    E = I(:,ind_zero);
    E = sparse(E); 
    % step2 is to compute the solution a_pol
    a_opt = [P,E;E',zeros(m,m)]\[q;zeros(m,1)];
    a_opt = a_opt(1:n);
    
    % Computing the residual and performance of the computed amplitudes
    res = real(a_opt'*P*a_opt) - 2*real(q'*a_opt) + s;
    performance = 100*sqrt(res/s);
    
    %Computing the resulting choosen modes
    lambda_chosen = lambda(nonzero_indices);
    
    end

    
