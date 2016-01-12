function [alpha,Xfeas]=MWLovasz(A,n,eps,beta)

% binary search bounds
alphamin = 1;
alphamax = n;

%R, delta
R = 2*n;
delta = eps/(2*R);

%iteration
niter = 0;

%Final Solution
Xfeas = zeros(n+1);

%BINARY SEARCH
while(((alphamax/alphamin) > 1 + eps) || (alphamax/alphamin) < 1 - eps)
    alpha = (alphamin + alphamax)/2;
    done = 0;
    feasible = 1;

    %initialize the weights
    w11 = 1;
    w1i_p = ones(1,n);
    w1i_n = ones(1,n);
    wii_p = ones(1,n);
    wii_n = ones(1,n);
    wij = zeros(n); %(i,j) corresponds to >= ; (j,i) corresponds to <=
     for i=1:n
        for j=i:n
            if (A(i,j) == 0)&&(i ~= j)
                wij(i,j) = 1;
                wij(j,i) = 1;
            end
        end
     end
    wnorm = w11 + sum(w1i_p) + sum(w1i_n) + sum(wii_p) + sum(wii_n) + sum(sum(wij));
    w11  = w11/wnorm;
    w1i_p = w1i_p/wnorm;
    w1i_n = w1i_n/wnorm;
    wii_p = wii_p/wnorm;
    wii_n = wii_n/wnorm;
    wij  = wij/wnorm;
    
    %Initialize result X
    X = zeros(n+1);
    time = 0;

    % RUN MW
    while(done == 0)
        %M is computed
        e0 = zeros(n+1,1);
        e0(1)  = 1;
        M = w11*(-(1/alpha)*(e0*e0') + (1/R)*eye(n+1));
        %add up the weights: equality to 1 constraints
        for i=1:n
            ei = zeros(n+1,1);
            ei(i+1) = 1;    
            M = M + (w1i_p(i) - w1i_n(i))*(ei*e0' - (1/R)*eye(n+1)) + (wii_p(i) - wii_n(i))*(ei*ei' - (1/R)*eye(n+1));
        end
        %add up the weights: equality to 0 constraints
        for i=1:n
            for j=i:n
                if (A(i,j) == 0)&&(i ~= j)
                    ei = zeros(n+1,1);
                    ej = zeros(n+1,1);
                    ei(i+1) = 1;
                    ej(j+1) = 1;
                    M = M + (wij(i,j) - wij(j,i))*(ej*ei');
                end
            end
        end
        
        [vec,eigv_max] = max_eigen(M,1);
        if (eigv_max < -delta) %infeasible
           done = 1;
           feasible = 0;
        end

        %update weights
        Xt = R*vec*vec';   
        w11 = w11*(1 - beta*(-(1/alpha)*Xt(1,1) + 1));
        for i=1:n
            w1i_p(i) = w1i_p(i)*(1 - beta*(Xt(1,i+1) - 1));
            w1i_n(i) = w1i_n(i)*(1 - beta*(-Xt(1,i+1) + 1));
            wii_p(i) = wii_p(i)*(1 - beta*(Xt(i+1,i+1) - 1));
            wii_n(i) = wii_n(i)*(1 - beta*(-Xt(i+1,i+1) + 1));
        end
        for i=1:n
            for j=i:n
                if (A(i,j) == 0)&&(i ~= j)
                    wij(i,j) = wij(i,j)*(1 - beta*( Xt(i+1,j+1) ));
                    wij(j,i) = wij(j,i)*(1 + beta*( Xt(i+1,j+1) ));
                end
            end
        end
        %normalisation
        wnorm = w11 + sum(w1i_p) + sum(w1i_n) + sum(wii_p) + sum(wii_n) + sum(sum(wij));
        w11 = w11/wnorm;
        w1i_p = w1i_p/wnorm;
        w1i_n = w1i_n/wnorm;
        wii_p = wii_p/wnorm;
        wii_n = wii_n/wnorm;
        wij = wij/wnorm;

        % Sanity check on updated weights
        wmin = min([w11,min(min(wij)),min(w1i_p),min(w1i_n),min(wii_p),min(wii_n)]);

        if wmin < 0
            warning('negative weight');
        end

        % update X, time
        time = time + 1;
        X = X + Xt;
        X_T = (1/time) * X;

        % check if we have an epsilon-feasible X
        ef = 1;
        %constraints
        for i=1:n
            for j=i:n
                if (A(i,j) == 0)&&(i ~= j)
                    if (X_T(i+1,j+1) < -eps)||(X_T(i+1,j+1) > eps)
                        ef = 0;
                    end
                end
            end
        end
        for i=1:n
            if (X_T(i+1,i+1) < 1-eps)||(X_T(i+1,i+1) > 1+eps)
                ef = 0;
            end
            if (X_T(1,i+1) < 1-eps)||(X_T(1,i+1) > 1+eps)
                ef = 0;
            end
        end
        if X_T(1,1) > alpha + eps
            ef = 0;
        end
        
        %are we done?
        if ef == 1
            done = 1;
            feasible = 1;
        end

    end
    
    % Print Results of this iteration, update alpha
    if (feasible) %feasible
       fprintf('Iteration %d is feasible (%d MW updates, alpha = %f4.4)\n',niter,time,alpha)
       alphamax = alpha;
       Xfeas = X_T;
    else %not feasible
       fprintf('Iteration %d is NOT feasible (%d MW updates, alpha = %f4.4)\n',niter,time,alpha)
       alphamin = alpha;       
    end

    niter = niter + 1;

end