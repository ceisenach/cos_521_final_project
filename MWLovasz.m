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
Xfeas = zeros(n);

%BINARY SEARCH
while(((alphamax/alphamin) > 1 + eps) || (alphamax/alphamin) < 1 - eps)
    alpha = (alphamin + alphamax)/2;
    done = 0;
    feasible = 1;

    %initialize the weights
    w0 = 1;
    w1p = ones(1,n);
    w1n = ones(1,n);
    w2p = ones(1,n);
    w2n = ones(1,n);
    w3 = zeros(n);
     for i=1:n
        for j=i:n
            if (A(i,j) == 0)&&(i ~= j)
                w3(i,j) = 1;
                w3(j,i) = 1;
            end
        end
     end
    wnorm = w0 + sum(w1p) + sum(w1n) + sum(w2p) + sum(w2n) + sum(sum(w3));
    w0  = w0/wnorm;
    w1p = w1p/wnorm;
    w1n = w1n/wnorm;
    w2p = w1p/wnorm;
    w2n = w1n/wnorm;
    w3  = w3/wnorm;
    
    %Initialize result X
    X = zeros(n);
    time = 0;

    % RUN MW
    while(done == 0)
        %M is computed
        e0 = zeros(n+1,1);
        e0(1)  = 1;
        M = w0*(-(1/alpha)*e0*e0' + (1/R)*eye(n+1));
        %add up the weights: equality to 1 constraints
        for i=1:n
            e  = zeros(n+1,1);
            e(i+1) = 1;    
            M = M + (w1p(i) - w1n(i))*(e*e0' - (1/R)*eye(n+1)) + (w2p(i) - w2n(i))*(e*e' - (1/R)*eye(n+1));
        end
        %add up the weights: equality to 0 constraints
        for i=1:n
            for j=i:n
                if (A(i,j) == 0)&&(i ~= j)
                    ei  = zeros(n+1,1);               
                    ej  = zeros(n+1,1);
                    ei(i+1) = 1;
                    ej(j+1) = 1;
                    M = M + (w3(i,j) - w3(j,i))*(ej*ei');
                end
            end
        end
        
        [vec,eigv_max] = max_eigen(M,2);
        if (eigv_max < -delta) %infeasible
           done = 1;
           feasible = 0;
        end

        %update weights
        Xt = R*vec*vec';   
        w0 = w0*(1 - beta*(-(1/alpha)*Xt(1,1) + 1));
        for i=1:n
            w1p(i) = w1p(i)*(1 - beta*(Xt(1,i+1) - 1));
            w1n(i) = w1n(i)*(1 - beta*(-Xt(1,i+1) + 1));
            w2p(i) = w1p(i)*(1 - beta*(Xt(i+1,i+1) - 1));
            w2n(i) = w1n(i)*(1 - beta*(-Xt(i+1,i+1) + 1));
        end
        for i=1:n
            for j=i:n
                if (A(i,j) == 0)&&(i ~= j)
                    w3(i,j) = w3(i,j)*(1 - beta*( Xt(i,j) ));
                    w3(j,i) = w3(j,i)*(1 + beta*( Xt(i,j) ));
                end
            end
        end
        %normalisation
        wnorm = w0 + sum(w1p) + sum(w1n) + sum(w2p) + sum(w2n) + sum(sum(w3));
        w0  = w0/wnorm;
        w1p = w1p/wnorm;
        w1n = w1n/wnorm;
        w2p = w1p/wnorm;
        w2n = w1n/wnorm;
        w3  = w3/wnorm;

        % Sanity check on updated weights
        wmin = min([w0,min(min(w3)),min(w1p),min(w1n),min(w2p),min(w2n)]);

        if wmin < 0
            warning('negative weight');
        end

        % update X, time
        time = time + 1;
        X = X + Xt;
        X_T = (1/time) * X;
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FIX THIS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check if we have an epsilon-feasible X
        ef = 1;
        for i=1:n
            if 1 - X_T(i,i) < -eps
                ef = 0;
            end
        end
        if (1/alpha)*trace(W*X_T) - 1 < - eps
            ef = 0;
        end
        if ef == 1
            done = 1;
            feasible = 1;
        end

    end
    
    % Print Results of this iteration, update alpha
    if (feasible) %feasible
       fprintf('Iteration %d is feasible (%d MW updates)\n',niter,time)
       alphamin = alpha;
       Xfeas = X_T;
    else %not feasible
       fprintf('Iteration %d is NOT feasible (%d MW updates)\n',niter,time)
       alphamax = alpha;       
    end

    niter = niter + 1;

end