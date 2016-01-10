function [alpha,Xfeas]=MWMaxCut(W,n,eps, beta)

% binary search bounds
alphamin = 1/n;
alphamax = 1;

%R, delta
R = n;
delta = eps/(2*R);

%iteration
niter = 0;

%Final Solution
Xfeas = zeros(n);

%BINARY SEARCH
while(alphamax - alphamin >eps/2)
    alpha = (alphamin + alphamax)/2;
    done = 0;
    feasible = 1;

    %initialize the weights
    w0 = 1;
    w = ones(1,n);
    Nw = sum(w) + w0;
    w0 = w0/Nw;
    w = w/Nw;
    
    %Initialize result X
    X = zeros(n);
    time = 0;

    % RUN MW
    while(done == 0)
        %M is computed
        M = w0*((1/alpha)*W - (1/R)*eye(n));
        for i=1:n
            e = zeros(n,1);
            e(i) = 1;    
            M = M + w(i)*(-e*e' + (1/R)*eye(n));
        end
        
        [V,D] = eig(M);
        eigv_max = max(diag(D));
        idx = find(diag(D) == eigv_max);
        vec = V(:,idx);
        if (eigv_max < -delta) %infeasible
           done = 1;
           feasible = 0;
        end

        %update weights
        Xt = R*vec*vec';
        w0 = w0*(1 - beta*((1/alpha)*trace(W*Xt) - 1));
        for i=1:n
            e = zeros(n,1);
            e(i) = 1;  
            w(i) = w(i)*(1 - beta*(trace(-e*e'*Xt) + 1));    
        end
        wn = sum(w) + w0;
        w0 = w0/wn;
        w = w/wn;
        
        % Sanity check on updated weights
        wmin = min(w0,min(w));
        if wmin < 0
            warning('negative weight');
        end

        % update X, time
        time = time + 1;
        X = X + Xt;
        X_T = (1/time) * X;
    
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