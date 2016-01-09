function [alpha,niter]=MWLovasz(A,n,eps)

%%
%parameter for the update
beta = 1e-1;
%bound on the optimum
alphamin = 1;
alphamax = n;
%R
R = 2*n;
%eps
delta = eps/(2*R);

%%
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

%%
%iteration
niter = 0;
while(alphamax - alphamin >eps/2)
   
    %alpha
    alpha = (alphamin + alphamax)/2;

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


    
    %matlab method... to be changed
    [V,D] = eig(M);
    value = max(diag(D));
    vec = V(:,n);


    if (value > -delta) %feasible
        alphamax = alpha;
    else %not feasible
        alphamin = alpha;
    end

    %update weights
    Xt = vec*vec';   
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

niter = niter + 1;

end
