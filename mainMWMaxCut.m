clear all
close all
clc

%generate the graph (uniform weights on the graph)
n = 50;
p = 1/2;
[ A ] = GraphGen( n, p );

C = -A;
%renormalize
Cn = C/sum(sum(abs(C)));

cvx_begin sdp
    variable X(n,n) symmetric
    maximize( trace(Cn*X) )
    for i=1:n
        X(i,i) == 1;
    end
    X >= 0;  
cvx_end

opt = trace(Cn*X);

% cvx_begin sdp
%     variable X(n,n) symmetric
%     maximize( trace(C*X) )
%     for i=1:n
%         X(i,i) <= 1;
%     end
%     trace(X) == n;
%     X >= 0;  
% cvx_end





alphamin = 1/n;
alphamax = 1;

eps = 1e-6;

%R
R = n;

%initialize the weights
w0 = 1;
w = ones(1,n);
Nw = sum(w) + w0;
w0 = w0/Nw;
w = w/Nw;


%iteration
niter = 0;
tic
while(alphamax - alphamin >eps/2)

    
    
    %alpha
    alpha = (alphamin + alphamax)/2;

    %M is computed
    M = w0*((1/alpha)*Cn - (1/R)*eye(n));
    for i=1:n
        e = zeros(n,1);
        e(i) = 1;    
        M = M + w(i)*(-e*e' + (1/R)*eye(n));
    end

    
    delta = eps/(2*R);
    [V,D] = eig(M);
    value = max(diag(D));
    vec = V(:,n);
%     %largest eigenvalue
%     delta = eps/(2*R);
%     [vec,value] = PowerMethod(rand(n,1),M,delta);

    if (value > -delta) %feasible
       alphamin = alpha;      
    else %not feasible
       alphamax = alpha;       
    end

    %update weights
    Xt = R*vec*vec';
    beta = 1e-1;
    w0 = w0*(1 - beta*((1/alpha)*trace(Cn*Xt) - 1));
    for i=1:n
        e = zeros(n,1);
        e(i) = 1;  
        w(i) = w(i)*(1 - beta*(trace(-e*e'*Xt)) + 1);    
    end

    wn = sum(w) + w0;
    w0 = w0/wn;
    w = w/wn;
    
    wmin = min(w0,min(w));

    if wmin < 0
        warning('negative weight');
    end

niter = niter + 1;


end
toc


opt/alpha
