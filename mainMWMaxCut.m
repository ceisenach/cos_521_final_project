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
        X(i,i) <= 1;
    end
    X >= 0;  
cvx_end

beta = 0.01;
tic;
[opt_mw,X_mw] = MWMaxCut(Cn,n,0.05, beta);
toc;

fprintf('Optimum From CVX: %5.5f\n',cvx_optval)
fprintf('Optimum From MW Algorithm: %5.5f\n',opt_mw)
fprintf('Ratio of CVX/MW: %5.5f\n',cvx_optval/opt_mw)