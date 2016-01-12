clear all
close all
clc

%generate the graph
n = 10;
p = 1/2;
[ A ] = GraphGen( n, p );

% %%
% %original lovasz
% cvx_begin sdp
%     variable Z(n,n) symmetric
%     maximize( trace(ones(n)*Z) )
%     trace( Z ) == 1;
%     for i=1:n
%         for j=1:n
%             if (A(i,j) == 1)
%                 Z(i,j) == 0;
%             end
%         end
%     end
%     Z >= 0;  
% cvx_end

%%
%dual
cvx_begin sdp
    variable Z(n+1,n+1) symmetric
    minimize( Z(1,1) )
    for i=2:n+1
        Z(i,i) == 1;
        Z(1,i) == 1;
    end    
    for i=1:n
        for j=1:n
            if (A(i,j) == 0) && (i ~= j)
                Z(i+1,j+1) == 0;
            end
        end
    end    
    Z >= 0;
cvx_end

% MW Algorithm
beta = 0.01;
tic;
[opt_mw,X_mw] = MWLovasz(A,n,0.01,beta);
toc;

fprintf('Optimum From CVX: %5.5f\n',cvx_optval)
fprintf('Optimum From MW Algorithm: %5.5f\n',opt_mw)
fprintf('Ratio of CVX/MW: %5.5f\n',cvx_optval/opt_mw)

