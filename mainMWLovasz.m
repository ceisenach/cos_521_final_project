clear all
close all
clc

%generate the graph
n = 50;
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

%%%%%%%%%%%
% study the dependence on beta
%precision
eps = 1e-6;
%parameter for the update
Alpha = zeros(1,20); 
for k=1:length(Alpha)
    
    beta = (2/3)^k;
    [alpha,niter]=MWLovasz(A,n,eps,beta);
    Alpha(k) = alpha;
end





