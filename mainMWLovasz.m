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



% %%
% %MW method
% alphamin = 1;
% alphamax = n;
% %R
% r = 2*n;
% 
% %define weights
% w0 = 1;
% wp1 = ones(1,n);
% wp2 = ones(1,n);
% wp3 = ones(n,n);
% wn1 = ones(1,n);
% wn2 = ones(1,n);
% wn3 = ones(n,n);
% 
% 
% %LOOP
% 
% %define alpha
% alpha = (alphamin + alphamax)/2;
% 
% %builds matrix C
% 
% 
% %compute largest eigenvector
% C = rand(n+1,n+1);
% toler = 1e-5;
% start = rand(n+1,1);
% [vec,value]=PowerMethod(start,C,toler);
% x = vec/value;
% 
% %update weights
% 
% %if feasible
% 
% %if not feasible
% 
% %end LOOP
% 





