function [alpha]=MWLovasz(A,n,eps)


alphamin = 1;
alphamax = n;

%R
R = 2n;

%initialize the weights
w0 = 1;
w1p = ones(1,n);
w1n = ones(1,n);
w2p = ones(1,n);
w2n = ones(1,n);
w3p = ones(1,n*(n-1)/2);
w3n = ones(1,n*(n-1)/2);


%iteration
niter = 0;
tic
while(alphamax - alphamin >eps/2)

    
    
    %alpha
    alpha = (alphamin + alphamax)/2;

%     %M is computed
%     M = w0*((1/alpha)*W - (1/R)*eye(n));
%     for i=1:n
%         e = zeros(n,1);
%         e(i) = 1;    
%         M = M + w(i)*(-e*e' + (1/R)*eye(n));
%     end
% 
%     
%     delta = eps/(2*R);
%     [V,D] = eig(M);
%     value = max(diag(D));
%     vec = V(:,n);
% %     %largest eigenvalue
% %     delta = eps/(2*R);
% %     [vec,value] = PowerMethod(rand(n,1),M,delta);
% 
%     if (value > -delta) %feasible
%        fprintf('Iteration %d is feasible\n',niter)
%        alphamin = alpha;      
%     else %not feasible
%        fprintf('Iteration %d is NOT feasible\n',niter)
%        alphamax = alpha;       
%     end
% 
%     %update weights
%     Xt = vec*vec';
%     beta = 1e-1;
%     w0 = w0*(1 - beta*((1/alpha)*trace(W*Xt) - 1));
%     for i=1:n
%         e = zeros(n,1);
%         e(i) = 1;  
%         w(i) = w(i)*(1 - beta*(trace(-e*e'*Xt) + 1));    
%     end
% 
%     wn = sum(w) + w0;
%     w0 = w0/wn;
%     w = w/wn;
%     
%     % Sanity check on updated weights
%     wmin = min(w0,min(w));
% 
%     if wmin < 0
%         warning('negative weight');
%     end

niter = niter + 1;

end
toc