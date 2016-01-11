% Adapted from  code Written by Danny Bickson, CMU
% Matlab code for running Lanczos algorithm
% Original Code available from: http://www.cs.cmu.edu/~bickson/gabp/

% m is number of eigenvalues desired.
% U,V can be used to reconstruct the Ritz eigenvectors
% U are the eigenvectors of the tridiagonal matrix and
% V is the transformation matrix
% D are the eigenvalues

function [U,V,D] = lanczos(A, m)

[n,k] = size(A);
V = zeros(k,m+1);
V(:,2) = rand(k,1);
V(:,2)=V(:,2)/norm(V(:,2),2);
beta(2)=0;

for j=2:m+2

    w = A*V(:,j);
    alpha(j) = w'*V(:,j);
    w = w - beta(j)*V(:,j-1);
    w = w - alpha(j)*V(:,j);

    %orthogonalize
%    for k=2:j-1
%      tmpalpha = w'*V(:,k);
%      w = w -tmpalpha*V(:,k);
%    end

    beta(j+1) = norm(w,2);
    V(:,j+1) = w/beta(j+1);
end

% prepare the tridiagonal matrix
T = diag(alpha(2:end)) + diag(beta(3:end-1),1) + diag(beta(3:end-1),-1);
% throw out extraneous part of V
V = V(:,2:end-1);
[U,D]=eig(full(T));