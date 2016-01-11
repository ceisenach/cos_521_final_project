% Return maximum eigenvector and eigenvalue
% mode == 1 means it uses naive eigendecomposition
% mode == 2 means it uses lanczsos method

function [eivec,eival] = max_eigen(A,mode)

if mode == 1
    [V,D] = eig(A);
    eival = max(diag(D));
    idx = find(diag(D) == eival);
    eivec = V(:,idx);
elseif mode == 2
    [U,V,D] = lanczos(A,5);
    eival = max(diag(D));
    idx = find(diag(D) == eival);
    eivec = V*U(:,idx);
end