function [V] = ordered_eigvecs(M,n_eigs)
% Function to evaluate eigenvectors and order them from largest to smallest
% eigenvalue in magnitude.
    [V,D]=eig(M');
    [~,I]=sort(real(diag(D)),'desc');
    V=V(:,I(1:n_eigs)); % top n_eigs eigenvectors according to real eigenvalue magnitude
    V(:,1) = abs(V(:,1)); % make V positive
end
    