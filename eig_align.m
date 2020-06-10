function [EA] = eig_align(V,n)
% [EA] = eig_align(V,n) assess eigenvector alignment between a selected
% node n and all other network nodes, using the system's dominant
% eigenvectors V.
    S1=norm(V(n,:));
    S2=vecnorm(V,2,2);
    Stemp=dot(repmat(V(n,:),size(V,1),1),V,2);
    EA = real(acos(Stemp./(S1*S2)));
end
    
 