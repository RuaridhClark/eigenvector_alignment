function [EA] = Eig_Align(V,n)
% [EA] = eig_align(V,n) assess eigenvector alignment between a selected
% node n and all other network nodes, using the system's dominant
% eigenvectors V.
    V(:,1) = ~any(V,2)*1e-15 + V(:,1);  % give nodes at origin small v_1 value
    S1=norm(V(n,:));
    S2=vecnorm(V,2,2);
    Stemp=dot(repmat(V(n,:),size(V,1),1),V,2);
    EA = real(acos(Stemp./(S1*S2)));
end
    
 