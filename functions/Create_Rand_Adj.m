% create (n) matrices of uniformly distributed random values between 0 and 1
function [Z] = Create_Rand_Adj(n)
    if ~exist('n','var')
        n=29;
    end
    sz=132;
    Z=zeros(sz,sz,n);
    for i = 1 : n
        adj=rand(sz,sz);            % populate matrix with uniformly distributed random values
        adj=adj-diag(diag(adj));    % set diagonal values to 0
        adj=(adj+adj');             % symmetrise
        adj=rescale(adj,0,1);       % ensure range between 0 and 1
        Z(:,:,i)=adj;
    end
end
